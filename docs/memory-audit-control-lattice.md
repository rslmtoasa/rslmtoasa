# Memory audit: control and lattice (first pass)

Base: 49c0031b6406a93cabc768a9d698796c1ed90fd9.
Branch: fix/control-lattice-memory-cleanup.
Scope: static inspection of pre_processing_bravais, control construction and
lattice construction/build_data/bravais/structb/atomlist. This is the first
cleanup pass, not a certification of every calculation mode.

## Workflow and ownership

| Stage | Allocation and lifetime |
| --- | --- |
| control(fname) | No active allocatable or pointer components. Defaults feed namelist input and later array sizes. |
| lattice(control_obj) | Borrows control through a pointer; the target in calculation.f90 outlives lattice operations. Owns its allocatable components. |
| restore_to_default | Allocates izp/no(ndim), crd(3,ndim), plus initially zero-sized inclu, ijpair, ijktrio, chargetrf_type and ct. ndim defaults to 9,900,000. |
| build_from_file | Moves arrays to namelist locals, reads, resizes, rereads and moves them back. Temporary capacity for ct is ndim, despite the comment mentioning 1000. |
| build_data | Populates basis data and allocates ib/iu/irec for built-in structures, or reads lattice.nml. |
| bravais | Allocates iz/num(ndim), cr/crbravais(3,ndim). Moves cr/iz/num into the object; crbravais is local. Retains capacity ndim after cutting to kk atoms. |
| structb | Allocates local nn(kk,5250), then the retained neighbour map and sbar. Calls nncal, remd, outmap and dbar1/clusba. |
| atomlist | Allocates atlist(ntype), ham_i(kk) and symbolic_atoms(ntype), whose construction delegates to symbolic_atom. |
| Scope exit | Allocatable components have automatic lifetime management. Absence from the handwritten destructor alone does not prove a heap leak. Allocation-tracker accounting is a separate concern. |

## Fixes in this pass

1. Control initializes nlim, nsp, nmdir and asd_atom. An omitted nsp now has an
   explicit invalid sentinel and is rejected; valid values remain 1 through 4.
   Negative nlim is rejected.
2. Control's optional filename fallback is actually used by OPEN and its error
   message. print_state initializes all namelist fields, obtains a NEWUNIT and
   closes only the file it opened.
3. Lattice's borrowed control pointer starts disassociated. Defaults now define
   alat, rc, r2, a, crystal_sym and no. build_from_file preserves alat before input.
4. The ijktrio resize condition checks its actual (njijk,6) shape rather than
   comparing its total size with 2*njijk.
5. build_from_lattice moves ct into the local namelist variable before reading
   and back afterwards. Previously a file containing ct could target an
   unallocated local.
6. fcc2/hcp initialize only izp(1:2) and no(1:2). The former full-array-section
   assignments had incompatible shapes for ndim other than two.
7. nncal checks both rows BEFORE inserting a neighbour. The former overflow
   check occurred after the invalid write. At least the count column is retained
   when no pairs are found.
8. nncal/mapa accept the actual CT storage through assumed-size dummies rather
   than declaring 50 entries when the caller allocates ntype entries.
   The existing cutoff rule using CT(1) is unchanged.
9. structb copies only nm valid columns and zeroes its spare column, avoiding
   reading column 5251 when nm reaches 5250. It releases the large temporary nn
   immediately, sizes remd workspace to the retained map, and supplies the actual
   width to remd/outmap. Previously their explicit-shape dummy arrays claimed
   5250 columns even when the actual array was smaller.
10. structb releases set/idnn before dbar1 starts its dense workspace allocation.

## Memory effect

With four-byte default integers, releasing nn removes kk*5250*4 bytes from
the live allocations during later structure-constant construction (about
200 MiB at kk=10000). This is a lifetime calculation, not a measured RSS
reduction; allocator behaviour and other workspaces determine process RSS.

set now has 3*ntot*(nm+1) real elements instead of 3*ntot*5250 and is released
before dbar1. The retained map's shape and sbar dimensions are preserved.
The initial nn search capacity is still 5250.

## Open findings and next steps

- **Constructor memory is still large.** At four-byte integers/eight-byte reals,
  izp+no+crd alone occupy about 302 MiB for default ndim; temporary ct adds about
  75.5 MiB. bravais adds two integer and two coordinate arrays of that capacity.
  Replace read/resize/reread with a robust sizing strategy before changing these
  capacities. A smaller guessed capacity would reproduce namelist bounds errors.
- **The first namelist read is not a safe sizing pass.** Zero-length pair/trio/
  impurity arrays and insufficient user-dependent capacity can fail before all
  dimensions are read. The first IOSTAT is not handled. Resized arrays also need
  deterministic defaults for omitted entries. More lattice namelist scalars are
  declared than are initialized or transferred by the reader.
- **build_from_lattice still has dimension assumptions.** Its ib/iu/irec sizing
  read needs the same systematic solution; this pass fixes only CT ownership.
- **Reinitialization is not generally safe.** restore_to_default, build_data,
  structb and atomlist contain unconditional allocations. Establish the public
  lifecycle contract before adding blanket deallocations that might invalidate
  downstream state.
- **Allocation accounting is mixed.** USE_SAFE_ALLOC coexists with raw allocate,
  deallocate, intrinsic assignment and move_alloc. Some tracked components
  (e.g. ham_i and reduced_acr) are absent from the manual destructor. Audit
  tracker balances separately from Fortran's automatic component cleanup.
- **bravais retains oversized arrays.** Trace surface/impurity consumers before
  shrinking cr/iz/num or releasing crd/izp/no. PBC generation and odd kk handling
  also deserve explicit numerical cases.
- **structb coordinate expressions create candidates for large temporaries.**
  cr*alat operates on the whole allocated array, not just kk. Review packing and
  scaled-coordinate reuse with a compiler/runtime measurement.
- **dbar1 uses a fixed 5250 capacity for local sbar before clusba returns the
  actual cluster size.** Allocate from that returned size and verify downstream
  dense-matrix extents in a dedicated follow-up.
- **symbolic_atom's nested allocation/finalization and downstream last uses**
  remain for the next object-level audit. This pass inspects the atomlist entry
  point, not every potential/element constructor.
- **Input validation remains incomplete.** Positive lattice dimensions, valid
  atom/type indices, cutoff availability and allocation failures need consistent
  validation rather than relying on default zeros.

## Validation status

The patch was reviewed against the pinned source, both preprocessor allocation
branches, and explicit dimensions of nncal/remd/outmap/mapa. No Fortran compiler
or shell is exposed in this session, so no compilation, runtime bounds checks,
SCF regression, or memory profile was executed locally. Existing GitHub build
and test workflows are triggered by the branch push; their result must be checked
before merging.

Required follow-up validation: bulk bccFe, file-based lattice with ct, fcc2/hcp,
a no-neighbour map, map-capacity exhaustion under bounds checking, and both
USE_SAFE_ALLOC configurations. Compare numerical outputs with the base commit;
measure peak memory separately.
