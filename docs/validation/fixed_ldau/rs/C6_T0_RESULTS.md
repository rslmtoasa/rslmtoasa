# C6 T0 — surface fixture and host setup

No SCF runs were performed.

The only fixture under `tests/scf/cases/surface/` is `fccCu001`. It uses
`pre_processing='buildsurf'`, `calctype='S'`, and `nlay=2`. The setup builds a
one-sided surface: `lattice_cluster` assigns `nbulk=nbulk_bulk` frozen host
types, `nrec=nlay` self-consistent recursion types, and maps sites through
`ib` (frozen) and `irec` (recursion). `is_frozen_hubbard_host` already accepts
`S` and classifies symbolic types from those maps. Here the frozen Cu host is
`Cu.nml`; the two active recursion types are `ES.nml` and `Cu-S.nml`.

The fixture explicitly sets `fermi=-0.088376` and `fix_fermi=.true.`. The
`energy` defaults for `S` also set `fix_fermi=.true.`, as they do for `I`, so
the host endpoint Fermi energy remains fixed when supplied in the deck.

Selected fixture: `fccCu001`, the simplest existing fcc surface. It is not a
valid T2/T3 candidate under the same-material host rule yet: its Cu bulk test
potential has U=J=0, and the accepted local U-carrying bulk endpoints found in
this checkout are Fe and Ni. No Fe or Ni surface fixture exists here. Under
C5.3 item 2, this material/host mismatch is recorded as uncovered; no deck is
launched or synthesized by changing the fixture's material or geometry.
