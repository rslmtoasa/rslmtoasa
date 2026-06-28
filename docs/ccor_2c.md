# Two-Centre Combined Correction

This implementation adds an optional two-centre LMTO-ASA combined correction as
a separate additive Hamiltonian operator:

```text
H = h - h*o*h + E_nu + Hcc_2c + H_LS
```

`Hcc_2c` is stored in `hamiltonian%eecc` for bulk/host blocks and
`hamiltonian%hallcc` for local/impurity blocks. It is not accumulated into
`ee` or `hall`, so it does not enter the existing `h-o-h` sweep.

For workflows that apply a first-order/no-`hoh` Hamiltonian, the accelerated
drivers may use cached merged work arrays

```text
ee_ccor_work   = ee   + eecc
hall_ccor_work = hall + hallcc
```

These merged arrays are an implementation detail of the operator application.
They are not used as the bare `h` input to an `h-o-h` product.

## Formula

The implemented scalar kernel is

```text
Kcc_ij = Ddot_ij + ccd1_i*D_ij + D_ij*ccd1_j
Kcc_ii = Kcc_ii + ccd0_i*d_i^2
Hcc_ij = (VMT - E_lin) * Kcc_ij
```

with `d = sqrt(Delta)` from the screened potential parameters. The noncollinear
path uses the same `wx0/wx1` and local moment convention as `ham0m_nc`; no
CCOR-spin-orbit cross terms are generated.

## Code Mapping

```text
S              -> lattice%sbar
Sdot_raw       -> lattice%sdot
sdot_cc        -> -avw**2 * lattice%sdot
alpha          -> potential%screening_alpha, falling back to lattice%alpha
adot           -> lattice%alpha_dot
wsr            -> symbolic_atoms(:)%potential%ws_r
avw            -> lattice%wav converted to Bohr when possible
sqrt(Delta)    -> potential%wx0 / wx1 spin decomposition
E_lin          -> hamiltonian%ccor_elin
VMTZ           -> potential%vrmax through hamiltonian%ccor_vmt_mode strategy
```

`makdia` coefficients are stored internally after Questaal's scaling:

```text
ccd = avw**2 * makdia_raw
```

`ccd(:,2)` is computed only as a diagnostic. It is not included in `Hcc_2c`.

## VMT Convention

The default `ccor_vmt_mode='surface_scalar'` uses the Stuttgart `getmtz`
surface average. `potential%vrmax(1)` stores the spin average
`-2*Z/wsr + (V_up(r=wsr)+V_down(r=wsr))/2`; `potential%vrmax(2)` stores
`V_up(r=wsr)-V_down(r=wsr)`. The scalar used by the currently active CCOR
Hamiltonian is the spin average of

```text
VMTZ_sigma = sum_type(n_type * wsr_type**2 * V_surface_type_sigma)
           / sum_type(n_type * wsr_type**2)
```

The NCOL/Stuttgart pair endpoint rule is also isolated in the implementation:

```text
VMT_pair_ij_sigma =
  (wsr_i**2 * V_surface_i_sigma + wsr_j**2 * V_surface_j_sigma)
  / (wsr_i**2 + wsr_j**2)
```

`ccor_vmt_mode='pair_surface'` uses this endpoint value in the ordered
noncollinear expression

```text
Hcc_ij =
    0.5*(Lambda_L*Ddot_ij + Ddot_ij*Lambda_R)
  + Lambda_L*A1_i*D_ij
  + D_ij*A1_j*Lambda_R
```

and adds `Lambda_L*A0_i` on site. The legacy debug fallback
`ccor_vmt_mode='vmad_scalar'` uses

```text
VMT = sum_type(wsr_type**2 * vmad_type) / sum_type(wsr_type**2)
```

This is not a verified `VMTZ` and should not be used for production CCOR
validation.

Older `*_out.nml` potential files do not contain `vrmax`. In non-strict mode,
surface VMT modes fall back endpoint-by-endpoint to `potential%vmad` with a
warning so old post-processing workflows do not silently produce `Hcc=0`.
Set `ccor_strict=.true.` to require persisted `vrmax` and fail instead.

## Controls

Add to `&hamiltonian`:

```text
ccor_2c = .true.
ccor_elin = 0.0
ccor_vmt_mode = 'surface_scalar'
ccor_debug = .false.
ccor_strict = .false.
```

`ccor_2c` requires `strux_want_sdot=.true.` and screened `sigma` or `fitted`
structure constants for meaningful results.

For reciprocal workflows, `Hcc(k)` is Bloch-summed from `eecc` in both first-
and second-order paths. The deprecated `kspace_ham_order='proper'` input is
treated as `second`.

## Backend Status

The legacy CPU recursion paths apply `Hcc_2c` exactly as a separate additive
operator. In `hoh` mode this preserves

```text
H = h - h*o*h + E_nu + Hcc_2c + H_LS
```

and does not feed `Hcc_2c` into either `h` factor.

For no-`hoh` runs, the fast CPU Chebyshev, batched/MKL Chebyshev,
block-Lanczos, stochastic, orbital, and CUDA plugin paths consume the merged
`ee_ccor_work/hall_ccor_work` blocks. The CUDA BSR export also builds the BSR
matrix from these merged blocks when `ccor_2c=.true.`, so the sparse GPU path
does not drop CCOR in no-`hoh` mode.

Accelerated `hoh+ccor_2c` remains routed to the exact legacy CPU path. It
cannot be represented by simply merging `ee+eecc`, because that would generate
unwanted `Hcc*o*h`, `h*o*Hcc`, and `Hcc*o*Hcc` terms. A future accelerated
implementation should keep the existing two-sweep `h-o-h` kernel and add a
separate `Hcc` matvec after the orthogonalization correction.

## Omissions

The implementation deliberately omits three-centre CCOR, active `ccd2`
contributions, mixed CCOR-spin-orbit terms, and a dedicated accelerated
`hoh+ccor_2c` kernel.
