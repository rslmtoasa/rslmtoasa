# TDDFT-02: Generalized response vertices

`response_vertices_mod` defines a `response_channel` as a site projector and
a charge/spin component, with optional `l_sector` and site-local `orbital`
selectors.  The first production choice is made by leaving both selectors at
`RESPONSE_UNRESOLVED`, which sums all local orbitals on that site.

The reciprocal spinor layout is site-major, and each site block is spin-major:

```
(site 1: orbital up..., orbital down..., site 2: orbital up..., orbital down..., ...).
```

`site_projected_operator` provides an explicit dense reference operator.
`response_transition_vertex` instead directly contracts arbitrary complex
spinors, and `response_transition_vectors` accepts a batch of transitions and
returns only `vertices(nchannels,npairs)`.  Later Lehmann accumulation can
immediately form `w*vertices(:,p)*conjg(vertices(:,p))`; this layer does not
store a `(k,n,m,A,B)` tensor.

Circular channels use the frozen TDDFT-00 measurement convention:
`PLUS = MX + i*MY` and `MINUS = MX - i*MY`.  They are therefore twice the
half-normalized ladder operators.

## Spin frame

Response projectors use the global spin frame at this stage.  This is the
frame in which the spinor Hamiltonian, `sigma_x/y/z`, and LR-TDDFT response
components are defined.  Ground-state local axes may still be useful for
reporting projected DOS, but inserting them into a response projector would
silently give each site a different response basis and obscure covariance.
Local-frame response can be added later as an explicit channel-basis
transformation, with the corresponding kernel transformation, rather than as
an implicit projection choice.

`UnitResponseVertices` verifies dense-operator agreement for complex,
non-collinear synthetic spinors, collinear selection/decoupling rules, and
global SU(2) rotation covariance.
