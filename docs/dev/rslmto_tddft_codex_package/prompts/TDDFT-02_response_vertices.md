# Codex Prompt — TDDFT-02: Generalized Response Basis and Transition Vertices

## Mission

Implement response channels and spinor transition matrix elements reusable by transverse, longitudinal, SOC, and non-collinear LR-TDDFT.

## Formal object

For channel `A=(a,mu)`,

\[
M^A_{nm}(\mathbf k,\mathbf q)
=
\langle n\mathbf k|P_a\sigma^\mu|m,\mathbf k+\mathbf q\rangle.
\]

The first production projection is site-summed/orbital-summed, but do not make the data structure incapable of later l- or orbital-resolved response.

## Tasks

- [x] Define a response-channel descriptor containing at least:
  - [x] site;
  - [x] density/spin component;
  - [x] optional l sector;
  - [x] optional orbital.
- [x] Implement site-projected operators using the basis conventions from TDDFT-00.
- [x] Implement Cartesian charge/spin vertices.
- [x] Implement circular `+/-` vertices as consistent transforms.
- [x] Compute matrix elements from full spinors; do not require pure up/down eigenstates.
- [x] Verify whether response projectors should use global or local spin frames at this stage.
- [x] For collinear no-SOC systems, demonstrate:
  - [x] charge and z vertices are spin conserving;
  - [x] `+/-` vertices are spin flipping;
  - [x] transverse and charge-longitudinal sectors decouple.
- [x] Add global spin-rotation covariance tests.
- [x] Provide a simple reference implementation and a batched production-friendly interface.
- [x] Avoid storing a full `(k,n,m,A,B)` tensor.

## Recommended computational shape

Build a response-space transition vector `v_A` and let later code accumulate

\[
w(\omega)\,v_Av_B^*.
\]

## Acceptance

For tiny synthetic spinors, explicit dense matrix multiplication must match the response-vertex helper to machine precision.

## Completion protocol

Tick all boxes and commit once:

`TDDFT-02: implement generalized response vertices`
