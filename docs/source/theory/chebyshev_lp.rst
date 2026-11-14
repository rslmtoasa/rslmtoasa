.. _theory/chebyshev_lp:

===============================================
Linear Prediction of Chebyshev Moments
===============================================

Overview
========

The real-space Chebyshev recursion (``recur = 'chebyshev'``) computes
``N = 2*lld + 2`` moments of the local Green function of each recursion atom
(and, for the exchange route, of each atom pair). The Jackson-damped
reconstruction of the DOS has a resolution of about :math:`\pi a / N` in energy
(``a`` is the scaling half-width of the Chebyshev window).

*Linear prediction* (LP) extends the computed moments to ``cheb_lp_factor * N``
moments before the reconstruction, at a small fraction of the cost of computing
them. The idea is that of Ganahl et al., PRB 90, 045144 (2014), applied to
the moments of the real-space cluster Green functions. LP approximates a run
with more computed moments **on the same cluster**; it does not add information
about the infinite solid (see :ref:`theory/chebyshev_lp_limits`).

The feature is off by default. With ``cheb_lp = .false.`` all results are
unchanged (checked bit for bit on one thread against the code before LP).

Method
======

For one moment block (one atom, or one pair slot) the moments of order
``n = 0 ... N-1`` are assumed to obey a linear recurrence

.. math::

   x(n) = -\sum_{i=1}^{p} a_i\, x(n-i).

* **One real coefficient vector per block.** The real and imaginary parts of all
  ``nb*nb`` matrix elements are channels of one fit, so the extension commutes
  with conjugation and transposition and keeps the block Hermitian. Hermitian
  blocks are fitted on the ``nb*nb`` unique real channels with weights 1 and 2;
  a block that is not Hermitian (relative deviation above 1e-8) is fitted on all
  ``2*nb*nb`` channels.
* **Order and window.** The training window is all ``N`` known moments and the
  order is ``p = N/2`` (``L = T/2`` in the notation of the paper).
* **Fit.** The normal equations ``R a = -r`` are built from the lag Gram matrix
  of the channel data (one ``dsyrk``) and solved with a pseudo-inverse from the
  eigendecomposition of ``R``, keeping eigenvalues above ``1e-6`` of the largest.
* **Extrapolation.** The companion matrix is diagonalised. Roots with
  :math:`|\lambda| \ge 1` are set to zero (the paper allows this or
  renormalising them to unit modulus, "as long as the corresponding coefficient
  is small"). The extension is the real part of the propagated state, evaluated
  as one real matrix product.
* **Reconstruction is unchanged.** All reconstruction routes use the length of
  the stored array, so the damping kernel (Jackson by default) is applied at the
  extended length.
* **Positivity.** The existing local spectral positivity projection still runs.
  The extended moments fix the zeroth moment, so clipping negative DOS adds
  charge and moves the Fermi level. With LP each atom, orbital and spin block
  is therefore rescaled after the projection so that its integrated weight equals
  the weight before clipping (the spin coherence by the geometric mean of the
  two factors, which keeps the 2x2 block positive semidefinite).

The constants (order, window, cutoff, root treatment) are module parameters in
``source/cheb_lp.f90``, not input keys.

Where it runs
-------------

* On-site moments: ``chebyshev_recur`` (CPU and CUDA paths), for every local atom.
* Pair moments: ``chebyshev_recur_ij`` (four slots per pair; self pairs fill the
  first slot only), used by the exchange route.
* The extension of all local atoms or slots runs in one OpenMP loop after the
  moment loops, so under MPI each rank extends its own blocks. On the CUDA path
  the extension runs on the host after the moments are back from the device.

Usage
=====

.. code-block:: fortran

   &energy
   energy_min = -1.2          ! the window must cover the whole spectrum of H
   energy_max =  1.2
   channels_ldos = 2500       ! the mesh must resolve pi*a/N_ext
   /
   &control
   recur = 'chebyshev'
   lld = 150                  ! 302 computed moments
   cheb_lp = .true.           ! default .false.
   cheb_lp_factor = 2         ! default 2; 604 moments
   /

``cheb_kernel`` (``'jackson'`` default, ``'lorentz'``, ``'none'``) and
``cheb_kernel_lambda`` select the damping kernel; see :ref:`keywords/control_parameters`.
Example decks: ``example/bulk/bccFe_chebyshev_lp`` (SCF) and
``example/exchange/bccFe_chebyshev_lp`` (J_ij).

The log shows the memory of the moments and one line per atom or pair slot::

   Chebyshev moments: 604 orders x 1 slots; mu_n = 3.1 MB
   Chebyshev LP atom 1: p=151 rank=76 unstable=4 (weight 0.032) lam_max=1.007467 fit_res=4.24E-02 imag=...
   Chebyshev LP charge-conserving rescale: max|1-s|=...

``weight`` is the share of the mode amplitude carried by the dropped roots. A
large value means the prediction is mostly discarded and should not be trusted.

Cost, memory, threading
-----------------------

* One extension call (``nb = 18``, 302 known moments, factor 2) takes about
  25 ms with one BLAS thread (about 105 ms before the Hermitian and lag-Gram
  formulations). The moment recursion takes about 22 s per atom on the
  bulk Fe example, so LP is about 0.1% of an on-site run.
* ``mu_n`` takes 16 bytes per element, ``nb*nb*N_ext`` per slot (3.1 MB at
  ``nb = 18``, 604 orders). The pair route has ``4*njij`` slots, allocated on
  every MPI rank.
* The extension is faster when the BLAS library uses one thread: on the 12-slot
  exchange example the LP step took 0.25 s with 1 OpenMP and 1 BLAS thread,
  0.08 s with 4 OpenMP and 1 BLAS thread, and 0.65 s with 4 and 4. Set
  ``OPENBLAS_NUM_THREADS=1`` (or ``MKL_NUM_THREADS=1``) when running with ``cheb_lp``.

What was measured
=================

All numbers are states/Ry (DOS) and are differences against a **computed** run
with the same number of moments on the same cluster, unless stated otherwise.
Machine: 4 shared cores, gfortran 13, OpenBLAS.

.. list-table::
   :header-rows: 1
   :widths: 24 46 30

   * - Case
     - Measured
     - Reference
   * - Unit oracles ``UnitChebyshevLP`` U1-U6
     - Closed-form damped poles, Hermitian and non-Hermitian blocks, growing
       pole (dropped-root weight 0.9209 against a closed form of 0.9209),
       Poisson-kernel density, inexact data against an independent
       least-squares reference: errors 1e-14 to 1e-13 against tolerances
       1e-11 / 1e-8
     - closed forms written in the tests
   * - bcc Fe, frozen potential, 5984-site cluster, lld = 150, factor 2
     - RMS DOS difference, occupied band and within 0.1 Ry of E_F, up/down
       spin: 0.13/0.13/0.14/0.14 (LP) against 0.89/1.45/0.47/0.75 (302
       moments); E_F -0.069642 (LP), -0.068259 (302), -0.069078 (reference)
       Ry; integrated charge 18.0027 against 18.0026
     - computed lld = 301 (604 moments)
   * - Ni(001) surface, LDA+U, E_F fixed, lld = 200, factor 2
     - Occupied-band RMS: Ni-S up/down 0.104/0.120 (LP) against 0.746/0.866;
       Ni-S-1 0.120/0.143 against 0.614/0.475; empty sphere 0.002/0.002
       against 0.014/0.025; DOS(E_F) of Ni-S 17.017 against 16.983 (reference)
       and 17.354 (plain)
     - computed lld = 401
   * - Si example deck, lld = 200, factor 2
     - RMS 0.31 against 0.85 (plain lld = 200); E_F inside the gap
     - computed lld = 400
   * - J_ij, triad Fe deck (1000-site periodic cluster), factor 4
     - shell 1/2/3 error +0.0206/+0.0074/-0.0074 against +0.229/+0.013/-0.010
       (302 moments)
     - computed lld = 603
   * - bcc Fe SCF, factor 2
     - Diff 3.1e-6 after 24 steps; magnetic moment 2.1105 muB (LP) against
       2.1101 (plain lld = 301); E_F -0.07212 against -0.06970 Ry
     - plain Chebyshev lld = 301 (8 steps)
   * - MPI, 2 ranks
     - exchange example (12 slots) J_ij identical to serial; Si deck E_F and
       LP fit identical to serial
     - serial run
   * - fp32 (``cheb_backend = 'fast'``), factor 4, earlier projection
     - maximum DOS difference to fp64 0.05 (302 moments) and 0.49 (LP) against
       a peak of about 40
     - fp64 run

.. _theory/chebyshev_lp_limits:

Limits
======

* **Finite cluster.** The moments of a finite cluster agree with those of a
  larger cluster only up to the order at which the signal reaches the cluster
  boundary: for bcc Fe, comparison to 1e-3 gave n = 88 (5984 sites), 110 (9720),
  136 (16888); for periodic supercells 24, 33, 50, 79, 127 for 8^3, 10^3, 14^3,
  20^3, 30^3 cells. About 70% of the 302 moments of the 5984-site cluster are
  boundary-contaminated, and LP extends them as they are. LP trained only on
  boundary-free moments did not predict a larger system's moments better than
  the smaller system's own moments (0 of 288 configurations). LP is therefore a
  cheaper way to obtain what more computed moments on the same cluster would
  give, not a route to the infinite solid.
* **Not equal to the computed reference.** On bulk Fe with factor 2 a
  residual RMS of about 0.13 states/Ry against the computed 604-moment run
  remains. Factor 4 is worse than factor 2 on bulk Fe (0.69/0.83/0.54/0.68
  against its computed 1208-moment run) and on Si (RMS 1.28 against a computed
  lld = 803 run, against 0.31 at factor 2).
* **Fit reliability.** The dropped-root weight was 0.03 on bulk Fe, 0.08 to 0.12
  on the triad J_ij deck and 0.865 on the Si deck, where most of the predicted
  signal is discarded although the DOS curve improves. Check the log line.
* **Moments do not decay.** The diagonal moments of the cluster Green functions
  stay near 0.1 up to n = 1200, so the paper's setting (exponentially decaying
  moments, no damping kernel) does not apply: on bcc Fe a reconstruction without
  a kernel had 473 negative-DOS points and RMS 1.0/2.1/0.4/0.6 against
  0.14/0.13/0.14/0.14 with Jackson; on Si the RMS was 1.00 without a kernel,
  0.59 with Lorentz and 0.31 with Jackson. A small imaginary shift (1e-5) did not
  change this.
* **Positivity projection.** It acts on thousands of local blocks per DOS
  evaluation with LP (about 5000 to 6500 on bcc Fe, 1800 on Si), against none
  for plain Jackson. The rescale removes the charge this adds.
* **Fermi level.** A gapped system is sensitive to small charge errors because
  the integrated DOS is flat in the gap: on Si, E_F moved from 0.1727 Ry (plain,
  inside the gap) to 0.0656 Ry (below the gap) with LP before the rescale and to
  0.205 Ry (inside the gap) after it. In the SCF on bulk Fe the converged E_F
  with LP is 2.4 mRy lower than plain Chebyshev at the same moment count; the
  moments and the DOS shape relative to E_F agree (RMS 0.32 states/Ry), and the
  offset was not diagnosed.
* **Cost grows with moments.** The reconstruction is linear in the stored
  length, so it costs 2 to 4 times more at factor 2 to 4; memory of the moments
  grows by the same factor.

Not tested
==========

* **GPU.** The LP hooks on the CUDA paths and the kernel selection
  (``cheb_kernel``) through the CUDA reconstruction were never compiled with
  ``nvcc`` or run on a GPU. Only a single GPU has been available. The CUDA
  reconstruction call in ``chebyshev_green_core`` indexes with the global atom
  range (``tests/KNOWN_ISSUES.md``); it was never run with more than one MPI
  rank on a GPU.
* **MPI** with more than 2 ranks.
* **Materials and settings.** Bulk Fe, Ni(001) and Si only. Not run: bases
  other than those of these decks (``nb = 18`` on Fe and Ni, ``nb = 8`` on Si;
  the unit tests use ``nb = 1`` to 3), ``nsp`` other than 2 (the Si deck uses
  ``nsp = 1``), long SCF runs on surfaces and with LDA+U (the Ni deck was run
  for one step on a frozen potential), and the conductivity and orbital-moment
  routes, which do not use LP.
* **Comparison with the Haydock recursion** (``recur = 'block'`` with a
  terminator) and with the k-space Lehmann route at matched broadening.
* **Memory under MPI.** ``mu_n`` is allocated for all slots on every rank.

References
==========

* M. Ganahl et al., Phys. Rev. B 90, 045144 (2014), arXiv:1403.1209.
* T. Barthel, U. Schollwoeck, S. R. White, Phys. Rev. B 79, 245101 (2009).
* A. Weisse et al., Rev. Mod. Phys. 78, 275 (2006), for the damping kernels.

Related code and tests
======================

``source/cheb_lp.f90`` (kernel and ``type cheb_lp``), ``source/recursion_chebyshev.f90``,
``source/recursion_transport.f90``, ``source/green_chebyshev.f90`` (reconstruction,
projection and rescale, kernel selection), ``source/chebyshev_fast.f90``;
tests ``UnitChebyshevLP``, ``UnitChebyshevKernel``, ``UnitChebyshevGreenPrecision``,
``UnitChebyshevDamping``. The campaign scripts (toy study and production check on
frozen potentials) are archived at the tag ``cheb-lp-campaign-2026-11``
(``tests/validation/val_cheb_lp_toy.py``, ``tests/validation/val_cheb_lp.py``) and are not in the tree.
