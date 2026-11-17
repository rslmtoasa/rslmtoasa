.. _keywords/spin_response:

=============================================
Transverse Spin Response (&spin_response)
=============================================

Overview
========

``post_processing = 'spin_response'`` computes the transverse spin
susceptibility :math:`\chi(\mathbf{q}, \omega)` of a converged collinear ground
state on the k-space route: the bare response :math:`\chi_0` from the
eigenpairs at :math:`\mathbf{k}` and :math:`\mathbf{k}+\mathbf{q}`, then
:math:`(1 + \chi_0 U)\chi = \chi_0`. The k mesh and Hamiltonian come from
``&reciprocal`` (:ref:`keywords/reciprocal`). The namelist is optional; every key
has a default, and a read error stops the run.

.. code-block:: fortran

   &calculation
      pre_processing = 'none'
      post_processing = 'spin_response'
   /

   &spin_response
      n_q = 2
      q_list(:, 1) = 0.0, 0.0, 0.0
      q_list(:, 2) = -0.0833333, 0.0833333, 0.0833333
      omega_min = -0.004
      omega_max = 0.02
      n_omega = 193
      eta = 5.0e-4
   /

Keys
====

.. list-table::
   :header-rows: 1
   :widths: 34 10 18 38

   * - Key
     - Type
     - Default
     - Meaning
   * - ``method``
     - character
     - ``'juelich'``
     - Transverse amplitudes: ``'juelich'`` or ``'mills'``.
   * - ``kernel``
     - character
     - ``'goldstone'``
     - Exchange kernel :math:`U`: ``'goldstone'`` (scalar U on the d channel, fixed by the
       Goldstone condition :math:`\chi_0(0,0)\,U M = -M`), ``'alsda'`` or ``'alsda_scaled'``
       (ALSDA radial kernel, unscaled or times the Goldstone factor). The ALSDA kernels are
       experimental and need ``txc = 1`` (LDA).
   * - ``channels``
     - character
     - ``'d'``
     - Orbital channels of the response space: letters of ``spdf``, must contain ``d``,
       limited by ``lmax``.
   * - ``n_q``
     - integer
     - 1
     - Number of q vectors read from ``q_list`` (1 to 64).
   * - ``q_list(:, i)``
     - real
     - 0.0
     - q vectors in fractional coordinates of the reciprocal vectors :math:`\mathbf{b}_1,
       \mathbf{b}_2, \mathbf{b}_3`.
   * - ``q_file``
     - character
     - ``''``
     - File with one q vector (three fractional coordinates) per line; when set it replaces
       ``n_q`` and ``q_list``.
   * - ``omega_min``, ``omega_max``
     - real (Ry)
     - 0.0, 0.05
     - Frequency window; ``omega_max > omega_min`` when ``n_omega > 1``.
   * - ``n_omega``
     - integer
     - 201
     - Number of frequencies, equally spaced; ``n_omega = 1`` gives the static response at
       ``omega_min``.
   * - ``eta``
     - real (Ry)
     - 0.002
     - Broadening :math:`\omega + i\eta`; must be :math:`\geq 0`.
   * - ``output_prefix``
     - character
     - ``'spin_response'``
     - Prefix of the output files.

With spin-orbit coupling (``nsp = 2`` or ``4`` with SOC) only ``method = 'juelich'``
with ``kernel = 'goldstone'`` is accepted.

Output files
============

- ``<prefix>_q<NNN>.dat`` — the response for each q.
- ``<prefix>_dispersion.dat`` — per q: q (fractional and Cartesian), peak and
  crossing frequencies (Ry).
- ``<prefix>_adiabatic.dat`` — per q: adiabatic frequencies from the static
  :math:`\chi_0`, and the force-theorem frequencies from the same :math:`\chi_0`.
- ``<prefix>_state.dat`` — the ground state the response used: method, Fermi level,
  electron count, k mesh, site moments and d occupations, and the kernel U per site.

Each file starts with a one-line header naming its columns and units.

**Related code:** ``source/spin_response.f90``, ``source/calculation_reciprocal.f90``.
