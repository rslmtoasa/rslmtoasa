.. _keywords/bdg:

=============================================
Bogoliubov-de Gennes pairing (&bdg)
=============================================

Overview
========

The ``&bdg`` namelist is read by ``post_processing = 'bdg'``. The route works on the frozen
k-space state: the k mesh (full Brillouin zone) and the temperature come from ``&reciprocal``
(:ref:`keywords/reciprocal`), the chemical potential :math:`\mu` from ``&energy``. The
on-site singlet gap :math:`\Delta_{iL}` is iterated with
:math:`\Delta_{iL} = \sum_{L'} K_{LL'}\,\chi_{iL'}`. The namelist is optional; every key has a
default, and a read error stops the run.

.. code-block:: fortran

   &calculation
      pre_processing = 'none'
      post_processing = 'bdg'
   /

   &bdg
      lambda(1, 1, 1) = 0.2
      lambda(3, 3, 1) = 0.5
      delta_seed = 1.0e-3
      mix = 0.3
      tol = 1.0e-6
      max_iter = 50
   /

Keys
====

.. list-table::
   :header-rows: 1
   :widths: 30 10 18 42

   * - Key
     - Type
     - Default
     - Meaning
   * - ``lambda(l+1, l'+1, t)``
     - real
     - 0.0 (eV)
     - Pairing coupling between shells l and l' of atom type t, array of size
       (lmax+1, lmax+1, ntype).
   * - ``delta_seed``
     - real
     - 1.0e-3 (eV)
     - Starting gap, the same on every orbital whose row of K is nonzero.
   * - ``delta_file``
     - character
     - ``''``
     - Restart from a ``bdg_delta.dat`` written by an earlier run; takes precedence over
       ``delta_seed``.
   * - ``mix``
     - real
     - 0.3
     - Linear mixing of the new gap into the old one; must be above 0.
   * - ``tol``
     - real
     - 1.0e-6 (eV)
     - Convergence threshold on :math:`\max|\Delta_\mathrm{new} - \Delta_\mathrm{old}|`; must be above 0.
   * - ``max_iter``
     - integer
     - 50
     - Maximum number of gap iterations; must not be negative. 0 does one pass at the given
       gap, for a DOS on a denser mesh.
   * - ``dos_window``
     - real
     - 0.01 (eV)
     - Half width of the energy window of the DOS around :math:`\mu`; must be above 0.
   * - ``n_dos``
     - integer
     - 401
     - Number of DOS energy points; at least 1.
   * - ``dos_sigma``
     - real
     - 1.0e-4 (eV)
     - Gaussian broadening of the DOS; must be above 0.

Pairing kernel
==============

The kernel is built from ``lambda``: :math:`K_{LL'} = \lambda(l(L)+1,\, l(L')+1,\, t)`, where
:math:`l(L)` is the shell of orbital L (s: 1, p: 2 to 4, d: 5 to 9, f: 10 to 16). It is the
same for all m of a shell. Values are read in eV and used in Ry.
