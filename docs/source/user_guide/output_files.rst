.. _user_guide/output_files:

=======================================
Output Files
=======================================

Overview
========

After running RS-LMTO-ASA, various output files are generated depending on your calculation 
settings. This section describes the main output files and how to interpret them.

Standard Output Files
=====================

**<label>_out.nml**

One file per atom type (for example ``Fe_out.nml``). It is a Fortran namelist with the
``&element`` and ``&par`` keys of that atom type (``mom``, ``lmom``, ``ql``, ``etot``, ...) at the
end of the run. A later run reads it instead of ``<label>.nml`` (see ``fresh_start`` in ``&control``).

Format: Same as the per-label input file (Fortran namelist)

Example:

.. code-block:: fortran

   &element
      symbol='Fe',
      ...
   /

**report.out**

Text report written at the end of the SCF run. It contains the total energy
(``Total energy of system:``), the band energy, the spin and orbital moments per atom
(``Spin moment of atom``, ``Orbital moment of atom``), the occupations and charge transfer per
atom (``Occupation at atom``, ``Charge transfer at atom``), the Fermi energy
(``Fermi energy:``) and the final mixing difference (``Total RMS Diff:``).

**Standard output / log**

The SCF progress is written to standard output (redirect it with ``> output.log``). After each
SCF iteration the code prints ``Not converged! Diff=`` followed by the mixing difference; the
last iteration prints ``Converged!`` followed by the mixing difference.

Density of States Output
========================

Real-space (recursion) calculations write, at the end of the DOS step of the SCF run:

- ``totaldos.out`` - two columns: energy relative to the Fermi level (Ry), total DOS
- ``magneticdos.out`` - columns: energy, ``dos_up``, ``dos_dw``, ``dos_mx``, ``dos_my``, ``dos_mz``, ``dos_nmag``
- ``<symbol>_dos.out`` - columns: energy, site DOS total, up, down (one file per atom type)
- ``<symbol>_orbital_dos.out`` - energy and orbital-resolved site DOS (one file per atom type)

Columns of ``totaldos.out``:

.. code-block:: text

   E-E_F (Ry)   Total DOS

With ``post_processing = 'density_of_states'`` the k-space DOS is written to ``dos_kspace.dat``
(header lines starting with ``#``, then ``Energy-E_F (Ry)`` and ``Total DOS``). A k-space SCF
(``use_kspace = .true.``) writes ``density_of_states.dat``. If projected DOS is available,
``projected_dos.dat`` is written as well.

**Usage:**

- Plot with gnuplot, matplotlib, etc.
- Integrate to verify electron count
- Identify band structure features

**Example plot (matplotlib):**

.. code-block:: python

   import numpy as np
   import matplotlib.pyplot as plt
   
   data = np.loadtxt('totaldos.out')
   E = data[:, 0]
   dos = data[:, 1]
   
   plt.plot(E, dos)
   plt.xlabel('E - E_F (Ry)')
   plt.ylabel('DOS')
   plt.show()

Band Structure Output
=====================

If ``post_processing = 'band_structure'``:

**band_structure.dat**

Energy bands along a k-path. Header lines start with ``#``; each data row has the k-path distance
followed by all eigenvalues at that k-point:

.. code-block:: text

   # Format: k_distance, eigenvalue_1, eigenvalue_2, ...

**band_structure_kpath.dat**

k-path information: k-path distance, ``kx, ky, kz`` and, at high-symmetry points, the label.

With ``post_processing = 'bsf'`` the Bloch spectral function files are ``bsf.dat`` and
``bsf_bands.dat``.

**Plotting:**

Use standard band structure visualization tools.

Magnetic Properties
===================

Spin and orbital moments per atom are printed in ``report.out`` (they are not written to a
separate file):

.. code-block:: text

   Spin moment of atom   1:  2.234000
   Spin moment projections of atom   1:  0.000000  0.000000  2.234000
   Orbital moment of atom   1:  0.045000

**Interpretation:**

- ``Spin moment of atom`` is the norm of the moment; the projections are its :math:`m_x, m_y, m_z` components
- Negative :math:`m_z` indicates opposite spin direction

**For non-collinear magnetism (nsp=3,4):**

All three components (:math:`m_x, m_y, m_z`) have meaning; they define the spin direction.

Exchange Interactions
=====================

If ``post_processing = 'exchange'``:

**jij.out**

Magnetic exchange interactions, one row per pair. Columns: ``iz(i)``, ``iz(j)``, the three
components of :math:`\mathbf{r}_j - \mathbf{r}_i`, ``jij``, and the distance :math:`|\mathbf{r}_i - \mathbf{r}_j|`:

.. code-block:: text

   iz_i  iz_j   dx  dy  dz   J_ij   distance

The file carries no header and does not state the unit of ``jij``. ``dij.out`` holds the
Dzyaloshinskii-Moriya vector (same first five columns, three components, distance), ``aij.out`` the
3x3 tensor ``aij`` (nine values) and ``jtens.out`` the exchange tensor. Further files named
``jijso.out``, ``jijfo.out``, ``dijso.out``, ``dijfo.out``, ``aijso.out``, ``aijfo.out``,
``jtensso.out``, ``jtensfo.out`` and ``*parts.out`` are written by other routines in
``source/exchange.f90``.

**Interpretation:**

- :math:`J > 0` -> Ferromagnetic coupling
- :math:`J < 0` -> Antiferromagnetic coupling
- Magnitude indicates strength

**Use for:**

- Mapping of exchange networks
- Estimation of ordering temperature (Curie/Néel)
- Comparison with experiment

**chi0_ij.dat** (with ``do_chi0 = .true.``)

Static transverse bare susceptibility of the d channel, one row per pair, after a one-line header:
atoms ``i``, ``j``, their types, :math:`\mathbf{r}_j - \mathbf{r}_i` (three components, cluster
coordinates), its length, and the real and imaginary parts of :math:`\chi^0_{ij}` (1/Ry).

**chi0_types.dat** (with ``do_chi0 = .true.``)

One row per atom type with a self pair, after a one-line header: type, seed atom, number of its pairs,
d moment (:math:`\mu_B`), Goldstone :math:`U` (Ry, by construction), Mills :math:`U` (Ry), Mills
residual, and :math:`\chi^0(q = 0)` towards every type (1/Ry), summed over the listed pairs. See
``do_chi0`` in :ref:`keywords/output_options`.

Transport Properties
====================

If ``post_processing = 'conductivity'``:

**cond_total.out**

Raw (bare-kernel) response of one tensor component, selected by ``v_alpha``, ``v_beta`` and
``pol_alpha``, ``pol_beta``. No header; three columns (``'(3es16.6)'``):

.. code-block:: text

   E-E_F (Ry)   Re   Im

The energy is in Ry and the values are not converted to physical units.
``cond_total_orb_real.out`` and ``cond_total_orb_im.out`` hold the orbital-resolved real and
imaginary parts, and ``<symbol>_cond.out``, ``<symbol>_cond_orb_real.out`` and
``<symbol>_cond_orb_im.out`` the per-atom-type contributions (``cond_calctype = 'per_type'``).

**response_total.out**

The same response converted to physical units. Header
``# E-Ef [eV], Re, Im; units: ...`` followed by three columns, energy in eV.
``<symbol>_response.out`` holds the per-type contributions, ``response_orbitals_real.out`` and
``response_orbitals_imag.out`` the orbital contributions, and ``response_units.txt`` the units,
directions and conversion factors used.

**Interpretation:**

- Which observable is computed (charge, spin or orbital current) is set by ``linear_out`` and ``linear_in`` in ``&control``
- Off-diagonal components (Hall response) are selected with ``v_alpha`` and ``v_beta``
- Sum over energies (with Fermi distribution) gives transport properties at temperature

Charge Transfer and Occupations
===============================

``report.out`` lists, for each atom, ``Occupation at atom``, the up and down orbital occupations
(``Up orbital occupation at atom``, ``Down orbital occupation at atom``) and
``Charge transfer at atom``.

Convergence Diagnostics
=======================

**Typical signs of good convergence:**

- The mixing difference (``Not converged! Diff=`` in the log, ``Total RMS Diff:`` in ``report.out``) decreases
- Fermi energy stabilizes
- Energy change becomes very small (< :math:`10^{-6}` Ry)
- Usually achieved in 10-100 iterations

**Warning signs (non-convergence):**

- Charge density oscillates
- Energy increases in some iterations
- The mixing difference plateaus at a high value
- Forces/torques don't decrease

**Action if non-convergent:**

- Reduce mixing parameter ``beta`` (``&mix``)
- Increase recursion cutoff (``llsp``, ``lld``)
- Use coarser energy mesh initially, then refine
- Check Hamiltonian parameters are reasonable

Reading and Plotting Output
============================

**Python script example:**

.. code-block:: python

   import numpy as np
   import matplotlib.pyplot as plt
   
   # Read DOS
   dos_data = np.loadtxt('totaldos.out')
   energy = dos_data[:, 0]
   dos = dos_data[:, 1]
   
   # Plot
   plt.plot(energy, dos)
   plt.xlabel('E - E_F (Ry)')
   plt.ylabel('DOS')
   plt.title('Density of States')
   plt.show()

Analyzing Magnetic Structure
============================

**Extract the spin moments from report.out:**

.. code-block:: python

   import re
   
   m_mag = []
   with open('report.out') as f:
       for line in f:
           m = re.match(r'Spin moment of atom\s*\d+:\s*(\S+)', line)
           if m:
               m_mag.append(float(m.group(1)))
   
   print("Spin moments per atom:", m_mag)

File Size and Storage
=====================

The size of the DOS, band-structure and response files grows with ``channels_ldos``, the number of
k-points and the number of energy points. For multiple energy windows or temperatures, output
files can be significantly larger; plan disk space accordingly.

Provenance
==========

Output is generated by:

- **Main output:** ``source/calculation.f90`` (overall workflow)
- **SCF log and report.out:** ``source/self.f90``
- **Real-space DOS:** ``source/bands.f90``
- **k-space DOS and band structure:** ``source/reciprocal_dos.f90``, ``source/reciprocal_bands.f90``
- **Exchange:** ``source/exchange.f90``
- **Conductivity:** ``source/conductivity.f90``

See Also
========

- :ref:`user_guide/input_files` - Input file format
- :doc:`examples` - Worked examples with sample output
- :ref:`keywords/output_options` - Parameters controlling output
