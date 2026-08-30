.. _theory/kspace_modes:

================================
K-Space Hamiltonian Options
================================

Overview
========

The reciprocal-space workflow is controlled from ``&reciprocal`` with
``kspace_ham_order``. All active paths build a standard Hermitian
Hamiltonian and solve with identity overlap.

Available options:

- ``kspace_ham_order='second'``: second-order ASA Hamiltonian.
- ``kspace_ham_order='first'``: first-order ``h(k) + H_LS`` Hamiltonian.
- ``kspace_ham_order='auto'``: use the default production path.
- ``kspace_ham_order='proper'``: deprecated alias for ``second``.

Current Functionality Inventory
===============================

Implemented now
---------------

- Reciprocal lattice generation and Monkhorst-Pack meshes.
- Optional symmetry-reduced k-point mesh generation when symmetry support is available.
- High-symmetry k-path generation and band-structure output.
- Real-space to reciprocal-space Fourier transform for multi-site ``H(k)`` assembly.
- Dense k-point diagonalization with LAPACK ``ZHEEV``.
- DOS workflows:
  
  - tetrahedron
  - Blochl-improved tetrahedron
  - Gaussian broadening

- Orbital/site projected DOS from eigenvectors.
- Fermi-level search from DOS and band-moment integration.
- Optional diagnostics:
  
  - k-space mode summary
  - ``H(Gamma)`` bounds diagnostics
  - experimental finite ``HALL`` diagonalization check.

Dense Fermi-surface post-processing
-----------------------------------

Fermi-surface data should normally be generated as a separate post-processing
calculation, since it needs a considerably denser k-mesh than SCF. Select the
route in ``&calculation`` and configure its mesh in ``&reciprocal``::

   &calculation
      post_processing = 'fermi_surface'
   /
   &reciprocal
      fs_nk1 = 64
      fs_nk2 = 64
      fs_nk3 = 64
      fs_use_symmetry_reduction = .false.
      write_eigenpair_projections = .true.  ! optional site/orbital/local-spin weights
      eigenpair_output_file = 'fermi_surface'
   /

Non-positive ``fs_nk1/2/3`` values use the regular ``nk1/2/3`` mesh. The route
rebuilds the converged bulk post-processing Hamiltonian and always writes the
complete Brillouin-zone mesh: ``fs_use_symmetry_reduction`` is accepted for
input compatibility but ignored because Fermi-surface colouring needs the
eigenvectors at every k-point.

This writes ``fermi_surface.bin`` and ``fermi_surface.meta``. The binary stream
contains fractional reciprocal coordinates, normalized k-point weights,
eigenvalues in Ry, and complex eigenvectors in Fortran column-major order. If
projections are enabled it also contains
``projection_weights(site, orbital, spin, band, kpoint)`` with orbital channels
``s,p,d,f`` and spin channels along each site's converged local moment axis.
The sidecar contains the basis-index map, so raw coefficient projections can
be constructed independently when needed. It also records ``site_N_species``
labels used by the viewer's species mode. For a Python reader, use
``tools/read_kspace_eigenpairs.py``::

   from tools.read_kspace_eigenpairs import read_kspace_eigenpairs
   data = read_kspace_eigenpairs('fermi_surface')
   energies = data['eigenvalues']       # (band, kpoint), Ry
   weights = data['projection_weights']  # optional (site, l, spin, band, kpoint)

Interactive PyVista viewer
--------------------------

Install the optional viewer dependency with ``python -m pip install pyvista``.
The standalone viewer ``tools/plot_fermi_surface.py`` extracts the selected
band's isosurface and provides sliders for band, Fermi level, and opacity. To
show every sampled band that crosses the current Fermi level in one scene, use
``--all-bands``::

   python tools/plot_fermi_surface.py fermi_surface --all-bands --color-by band

The ``band`` colour mode gives each Fermi-surface sheet a colour according to
its band index and works even when projections were not exported. Projection
colouring can instead show a selected spin-channel weight or the weight of a
selected species (the species slider or ``--species Fe`` chooses the species)::

   python tools/plot_fermi_surface.py fermi_surface --all-bands --color-by spin
   python tools/plot_fermi_surface.py fermi_surface --all-bands --color-by species --species Fe

When projections were exported, ``P`` toggles projection colouring and keys
``1``--``5`` select site, orbital, spin, site/orbital/spin channel, or species
mode. The site, orbital, and species sliders select the corresponding channel;
the spin slider applies to the site/orbital/spin channel mode. Aggregate spin
colouring uses the signed polarization
``P=(w_up-w_down)/(w_up+w_down)``, fixed to ``[-1,1]``: ``-1`` is fully down,
``+1`` is fully up, and zero is unpolarized (or has zero weight). It uses the
``Reds`` colour map, from white through dark red. Species colouring sums all
orbitals and spin channels on sites with the selected element symbol. Press
``B`` to switch to band-index colouring.
Projection colours are computed directly from the exported eigenvectors using
the recorded site-major spin-blocked basis: ``s,p,d,f`` correspond to the
``1,3,5,7`` orbital blocks, and ``up/down`` to the two consecutive spin
blocks. This is the default ``--spin-frame basis`` behaviour and works even
when ``write_eigenpair_projections`` was not enabled. Use
``--spin-frame local`` to use the optional spin projections resolved along
each site's converged local moment axis.
The viewer defaults to Cartesian reciprocal coordinates, transforming the
uniform Monkhorst--Pack grid from fractional/direct coordinates with the
stored reciprocal basis vectors. Use ``--coordinate-system fractional`` to
inspect the direct-coordinate cell. BZ seam closure is enabled by default;
use ``--no-close-bz-seams`` for an open diagnostic plot. The old
``--periodic`` spelling remains accepted as an alias. A non-interactive image
can be produced with ``--screenshot surface.png``. Older sidecars without
stored reciprocal vectors can still be plotted in fractional mode, or in
Cartesian mode by supplying ``--reciprocal-vectors``.

For a collinear, no-spin-orbit calculation, pure global spin sheets can also
be exported. Enable the independent sector solve with
``write_spin_resolved_eigenpairs = .true.``::

   &control
      nsp = 1
   /
   &reciprocal
      write_spin_resolved_eigenpairs = .true.
   /

The writer checks the complete assembled ``H(k)`` (and ``S(k)`` in generalized
mode) before independently diagonalizing the up and down blocks. It writes
those eigenpairs to ``fermi_surface.spin.bin`` and records the extra payload
in the main sidecar. In the normal site-major basis, an ``N``-site ``spd``
calculation produces two ``9*N`` by ``9*N`` sector problems; the sector basis
is reconstructed from the interleaved per-site full basis. The option fails
closed if transverse terms or spin-orbit coupling are present, and is not
defined for a finite-q GBT spiral, noncollinear state, or BdG/Nambu problem.
The full combined export remains available and is unchanged.

Select the payload in the viewer with ``--spin-sheets combined`` (the default),
``--spin-sheets up``, ``--spin-sheets down``, or ``--spin-sheets both``::

   python tools/plot_fermi_surface.py fermi_surface --spin-sheets up --all-bands \
       --color-by band

Sector eigenvectors also support site, orbital, species, and global-spin
colouring. The local-moment frame is intentionally unavailable for these
independently solved global sectors.

For a clean, publication-style still image, enable the high-quality renderer
and choose the base resolution explicitly::

   python tools/plot_fermi_surface.py fermi_surface --all-bands \
       --color-by band --high-quality --window-size 1600 1200 \
       --screenshot fermi_surface.png

This uses a solid white background with black annotations, restrained Phong
diffuse lighting, soft screen-space ambient occlusion, FXAA, depth peeling
for translucent sheets, and a 2x output scale by default (3200x2400 for the
command above). Use ``--render-scale`` to change that multiplier,
``--background`` to override the background, and ``--ao-radius``,
``--ao-bias``, or ``--ao-kernel-size`` to tune the AO. The interactive widgets
and status text are omitted from high-quality screenshots, while the
coordinate frame and colour bar are retained. Some headless VTK builds
substitute SSAA for FXAA; the resulting image remains anti-aliased.

In the interactive viewer, the ``Axes``, ``Labels``, and ``Grid`` checkbox
buttons independently control the orientation marker, coordinate annotations,
and cube grid lines. If the marching-cubes surface is still visibly faceted,
enable scalar-field refinement, for example
``--smooth-interpolation 2``. This evaluates a cubic interpolation on twice
as many cells per direction before contouring; it is optional because it
increases memory and contouring time, and requires SciPy when enabled.

For a quick approximate nesting scan on the same export, use
``--nesting``::

   python tools/plot_fermi_surface.py fermi_surface --all-bands --nesting \
       --nesting-width 0.01 --nesting-top 10

This reports the strongest nonzero, non-equivalent ``q`` vectors from a
periodic FFT Fermi-surface joint-density-of-states (FS-JDOS) correlation,
both in direct and Cartesian reciprocal coordinates. By default all pairs of
Fermi-level bands contribute. Restrict the pair channel with
``--nesting-pairs intraband`` or ``--nesting-pairs interband``. These two
labels refer to equal or unequal *sorted eigenvalue indices* at each k-point;
they are not connectivity-tracked bands through avoided crossings. For an
export containing independently solved collinear sectors, ``--nesting-sheets
up``, ``down``, or ``both`` selects the spin sheets used in the correlation.
Use ``--nesting-pairs cross-spin --nesting-sheets both`` for the ordered
``N_up,down(q)`` channel. Unlike the other channels, its q and -q values are
not generally equal and are therefore reported separately.

The raw mesh sets the q resolution. To reduce mesh locking in the Gaussian
shells, ``--nesting-interpolation 4`` periodically zero-pads the Fourier
coefficients of each sampled band energy before evaluating the nesting map.
This is a spectral interpolant on the k-space torus; it adds q sampling but
does not add ab-initio information and can ring near sorted-band crossings.
The code warns when ``--nesting-width`` is substantially narrower than the
typical energy change between neighbouring raw k-points.

Add ``--nesting-only`` for a headless analysis that exits without opening the
PyVista viewer. Use ``--nesting-output nesting.npz`` to save the normalized map, fractional
and (when reciprocal vectors are available) Cartesian q-grids, metadata, and
the labels of bands contributing within five Gaussian widths of the selected
Fermi level. The reported map is
normalized to ``N(q=0)=1`` when possible. This is a useful structure-factor-
like screening diagnostic, not a full Lindhard susceptibility calculation:
the latter additionally requires occupations, an energy denominator and its
temperature/broadening prescription, and potentially matrix elements.

For a one-dimensional summary, ``--nesting-radial-output nesting_vs_q.png``
saves a high-resolution plot of the shell-averaged ``S(|q|)``. The horizontal
axis is the Cartesian reciprocal-vector length, so this option requires the
reciprocal vectors in the sidecar or via ``--reciprocal-vectors``. Adjust the
number of shell bins with ``--nesting-radial-bins``. This radial average is
useful for a quick overview, but directional peaks should be read from the
full three-dimensional ``nesting`` map or its angular-baseline-normalized
strongest-vector report.

The export itself does not apply Gaussian or tetrahedron smearing: a Fermi
surface is the ``E_n(k) = E_F`` isosurface of the sampled bands. If
``auto_find_fermi = .true.`` is selected, the Fermi level is determined from
the eigenvalue occupations on this dense mesh; otherwise the value from
``&energy`` is recorded. Tetrahedron (preferably Blöchl-improved) is the
better choice for a separate DOS or electron-count/Fermi-level calculation,
whereas Gaussian broadening is mainly useful for a smooth visual DOS.

Wired entry points
------------------

- ``post_processing='band_structure'`` routes through reciprocal-space workflow.
- ``post_processing='density_of_states'`` routes through reciprocal-space workflow.
- ``post_processing='fermi_surface'`` writes a dense full-BZ eigenpair payload.
- ``&self use_kspace=.true.`` enables an optional k-space SCF DOS/moment branch.

Hamiltonian Semantics
=====================

``second``
----------

- Builds the second-order base Hamiltonian:

  .. math::

     H(k) = h(k) - [hoh](k) + E_\nu + H_{LS}

- Adds the Bloch-summed ``Hcc_2c(k)`` when ``ccor_2c=.true.``.
- Solves a standard Hermitian eigenproblem with identity overlap.

``first``
---------

- Builds the first-order Hamiltonian ``h(k) + H_LS``.
- Adds the Bloch-summed ``Hcc_2c(k)`` when ``ccor_2c=.true.``.
- Useful for comparisons against the lower-order reciprocal-space path.

``proper``
----------

- Deprecated compatibility spelling for ``second``.

Minimal Example
===============

.. code-block:: fortran

   &calculation
      post_processing = 'band_structure'
   /

   &reciprocal
      kspace_ham_order = 'second'
      nk1 = 16, nk2 = 16, nk3 = 16
      use_symmetry_reduction = .true.
      use_time_reversal = .true.
      kspace_diagnostics = .true.
      suppress_internal_logs = .true.
   /

Optional SCF K-Space Branch
===========================

The SCF module has an optional k-space DOS/moment branch through ``&self``:

.. code-block:: fortran

   &self
      use_kspace = .true.
   /

Notes:

- Recursion remains the default production SCF path.
- With ``use_kspace=.true.``, DOS, band moments, Fermi level, and LDA+U LDM
  updates are taken from the reciprocal-space workflow.
- This branch is best suited for bulk validation runs.

For Chebyshev bounds used by recursion/KPM, see :ref:`theory/recursion_method`.
