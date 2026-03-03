.. _user_guide/examples:

========
Examples
========

This section provides comprehensive, step-by-step examples of RS-LMTO-ASA calculations. Each example includes complete input files with detailed explanations to enable full reproduction of the results.

Overview
========

All examples are located in the ``example/`` directory of the RS-LMTO-ASA distribution. The examples cover:

**Self-Consistent Field (SCF) Calculations:**

- :doc:`examples/bulk_si` - Diamond cubic silicon (non-magnetic bulk)
- :doc:`examples/bulk_bccfe` - Body-centered cubic iron (magnetic bulk)
- :doc:`examples/bulk_mn3sn` - Mn₃Sn kagome antiferromagnet (non-collinear magnetism)
- :doc:`examples/surface_fcccu001` - Cu(001) surface (semi-infinite system)
- :doc:`examples/impurity_b2feco` - Fe impurity in B2 FeCo (embedded cluster)

**Post-Processing Calculations:**

- :doc:`examples/exchange_bccfe` - Magnetic exchange interactions J_ij (LKAG method)
- :doc:`examples/conductivity_fcccu` - Electrical conductivity of copper
- :doc:`examples/conductivity_fccpt` - Spin Hall effect in platinum
- :doc:`examples/conductivity_bccfe` - Anomalous Hall effect in iron

Running the Examples
====================

All examples assume you have compiled the code as described in :doc:`../getting_started`. The executable is located at ``build/bin/rslmto.x``.

**Basic workflow:**

.. code-block:: bash

   cd example/<category>/<system>
   ../../../build/bin/rslmto.x    # Reads input.nml by default

The code reads ``input.nml`` from the current directory and writes output files with ``_out.nml`` suffix.

Input File Structure
====================

Each example uses three types of input files:

1. **input.nml** - Main calculation parameters (lattice, energy, convergence, etc.)
2. **lattice.nml** - Crystal structure (auto-generated or user-provided)
3. **atom.nml** - Tight-binding parameters for each atom type (e.g., ``Fe.nml``, ``Si1.nml``)

The following sections explain each example in detail.

.. toctree::
   :maxdepth: 1
   :caption: SCF Examples:

   examples/bulk_si
   examples/bulk_bccfe
   examples/bulk_mn3sn
   examples/surface_fcccu001
   examples/impurity_b2feco

.. toctree::
   :maxdepth: 1
   :caption: Post-Processing Examples:

   examples/exchange_bccfe
   examples/conductivity_fcccu
   examples/conductivity_fccpt
   examples/conductivity_bccfe

Best Practices
==============

**Parameter Selection Guidelines:**

- **Cluster radius (rc):** 40-60 a.u. for bulk, 80-120 a.u. for surfaces
- **Mixing parameter (beta):** 0.1-0.5 for stable systems, 0.001-0.01 for surfaces
- **SCF iterations (nstep):** 25-50 for non-magnetic, 100-500 for magnetic/surface
- **Recursion levels (lld):** 21-50 for simple systems, 150-300 for conductivity

**Convergence Criteria:**

Check that the final iterations show:

.. code-block:: bash

   grep "Total energy" input_out.nml | tail -5   # Should vary by < 0.001 Ry
   grep "Total charge" *_out.nml                  # Should match valence electrons
   grep "Fermi" input_out.nml | tail -5           # Should be stable to 0.001 Ry

**Common Issues:**

- **SCF not converging:** Reduce ``beta``, increase ``nstep``
- **Wrong total charge:** Adjust ``ct`` parameter or increase ``rc``
- **Slow performance:** Reduce ``lld`` or ``channels_ldos``, use OpenMP

See Also
========

- :doc:`input_files` - Complete input file format reference
- :doc:`output_files` - Output file descriptions
- :doc:`../keywords/index` - All input keywords
- :doc:`../theory/index` - Theoretical background
