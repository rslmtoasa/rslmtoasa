.. _user_guide_index:

=======================================
User Guide
=======================================

This section provides practical guidance for running RS-LMTO-ASA calculations.

.. toctree::
   :maxdepth: 2

   input_files
   output_files
   examples

Quick Start
===========

**1. Prepare input file** (namelist format)

Create ``input.nml`` (the default input filename) with structure, control parameters, and settings.

**2. Run calculation**

.. code-block:: bash

   ./rslmto.x

**3. Analyze results**

Output files contain energies, charge densities, magnetic moments, and requested properties.

Key Sections in User Guide
===========================

- **Input Files** (:doc:`input_files`) - Complete namelist syntax and all parameters
- **Output Files** (:doc:`output_files`) - Understanding calculation output
- **Examples** (:doc:`examples`) - Worked examples and common calculations

Before You Start
================

Ensure you have:

✓ Compiled RS-LMTO-ASA (see :doc:`../getting_started`)
✓ Basic understanding of LMTO-ASA (see :ref:`theory/lmto_asa_overview`)
✓ Example files (in ``example/`` directory)
✓ Element/potential database files (required for your system)

Common Tasks
============

**Run a bulk crystal calculation:**

See ``example/bulk/Si/Si1.nml`` for an example (copy to ``input.nml`` or specify as argument).

**Calculate magnetic properties:**

Set ``nsp=1`` (scalar relativistic) or ``nsp=2`` (collinear magnetism).
Provide exchange field or constrain moments.

**Compute band structure:**

Include ``post_processing='dos'`` and appropriate energy mesh.

**Surface calculation:**

Use ``pre_processing='buildsurf'`` to generate surface geometry.

**Magnetic impurity:**

Use ``pre_processing='newclusurf'`` or ``newclubulk`` to insert impurity.

Getting Help
============

- Check keyword documentation: :ref:`keywords/index`
- Review example input files in ``example/``
- Enable verbose output: ``verbose = .true.``
- Check log output for error messages

See Also
========

- :doc:`../getting_started` - Installation and running
- :doc:`../code_structure` - Code organization
- :ref:`theory/lmto_asa_overview` - Theoretical background
