.. RS-LMTO-ASA Documentation master file

=======================================
RS-LMTO-ASA Documentation
=======================================

**RS-LMTO-ASA** is a massively parallel real-space linear muffin-tin orbital 
(LMTO) code based on the atomic sphere approximation (ASA). The code implements 
density functional theory (DFT) within the LMTO formalism to perform first-principles 
electronic structure calculations of bulk materials, surfaces, and clusters with support 
for magnetic, non-collinear, and relativistic calculations.

This documentation provides:

- **Getting Started**: Installation, compilation, and basic usage
- **Theory**: Overview of the LMTO-ASA method, Green's function formalism, and key algorithms
- **User Guide**: Input file syntax, running calculations, and interpreting outputs
- **Reference**: Module and routine documentation, code structure overview
- **Keywords**: Comprehensive listing of all input parameters and their meanings

Table of Contents
=================

.. toctree::
   :maxdepth: 2
   :caption: Introduction

   getting_started
   code_structure

.. toctree::
   :maxdepth: 2
   :caption: Theory and Methods

   theory/index

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   user_guide/index

.. toctree::
   :maxdepth: 2
   :caption: Input Keywords

   keywords/index

.. toctree::
   :maxdepth: 2
   :caption: Reference

   reference/index

Indices and Tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Quick Links
===========

- `GitHub Repository <https://github.com/rslmto/rslmto_devel>`_
- `DFT Fundamentals <https://theory.physics.umass.edu/kieron/>`_
- `LMTO Reviews <https://doi.org/10.1103/RevModPhys.84.1419>`_

Citation
=========

If you use RS-LMTO-ASA in your research, please cite:

.. code-block:: bibtex

   @article{rslmto2024,
     title={RS-LMTO-ASA: A Massively Parallel Real-Space LMTO Code},
     author={Development Team},
     journal={To be published},
     year={2024}
   }

.. note::

   This documentation is automatically generated from the source code.
   For the latest development version, visit the 
   `GitHub repository <https://github.com/rslmto/rslmto_devel>`_.
