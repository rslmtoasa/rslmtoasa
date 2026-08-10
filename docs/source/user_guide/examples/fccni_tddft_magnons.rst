.. _examples/fccni_tddft_magnons:

======================================
Lightweight fcc Nickel TDDFT Magnons
======================================

**Location:** ``example/susceptibility/fccNi/``

This example uses the one-atom fcc-Ni block-recursion SCF restart in
``SCF/`` and a separate transverse LR-TDDFT post-processing deck in
``TDDFT/``.  The TDDFT input reads ``SCF/Ni_out.nml`` without modifying it,
then writes bare χ_KS, Dyson-enhanced response, Goldstone diagnostics, and
mode analysis in its own directory.

The deliberately light settings (4×4×4 k mesh, three q points, and 81
frequencies) make it useful as an interface demonstration.  They are not
converged magnon parameters.  The stored SCF state has a 0.594512 μB
spin moment per Ni site.

Run the SCF deck first and then the response deck:

.. code-block:: bash

   cd example/susceptibility/fccNi/SCF
   ../../../../build_13/bin/rslmto.x
   cd ../TDDFT
   ../../../../build_13/bin/rslmto.x

For material results, converge the SCF moment and Fermi level, then increase
the k mesh, response band window, q path, and frequency resolution.  Reduce
and extrapolate the numerical η broadening before assigning linewidths.  See
the ``TDDFT/README.md`` for output names and run details.
