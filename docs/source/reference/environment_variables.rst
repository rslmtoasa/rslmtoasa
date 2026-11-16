.. _reference/environment_variables:

=====================
Environment Variables
=====================

The code reads the following environment variables. None is needed for a normal run;
all are optional.

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Variable
     - Effect
   * - ``RSLMTO_CUDA_DEVICE``
     - CUDA device index (non-negative, in the ``CUDA_VISIBLE_DEVICES`` numbering) for the
       CUDA recursion plugin and the CUDA reciprocal backend. Overrides the per-rank device
       mapping; invalid text stops the plugin run and is reported by the reciprocal backend.
   * - ``RSLMTO_KPM_RANDOM_SEED``
     - Integer seed of the random vectors of the stochastic transport
       (``cond_calctype = 'random_vec'``), set once per transport request so that runs can
       be repeated. Unset keeps the compiler's default random seed.
   * - ``RSLMTO_SPINOR_TRACE``
     - File name; when set, rank 0 writes the spinor d-shell potential and occupation of every
       recursion site to that file.
   * - ``RSLMTO_FIXED_U_TRACE``
     - File name of a CSV trace of the fixed-U SCF (occupation matrices, potentials, moments per
       iteration). Active only for ``nsp = 1``, ``lmax = 2``, ``hoh = .false.`` and fixed
       Liechtenstein U; with no U it records the U = 0 reference.
   * - ``RSLMTO_FIXED_U_RUN_ID``
     - Run identifier written into every row of the fixed-U trace.
   * - ``RSLMTO_FIXED_U_LAMBDA``
     - Real value written into the ``lambda`` column of the fixed-U trace.
   * - ``RSLMTO_FIXED_U_TRACE_ITERATION``
     - SCF iteration to trace; 0 keeps only the last completed iteration. A negative or
       unreadable value switches the trace off.
   * - ``RSLMTO_FIXED_U_TRACE_CAPTURE_ONLY``
     - ``1`` or ``true``: capture the radial state for the trace without the extra replay
       check.
   * - ``RSLMTO_FIXED_U_TRACE_MOMENTS_ONLY``
     - ``1`` or ``true``: leave the Hamiltonian rows out of the trace.
   * - ``RSLMTO_FIXED_U_TRACE_QL_ONLY``
     - ``1`` or ``true``: write only the ``ql`` rows of the output stage.

**Related code:** ``source/mpi.f90`` (``get_cuda_device_override``),
``source/recursion_transport.f90``, ``source/ldau_diagnostics.f90``, ``source/bands.f90``.
