.. _theory/gpu_acceleration:

================================
GPU Acceleration (CUDA plugin)
================================

Overview
========

The recursion kernels — the computational hot path of a real-space RS-LMTO
calculation — can optionally be offloaded to an NVIDIA GPU through a CUDA
*plugin*. The plugin is a self-contained backend: a Fortran ``iso_c_binding``
wrapper (``source/rsrec_cuda_plugin.f90``) over a C/CUDA library
(``source/cuda/rsrec_gpu.cu`` + ``rsrec_cuda.h``, dispatched through
``rsrec_cuda.cpp``). When enabled and available, it replaces the CPU recursion
kernels; when unavailable it falls back transparently to the CPU path, so the
same input runs on both GPU and non-GPU builds.

.. note::

   The CUDA plugin is **off by default** and must be compiled in explicitly
   (see below). A build without it behaves exactly as before — the ``gpu_*``
   keywords are simply ignored on the CPU path.

Building with the CUDA plugin
=============================

Configure with ``ENABLE_CUDA_PLUGIN=ON`` (requires the CUDA toolkit —
``cudart``, ``cusparse``, ``cufft``, ``cublas``):

.. code-block:: bash

   cmake -DENABLE_CUDA_PLUGIN=ON ..
   make -j 4

This defines ``USE_CUDA_PLUGIN`` and links the ``rsrec_cuda`` backend into
``rslmto.x``. See :ref:`getting_started` for the full build-options table.

Enabling the GPU at runtime
===========================

Two ``&control`` keywords select the GPU path:

.. list-table::
   :widths: 25 20 55
   :header-rows: 1

   * - Keyword
     - Default
     - Meaning
   * - ``gpu_plugin``
     - ``.false.``
     - Master switch. When ``.true.`` (and the plugin is compiled in and a
       device is available), the supported recursion kernels dispatch to the
       GPU. When ``.false.``, or when no device is available, the CPU path is
       used.
   * - ``gpu_backend``
     - ``'csr'``
     - Sparse storage / matvec strategy on the device. One of ``'csr'``,
       ``'bsr'``, ``'fft'``, ``'conv'`` (see below).

Example:

.. code-block:: fortran

   &control
      recur       = 'chebyshev'
      gpu_plugin  = .true.
      gpu_backend = 'bsr'
   /

``gpu_backend`` values
----------------------

- ``'csr'`` — compressed sparse row matvec (default, general).
- ``'bsr'`` — block sparse row; groups per-site orbital blocks, typically
  faster for the dense per-atom blocks of the LMTO Hamiltonian.
- ``'fft'`` — FFT-based periodic operator apply. **Requires a fully periodic
  lattice** (``pbc=.true.`` with all of ``b1``/``b2``/``b3``); otherwise the
  run warns and falls back to the CPU path.
- ``'conv'`` — real-space convolution periodic operator apply; same periodicity
  requirement as ``'fft'``.

An unrecognized ``gpu_backend`` is not fatal: the run logs a warning and falls
back to the current CPU recursion path.

Which kernels are offloaded
===========================

The plugin covers the recursion features that dominate runtime. Each is
guarded by ``gpu_plugin_ready`` (``recursion_core.f90``), which checks device
availability, backend compatibility, and per-feature constraints before
choosing the GPU path:

- Chebyshev moment recursion (``chebyshev_recur``)
- Block-Lanczos recursion (``recur_b``)
- Scalar Lanczos recursion (``recur``) — scalar-relativistic only (``nsp=1``)
- Stochastic Chebyshev moments (``compute_moments_stochastic``)
- Intersite / transport variants (``recur_b_ij``, ``chebyshev_recur_ij``)
- Orbital moments (``chebyshev_orbital_mod``)

When ``gpu_plugin=.true.``, these bypass the CPU Chebyshev/Lanczos backends
(``chebyshev_fast.f90`` / ``haydock_fast.f90``) entirely — the GPU backend is
not layered on top of them.

Fallback and correctness
========================

The GPU path is designed to be **numerically consistent** with the CPU path,
not a separate approximation: it computes the same recursion coefficients /
moments to the same tolerances. If any precondition fails (plugin not compiled
in, no device, incompatible backend, feature guard such as the ``nsp=1``
restriction on scalar Lanczos or the periodicity requirement on ``fft``/``conv``),
the run emits a warning and continues on the CPU path rather than aborting.

Device data is fingerprinted and cached (context creation, BSR upload) so that
repeated calls with the same Hamiltonian — e.g. across energy points — do not
re-upload or re-pack, mirroring the CPU ``cheb_cache_t`` lifecycle.

Testing
=======

Real-GPU consistency is validated out-of-CI via ``tests/run_gpu_matrix.sh``
(compares GPU output against the CPU reference across the recursion matrix).
CI itself only *compiles* the plugin (the ``cuda_compile`` job) since the
runners have no GPU. See ``docs/DEVELOPER_MAP.md`` for the kernel-level plugin
surface and ``tests/README.md`` for the GPU coverage tiers.

See Also
========

- :ref:`theory/recursion_method` — the recursion kernels the plugin accelerates
- :ref:`code_structure` — where the plugin sits in the source tree
- :ref:`getting_started` — build-options table
