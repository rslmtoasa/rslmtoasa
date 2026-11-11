# Bundled RS-utils release-prep wheels

These cibuildwheel archives were supplied by the package author for inclusion
in the Codespaces tutorial. They distribute `rslmto-utils` 0.11.0rc1 from the
`release-prep` branch; the CLI and Python import are named `rsutil`.

- `wheels-linux-x86_64.zip`: CPython 3.11, 3.12, 3.13, and 3.14 Linux x86_64 wheels.
- `wheels-linux-aarch64.zip`: the corresponding ARM64 wheels.

Codespaces setup installs the matching wheel from these local archives,
checks its version and native library, and adds the venv CLI to terminal PATH.
The repository clone includes the package; students do not upload it or need
access to its source repository. Python dependencies are resolved from PyPI.
Package license/notice files are contained in the supplied wheels.

Pinned SHA-256 values, also checked by `.devcontainer/install-rsutils.sh`:

```text
63f1c9b13462c43127d0c375ddda6c1f0749cbd2df96669b7260576250efc826  wheels-linux-x86_64.zip
e96353595e8aa767b947e439d89a105351bb45703c0e3974009f64d385b3492b  wheels-linux-aarch64.zip
```

To update, replace the author-provided archives and update the installer's
checksums and expected version together. An instructor can test an individual
wheel or the pinned release-prep source ZIP by passing its path explicitly.
