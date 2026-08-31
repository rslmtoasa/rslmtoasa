# RS-LMTO-ASA

Real-space Linearized Muffin-Tin Orbital method in the Atomic Sphere Approximation.

## New Features

### libXC Integration

RS-LMTO-ASA provides an optional, validated [libXC](https://libxc.gitlab.io/)
interface for spherical-ASA LDA and GGA functionals while preserving the
historical XC implementation.  Meta-GGA, hybrid, orbital-dependent, and
kinetic-density-dependent families are rejected explicitly because the ASA
interface does not provide all of their required ingredients.

```bash
cmake .. -DENABLE_LIBXC=ON
make

# Predefined PBE alias, or a direct native libXC ID:
echo "txc = 108" > input.nml
echo "txc = 1101" > input.nml  # native libXC ID 101, exchange-only PBE GGA
```

See [the production contract](docs/LIBXC_PRODUCTION_CONTRACT.md) for selector
namespaces, spin and unit conventions, radial GGA evaluation, aliases, and
validation evidence.

### TDB

<!--
**rslmtoasa/rslmtoasa** is a ✨ _special_ ✨ repository because its `README.md` (this file) appears on your GitHub profile.

Here are some ideas to get you started:

- 🔭 I’m currently working on ...
- 🌱 I’m currently learning ...
- 👯 I’m looking to collaborate on ...
- 🤔 I’m looking for help with ...
- 💬 Ask me about ...
- 📫 How to reach me: ...
- 😄 Pronouns: ...
- ⚡ Fun fact: ...
-->
