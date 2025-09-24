# RS-LMTO-ASA

Real-space Linearized Muffin-Tin Orbital method in the Atomic Sphere Approximation.

## New Features

### libXC Integration 🎉

RS-LMTO-ASA now supports integration with the [libXC library](https://libxc.gitlab.io/), providing access to over 600 exchange-correlation functionals while maintaining full backward compatibility.

**Key Benefits:**
- Access to modern functionals (SCAN, B3LYP, PBE0, etc.)
- Automatic updates as new functionals are added to libXC
- Validation against standard implementations
- Optional dependency - existing calculations remain unchanged

**Quick Start:**
```bash
# Build with libXC support
cmake .. -DENABLE_LIBXC=ON
make

# Use libXC functionals (txc = 1000 + libXC_ID)
echo "txc = 1135" > input.nml  # PBE via libXC
echo "txc = 1402" > input.nml  # B3LYP hybrid
```

📖 **[Complete Integration Guide](LIBXC_INTEGRATION_GUIDE.md)**

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
