# Schrodinger's Sandbox

A native desktop application for interactive quantum mechanical simulation and 3D visualisation of atomic orbitals, molecular orbitals, and coordination compounds. Built from first principles in C++17 with real-time GPU rendering via OpenGL.

Unlike existing educational chemistry software that either provides toy-level approximations or exposes impenetrable research-grade interfaces, Schrodinger's Sandbox occupies the gap between: it solves real quantum mechanics and renders the results in real time, while remaining accessible to a first-year chemistry student.

---

## Table of Contents

- [Features](#features)
- [Screenshots](#screenshots)
- [Theory](#theory)
- [Architecture](#architecture)
- [Prerequisites](#prerequisites)
- [Building from Source](#building-from-source)
- [Usage](#usage)
- [Project Structure](#project-structure)
- [Computational Methods](#computational-methods)
- [Rendering Pipeline](#rendering-pipeline)
- [Testing](#testing)
- [Roadmap](#roadmap)
- [Known Limitations](#known-limitations)
- [References](#references)
- [License](#license)

---

## Features

### Atomic Orbital Viewer
- Interactive periodic table with all 118 elements and colour-coded categories
- Hydrogen-like wavefunctions with Slater screening for effective nuclear charges
- Real-time 3D orbital rendering via GPU volumetric ray marching
- Three visualisation modes: volume rendering, isosurface, and phase-coloured isosurface
- Adjustable isosurface threshold and volume rendering gamma
- Full arcball camera with rotate, zoom, and pan
- Radial probability density plots for every occupied subshell
- Energy level diagrams with electron occupancy
- Spectroscopic notation for all electron configurations, including Aufbau exceptions

### Molecular Orbital Viewer (Planned)
- Interactive molecule builder with drag-and-drop construction
- Hartree-Fock self-consistent field solver
- Kohn-Sham density functional theory (B3LYP, PBE, PBE0)
- Molecular orbital correlation diagrams
- Mulliken and Lowdin population analysis

### Coordination Chemistry (Planned)
- GFN2-xTB tight-binding for large systems (50-500+ atoms)
- d-orbital splitting diagrams
- Spectrochemical series verification from computed splittings
- Ligand library with common mono- and polydentate ligands
- Coordination geometry templates (octahedral, tetrahedral, square planar)

### Reactions (Planned)
- Potential energy surface scans along reaction coordinates
- Transition state finding via nudged elastic band
- Orbital evolution animation along reaction paths

---

## Theory

### Single-Atom Orbitals

For isolated atoms, the application uses analytical hydrogen-like wavefunctions:

```
psi_nlm(r, theta, phi) = R_nl(r) * Y_lm(theta, phi)
```

where R_nl are radial wavefunctions constructed from associated Laguerre polynomials and Y_lm are real spherical harmonics constructed from associated Legendre polynomials.

Multi-electron atoms are treated in the independent-particle approximation with effective nuclear charges computed via Slater's rules. Each electron experiences a screened nuclear charge:

```
Z_eff = Z - sigma
```

where sigma is the total screening constant from all other electrons, computed according to Slater's grouping and shielding rules. Orbital energies are then:

```
E_n = -13.6 * (Z_eff / n)^2 eV
```

All 118 elements include their correct NIST ground-state electron configurations, with hardcoded Aufbau exceptions for Cr, Cu, Nb, Mo, Ru, Rh, Pd, Ag, Pt, Au, and the relevant lanthanides and actinides.

### Molecular Orbitals (Planned)

For molecules, the application will solve the Roothaan-Hall equations iteratively:

```
FC = SCe
```

using Gaussian basis sets (STO-3G, 6-31G*, cc-pVDZ) with integrals computed via libint2. Convergence acceleration via DIIS, level shifting, and density damping.

### Large Systems (Planned)

For coordination compounds and large molecular systems, GFN2-xTB tight-binding (via the tblite library) provides qualitatively correct orbital pictures at a fraction of the ab initio cost.

---

## Architecture

```
                    +---------------------------+
                    |     UI / Application      |
                    |  ImGui panels, periodic   |
                    |  table, orbital browser   |
                    +------------+--------------+
                                 |
                    +------------+--------------+
                    |       Renderer            |
                    |  OpenGL 4.1, GLSL ray     |
                    |  marching, isosurfaces,   |
                    |  arcball camera            |
                    +------------+--------------+
                                 |
          +----------------------+----------------------+
          |                      |                      |
+---------+--------+  +----------+---------+  +---------+---------+
| Analytical Atom  |  |   HF / DFT Engine  |  | Semi-Empirical    |
| Solver           |  |   (Planned)        |  | xTB (Planned)     |
| Slater screening |  |   libint2 + libxc  |  | tblite            |
+------------------+  +--------------------+  +-------------------+
          |                      |                      |
          +----------------------+----------------------+
                                 |
                    +------------+--------------+
                    |     Core / Math           |
                    |  Special functions,       |
                    |  basis set types,         |
                    |  physical constants       |
                    +---------------------------+
```

The application is structured as a set of static libraries linked into a single executable. Dependencies flow strictly downward. The UI layer never directly calls integral routines; it communicates through an orchestrator that selects the appropriate computational method based on system size.

---

## Prerequisites

### macOS

```
brew install cmake eigen python3
```

Xcode Command Line Tools must be installed:

```
xcode-select --install
```

### Linux (Ubuntu/Debian)

```
sudo apt install cmake build-essential pkg-config libeigen3-dev \
  libxrandr-dev libxinerama-dev libxcursor-dev libxi-dev libgl-dev python3-pip
```

### Windows

- Visual Studio 2022 with the C++ desktop workload
- CMake 3.20+
- Git
- Python 3.8+
- Eigen via vcpkg: `vcpkg install eigen3:x64-windows`

---

## Building from Source

### 1. Clone the Repository

```
git clone https://github.com/<your-username>/schrodingerssandbox.git
cd schrodingerssandbox
```

### 2. Set Up External Dependencies

Run the setup script to clone GLFW, Dear ImGui, GoogleTest, and generate the glad loader:

```
./setup.sh
```

On Windows, use:

```
setup_windows.bat
```

If the setup script fails (for instance if glad2 is not installed), you can set up dependencies manually:

```
git clone --depth 1 --branch 3.4 https://github.com/glfw/glfw.git external/glfw
git clone --depth 1 --branch docking https://github.com/ocornut/imgui.git external/imgui
git clone --depth 1 --branch v1.14.0 https://github.com/google/googletest.git external/googletest

pip3 install glad2
python3 -m glad --api gl:core=4.6 --out-path external/glad c
```

### 3. Configure and Build

```
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j8
```

On Windows with Visual Studio:

```
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

### 4. Run Tests

```
ctest --output-on-failure
```

### 5. Launch

```
./schrodingers_sandbox
```

---

## Usage

### Exploring Atomic Orbitals

1. Open the application. The periodic table is displayed in the bottom panel.
2. Click any element to select it. The orbital browser on the left updates to show all occupied subshells.
3. Click a subshell in the orbital browser (e.g., "3d") to visualise it in the 3D viewport.
4. Use the magnetic quantum number buttons to switch between m values (e.g., m = -2, -1, 0, +1, +2 for d orbitals).
5. Use the render mode dropdown in the viewport to switch between volume, isosurface, and phase-coloured isosurface rendering.
6. Adjust the gamma slider (volume mode) or iso value slider (isosurface mode) to control the visualisation.

### Camera Controls

- Left-click and drag in the viewport to rotate
- Scroll to zoom in/out
- Middle-click and drag to pan

### Search

Use the search bar at the top of the periodic table to filter elements by name, symbol, or atomic number.

---

## Project Structure

```
schrodingerssandbox/
|-- CMakeLists.txt
|-- setup.sh                        # Dependency setup (Linux/macOS)
|-- setup_windows.bat               # Dependency setup (Windows)
|-- README.md
|-- external/
|   |-- glfw/                       # GLFW 3.4 (windowing + OpenGL context)
|   |-- imgui/                      # Dear ImGui docking branch (UI)
|   |-- glad/                       # glad2 OpenGL 4.6 core loader
|   +-- googletest/                 # GoogleTest v1.14.0
|-- data/
|   +-- shaders/
|       |-- fullscreen_quad.vert    # Fullscreen triangle vertex shader
|       |-- test_gradient.frag      # Test gradient (fallback)
|       |-- orbital_raymarch.vert   # Orbital rendering vertex shader
|       +-- orbital_raymarch.frag   # Orbital rendering fragment shader
|-- src/
|   |-- core/
|   |   |-- constants.h             # Physical constants (Bohr radius, Hartree, etc.)
|   |   |-- types.h                 # Type aliases and forward declarations
|   |   |-- special_functions.h/cpp # Laguerre, Legendre, spherical harmonics
|   |   |-- hydrogen.h/cpp          # Hydrogen-like radial wavefunctions
|   |   |-- slater.h/cpp            # Slater's rules for Z_eff
|   |   +-- elements.h/cpp          # Element database (all 118 elements)
|   |-- renderer/
|   |   |-- window.h/cpp            # GLFW window wrapper
|   |   |-- shader.h/cpp            # Shader compilation and uniform management
|   |   +-- camera.h/cpp            # Arcball camera
|   |-- ui/
|   |   |-- app.h/cpp               # Main application loop and orchestration
|   |   |-- app_state.h             # Shared application state
|   |   |-- periodic_table.h/cpp    # Interactive periodic table widget
|   |   |-- orbital_browser.h/cpp   # Orbital selection and quantum numbers
|   |   +-- properties_panel.h/cpp  # Element info and radial probability plot
|   +-- main.cpp
+-- tests/
    |-- test_constants.cpp           # Physical constant validation
    |-- test_hydrogen.cpp            # Radial wavefunctions and normalisation
    +-- test_slater.cpp              # Slater Z_eff values and element database
```

---

## Computational Methods

### Currently Implemented

| Method | Scope | Basis | Speed |
|--------|-------|-------|-------|
| Analytical hydrogen-like | Single atoms (Z = 1-118) | Exact (with Slater screening) | Instant (GPU) |

### Planned

| Method | Scope | Basis | Speed |
|--------|-------|-------|-------|
| Restricted Hartree-Fock | Small molecules (2-15 heavy atoms) | STO-3G, 6-31G*, cc-pVDZ | < 5 seconds |
| Unrestricted Hartree-Fock | Open-shell molecules, radicals | STO-3G, 6-31G*, cc-pVDZ | < 10 seconds |
| DFT (B3LYP, PBE, PBE0) | Medium molecules (2-20 heavy atoms) | 6-31G*, cc-pVDZ | < 30 seconds |
| GFN2-xTB | Large systems (up to 500 atoms) | Parameterised tight-binding | < 5 seconds |

---

## Rendering Pipeline

### Volume Rendering

For each pixel, a ray is cast from the camera through the scene. The ray is marched through a bounding sphere centred on the nucleus (or molecular centre of mass) in 192 steps. At each sample point, the wavefunction psi(r, theta, phi) is evaluated directly on the GPU using GLSL implementations of the associated Laguerre and Legendre recurrences. The probability density |psi|^2 is mapped to colour and opacity via a transfer function with adjustable gamma, and accumulated using front-to-back compositing.

### Isosurface Rendering

The ray march detects where |psi|^2 crosses a user-defined threshold. A bisection refinement (5 iterations) locates the surface precisely. The surface normal is computed via central differences and lit with Phong shading (ambient + diffuse + specular + rim light). In phase mode, the surface is coloured by the sign of psi (blue for positive, red for negative), allowing students to identify bonding and antibonding character.

### Performance

All orbital rendering targets 60 fps at 1080p resolution for orbitals up to n = 7 on any GPU supporting OpenGL 4.1. The fragment shader is the sole bottleneck; adaptive step sizing and early ray termination keep performance consistent.

---

## Testing

The test suite uses GoogleTest and covers four areas:

```
ctest --output-on-failure
```

| Test Suite | What It Validates |
|------------|-------------------|
| test_constants | Physical constants (Bohr radius, Hartree energy) match NIST values |
| test_special_functions | Associated Laguerre and Legendre polynomials, real spherical harmonics against analytical values |
| test_hydrogen | Radial wavefunctions (values at specific points, nodes, normalisation integrals via Simpson's rule) |
| test_slater | Slater Z_eff for H, He, Li, C, Na, Fe against published values; element database integrity; Aufbau exceptions |

All floating-point comparisons use EXPECT_NEAR with tolerances appropriate to each quantity (1e-10 for special functions, 1e-6 for wavefunctions, 1e-2 for Slater values which are inherently approximate).

---

## Roadmap

The project is developed in phases. Each phase produces a working application with strictly increasing capability.

| Phase | Description | Status |
|-------|-------------|--------|
| 0 | Project skeleton, build system, OpenGL window, ImGui docking | Complete |
| 1 | Atomic orbital viewer (analytical wavefunctions, GPU ray marching, periodic table) | Complete |
| 2 | Gaussian basis set infrastructure and integral engine (libint2 wrapper) | Planned |
| 3 | Hartree-Fock SCF solver (RHF, UHF, DIIS convergence) | Planned |
| 4 | Molecule builder and molecular orbital visualisation on GPU | Planned |
| 5 | Density functional theory (numerical grid, libxc, B3LYP/PBE/PBE0) | Planned |
| 6 | Semi-empirical methods (GFN2-xTB via tblite) and coordination chemistry tools | Planned |
| 7 | Geometry optimisation, vibrational analysis, file I/O, packaging | Planned |
| 8 | Potential energy surfaces and transition state finding | Planned |
| 9 | Ab initio molecular dynamics and reaction trajectory playback | Planned |

A detailed implementation roadmap with algorithms, data structures, milestones, and risk analysis is available in `docs/qcw-roadmap.docx`.

---

## Known Limitations

- macOS limits OpenGL to version 4.1, so compute shaders and shader storage buffer objects are not available. All GPU computation uses fragment shaders. A future version may transition to MoltenVK (Vulkan on Metal) to lift this restriction.
- The Slater screening model is approximate. Effective nuclear charges and orbital energies for heavy elements will differ from Hartree-Fock or DFT values. This is by design for the atomic viewer; the planned HF/DFT engines provide quantitative accuracy for molecules.
- Electron configuration display uses simple spectroscopic notation. Term symbols (Russell-Saunders coupling) and relativistic effects (spin-orbit splitting) are not shown.
- The current atomic viewer treats each electron independently. Electron-electron repulsion, exchange, and correlation are accounted for only through Z_eff screening, not through self-consistent solution of the many-electron problem. This limitation is addressed in Phases 3-6.

---

## References

### Textbooks

- A. Szabo and N. S. Ostlund, "Modern Quantum Chemistry: Introduction to Advanced Electronic Structure Theory" (Dover, 1996). The primary reference for Hartree-Fock implementation.
- T. Helgaker, P. Jorgensen, and J. Olsen, "Molecular Electronic-Structure Theory" (Wiley, 2000). Comprehensive reference for integral evaluation and DFT.
- D. J. Griffiths and D. F. Schroeter, "Introduction to Quantum Mechanics" (Cambridge, 3rd ed., 2018). For hydrogen-like wavefunctions and angular momentum theory.

### Libraries

- [libint2](https://github.com/evaleev/libint) -- Gaussian integral evaluation
- [libxc](https://gitlab.com/libxc/libxc) -- Exchange-correlation functionals
- [tblite](https://github.com/tblite/tblite) -- GFN2-xTB tight binding
- [Eigen](https://eigen.tuxfamily.org) -- Linear algebra
- [Open Babel](https://openbabel.org) -- Molecular file formats and force fields
- [Dear ImGui](https://github.com/ocornut/imgui) -- Immediate-mode GUI
- [GLFW](https://www.glfw.org) -- Window and OpenGL context management
- [Basis Set Exchange](https://www.basissetexchange.org) -- Gaussian basis set data

### Reference Codes

- [PySCF](https://github.com/pyscf/pyscf) -- Python quantum chemistry, excellent for verification
- [ERKALE](https://github.com/susilehtola/erkale) -- Clean, readable HF/DFT implementation
- [Crawford Group Programming Projects](https://github.com/CrawfordGroup/ProgrammingProjects) -- Step-by-step SCF tutorials with test data

---

## License

[Specify your license here]
