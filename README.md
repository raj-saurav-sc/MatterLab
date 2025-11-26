# MatterLab: Universal Material & Fluid Simulator

**MatterLab** is a multi-tiered physics simulator designed for both educational exploration and scientific research. It allows users to interact with various states of matter and physical phenomena in a real-time 3D environment with advanced visualization and analysis tools.

## 🌟 Features

### 🏗️ Solid Mechanics
*   **Elastic Deformation**: Visualize stress and strain using Hooke's Law
*   **Plasticity & Fracture**: Simulate permanent deformation and material failure
*   **FEM Solver**: Finite Element Method with truss elements and dynamic relaxation
*   **Non-linear Materials**: Yield strength, strain hardening, and fracture mechanics

### 💧 Fluid Mechanics
*   **Buoyancy**: Archimedes' Principle simulation
*   **Laminar Pipe Flow**: Poiseuille's Law with parabolic velocity profiles
*   **Advanced CFD**: Navier-Stokes solver using Stable Fluids algorithm
*   **Turbulence Modeling**: k-ε turbulence model with Reynolds number calculation
*   **Flow Regime Detection**: Automatic laminar/transitional/turbulent classification

### 💨 Gas Dynamics
*   **Ideal Gas Law**: Interactive PV = nRT simulation
*   **Real Gas Equations**: Van der Waals equation for high-pressure accuracy
*   **Compressibility Factor**: Z-factor calculation and visualization
*   **Gas-Specific Constants**: Air, CO₂, Helium with accurate properties

### 🔥 Thermodynamics & Phase Changes
*   **Phase Transitions**: Ice → Water → Steam with latent heat
*   **Thermal Conductivity**: Heat equation solver (1D)
*   **Convection**: Natural and forced convection simulation

### 🌌 Astrophysics & Motion
*   **Orbital Mechanics**: N-body gravity simulation
*   **Projectile Motion**: Ballistic trajectories with air resistance

### 🎨 Visualization & Analysis
*   **Vector Fields**: Professional arrow glyphs with magnitude-based coloring
*   **Real-time Graphs**: Stress-strain curves, displacement history, flow analysis
*   **VTK Export**: Export to Paraview for advanced visualization (.vts, .vtu)
*   **Image Export**: Save plots and visualizations as PNG
*   **3D Rendering**: OpenGL with camera controls (pan, zoom, orbit, reset)

### 🛠️ User Interface
*   **Material Database**: 10 preset materials (Steel, Aluminum, Water, Air, etc.)
*   **Custom Material Editor**: Create materials without coding
*   **Preset Templates**: Metal, Fluid, Gas quick-start
*   **Scientific Mode**: Detailed equations and parameters
*   **Project Save/Load**: JSON-based session persistence

## 🛠️ Technology Stack
*   **Language**: C++17
*   **Rendering**: OpenGL 3.0
*   **UI**: [Dear ImGui](https://github.com/ocornut/imgui)
*   **Plotting**: [ImPlot](https://github.com/epezent/implot)
*   **Math**: GLM
*   **Windowing**: GLFW, GLEW
*   **Export**: VTK XML format, stb_image_write

## 🚀 Building and Running

### Prerequisites
*   C++ Compiler (GCC/Clang/MSVC with C++17 support)
*   CMake (3.10+)
*   OpenGL Libraries
*   GLFW3
*   GLEW

### Linux (Ubuntu/Debian)
```bash
sudo apt-get install build-essential cmake libglfw3-dev libglew-dev
```

### macOS (Homebrew)
```bash
brew install cmake glfw glew glm
```

### Windows (vcpkg)
1.  Install [vcpkg](https://github.com/microsoft/vcpkg)
2.  Install dependencies:
    ```powershell
    vcpkg install glfw3 glew glm opengl
    ```
3.  When running CMake, specify the vcpkg toolchain file:
    ```powershell
    cmake -B build -S . -DCMAKE_TOOLCHAIN_FILE=[path to vcpkg]/scripts/buildsystems/vcpkg.cmake
    ```

### Build Steps
```bash
# Clone the repository
git clone https://github.com/raj-saurav-sc/MatterLab.git
cd MatterLab

# Initialize submodules (Critical for third-party dependencies)
git submodule update --init --recursive

# Create build directory
cmake -B build -S .

# Build the project
cmake --build build

# Run the simulator
./build/UniversalSimulator
```

## 🎮 Controls

### Camera
*   **Orbit**: Right Mouse Button Drag
*   **Pan**: Left/Middle Mouse Button Drag or Arrow Keys
*   **Zoom**: Mouse Scroll Wheel or W/S keys
*   **Reset**: R key or "Reset View" button

### Simulation
*   **Select Scenario**: 12 scenarios including FEM, CFD, Gas Dynamics
*   **Material Selection**: Choose from 10 presets or create custom materials
*   **Parameters**: Real-time sliders for all physics properties
*   **Export**: VTK button for Paraview, Save Image for plots

### Advanced Features
*   **Turbulence Toggle**: Enable k-ε model for CFD
*   **Real Gas Mode**: Switch between ideal and Van der Waals
*   **Vector Visualization**: Adjust scale and density
*   **Material Editor**: Create custom materials with all properties

## 🏗️ Architecture
Modular design with clean separation of concerns:
*   **Core** (`src/core`): Application lifecycle, windowing, camera
*   **Physics** (`src/physics`): ISolver interface, Material system
*   **Solvers** (`src/physics/solvers`): 12 independent simulation modules
*   **Visualization** (`src/visualization`): Vector fields, rendering utilities
*   **IO** (`src/io`): VTK export, image capture

## 📊 Simulation Scenarios

1. **Solid Mechanics** - Hooke's Law deformation
2. **Fluid Mechanics** - Buoyancy and Archimedes
3. **Gas Dynamics** - Ideal/Real gas with Van der Waals
4. **Phase Change** - Melting and boiling
5. **Thermal Conductivity** - Heat diffusion
6. **Projectile Motion** - Ballistic trajectories
7. **Orbital Mechanics** - N-body gravity
8. **Pipe Flow** - Poiseuille flow
9. **Fluid Convection** - Natural convection
10. **Coupled Heat Transfer** - Thermal-fluid interaction
11. **Advanced FEM** - Truss solver with plasticity
12. **Advanced CFD** - Navier-Stokes with turbulence

## 🗺️ Roadmap

### Phase 1 (Completed) ✅
*   Core Physics Engines (Solids, Fluids, Gases)
*   Basic Visualization
*   Interactive UI
*   12 Simulation Scenarios

### Phase 2 (16% Complete) 🚧

#### Completed ✅
*   **Architecture Refactor**: Modular ISolver interface
*   **FEM Solver**: Truss elements with plasticity
*   **CFD Solver**: Navier-Stokes (Stable Fluids)
*   **Turbulence Models**: k-ε with Reynolds number
*   **Real Gas Equations**: Van der Waals, compressibility
*   **Vector Field Visualization**: Arrow glyphs, magnitude coloring
*   **VTK Export**: Paraview integration (.vts, .vtu)
*   **Custom Material Editor**: UI-based material creation
*   **Enhanced Controls**: Camera and plot interaction
*   **Data Persistence**: JSON save/load

#### In Progress 🚧
*   Thermal-Structural Coupling
*   Shock Wave Simulation
*   Mesh Visualization
*   Contour Plots
*   Parametric Studies

### Phase 3 (Planned) 📋
*   GPU Acceleration (CUDA/OpenCL)
*   Parallel Computing (OpenMP)
*   Advanced Mesh Generation
*   Topology Optimization
*   AI Material Prediction

## 📄 License
This project is open source.

## 🙏 Acknowledgments
*   Dear ImGui for the UI framework
*   ImPlot for scientific plotting
*   Jos Stam for the Stable Fluids algorithm
*   VTK community for visualization standards
