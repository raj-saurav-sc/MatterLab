# MatterLab: Universal Material & Fluid Simulator

**MatterLab** is a multi-tiered physics simulator designed for both educational exploration and scientific visualization. It allows users to interact with various states of matter and physical phenomena in a real-time 3D environment.

## 🌟 Features

### 🏗️ Solid Mechanics
*   **Elastic Deformation**: Visualize stress and strain on materials using Hooke's Law.
*   **Plasticity & Fracture**: Simulate permanent deformation and material failure under high stress.

### 💧 Fluid Mechanics
*   **Buoyancy**: Drop objects into fluids and watch them float or sink based on density (Archimedes' Principle).
*   **Laminar Pipe Flow**: Visualize Poiseuille's Law with a parabolic velocity profile.

### 💨 Gas Dynamics
*   **Ideal Gas Law**: Interact with pressure, volume, and temperature ($PV = nRT$) in a confined container.

### 🔥 Thermodynamics & Phase Changes
*   **Phase Transitions**: Watch ice melt to water and boil to steam as you add heat.
*   **Thermal Conductivity**: Visualize heat propagation through a 1D rod (Heat Equation).

### 🌌 Astrophysics & Motion
*   **Orbital Mechanics**: Create solar systems with N-body gravity simulation.
*   **Projectile Motion**: Launch objects and observe the effects of gravity and air resistance.

## 🛠️ Technology Stack
*   **Language**: C++
*   **Rendering**: OpenGL (Legacy/Fixed Function for simplicity)
*   **UI**: [Dear ImGui](https://github.com/ocornut/imgui)
*   **Plotting**: [ImPlot](https://github.com/epezent/implot)
*   **Math**: GLM
*   **Windowing**: GLFW, GLEW

## 🚀 Building and Running

### Prerequisites
*   C++ Compiler (GCC/Clang/MSVC)
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
1.  Install [vcpkg](https://github.com/microsoft/vcpkg).
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

# Build the project (Assets are automatically copied to the build directory)
cmake --build build

# Run the simulator
./build/UniversalSimulator
```

## 🎮 Controls

### Camera
*   **Orbit**: Right Mouse Button Drag
*   **Pan**: Middle Mouse Button Drag
*   **Zoom**: Mouse Scroll Wheel
*   **UI Buttons**: Use the on-screen controls in the "Camera Controls" panel.

### Simulation
*   **Select Scenario**: Use the dropdown menu to switch between simulations (Solids, Fluids, Orbit, etc.).
*   **Parameters**: Adjust sliders to change mass, gravity, viscosity, temperature, etc.
*   **Modes**: Toggle "Scientific Mode" for detailed equations or "Educational Mode" for simple explanations.

## 🏗️ Architecture
The project has been refactored into a modular architecture to support advanced features:
*   **Core**: Handles application lifecycle, windowing, and input (`src/core`).
*   **Physics**: Defines the `ISolver` interface and material system (`src/physics`).
*   **Solvers**: Individual simulation modules (`src/physics/solvers`) that can be easily extended.

## 🗺️ Roadmap

### Phase 1 (Completed) ✅
*   Core Physics Engines (Solids, Fluids, Gases)
*   Basic Visualization
*   Interactive UI

### Phase 2 (In Progress) 🚧
*   **Architecture Refactor** (Completed) ✅
    *   Modular `ISolver` interface
    *   Core/Physics separation
*   Advanced Fluid Dynamics (Convection, Turbulence)
*   Coupled Heat Transfer (Solid-Fluid interaction)
*   Save/Load Simulation States
*   Data Export (CSV, Images)

## 📄 License
This project is open source.
