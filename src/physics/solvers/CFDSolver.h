#pragma once

#include "../ISolver.h"
#include "../../core/RandomGenerator.h"
#include "../../visualization/VectorField.h"
#include "../../io/VTKExporter.h"
#include "../Material.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <imgui.h>
#include "implot.h"
#include <GL/glew.h>

#define IX(x, y) ((x) + (N + 2) * (y))

class CFDSolver : public ISolver {
private:
    static const int N = 64; // Grid size (N x N)
    static const int SIZE = (N + 2) * (N + 2);
    
    // Fluid Data
    std::vector<float> u, v, u_prev, v_prev;
    std::vector<float> dens, dens_prev;
    
    // Turbulence Data (k-ε model)
    std::vector<float> k, k_prev;           // Turbulent kinetic energy
    std::vector<float> epsilon, eps_prev;   // Dissipation rate
    std::vector<float> nu_t;                // Turbulent viscosity
    
    // Simulation Parameters
    float dt = 0.1f;
    float diff = 0.0001f; // Diffusion rate
    float visc = 0.0001f; // Viscosity
    float force = 5.0f;
    float source = 100.0f;
    
    // Turbulence Model Constants (standard k-ε)
    float C_mu = 0.09f;
    float C_1 = 1.44f;
    float C_2 = 1.92f;
    float sigma_k = 1.0f;
    float sigma_e = 1.3f;
    
    // Flow Parameters
    float reynoldsNumber = 0.0f;
    bool useTurbulence = false;
    
    // Visualization
    bool showVelocity = false;
    bool showGrid = false;
    bool showTurbulence = false;
    float vectorScale = 10.0f;
    int vectorDensity = 2;  // Show every Nth vector
    bool colorByMagnitude = true;
    
    // Interaction
    int lastMouseX = -1, lastMouseY = -1;

public:
    CFDSolver() {
        u.resize(SIZE, 0.0f); v.resize(SIZE, 0.0f);
        u_prev.resize(SIZE, 0.0f); v_prev.resize(SIZE, 0.0f);
        dens.resize(SIZE, 0.0f); dens_prev.resize(SIZE, 0.0f);
        
        // Initialize turbulence fields
        k.resize(SIZE, 0.001f);
        k_prev.resize(SIZE, 0.001f);
        epsilon.resize(SIZE, 0.001f);
        eps_prev.resize(SIZE, 0.001f);
        nu_t.resize(SIZE, 0.0f);
    }

    std::string getName() const override { return "Advanced CFD (Navier-Stokes)"; }

    void initialize() override {
        reset();
    }

    void reset() override {
        std::fill(u.begin(), u.end(), 0.0f);
        std::fill(v.begin(), v.end(), 0.0f);
        std::fill(u_prev.begin(), u_prev.end(), 0.0f);
        std::fill(v_prev.begin(), v_prev.end(), 0.0f);
        std::fill(dens.begin(), dens.end(), 0.0f);
        std::fill(dens_prev.begin(), dens_prev.end(), 0.0f);
        std::fill(k.begin(), k.end(), 0.001f);
        std::fill(k_prev.begin(), k_prev.end(), 0.001f);
        std::fill(epsilon.begin(), epsilon.end(), 0.001f);
        std::fill(eps_prev.begin(), eps_prev.end(), 0.001f);
        std::fill(nu_t.begin(), nu_t.end(), 0.0f);
    }

    void update(float deltaTime) override {
        // Handle Input (Add Source)
        handleInput();
        
        // Calculate Reynolds number
        reynoldsNumber = calculateReynoldsNumber();
        
        // Solve turbulence (if enabled)
        if (useTurbulence) {
            solve_k_epsilon(dt);
        }
        
        // Velocity Step
        vel_step(N, u.data(), v.data(), u_prev.data(), v_prev.data(), visc, dt);
        
        // Density Step
        dens_step(N, dens.data(), dens_prev.data(), u.data(), v.data(), diff, dt);
    }

    void render3D() override {
        // Render Density Field
        glDisable(GL_LIGHTING);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        
        float cellWidth = 10.0f / N; // Assuming world size is roughly 10x10
        float startX = -5.0f;
        float startY = -5.0f;
        
        glBegin(GL_QUADS);
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                float d = dens[IX(i, j)];
                if (d > 0.01f) {
                    glColor4f(1.0f, 1.0f, 1.0f, std::min(d, 1.0f)); // White smoke
                    
                    float x = startX + (i - 1) * cellWidth;
                    float y = startY + (j - 1) * cellWidth;
                    
                    glVertex3f(x, y, 0.0f);
                    glVertex3f(x + cellWidth, y, 0.0f);
                    glVertex3f(x + cellWidth, y + cellWidth, 0.0f);
                    glVertex3f(x, y + cellWidth, 0.0f);
                }
            }
        }
        glEnd();
        
        // Render Velocity Vectors with improved arrows
        if (showVelocity) {
            std::vector<glm::vec3> positions;
            std::vector<glm::vec3> vectors;
            
            for (int i = 1; i <= N; i += vectorDensity) {
                for (int j = 1; j <= N; j += vectorDensity) {
                    float x = startX + (i - 1) * cellWidth + cellWidth * 0.5f;
                    float y = startY + (j - 1) * cellWidth + cellWidth * 0.5f;
                    
                    float vx = u[IX(i, j)];
                    float vy = v[IX(i, j)];
                    
                    positions.push_back(glm::vec3(x, y, 0.1f));
                    vectors.push_back(glm::vec3(vx, vy, 0.0f));
                }
            }
            
            VectorField::renderVectorField(positions, vectors, vectorScale, colorByMagnitude, 1);
        }
        
        // Render Grid
        if (showGrid) {
            glColor4f(0.3f, 0.3f, 0.3f, 0.5f);
            glBegin(GL_LINES);
            for (int i = 0; i <= N; i++) {
                float x = startX + i * cellWidth;
                glVertex3f(x, startY, 0.0f);
                glVertex3f(x, startY + N * cellWidth, 0.0f);
                
                float y = startY + i * cellWidth;
                glVertex3f(startX, y, 0.0f);
                glVertex3f(startX + N * cellWidth, y, 0.0f);
            }
            glEnd();
        }
        
        glDisable(GL_BLEND);
    }

    void renderUI() override {
        ImGui::Text("Navier-Stokes Fluid Solver");
        ImGui::Separator();
        
        // Turbulence Model
        if (ImGui::Checkbox("Enable Turbulence (k-\u03b5)", &useTurbulence)) {
            if (!useTurbulence) {
                // Reset turbulence fields when disabled
                std::fill(k.begin(), k.end(), 0.001f);
                std::fill(epsilon.begin(), epsilon.end(), 0.001f);
                std::fill(nu_t.begin(), nu_t.end(), 0.0f);
            }
        }
        
        ImGui::Separator();
        
        ImGui::SliderFloat("Viscosity", &visc, 0.0f, 0.01f, "%.5f");
        ImGui::SliderFloat("Diffusion", &diff, 0.0f, 0.01f, "%.5f");
        ImGui::SliderFloat("Input Force", &force, 1.0f, 20.0f);
        ImGui::SliderFloat("Input Density", &source, 10.0f, 500.0f);
        
        ImGui::Separator();
        
        // Flow Information
        ImGui::Text("Flow Analysis:");
        ImGui::Text("Reynolds Number: %.1f", reynoldsNumber);
        ImGui::Text("Flow Regime: %s", getFlowRegime().c_str());
        
        if (useTurbulence) {
            // Calculate average turbulent viscosity
            float avg_nu_t = 0.0f;
            for (int i = 1; i <= N; i++) {
                for (int j = 1; j <= N; j++) {
                    avg_nu_t += nu_t[IX(i, j)];
                }
            }
            avg_nu_t /= (N * N);
            
            float ratio = avg_nu_t / (visc + 1e-10f);
            ImGui::Text("Turb. Viscosity Ratio: %.2f", ratio);
            
            if (ratio > 10.0f) {
                ImGui::TextColored(ImVec4(1.0f, 0.5f, 0.0f, 1.0f), "High turbulence!");
            }
        }
        
        ImGui::Separator();
        
        
        ImGui::Checkbox("Show Velocity", &showVelocity);
        if (showVelocity) {
            ImGui::SliderFloat("Vector Scale", &vectorScale, 1.0f, 50.0f);
            ImGui::SliderInt("Vector Density", &vectorDensity, 1, 8);
            ImGui::Checkbox("Color by Magnitude", &colorByMagnitude);
        }
        
        ImGui::Checkbox("Show Grid", &showGrid);
        if (useTurbulence) {
            ImGui::Checkbox("Show Turbulence", &showTurbulence);
        }
        
        ImGui::Separator();
        
        if (ImGui::Button("Reset")) {
            reset();
        }
        ImGui::SameLine();
        if (ImGui::Button("Export VTK")) {
            exportToVTK("cfd_output.vts");
        }
        
        ImGui::TextWrapped("Use Left Mouse to add density and velocity.");
    }
    
    void exportToVTK(const std::string& filename) {
        // Prepare velocity vector data
        std::vector<glm::vec3> velocityVectors;
        for (int i = 0; i < SIZE; i++) {
            velocityVectors.push_back(glm::vec3(u[i], v[i], 0.0f));
        }
        
        // Export structured grid with density and velocity
        bool success = VTKExporter::exportStructuredGrid(
            filename,
            N, N, 0,  // 2D grid (nz = 0)
            dens,     // Scalar: density
            "density",
            velocityVectors,  // Vector: velocity
            "velocity"
        );
        
        if (success) {
            std::cout << "VTK export successful: " << filename << std::endl;
        } else {
            std::cerr << "VTK export failed!" << std::endl;
        }
    }
    
    void renderGraph() override {
        // No graph for now, maybe total mass?
    }
    
    json saveState() const override {
        json state;
        state["visc"] = visc;
        state["diff"] = diff;
        return state;
    }
    
    void loadState(const json& state) override {
        if (state.contains("visc")) visc = state["visc"];
        if (state.contains("diff")) diff = state["diff"];
    }

private:
    // --- Turbulence Methods ---
    
    float calculateReynoldsNumber() {
        // Re = ρ * U * L / μ
        // U = characteristic velocity (max velocity magnitude)
        // L = characteristic length (grid cell size)
        float maxVel = 0.0f;
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                float vel = std::sqrt(u[IX(i, j)] * u[IX(i, j)] + v[IX(i, j)] * v[IX(i, j)]);
                maxVel = std::max(maxVel, vel);
            }
        }
        
        float rho = 1.0f;  // Fluid density (kg/m³)
        float L = 10.0f / N;  // Grid cell size (m)
        return (rho * maxVel * L) / (visc + 1e-10f);  // Avoid division by zero
    }
    
    std::string getFlowRegime() {
        if (reynoldsNumber < 2300) return "Laminar";
        else if (reynoldsNumber < 4000) return "Transitional";
        else return "Turbulent";
    }
    
    void solve_k_epsilon(float dt) {
        if (!useTurbulence) return;
        
        // Calculate production term (simplified)
        std::vector<float> P_k(SIZE, 0.0f);
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                // Velocity gradients (simplified)
                float du_dx = (u[IX(i+1, j)] - u[IX(i-1, j)]) / 2.0f;
                float dv_dy = (v[IX(i, j+1)] - v[IX(i, j-1)]) / 2.0f;
                float du_dy = (u[IX(i, j+1)] - u[IX(i, j-1)]) / 2.0f;
                float dv_dx = (v[IX(i+1, j)] - v[IX(i-1, j)]) / 2.0f;
                
                // Production: P_k = ν_t * (∂u_i/∂x_j + ∂u_j/∂x_i) * ∂u_i/∂x_j
                float S = std::sqrt(2.0f * (du_dx*du_dx + dv_dy*dv_dy + 0.5f*(du_dy + dv_dx)*(du_dy + dv_dx)));
                P_k[IX(i, j)] = nu_t[IX(i, j)] * S * S;
            }
        }
        
        // k equation: ∂k/∂t + u·∇k = ∇·((ν + ν_t/σ_k)∇k) + P_k - ε
        std::copy(k.begin(), k.end(), k_prev.begin());
        
        // Diffusion with turbulent viscosity
        float eff_diff_k = visc + 0.001f; // Simplified: average ν_t
        diffuse(N, 0, k.data(), k_prev.data(), eff_diff_k / sigma_k, dt);
        
        // Advection
        std::copy(k.begin(), k.end(), k_prev.begin());
        advect(N, 0, k.data(), k_prev.data(), u.data(), v.data(), dt);
        
        // Source terms
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                k[IX(i, j)] += dt * (P_k[IX(i, j)] - epsilon[IX(i, j)]);
                if (k[IX(i, j)] < 0.001f) k[IX(i, j)] = 0.001f; // Prevent negative k
            }
        }
        
        // ε equation: ∂ε/∂t + u·∇ε = ∇·((ν + ν_t/σ_ε)∇ε) + C_1*ε/k*P_k - C_2*ε²/k
        std::copy(epsilon.begin(), epsilon.end(), eps_prev.begin());
        
        // Diffusion
        float eff_diff_e = visc + 0.001f;
        diffuse(N, 0, epsilon.data(), eps_prev.data(), eff_diff_e / sigma_e, dt);
        
        // Advection
        std::copy(epsilon.begin(), epsilon.end(), eps_prev.begin());
        advect(N, 0, epsilon.data(), eps_prev.data(), u.data(), v.data(), dt);
        
        // Source terms
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                int idx = IX(i, j);
                float k_val = std::max(k[idx], 0.001f);
                float prod = C_1 * (epsilon[idx] / k_val) * P_k[idx];
                float diss = C_2 * (epsilon[idx] * epsilon[idx]) / k_val;
                epsilon[idx] += dt * (prod - diss);
                if (epsilon[idx] < 0.001f) epsilon[idx] = 0.001f;
                
                // Clamp k and epsilon to prevent negative values and division by zero
                k[idx] = std::max(k[idx], 1e-10f);
                epsilon[idx] = std::max(epsilon[idx], 1e-10f);
                
                // Update turbulent viscosity: nu_t = C_mu * k^2 / epsilon
                nu_t[idx] = C_mu * k[idx] * k[idx] / epsilon[idx];
                
                // Limit turbulent viscosity to reasonable values
                nu_t[idx] = std::min(nu_t[idx], 1000.0f * visc);
            }
        }
    }
    
    void handleInput() {
        ImGuiIO& io = ImGui::GetIO();
        if (!io.WantCaptureMouse && ImGui::IsMouseDown(ImGuiMouseButton_Left)) {
            // Map mouse pos to grid
            // This is tricky without proper raycasting, so we'll do a simplified version
            // assuming the camera is looking at Z=0
            
            // For now, let's just add to the center if mouse is pressed, 
            // or implement a proper raycast later.
            // Add density/velocity at center + some noise
            
            int i = N/2 + RandomGenerator::randInt(-2, 2);
            int j = N/2 + RandomGenerator::randInt(-2, 2);
            
            float fx = RandomGenerator::randFloat(-5.0f, 5.0f);
            float fy = RandomGenerator::randFloat(-5.0f, 5.0f);
            
            add_source(N, dens.data(), &source, dt); // Add density everywhere? No.
            
            // Proper interaction:
            dens[IX(i, j)] += source * dt;
            u[IX(i, j)] += fx * dt;
            v[IX(i, j)] += fy * dt;
        }
    }

    // --- Stable Fluids Core ---
    
    void add_source(int N, float* x, float* s, float dt) {
        for (int i = 0; i < SIZE; i++) x[i] += dt * s[i];
    }

    void set_bnd(int N, int b, float* x) {
        for (int i = 1; i <= N; i++) {
            x[IX(0, i)] = b == 1 ? -x[IX(1, i)] : x[IX(1, i)];
            x[IX(N + 1, i)] = b == 1 ? -x[IX(N, i)] : x[IX(N, i)];
            x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
            x[IX(i, N + 1)] = b == 2 ? -x[IX(i, N)] : x[IX(i, N)];
        }
        x[IX(0, 0)] = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
        x[IX(0, N + 1)] = 0.5f * (x[IX(1, N + 1)] + x[IX(0, N)]);
        x[IX(N + 1, 0)] = 0.5f * (x[IX(N, 0)] + x[IX(N + 1, 1)]);
        x[IX(N + 1, N + 1)] = 0.5f * (x[IX(N, N + 1)] + x[IX(N + 1, N)]);
    }

    void lin_solve(int N, int b, float* x, float* x0, float a, float c) {
        for (int k = 0; k < 20; k++) {
            for (int i = 1; i <= N; i++) {
                for (int j = 1; j <= N; j++) {
                    x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i - 1, j)] + x[IX(i + 1, j)] +
                                                       x[IX(i, j - 1)] + x[IX(i, j + 1)])) / c;
                }
            }
            set_bnd(N, b, x);
        }
    }

    void diffuse(int N, int b, float* x, float* x0, float diff, float dt) {
        float a = dt * diff * N * N;
        lin_solve(N, b, x, x0, a, 1 + 4 * a);
    }

    void advect(int N, int b, float* d, float* d0, float* u, float* v, float dt) {
        float dt0 = dt * N;
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                float x = i - dt0 * u[IX(i, j)];
                float y = j - dt0 * v[IX(i, j)];
                
                if (x < 0.5f) x = 0.5f; if (x > N + 0.5f) x = N + 0.5f;
                int i0 = (int)x; int i1 = i0 + 1;
                
                if (y < 0.5f) y = 0.5f; if (y > N + 0.5f) y = N + 0.5f;
                int j0 = (int)y; int j1 = j0 + 1;
                
                float s1 = x - i0; float s0 = 1 - s1;
                float t1 = y - j0; float t0 = 1 - t1;
                
                d[IX(i, j)] = s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) +
                              s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
            }
        }
        set_bnd(N, b, d);
    }

    void project(int N, float* u, float* v, float* p, float* div) {
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                div[IX(i, j)] = -0.5f * (u[IX(i + 1, j)] - u[IX(i - 1, j)] +
                                         v[IX(i, j + 1)] - v[IX(i, j - 1)]) / N;
                p[IX(i, j)] = 0;
            }
        }
        set_bnd(N, 0, div); set_bnd(N, 0, p);
        
        lin_solve(N, 0, p, div, 1, 4);
        
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                u[IX(i, j)] -= 0.5f * N * (p[IX(i + 1, j)] - p[IX(i - 1, j)]);
                v[IX(i, j)] -= 0.5f * N * (p[IX(i, j + 1)] - p[IX(i, j - 1)]);
            }
        }
        set_bnd(N, 1, u); set_bnd(N, 2, v);
    }

    void dens_step(int N, float* x, float* x0, float* u, float* v, float diff, float dt) {
        // Add source handled externally for now
        std::swap(x0, x); diffuse(N, 0, x, x0, diff, dt);
        std::swap(x0, x); advect(N, 0, x, x0, u, v, dt);
    }

    void vel_step(int N, float* u, float* v, float* u0, float* v0, float visc, float dt) {
        // Add source handled externally
        std::swap(u0, u); diffuse(N, 1, u, u0, visc, dt);
        std::swap(v0, v); diffuse(N, 2, v, v0, visc, dt);
        project(N, u, v, u0, v0);
        std::swap(u0, u); std::swap(v0, v);
        advect(N, 1, u, u0, u0, v0, dt); advect(N, 2, v, v0, u0, v0, dt);
        project(N, u, v, u0, v0);
    }
};
