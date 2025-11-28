#pragma once

#include "../ISolver.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <imgui.h>
#include "implot.h"

class ShockWaveSolver : public ISolver {
private:
    // Grid parameters
    static const int NUM_CELLS = 200;
    float dx;
    float length = 1.0f; // Tube length in meters
    
    // Physics constants
    float gamma = 1.4f; // Heat capacity ratio (adiabatic index) for air
    float cfl = 0.5f;   // CFL condition number
    
    // State variables (Conservative variables: rho, rho*u, E)
    // U = [rho, momentum, energy]
    struct State {
        float rho;      // Density
        float momentum; // rho * u
        float energy;   // Total energy per unit volume
    };
    
    std::vector<State> U;
    std::vector<State> U_new;
    
    // Visualization data
    std::vector<float> xCoords;
    std::vector<float> densityData;
    std::vector<float> pressureData;
    std::vector<float> velocityData;
    std::vector<float> energyData;
    
    // Simulation control
    bool isRunning = false;
    float simulationTime = 0.0f;
    
public:
    ShockWaveSolver() {
        dx = length / NUM_CELLS;
        U.resize(NUM_CELLS);
        U_new.resize(NUM_CELLS);
        
        // Initialize x coordinates for plotting
        for (int i = 0; i < NUM_CELLS; ++i) {
            xCoords.push_back((i + 0.5f) * dx);
        }
        
        resetSimulation();
    }
    
    std::string getName() const override { return "Shock Wave (Gas Dynamics)"; }
    
    void resetSimulation() {
        simulationTime = 0.0f;
        isRunning = false;
        
        // Sod Shock Tube Initial Conditions
        // Left side (High pressure/density)
        // Right side (Low pressure/density)
        for (int i = 0; i < NUM_CELLS; ++i) {
            if (i < NUM_CELLS / 2) {
                // Left state
                float rho = 1.0f;
                float p = 1.0f;
                float u = 0.0f;
                
                U[i].rho = rho;
                U[i].momentum = rho * u;
                U[i].energy = p / (gamma - 1.0f) + 0.5f * rho * u * u;
            } else {
                // Right state
                float rho = 0.125f;
                float p = 0.1f;
                float u = 0.0f;
                
                U[i].rho = rho;
                U[i].momentum = rho * u;
                U[i].energy = p / (gamma - 1.0f) + 0.5f * rho * u * u;
            }
        }
        
        updateVisualizationData();
    }
    
    // Calculate flux vector F(U)
    State calculateFlux(const State& u) {
        float rho = u.rho;
        float velocity = u.momentum / rho;
        float pressure = (gamma - 1.0f) * (u.energy - 0.5f * rho * velocity * velocity);
        
        State flux;
        flux.rho = u.momentum;
        flux.momentum = u.momentum * velocity + pressure;
        flux.energy = (u.energy + pressure) * velocity;
        
        return flux;
    }
    
    void update(float deltaTime) override {
        if (!isRunning) return;
        
        // Calculate max wave speed for time step (CFL condition)
        float maxSpeed = 0.0f;
        for (const auto& u : U) {
            float rho = u.rho;
            float velocity = u.momentum / rho;
            float pressure = (gamma - 1.0f) * (u.energy - 0.5f * rho * velocity * velocity);
            pressure = std::max(pressure, 1e-6f); // Prevent negative pressure
            float soundSpeed = std::sqrt(gamma * pressure / rho);
            maxSpeed = std::max(maxSpeed, std::abs(velocity) + soundSpeed);
        }
        
        // Adaptive time step
        float dt = cfl * dx / maxSpeed;
        dt = std::min(dt, deltaTime); // Don't exceed frame time
        
        // Lax-Friedrichs Method (First order, stable but diffusive)
        // U_i^{n+1} = 0.5*(U_{i+1}^n + U_{i-1}^n) - 0.5*(dt/dx)*(F_{i+1}^n - F_{i-1}^n)
        
        for (int i = 1; i < NUM_CELLS - 1; ++i) {
            State U_left = U[i-1];
            State U_right = U[i+1];
            
            State F_left = calculateFlux(U_left);
            State F_right = calculateFlux(U_right);
            
            U_new[i].rho = 0.5f * (U_left.rho + U_right.rho) - 0.5f * (dt / dx) * (F_right.rho - F_left.rho);
            U_new[i].momentum = 0.5f * (U_left.momentum + U_right.momentum) - 0.5f * (dt / dx) * (F_right.momentum - F_left.momentum);
            U_new[i].energy = 0.5f * (U_left.energy + U_right.energy) - 0.5f * (dt / dx) * (F_right.energy - F_left.energy);
        }
        
        // Boundary conditions (Transmissive/Neumann zero gradient)
        U_new[0] = U_new[1];
        U_new[NUM_CELLS-1] = U_new[NUM_CELLS-2];
        
        // Update state
        U = U_new;
        simulationTime += dt;
        
        updateVisualizationData();
    }
    
    void updateVisualizationData() {
        densityData.clear();
        pressureData.clear();
        velocityData.clear();
        energyData.clear();
        
        for (const auto& u : U) {
            float rho = u.rho;
            float velocity = (rho > 1e-6f) ? (u.momentum / rho) : 0.0f;
            float pressure = (gamma - 1.0f) * (u.energy - 0.5f * rho * velocity * velocity);
            float internalEnergy = pressure / ((gamma - 1.0f) * rho);
            
            densityData.push_back(rho);
            pressureData.push_back(pressure);
            velocityData.push_back(velocity);
            energyData.push_back(internalEnergy);
        }
    }
    
    void renderUI() override {
        ImGui::Text("Shock Wave Simulation (1D Euler Equations)");
        ImGui::Text("Scenario: Sod Shock Tube");
        ImGui::Text("Time: %.4f s", simulationTime);
        
        if (ImGui::Button(isRunning ? "Pause" : "Run")) {
            isRunning = !isRunning;
        }
        ImGui::SameLine();
        if (ImGui::Button("Reset")) {
            resetSimulation();
        }
        
        ImGui::Separator();
        ImGui::Text("Visualization:");
        
        // Render 1D strip visualization
        ImDrawList* draw_list = ImGui::GetWindowDrawList();
        ImVec2 p = ImGui::GetCursorScreenPos();
        float width = ImGui::GetContentRegionAvail().x;
        float height = 30.0f;
        
        for (int i = 0; i < NUM_CELLS; ++i) {
            float x0 = p.x + (i * width / NUM_CELLS);
            float x1 = p.x + ((i + 1) * width / NUM_CELLS);
            float y0 = p.y;
            float y1 = p.y + height;
            
            // Color based on density (Blue -> Red)
            float rho = densityData[i];
            float t = std::clamp((rho - 0.125f) / (1.0f - 0.125f), 0.0f, 1.0f);
            ImU32 color = ImGui::ColorConvertFloat4ToU32(ImVec4(t, 0.0f, 1.0f - t, 1.0f));
            
            draw_list->AddRectFilled(ImVec2(x0, y0), ImVec2(x1, y1), color);
        }
        ImGui::Dummy(ImVec2(width, height));
        ImGui::Text("Density Map (Blue=Low, Red=High)");
        
        ImGui::Separator();
        
        // Render Graphs
        if (ImPlot::BeginPlot("Flow Properties", ImVec2(-1, 300))) {
            ImPlot::SetupAxes("Position (m)", "Value");
            ImPlot::SetupAxisLimits(ImAxis_X1, 0, length, ImPlotCond_Always);
            ImPlot::SetupAxisLimits(ImAxis_Y1, -0.2, 1.2, ImPlotCond_Once);
            
            ImPlot::PlotLine("Density", xCoords.data(), densityData.data(), NUM_CELLS);
            ImPlot::PlotLine("Pressure", xCoords.data(), pressureData.data(), NUM_CELLS);
            ImPlot::PlotLine("Velocity", xCoords.data(), velocityData.data(), NUM_CELLS);
            
            ImPlot::EndPlot();
        }
    }
    
    // Required but unused for 1D solver
    void setMaterial(const Material& mat) override {}
    json saveState() const override { return json(); }
    void loadState(const json& state) override {}
    void initialize() override { resetSimulation(); }
    void render3D() override {}
};
