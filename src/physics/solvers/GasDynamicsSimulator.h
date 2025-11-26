#pragma once

#include "../ISolver.h"
#include <imgui.h>
#include "implot.h"
#include <GL/glew.h>
#include <vector>
#include <cmath>

class GasDynamicsSimulator : public ISolver {
private:
    Material gas;
    float pressure = 101325.0f; // Pa
    float volume = 1.0f;        // m³
    float temperature = 300.0f; // K
    float moles = 1.0f;         // mol
    
    // Real gas mode
    bool useRealGas = false;
    float compressibilityFactor = 1.0f;
    float idealPressure = 101325.0f;
    
    // Visualization
    std::vector<glm::vec3> particles;
    std::vector<glm::vec3> velocities;
    
    // Graph data
    std::vector<glm::vec2> pvData;
    std::vector<glm::vec2> tvData;
    std::vector<glm::vec2> zData;  // Compressibility factor vs pressure
    
public:
    std::string getName() const override { return "Gas Dynamics"; }

    json saveState() const override {
        json state;
        state["temperature"] = temperature;
        state["volume"] = volume;
        state["moles"] = moles;
        state["gasName"] = gas.name;
        return state;
    }

    void loadState(const json& state) override {
        if (state.contains("temperature")) temperature = state["temperature"];
        if (state.contains("volume")) volume = state["volume"];
        if (state.contains("moles")) moles = state["moles"];
        calculate();
    }

    void setGas(const Material& g) { gas = g; }
    void setTemperature(float t) { temperature = t; calculate(); }
    
    void initialize() override {
        // Initialize particles
        for (int i = 0; i < 100; ++i) {
            particles.push_back(glm::vec3(
                (rand() % 100) / 50.0f - 1.0f,
                (rand() % 100) / 50.0f - 1.0f,
                (rand() % 100) / 50.0f - 1.0f
            ));
            velocities.push_back(glm::vec3(
                (rand() % 100) / 100.0f - 0.5f,
                (rand() % 100) / 100.0f - 0.5f,
                (rand() % 100) / 100.0f - 0.5f
            ));
        }
        calculate();
    }
    
    void calculate() {
        const float R = 8.314f;
        
        // Calculate ideal gas pressure first
        idealPressure = (moles * R * temperature) / volume;
        
        if (useRealGas && (gas.vanDerWaalsA > 0 || gas.vanDerWaalsB > 0)) {
            // Van der Waals equation: (P + a*n²/V²)(V - n*b) = nRT
            // Solving for P: P = nRT/(V - n*b) - a*n²/V²
            float a = gas.vanDerWaalsA;
            float b = gas.vanDerWaalsB;
            
            float term1 = (moles * R * temperature) / (volume - moles * b);
            float term2 = (a * moles * moles) / (volume * volume);
            pressure = term1 - term2;
            
            // Compressibility factor Z = PV/nRT
            compressibilityFactor = (pressure * volume) / (moles * R * temperature);
        } else {
            // Ideal gas
            pressure = idealPressure;
            compressibilityFactor = 1.0f;
        }
        
        // Update graph data
        pvData.push_back(glm::vec2(volume, pressure / 1000.0f)); // m³ vs kPa
        if (pvData.size() > 100) pvData.erase(pvData.begin());
        
        tvData.push_back(glm::vec2(temperature, volume));
        if (tvData.size() > 100) tvData.erase(tvData.begin());
        
        zData.push_back(glm::vec2(pressure / 1000.0f, compressibilityFactor));
        if (zData.size() > 100) zData.erase(zData.begin());
    }
    
    void update(float deltaTime) override {
        // Update particle positions based on temperature (kinetic energy)
        float speedMultiplier = sqrt(temperature / 300.0f);
        
        for (size_t i = 0; i < particles.size(); ++i) {
            particles[i] += velocities[i] * speedMultiplier * deltaTime;
            
            // Bounce off walls (container size based on volume)
            float containerSize = pow(volume, 1.0f/3.0f);
            float limit = containerSize / 2.0f;
            
            if (abs(particles[i].x) > limit) velocities[i].x *= -1;
            if (abs(particles[i].y) > limit) velocities[i].y *= -1;
            if (abs(particles[i].z) > limit) velocities[i].z *= -1;
            
            // Clamp positions
            particles[i].x = std::max(-limit, std::min(limit, particles[i].x));
            particles[i].y = std::max(-limit, std::min(limit, particles[i].y));
            particles[i].z = std::max(-limit, std::min(limit, particles[i].z));
        }
    }
    
    void renderUI() override {
        ImGui::Text("Gas Dynamics Simulation");
        ImGui::Separator();
        
        if (ImGui::Checkbox("Use Real Gas (Van der Waals)", &useRealGas)) {
            calculate();
        }
        
        ImGui::Separator();
        
        if (ImGui::SliderFloat("Temperature (K)", &temperature, 0.0f, 1000.0f)) {
            calculate();
        }
        if (ImGui::SliderFloat("Volume (m³)", &volume, 0.1f, 5.0f)) {
            calculate();
        }
        if (ImGui::SliderFloat("Moles", &moles, 0.1f, 10.0f)) {
            calculate();
        }
        
        ImGui::Separator();
        ImGui::Text("Results:");
        ImGui::Text("Pressure: %.2f kPa", pressure / 1000.0f);
        
        if (useRealGas) {
            ImGui::Text("Ideal Pressure: %.2f kPa", idealPressure / 1000.0f);
            ImGui::Text("Compressibility (Z): %.4f", compressibilityFactor);
            ImGui::Text("Deviation: %.2f%%", (compressibilityFactor - 1.0f) * 100.0f);
            
            if (compressibilityFactor < 0.95f || compressibilityFactor > 1.05f) {
                ImGui::TextColored(ImVec4(1.0f, 0.5f, 0.0f, 1.0f), "Significant real gas effects!");
            }
        }
        
        ImGui::Separator();
        ImGui::Text("Van der Waals Constants:");
        ImGui::Text("a = %.4f Pa·m⁶/mol²", gas.vanDerWaalsA);
        ImGui::Text("b = %.2e m³/mol", gas.vanDerWaalsB);
    }
    
    void renderGraph() override {
        const float plot_height = (ImGui::GetContentRegionAvail().y - ImGui::GetStyle().ItemSpacing.y) * 0.5f;
        const ImVec2 plot_size(-1, plot_height > 1.0f ? plot_height : 1.0f);
        
        if (ImPlot::BeginPlot("Pressure vs Volume", plot_size)) {
            if (!pvData.empty()) {
                std::vector<float> volVec, pressVec;
                for (const auto& point : pvData) {
                    volVec.push_back(point.x);
                    pressVec.push_back(point.y);
                }
                ImPlot::PlotLine("P vs V", volVec.data(), pressVec.data(), volVec.size());
            }
            ImPlot::EndPlot();
        }
        
        ImVec2 pvMin = ImGui::GetItemRectMin();
        ImVec2 pvMax = ImGui::GetItemRectMax();
        
        if (ImGui::Button("Save Image")) {
            requestPlotCapture("phase_change_plot.png", pvMin, pvMax); // Assuming plotMin/plotMax refers to the current plot's bounds
        }
        ImGui::SameLine();
        if (ImGui::Button("Save PV Plot")) {
            requestPlotCapture("gas_pv_plot.png", pvMin, pvMax);
        }
        
        if (ImPlot::BeginPlot("Temperature vs Volume", plot_size)) {
            if (!tvData.empty()) {
                std::vector<float> tempVec, volumeVec;
                for (const auto& point : tvData) {
                    tempVec.push_back(point.x);
                    volumeVec.push_back(point.y);
                }
                ImPlot::PlotLine("T vs V", tempVec.data(), volumeVec.data(), tempVec.size());
            }
            ImPlot::EndPlot();
        }
        ImVec2 tvMin = ImGui::GetItemRectMin();
        ImVec2 tvMax = ImGui::GetItemRectMax();
        
        if (ImGui::Button("Save TV Plot")) {
            requestPlotCapture("gas_tv_plot.png", tvMin, tvMax);
        }
        
        // Compressibility Factor graph (only if using real gas)
        if (useRealGas) {
            if (ImPlot::BeginPlot("Compressibility Factor (Z) vs Pressure", plot_size)) {
                if (!zData.empty()) {
                    std::vector<float> pressVec, zVec;
                    for (const auto& point : zData) {
                        pressVec.push_back(point.x);
                        zVec.push_back(point.y);
                    }
                    ImPlot::PlotLine("Z Factor", pressVec.data(), zVec.data(), pressVec.size());
                }
                ImPlot::EndPlot();
            }
            
            ImVec2 zMin = ImGui::GetItemRectMin();
            ImVec2 zMax = ImGui::GetItemRectMax();
            
            if (ImGui::Button("Save Z Plot")) {
                requestPlotCapture("gas_z_plot.png", zMin, zMax);
            }
        }
    }
    
    void render3D() override {
        // Render container
        float containerSize = pow(volume, 1.0f/3.0f);
        float s = containerSize / 2.0f;
        
        glColor3f(1.0f, 1.0f, 1.0f);
        glBegin(GL_LINE_LOOP);
        glVertex3f(-s, -s, -s); glVertex3f(-s, s, -s); glVertex3f(s, s, -s); glVertex3f(s, -s, -s);
        glEnd();
        glBegin(GL_LINE_LOOP);
        glVertex3f(-s, -s, s); glVertex3f(-s, s, s); glVertex3f(s, s, s); glVertex3f(s, -s, s);
        glEnd();
        glBegin(GL_LINES);
        glVertex3f(-s, -s, -s); glVertex3f(-s, -s, s);
        glVertex3f(-s, s, -s); glVertex3f(-s, s, s);
        glVertex3f(s, s, -s); glVertex3f(s, s, s);
        glVertex3f(s, -s, -s); glVertex3f(s, -s, s);
        glEnd();
        
        // Render particles
        glPointSize(5.0f);
        glBegin(GL_POINTS);
        glColor3f(gas.color.x, gas.color.y, gas.color.z);
        for (const auto& p : particles) {
            glVertex3f(p.x, p.y, p.z);
        }
        glEnd();
        glPointSize(1.0f);
    }
};
