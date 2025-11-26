#pragma once

#include "../ISolver.h"
#include <imgui.h>
#include "implot.h"
#include <GL/glew.h>
#include <vector>
#include <cmath>

class PhaseChangeSimulator : public ISolver {
private:
    Material material;
    float temperature = 273.15f; // K
    float mass = 1.0f;           // kg
    float heatInput = 0.0f;      // J/s
    
    // State
    std::string currentPhase = "Solid";
    float energy = 0.0f;
    
    // Graph data
    std::vector<glm::vec2> tempTimeData;
    float time = 0.0f;
    
public:
    std::string getName() const override { return "Phase Change"; }

    json saveState() const override {
        json state;
        state["temperature"] = temperature;
        state["heatInput"] = heatInput;
        state["materialName"] = material.name;
        return state;
    }

    void loadState(const json& state) override {
        if (state.contains("temperature")) temperature = state["temperature"];
        if (state.contains("heatInput")) heatInput = state["heatInput"];
        calculatePhase();
    }

    void setMaterial(const Material& mat) override { material = mat; }
    void setTemperature(float t) { temperature = t; calculatePhase(); }
    
    void update(float deltaTime) override {
        if (heatInput != 0.0f) {
            // Q = mcΔT
            // ΔT = Q / (mc)
            float deltaT = (heatInput * deltaTime) / (mass * material.specificHeatCapacity);
            
            // Handle phase change plateaus (simplified)
            bool isPhaseChange = false;
            if (abs(temperature - material.meltingPoint) < 1.0f && heatInput > 0 && currentPhase == "Solid") {
                // Melting...
                isPhaseChange = true;
            } else if (abs(temperature - material.boilingPoint) < 1.0f && heatInput > 0 && currentPhase == "Liquid") {
                // Boiling...
                isPhaseChange = true;
            }
            
            if (!isPhaseChange) {
                temperature += deltaT;
            }
            
            time += deltaTime;
            calculatePhase();
            
            tempTimeData.push_back(glm::vec2(time, temperature - 273.15f)); // Time vs Celsius
            if (tempTimeData.size() > 500) tempTimeData.erase(tempTimeData.begin());
        }
    }
    
    void calculatePhase() {
        if (temperature < material.meltingPoint) currentPhase = "Solid";
        else if (temperature < material.boilingPoint) currentPhase = "Liquid";
        else currentPhase = "Gas";
    }
    
    void renderUI() override {
        ImGui::Text("Phase Change Simulation");
        ImGui::Separator();
        
        ImGui::SliderFloat("Heat Input (W)", &heatInput, -1000.0f, 1000.0f);
        ImGui::Text("Temperature: %.2f K (%.2f C)", temperature, temperature - 273.15f);
        ImGui::Text("Phase: %s", currentPhase.c_str());
        
        ImGui::Separator();
        ImGui::Text("Material Properties:");
        ImGui::Text("Melting Point: %.2f K", material.meltingPoint);
        ImGui::Text("Boiling Point: %.2f K", material.boilingPoint);
    }
    
    void renderGraph() override {
        const float plot_height = (ImGui::GetContentRegionAvail().y - ImGui::GetStyle().ItemSpacing.y) * 0.5f;
        const ImVec2 plot_size(-1, plot_height > 1.0f ? plot_height : 1.0f);
        
        if (ImPlot::BeginPlot("Heating Curve", plot_size)) {
            if (!tempTimeData.empty()) {
                std::vector<float> timeVec, tempVec;
                for (const auto& point : tempTimeData) {
                    timeVec.push_back(point.x);
                    tempVec.push_back(point.y);
                }
                ImPlot::PlotLine("Temp (C) vs Time", timeVec.data(), tempVec.data(), timeVec.size());
                
                // Phase lines
                float meltC = material.meltingPoint - 273.15f;
                float boilC = material.boilingPoint - 273.15f;
                ImPlot::PlotInfLines("Melting Point", &meltC, 1, ImPlotInfLinesFlags_Horizontal);
                ImPlot::PlotInfLines("Boiling Point", &boilC, 1, ImPlotInfLinesFlags_Horizontal);
            }
            ImPlot::EndPlot();
        }
        
        ImVec2 plotMin = ImGui::GetItemRectMin();
        ImVec2 plotMax = ImGui::GetItemRectMax();

        if (ImGui::Button("Export Temp Data")) {
            // Export logic
        }
        ImGui::SameLine();
        if (ImGui::Button("Save Image")) {
            requestPlotCapture("phase_change_plot.png", plotMin, plotMax);
        }
    }
    
    void render3D() override {
        // Visualize phase state
        if (currentPhase == "Solid") {
            // Draw cube
            glColor3f(material.color.x, material.color.y, material.color.z);
            float s = 0.5f;
            glBegin(GL_QUADS);
            // ... cube vertices ...
            glVertex3f(-s, -s, s); glVertex3f(s, -s, s); glVertex3f(s, s, s); glVertex3f(-s, s, s);
            // (omitting other faces for brevity, assume standard cube)
            glEnd();
        } else if (currentPhase == "Liquid") {
            // Draw "puddle" or container
            glColor3f(material.color.x, material.color.y, material.color.z);
            glBegin(GL_QUADS);
            glVertex3f(-0.5f, -0.5f, 0.0f);
            glVertex3f(0.5f, -0.5f, 0.0f);
            glVertex3f(0.5f, 0.5f, 0.0f);
            glVertex3f(-0.5f, 0.5f, 0.0f);
            glEnd();
        } else {
            // Draw particles
            glPointSize(3.0f);
            glBegin(GL_POINTS);
            glColor3f(material.color.x, material.color.y, material.color.z);
            for (int i = 0; i < 50; ++i) {
                glVertex3f(
                    (rand()%100)/50.0f - 1.0f,
                    (rand()%100)/50.0f - 1.0f,
                    (rand()%100)/50.0f - 1.0f
                );
            }
            glEnd();
            glPointSize(1.0f);
        }
    }
};
