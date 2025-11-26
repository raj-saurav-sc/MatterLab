#pragma once

#include "../ISolver.h"
#include <imgui.h>
#include "implot.h"
#include <GL/glew.h>
#include <vector>
#include <cmath>

class ThermalConductivitySimulator : public ISolver {
private:
    struct Node {
        float temperature;
        float x;
    };
    
    std::vector<Node> nodes;
    int numNodes = 20;
    float rodLength = 1.0f;
    float leftTemp = 300.0f;
    float rightTemp = 300.0f;
    float thermalConductivity = 50.0f; // W/(m*K)
    
    std::vector<glm::vec2> tempProfileData;
    
public:
    std::string getName() const override { return "Thermal Conductivity"; }

    json saveState() const override {
        json state;
        state["leftTemp"] = leftTemp;
        state["rightTemp"] = rightTemp;
        state["thermalConductivity"] = thermalConductivity;
        return state;
    }

    void loadState(const json& state) override {
        if (state.contains("leftTemp")) leftTemp = state["leftTemp"];
        if (state.contains("rightTemp")) rightTemp = state["rightTemp"];
        if (state.contains("thermalConductivity")) thermalConductivity = state["thermalConductivity"];
        initialize(); // Reset with new params
    }

    void initialize() override {
        nodes.clear();
        for (int i = 0; i < numNodes; ++i) {
            nodes.push_back({300.0f, (float)i / (numNodes - 1) * rodLength});
        }
    }
    
    void update(float deltaTime) override {
        // 1D Heat Equation: dT/dt = alpha * d^2T/dx^2
        // alpha = k / (rho * cp) -> simplified here
        
        float dx = rodLength / (numNodes - 1);
        float alpha = thermalConductivity * 0.0001f; // Simplified diffusivity
        
        std::vector<float> newTemps(numNodes);
        
        // Boundary conditions
        nodes[0].temperature = leftTemp;
        nodes[numNodes-1].temperature = rightTemp;
        newTemps[0] = leftTemp;
        newTemps[numNodes-1] = rightTemp;
        
        // Explicit finite difference
        for (int i = 1; i < numNodes - 1; ++i) {
            float d2T = (nodes[i+1].temperature - 2*nodes[i].temperature + nodes[i-1].temperature) / (dx*dx);
            newTemps[i] = nodes[i].temperature + alpha * d2T * deltaTime * 1000.0f; // Speed up
        }
        
        for (int i = 0; i < numNodes; ++i) {
            nodes[i].temperature = newTemps[i];
        }
        
        // Update graph data
        tempProfileData.clear();
        for (const auto& node : nodes) {
            tempProfileData.push_back(glm::vec2(node.x, node.temperature - 273.15f));
        }
    }
    
    void renderUI() override {
        ImGui::Text("1D Heat Conduction");
        ImGui::Separator();
        
        ImGui::SliderFloat("Left Temp (K)", &leftTemp, 200.0f, 1000.0f);
        ImGui::SliderFloat("Right Temp (K)", &rightTemp, 200.0f, 1000.0f);
        ImGui::SliderFloat("Conductivity", &thermalConductivity, 1.0f, 400.0f);
        
        if (ImGui::Button("Reset")) {
            initialize();
        }
    }
    
    void renderGraph() override {
        const float plot_height = (ImGui::GetContentRegionAvail().y - ImGui::GetStyle().ItemSpacing.y) * 0.5f;
        const ImVec2 plot_size(-1, plot_height > 1.0f ? plot_height : 1.0f);
        
        if (ImPlot::BeginPlot("Temperature Profile", plot_size)) {
            if (!tempProfileData.empty()) {
                std::vector<float> xVec, tVec;
                for (const auto& point : tempProfileData) {
                    xVec.push_back(point.x);
                    tVec.push_back(point.y);
                }
                ImPlot::PlotLine("Temp (C) vs Position", xVec.data(), tVec.data(), xVec.size());
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
            requestPlotCapture("thermal_conductivity_plot.png", plotMin, plotMax);
        }
    }
    
    void render3D() override {
        // Render rod with temperature gradient
        glLineWidth(10.0f);
        glBegin(GL_LINE_STRIP);
        for (const auto& node : nodes) {
            float tNormalized = (node.temperature - 200.0f) / 800.0f;
            glColor3f(tNormalized, 0.0f, 1.0f - tNormalized);
            glVertex3f(node.x - rodLength/2.0f, 0.0f, 0.0f);
        }
        glEnd();
        glLineWidth(1.0f);
    }
};
