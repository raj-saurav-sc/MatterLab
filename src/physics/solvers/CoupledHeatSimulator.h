#pragma once

#include "../ISolver.h"
#include <imgui.h>
#include "implot.h"
#include <GL/glew.h>
#include <vector>
#include <cmath>

class CoupledHeatSimulator : public ISolver {
private:
    float solidTemp = 300.0f;
    float fluidTemp = 350.0f;
    float h_coeff = 10.0f; // Convection coefficient
    float surfaceArea = 1.0f;
    float solidMass = 1.0f;
    float solidCp = 500.0f;
    
    std::vector<glm::vec2> solidTempHistory;
    std::vector<glm::vec2> fluidTempHistory;
    float time = 0.0f;
    
public:
    std::string getName() const override { return "Coupled Heat Transfer"; }

    void initialize() override {
        resetSimulation();
    }
    
    void resetSimulation() {
        solidTemp = 300.0f;
        fluidTemp = 350.0f;
        time = 0.0f;
        solidTempHistory.clear();
        fluidTempHistory.clear();
    }
    
    void update(float deltaTime) override {
        // Newton's Law of Cooling: Q = h * A * (T_fluid - T_solid)
        float Q = h_coeff * surfaceArea * (fluidTemp - solidTemp);
        
        // Energy balance: Q = m * Cp * dT/dt
        float dT_solid = (Q * deltaTime) / (solidMass * solidCp);
        
        solidTemp += dT_solid;
        // Assume fluid temp is constant (infinite reservoir) or update it too
        
        time += deltaTime;
        solidTempHistory.push_back(glm::vec2(time, solidTemp - 273.15f));
        fluidTempHistory.push_back(glm::vec2(time, fluidTemp - 273.15f));
        
        if (solidTempHistory.size() > 500) solidTempHistory.erase(solidTempHistory.begin());
        if (fluidTempHistory.size() > 500) fluidTempHistory.erase(fluidTempHistory.begin());
    }
    
    void renderUI() override {
        ImGui::Text("Conjugate Heat Transfer");
        ImGui::Separator();
        
        ImGui::SliderFloat("Fluid Temp (K)", &fluidTemp, 200.0f, 500.0f);
        ImGui::SliderFloat("Convection Coeff (h)", &h_coeff, 1.0f, 100.0f);
        ImGui::Text("Solid Temp: %.2f K", solidTemp);
        
        if (ImGui::Button("Reset")) resetSimulation();
    }
    
    void renderGraph() override {
        const float plot_height = (ImGui::GetContentRegionAvail().y - ImGui::GetStyle().ItemSpacing.y) * 0.5f;
        const ImVec2 plot_size(-1, plot_height > 1.0f ? plot_height : 1.0f);
        
        if (ImPlot::BeginPlot("Temperature History", plot_size)) {
            if (!solidTempHistory.empty()) {
                std::vector<float> tVec, sVec, fVec;
                for (size_t i = 0; i < solidTempHistory.size(); ++i) {
                    tVec.push_back(solidTempHistory[i].x);
                    sVec.push_back(solidTempHistory[i].y);
                    fVec.push_back(fluidTempHistory[i].y);
                }
                ImPlot::PlotLine("Solid Temp", tVec.data(), sVec.data(), tVec.size());
                ImPlot::PlotLine("Fluid Temp", tVec.data(), fVec.data(), tVec.size());
            }
            ImPlot::EndPlot();
        }
        
        ImVec2 plotMin = ImGui::GetItemRectMin();
        ImVec2 plotMax = ImGui::GetItemRectMax();
        
        if (ImGui::Button("Export Heat Data")) {
            // Export logic
        }
        ImGui::SameLine();
        if (ImGui::Button("Save Image")) {
            requestPlotCapture("coupled_heat_plot.png", plotMin, plotMax);
        }
    }
    
    void render3D() override {
        // Draw Solid
        float tNorm = (solidTemp - 300.0f) / 100.0f;
        glColor3f(tNorm, 0.0f, 1.0f - tNorm);
        glPushMatrix();
        glTranslatef(-0.5f, 0, 0);
        // Cube
        float s = 0.4f;
        glBegin(GL_QUADS);
        glVertex3f(-s,-s,s); glVertex3f(s,-s,s); glVertex3f(s,s,s); glVertex3f(-s,s,s);
        // ...
        glEnd();
        glPopMatrix();
        
        // Draw Fluid Flow arrows
        glColor3f(0.0f, 0.5f, 1.0f);
        glBegin(GL_LINES);
        for (float y = -0.5f; y <= 0.5f; y += 0.2f) {
            glVertex3f(0.0f, y, 0.0f);
            glVertex3f(1.0f, y, 0.0f);
        }
        glEnd();
    }
};
