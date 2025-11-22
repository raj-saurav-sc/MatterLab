#pragma once

#include "../ISolver.h"
#include <imgui.h>
#include "implot.h"
#include <GL/glew.h>
#include <vector>
#include <cmath>

class PipeFlowSimulator : public ISolver {
private:
    Material fluid;
    float pipeRadius = 0.1f; // m
    float pipeLength = 2.0f; // m
    float pressureDiff = 100.0f; // Pa
    
    // Visualization
    std::vector<glm::vec3> particles;
    
public:
    std::string getName() const override { return "Pipe Flow"; }

    void setFluid(const Material& f) { fluid = f; }
    
    void initialize() override {
        resetParticles();
    }
    
    void resetParticles() {
        particles.clear();
        for (int i = 0; i < 200; ++i) {
            float r = (float)sqrt((rand() % 100) / 100.0f) * pipeRadius;
            float theta = (rand() % 360) * M_PI / 180.0f;
            float z = (rand() % 100) / 100.0f * pipeLength - pipeLength/2;
            particles.push_back(glm::vec3(r * cos(theta), r * sin(theta), z));
        }
    }
    
    void update(float deltaTime) override {
        // Poiseuille Flow: v(r) = (ΔP / 4ηL) * (R² - r²)
        float viscosity = std::max(fluid.viscosity, 0.0001); // Avoid div by zero
        float factor = pressureDiff / (4.0f * viscosity * pipeLength);
        
        for (auto& p : particles) {
            float r2 = p.x*p.x + p.y*p.y;
            float v = factor * (pipeRadius*pipeRadius - r2);
            
            p.z += v * deltaTime;
            
            // Wrap around
            if (p.z > pipeLength/2) p.z -= pipeLength;
        }
    }
    
    void renderUI() override {
        ImGui::Text("Pipe Flow (Poiseuille)");
        ImGui::Separator();
        
        ImGui::SliderFloat("Pipe Radius (m)", &pipeRadius, 0.01f, 0.5f);
        ImGui::SliderFloat("Pressure Diff (Pa)", &pressureDiff, 0.0f, 1000.0f);
        ImGui::Text("Fluid Viscosity: %.4f Pa·s", fluid.viscosity);
        
        if (ImGui::Button("Reset Particles")) resetParticles();
    }
    
    void renderGraph() override {
        const float plot_height = (ImGui::GetContentRegionAvail().y - ImGui::GetStyle().ItemSpacing.y) * 0.5f;
        const ImVec2 plot_size(-1, plot_height > 1.0f ? plot_height : 1.0f);
        
        if (ImPlot::BeginPlot("Velocity Profile", plot_size)) {
            std::vector<float> rVec, vVec;
            float viscosity = std::max(fluid.viscosity, 0.0001);
            float factor = pressureDiff / (4.0f * viscosity * pipeLength);
            
            for (float r = -pipeRadius; r <= pipeRadius; r += pipeRadius/20.0f) {
                rVec.push_back(r);
                vVec.push_back(factor * (pipeRadius*pipeRadius - r*r));
            }
            
            ImPlot::PlotLine("Velocity vs Radius", rVec.data(), vVec.data(), rVec.size());
            ImPlot::EndPlot();
        }
        
        ImVec2 plotMin = ImGui::GetItemRectMin();
        ImVec2 plotMax = ImGui::GetItemRectMax();
        
        if (ImGui::Button("Export Velocity Data")) {
            // Export logic
        }
        ImGui::SameLine();
        if (ImGui::Button("Save Image")) {
            requestPlotCapture("pipe_flow_plot.png", plotMin, plotMax);
        }
    }
    
    void render3D() override {
        // Render pipe walls (wireframe cylinder)
        glColor3f(0.5f, 0.5f, 0.5f);
        const int segments = 20;
        for (int i = 0; i < segments; ++i) {
            float theta = 2 * M_PI * i / segments;
            float x = pipeRadius * cos(theta);
            float y = pipeRadius * sin(theta);
            
            glBegin(GL_LINES);
            glVertex3f(x, y, -pipeLength/2);
            glVertex3f(x, y, pipeLength/2);
            glEnd();
        }
        
        // Render particles
        glPointSize(3.0f);
        glBegin(GL_POINTS);
        for (const auto& p : particles) {
            // Color based on speed
            float r2 = p.x*p.x + p.y*p.y;
            float maxV = (pressureDiff / (4.0f * std::max(fluid.viscosity, 0.0001) * pipeLength)) * (pipeRadius*pipeRadius);
            float v = (pressureDiff / (4.0f * std::max(fluid.viscosity, 0.0001) * pipeLength)) * (pipeRadius*pipeRadius - r2);
            float intensity = v / maxV;
            
            glColor3f(intensity, 0.0f, 1.0f - intensity);
            glVertex3f(p.x, p.y, p.z);
        }
        glEnd();
        glPointSize(1.0f);
    }
};
