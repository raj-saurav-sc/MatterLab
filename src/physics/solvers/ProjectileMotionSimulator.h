#pragma once

#include "../ISolver.h"
#include <imgui.h>
#include "implot.h"
#include <GL/glew.h>
#include <vector>
#include <cmath>

class ProjectileMotionSimulator : public ISolver {
private:
    glm::vec3 position{0,0,0};
    glm::vec3 velocity{10,10,0};
    float gravity = 9.81f;
    float mass = 1.0f;
    bool isSimulating = false;
    
    std::vector<glm::vec3> trail;
    
public:
    std::string getName() const override { return "Projectile Motion"; }

    void update(float deltaTime) override {
        if (isSimulating && position.y >= 0) {
            velocity.y -= gravity * deltaTime;
            position += velocity * deltaTime;
            
            trail.push_back(position);
            
            if (position.y < 0) {
                position.y = 0;
                isSimulating = false;
            }
        }
    }
    
    void renderUI() override {
        ImGui::Text("Projectile Motion");
        ImGui::Separator();
        
        static float launchSpeed = 15.0f;
        static float launchAngle = 45.0f;
        
        ImGui::SliderFloat("Speed (m/s)", &launchSpeed, 1.0f, 100.0f);
        ImGui::SliderFloat("Angle (deg)", &launchAngle, 0.0f, 90.0f);
        ImGui::SliderFloat("Gravity (m/s²)", &gravity, 1.0f, 20.0f);
        
        if (ImGui::Button("Launch")) {
            position = glm::vec3(-4, 0, 0);
            float rad = launchAngle * M_PI / 180.0f;
            velocity = glm::vec3(launchSpeed * cos(rad), launchSpeed * sin(rad), 0);
            trail.clear();
            trail.push_back(position);
            isSimulating = true;
        }
        
        if (ImGui::Button("Reset")) {
            position = glm::vec3(-4, 0, 0);
            trail.clear();
            isSimulating = false;
        }
    }
    
    void renderGraph() override {
        const float plot_height = (ImGui::GetContentRegionAvail().y - ImGui::GetStyle().ItemSpacing.y) * 0.5f;
        const ImVec2 plot_size(-1, plot_height > 1.0f ? plot_height : 1.0f);
        
        if (ImPlot::BeginPlot("Trajectory", plot_size)) {
            if (!trail.empty()) {
                std::vector<float> xVec, yVec;
                float startX = -4.0f;
                for (const auto& p : trail) {
                    xVec.push_back(p.x - startX);
                    yVec.push_back(p.y);
                }
                ImPlot::PlotLine("Height vs Distance", xVec.data(), yVec.data(), xVec.size());
            }
            ImPlot::EndPlot();
        }
        
        ImVec2 plotMin = ImGui::GetItemRectMin();
        ImVec2 plotMax = ImGui::GetItemRectMax();
        
        if (ImGui::Button("Export Trajectory")) {
            // Export logic
        }
        ImGui::SameLine();
        if (ImGui::Button("Save Image")) {
            requestPlotCapture("projectile_plot.png", plotMin, plotMax);
        }
    }
    
    void render3D() override {
        // Render ground
        glColor3f(0.3f, 0.8f, 0.3f);
        glBegin(GL_QUADS);
        glVertex3f(-10, 0, -5);
        glVertex3f(10, 0, -5);
        glVertex3f(10, 0, 5);
        glVertex3f(-10, 0, 5);
        glEnd();
        
        // Render projectile
        glColor3f(1.0f, 0.0f, 0.0f);
        glPushMatrix();
        glTranslatef(position.x, position.y, position.z);
        // Sphere...
        glPointSize(10.0f);
        glBegin(GL_POINTS);
        glVertex3f(0,0,0);
        glEnd();
        glPopMatrix();
        
        // Render trail
        glColor3f(1.0f, 1.0f, 0.0f);
        glBegin(GL_LINE_STRIP);
        for (const auto& p : trail) {
            glVertex3f(p.x, p.y, p.z);
        }
        glEnd();
    }
};
