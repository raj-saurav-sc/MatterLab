#pragma once

#include "../ISolver.h"
#include <imgui.h>
#include <GL/glew.h>
#include <cmath>

class FluidMechanicsSolver : public ISolver {
private:
    Material object;
    Material fluid;
    float objectVolume = 0.001f; // m³
    float objectRadius = 0.1f;   // m
    
    // Results
    double buoyantForce = 0.0;
    double weight = 0.0;
    bool floats = false;
    double terminalVelocity = 0.0;
    
    // Animation
    float objectPosition = 0.0f;
    float targetPosition = 0.0f;
    
public:
    std::string getName() const override { return "Fluid Mechanics"; }

    void setObject(const Material& obj) { object = obj; }
    void setFluid(const Material& fl) { fluid = fl; }
    void setObjectProperties(float vol, float rad) { objectVolume = vol; objectRadius = rad; }
    
    // ISolver interface
    void setMaterial(const Material& mat) override {
        // Heuristic: if viscosity > 0, it's a fluid, else object
        if (mat.viscosity > 0) {
            setFluid(mat);
        } else {
            setObject(mat);
        }
    }
    
    void calculate() {
        weight = object.density * objectVolume * 9.81;
        buoyantForce = fluid.density * objectVolume * 9.81;
        floats = object.density < fluid.density;
        
        // Calculate terminal velocity (simplified)
        double netForce = buoyantForce - weight;
        terminalVelocity = std::sqrt(std::abs(netForce) / (0.5 * fluid.density * 3.14159 * objectRadius * objectRadius));
        
        // Set target position for animation
        targetPosition = floats ? 0.5f : -0.5f;
    }
    
    void update(float deltaTime) override {
        // Animate object movement
        float speed = 2.0f * deltaTime;
        if (std::abs(objectPosition - targetPosition) > 0.01f) {
            if (objectPosition < targetPosition) {
                objectPosition += speed;
            } else {
                objectPosition -= speed;
            }
        }
    }
    
    void renderUI() override {
        ImGui::Text("Fluid Mechanics Simulation");
        ImGui::Separator();
        
        ImGui::SliderFloat("Object Volume (m³)", &objectVolume, 0.0001f, 0.01f, "%.4f");
        ImGui::SliderFloat("Object Radius (m)", &objectRadius, 0.01f, 0.5f);
        
        if (ImGui::Button("Calculate")) {
            calculate();
        }
        
        ImGui::Separator();
        ImGui::Text("Results:");
        ImGui::Text("Weight: %.2f N", weight);
        ImGui::Text("Buoyant Force: %.2f N", buoyantForce);
        ImGui::Text("Net Force: %.2f N", buoyantForce - weight);
        ImGui::Text("Object %s", floats ? "FLOATS" : "SINKS");
        ImGui::Text("Terminal Velocity: %.2f m/s", terminalVelocity);
    }
    
    void render3D() override {
        // Render fluid container
        glColor3f(fluid.color.x * 0.5f, fluid.color.y * 0.5f, fluid.color.z * 0.5f);
        renderContainer();
        
        // Render object
        glPushMatrix();
        glTranslatef(0.0f, objectPosition, 0.0f);
        glColor3f(object.color.x, object.color.y, object.color.z);
        renderSphere(objectRadius);
        glPopMatrix();
    }
    
private:
    void renderContainer() {
        // Render transparent fluid container
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glColor4f(fluid.color.x, fluid.color.y, fluid.color.z, 0.3f);
        
        glBegin(GL_QUADS);
        // Render container walls (simplified box)
        float s = 0.5f;
        // Back
        glVertex3f(-s, -s, -s); glVertex3f(-s, s, -s); glVertex3f(s, s, -s); glVertex3f(s, -s, -s);
        // Left
        glVertex3f(-s, -s, -s); glVertex3f(-s, -s, s); glVertex3f(-s, s, s); glVertex3f(-s, s, -s);
        // Right
        glVertex3f(s, -s, -s); glVertex3f(s, s, -s); glVertex3f(s, s, s); glVertex3f(s, -s, s);
        // Bottom
        glVertex3f(-s, -s, -s); glVertex3f(s, -s, -s); glVertex3f(s, -s, s); glVertex3f(-s, -s, s);
        glEnd();
        
        glDisable(GL_BLEND);
    }
    
    void renderSphere(float radius) {
        // Simple sphere rendering using triangulation
        const int stacks = 20;
        const int slices = 20;
        
        for (int i = 0; i < stacks; ++i) {
            float phi1 = M_PI * i / stacks;
            float phi2 = M_PI * (i + 1) / stacks;
            
            glBegin(GL_TRIANGLE_STRIP);
            for (int j = 0; j <= slices; ++j) {
                float theta = 2 * M_PI * j / slices;
                
                float x1 = radius * sin(phi1) * cos(theta);
                float y1 = radius * cos(phi1);
                float z1 = radius * sin(phi1) * sin(theta);
                
                float x2 = radius * sin(phi2) * cos(theta);
                float y2 = radius * cos(phi2);
                float z2 = radius * sin(phi2) * sin(theta);
                
                glVertex3f(x1, y1, z1);
                glVertex3f(x2, y2, z2);
            }
            glEnd();
        }
    }
};
