#pragma once

#include "../ISolver.h"
#include <imgui.h>
#include "implot.h"
#include <GL/glew.h>
#include <vector>
#include <cmath>

class OrbitalMechanicsSimulator : public ISolver {
public:
    struct Body {
        glm::vec3 position;
        glm::vec3 velocity;
        float mass;
        float radius;
        glm::vec3 color;
        std::vector<glm::vec3> trail;
    };
    
private:
    std::vector<Body> bodies;
    float G = 1.0f; // Gravitational constant (scaled)
    float timeScale = 1.0f;
    
public:
    std::string getName() const override { return "Orbital Mechanics"; }

    json saveState() const override {
        json state;
        state["timeScale"] = timeScale;
        state["G"] = G;
        state["bodyCount"] = bodies.size();
        // Could save individual bodies too, but for now just params
        return state;
    }

    void loadState(const json& state) override {
        if (state.contains("timeScale")) timeScale = state["timeScale"];
        if (state.contains("G")) G = state["G"];
        // Reset system to apply defaults if needed, or load bodies if we implemented that
        if (state.contains("bodyCount") && state["bodyCount"] == 3) {
             resetSolarSystem();
        }
    }

    void initialize() override {
        resetSolarSystem();
    }
    
    void resetSolarSystem() {
        bodies.clear();
        // Sun
        bodies.push_back({
            glm::vec3(0,0,0),
            glm::vec3(0,0,0),
            1000.0f,
            1.0f,
            glm::vec3(1.0f, 1.0f, 0.0f),
            {}
        });
        // Planet
        bodies.push_back({
            glm::vec3(10,0,0),
            glm::vec3(0,0,10), // v = sqrt(GM/r) = sqrt(1000/10) = 10
            10.0f,
            0.3f,
            glm::vec3(0.2f, 0.4f, 1.0f),
            {}
        });
        // Moon
        bodies.push_back({
            glm::vec3(11,0,0),
            glm::vec3(0,0,13),
            1.0f,
            0.1f,
            glm::vec3(0.6f, 0.6f, 0.6f),
            {}
        });
    }
    
    void addRandomPlanet() {
        float r = 5.0f + (rand() % 150) / 10.0f;
        float theta = (rand() % 360) * M_PI / 180.0f;
        float x = r * cos(theta);
        float z = r * sin(theta);
        
        // Protect division by zero
        float r_safe = std::max(r, 0.1f);
        float vMag = sqrt(G * bodies[0].mass / r_safe);
        float vx = -vMag * sin(theta);
        float vz = vMag * cos(theta);
        
        bodies.push_back({
            glm::vec3(x, 0, z),
            glm::vec3(vx, 0, vz),
            1.0f + (rand() % 50) / 10.0f,
            0.2f + (rand() % 30) / 100.0f,
            glm::vec3((rand()%10)/10.0f, (rand()%10)/10.0f, (rand()%10)/10.0f),
            {}
        });
    }
    
    void update(float deltaTime) override {
        float dt = deltaTime * timeScale;
        
        // Calculate forces
        for (size_t i = 0; i < bodies.size(); ++i) {
            glm::vec3 force(0,0,0);
            for (size_t j = 0; j < bodies.size(); ++j) {
                if (i == j) continue;
                
                glm::vec3 diff = bodies[j].position - bodies[i].position;
                float dist = glm::length(diff);
                if (dist < 0.1f) continue; // Collision/Singularity avoidance
                
                float fMag = G * bodies[i].mass * bodies[j].mass / (dist * dist);
                force += glm::normalize(diff) * fMag;
            }
            
            // F = ma -> a = F/m (protect division by zero)
            if (bodies[i].mass > 1e-6f) {
                glm::vec3 acceleration = force / bodies[i].mass;
                bodies[i].velocity += acceleration * dt;
            }
        }
        
        // Update positions
        for (auto& body : bodies) {
            body.position += body.velocity * dt;
            
            // Update trail
            if (body.trail.size() > 500) body.trail.erase(body.trail.begin());
            body.trail.push_back(body.position);
        }
    }
    
    void renderUI() override {
        ImGui::Text("Orbital Mechanics");
        ImGui::Separator();
        
        ImGui::SliderFloat("Time Scale", &timeScale, 0.1f, 5.0f);
        ImGui::SliderFloat("G Constant", &G, 0.1f, 10.0f);
        
        if (ImGui::Button("Reset System")) resetSolarSystem();
        if (ImGui::Button("Add Random Planet")) addRandomPlanet();
        
        ImGui::Text("Bodies: %lu", bodies.size());
    }
    
    void renderGraph() override {
        const float plot_height = (ImGui::GetContentRegionAvail().y - ImGui::GetStyle().ItemSpacing.y) * 0.5f;
        const ImVec2 plot_size(-1, plot_height > 1.0f ? plot_height : 1.0f);
        
        if (ImPlot::BeginPlot("Orbit X-Z", plot_size)) {
            for (size_t i = 0; i < bodies.size(); ++i) {
                if (bodies[i].trail.empty()) continue;
                
                std::vector<float> xVec, zVec;
                for (const auto& p : bodies[i].trail) {
                    xVec.push_back(p.x);
                    zVec.push_back(p.z);
                }
                ImPlot::PlotLine(("Body " + std::to_string(i)).c_str(), xVec.data(), zVec.data(), xVec.size());
            }
            ImPlot::EndPlot();
        }
        
        ImVec2 plotMin = ImGui::GetItemRectMin();
        ImVec2 plotMax = ImGui::GetItemRectMax();
        
        if (ImGui::Button("Export Orbits")) {
            // Export logic
        }
        ImGui::SameLine();
        if (ImGui::Button("Save Image")) {
            requestPlotCapture("orbital_plot.png", plotMin, plotMax);
        }
    }
    
    void render3D() override {
        for (const auto& body : bodies) {
            glPushMatrix();
            glTranslatef(body.position.x, body.position.y, body.position.z);
            glColor3f(body.color.x, body.color.y, body.color.z);
            
            // Render sphere (point for now)
            glPointSize(body.radius * 10.0f);
            glBegin(GL_POINTS);
            glVertex3f(0,0,0);
            glEnd();
            glPopMatrix();
            
            // Render trail
            glColor3f(body.color.x, body.color.y, body.color.z);
            glBegin(GL_LINE_STRIP);
            for (const auto& p : body.trail) {
                glVertex3f(p.x, p.y, p.z);
            }
            glEnd();
        }
        glPointSize(1.0f);
    }
};
