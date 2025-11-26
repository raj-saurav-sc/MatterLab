#include "../ISolver.h"
#include "../../core/RandomGenerator.h"
#include <imgui.h>
#include <GL/glew.h>
#include <vector>
#include <cmath>

class FluidConvectionSimulator : public ISolver {
private:
    struct Particle {
        glm::vec3 position;
        glm::vec3 velocity;
        float temperature; // 0.0 to 1.0
        float life;
    };
    
    std::vector<Particle> particles;
    float heatSourceTemp = 100.0f; // C
    float ambientTemp = 20.0f;     // C
    float buoyancyFactor = 1.0f;
    
public:
    std::string getName() const override { return "Fluid Convection"; }

    json saveState() const override {
        json state;
        state["buoyancyFactor"] = buoyancyFactor;
        return state;
    }

    void loadState(const json& state) override {
        if (state.contains("buoyancyFactor")) buoyancyFactor = state["buoyancyFactor"];
    }

    void initialize() override {
        resetSimulation();
    }
    
    void resetSimulation() {
        particles.clear();
        for (int i = 0; i < 500; ++i) {
            spawnParticle();
        }
    }
    
    void spawnParticle() {
        Particle p;
        p.position = glm::vec3(
            RandomGenerator::randFloat(-1.0f, 1.0f), // x: -1 to 1
            -1.0f,                                    // y: bottom
            RandomGenerator::randFloat(-1.0f, 1.0f)  // z: -1 to 1
        );
        p.velocity = glm::vec3(0,0,0);
        p.temperature = 0.0f; // Start cold
        p.life = 1.0f;
        particles.push_back(p);
    }
    
    void update(float deltaTime) override {
        for (auto& p : particles) {
            // Heat up near bottom
            if (p.position.y < -0.8f) {
                p.temperature += deltaTime * 2.0f;
            } else {
                // Cool down
                p.temperature -= deltaTime * 0.5f;
            }
            p.temperature = std::max(0.0f, std::min(1.0f, p.temperature));
            
            // Buoyancy force: F = alpha * deltaT * g
            float buoyancy = p.temperature * buoyancyFactor;
            p.velocity.y += buoyancy * deltaTime;
            p.velocity.y -= 0.5f * deltaTime; // Gravity/Drag
            
            // Turbulence/Diffusion
            p.velocity.x += RandomGenerator::randFloat(-0.1f, 0.1f) * deltaTime;
            p.velocity.z += RandomGenerator::randFloat(-0.1f, 0.1f) * deltaTime;
            
            p.position += p.velocity * deltaTime;
            
            // Reset if too high or dead
            if (p.position.y > 2.0f) {
                p.position.y = -1.0f;
                p.velocity = glm::vec3(0,0,0);
                p.temperature = 0.0f;
            }
            // Walls
            if (abs(p.position.x) > 1.0f) p.velocity.x *= -1;
            if (abs(p.position.z) > 1.0f) p.velocity.z *= -1;
        }
    }
    
    void renderUI() override {
        ImGui::Text("Convection Simulation");
        ImGui::Separator();
        ImGui::SliderFloat("Buoyancy Factor", &buoyancyFactor, 0.1f, 5.0f);
        if (ImGui::Button("Reset")) resetSimulation();
    }
    
    void render3D() override {
        // Draw heat source
        glColor3f(1.0f, 0.2f, 0.0f);
        glBegin(GL_QUADS);
        glVertex3f(-1, -1.05f, -1);
        glVertex3f(1, -1.05f, -1);
        glVertex3f(1, -1.05f, 1);
        glVertex3f(-1, -1.05f, 1);
        glEnd();
        
        // Draw particles as quads (billboards would be better)
        for (const auto& p : particles) {
            // Color map: Blue (cold) -> Red (hot)
            glColor3f(p.temperature, 0.0f, 1.0f - p.temperature);
            
            float s = 0.05f;
            glPushMatrix();
            glTranslatef(p.position.x, p.position.y, p.position.z);
            glBegin(GL_QUADS);
            glVertex3f(-s, -s, 0);
            glVertex3f(s, -s, 0);
            glVertex3f(s, s, 0);
            glVertex3f(-s, s, 0);
            glEnd();
            glPopMatrix();
        }
    }
};
