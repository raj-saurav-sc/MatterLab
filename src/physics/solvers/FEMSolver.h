#pragma once

#include "../ISolver.h"
#include "../../io/VTKExporter.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <glm/glm.hpp>
#include <imgui.h>

struct Node {
    glm::vec3 position;
    glm::vec3 originalPosition;
    glm::vec3 velocity;
    glm::vec3 force;
    float mass = 1.0f;
    bool fixed = false;
};

struct Element {
    int nodeIndices[2]; // Truss element connects 2 nodes
    float restLength;
    float stiffness; // k = EA/L
    float currentStress;
    
    // Plasticity
    float plasticStrain = 0.0f;
    float yieldStrength = 250.0f; // Reduced for visibility (simulating a weaker material or high load)
    bool isYielding = false;
};

class FEMSolver : public ISolver {
private:
    std::vector<Node> nodes;
    std::vector<Element> elements;
    Material material;
    
    // Simulation parameters
    float damping = 0.95f;
    float gravity = -9.81f;
    int subSteps = 10;
    
    // UI State
    bool showNodes = true;
    bool showStress = true;
    float loadMultiplier = 1.0f;
    
    // Graph Data
    std::vector<glm::vec2> displacementHistory;
    float time = 0.0f;

public:
    std::string getName() const override { return "Advanced FEM (Truss)"; }

    void initialize() override {
        resetSimulation();
    }

    void resetSimulation() {
        nodes.clear();
        elements.clear();
        displacementHistory.clear();
        time = 0.0f;
        
        // Create a simple bridge/truss structure
        int spans = 4;
        float spanLength = 1.0f;
        float height = 1.0f;
        
        // Bottom chord
        for (int i = 0; i <= spans; ++i) {
            Node n;
            n.position = glm::vec3(i * spanLength - (spans * spanLength)/2.0f, 0, 0);
            n.originalPosition = n.position;
            n.velocity = glm::vec3(0);
            n.fixed = (i == 0 || i == spans); // Fix ends
            nodes.push_back(n);
        }
        
        // Top chord
        for (int i = 0; i < spans; ++i) {
            Node n;
            n.position = glm::vec3(i * spanLength + spanLength/2.0f - (spans * spanLength)/2.0f, height, 0);
            n.originalPosition = n.position;
            n.velocity = glm::vec3(0);
            nodes.push_back(n);
        }
        
        // Create elements (Simple Truss)
        auto addElement = [&](int idx1, int idx2) {
            Element e;
            e.nodeIndices[0] = idx1;
            e.nodeIndices[1] = idx2;
            e.restLength = glm::length(nodes[idx1].position - nodes[idx2].position);
            e.stiffness = 10000.0f; 
            e.yieldStrength = 500.0f; // Set yield strength
            elements.push_back(e);
        };
        
        // Connect bottom chord
        for (int i = 0; i < spans; ++i) addElement(i, i+1);
        
        // Connect top chord
        for (int i = 0; i < spans - 1; ++i) addElement(spans + 1 + i, spans + 1 + i + 1);
        
        // Connect verticals and diagonals
        for (int i = 0; i < spans; ++i) {
            addElement(i, spans + 1 + i); // Vertical left
            addElement(i + 1, spans + 1 + i); // Diagonal
            if (i < spans - 1) addElement(i + 1, spans + 1 + i + 1); // Vertical right (shared)
        }
        addElement(spans, spans * 2); // Last vertical
    }

    void update(float deltaTime) override {
        float dt = deltaTime / subSteps;
        time += deltaTime;
        
        for (int step = 0; step < subSteps; ++step) {
            // Clear forces
            for (auto& node : nodes) {
                node.force = glm::vec3(0, gravity * node.mass, 0);
                // Apply external load to center bottom node
                if (&node == &nodes[nodes.size()/2]) { // Rough approximation of center
                     node.force.y -= 50.0f * loadMultiplier;
                }
            }
            
            // Calculate elastic forces & Plasticity
            for (auto& element : elements) {
                Node& n1 = nodes[element.nodeIndices[0]];
                Node& n2 = nodes[element.nodeIndices[1]];
                
                glm::vec3 delta = n2.position - n1.position;
                float currentLength = glm::length(delta);
                if (currentLength == 0) continue;
                
                glm::vec3 dir = delta / currentLength;
                float displacement = currentLength - element.restLength;
                
                float forceMag = element.stiffness * displacement;
                
                // Plasticity Check
                element.isYielding = false;
                if (std::abs(forceMag) > element.yieldStrength) {
                    element.isYielding = true;
                    
                    // Calculate plastic deformation
                    // If force > yield, the element stretches permanently
                    // New rest length should be such that force == yieldStrength
                    
                    float sign = (forceMag > 0) ? 1.0f : -1.0f;
                    float maxElasticDisplacement = element.yieldStrength / element.stiffness;
                    float targetLength = currentLength - (sign * maxElasticDisplacement);
                    
                    // Update rest length (plastic flow)
                    // Rate dependent or instant? Let's do instant for stability in this simple model
                    element.restLength = targetLength;
                    
                    // Cap force at yield strength
                    forceMag = sign * element.yieldStrength;
                }
                
                glm::vec3 force = dir * forceMag;
                
                n1.force += force;
                n2.force -= force;
                
                element.currentStress = forceMag;
            }
            
            // Integrate
            for (auto& node : nodes) {
                if (node.fixed) continue;
                
                glm::vec3 acceleration = node.force / node.mass;
                node.velocity += acceleration * dt;
                node.velocity *= damping; // Simple damping
                node.position += node.velocity * dt;
            }
        }
        
        // Record data for graph (Max displacement)
        float maxDisp = 0.0f;
        for (const auto& node : nodes) {
            float disp = glm::length(node.position - node.originalPosition);
            if (disp > maxDisp) maxDisp = disp;
        }
        displacementHistory.push_back(glm::vec2(time, maxDisp));
        if (displacementHistory.size() > 500) displacementHistory.erase(displacementHistory.begin());
    }

    void render3D() override {
        // Render Elements
        glLineWidth(3.0f);
        glBegin(GL_LINES);
        for (const auto& element : elements) {
            // Color based on stress/yield
            if (element.isYielding) {
                glColor3f(1.0f, 1.0f, 0.0f); // Yellow for yielding
            } else {
                float stressRatio = std::abs(element.currentStress) / element.yieldStrength;
                if (element.currentStress > 0) glColor3f(stressRatio, 0.0f, 0.0f); // Red Tension
                else glColor3f(0.0f, 0.0f, stressRatio); // Blue Compression
            }
            
            glm::vec3 p1 = nodes[element.nodeIndices[0]].position;
            glm::vec3 p2 = nodes[element.nodeIndices[1]].position;
            
            glVertex3f(p1.x, p1.y, p1.z);
            glVertex3f(p2.x, p2.y, p2.z);
        }
        glEnd();
        glLineWidth(1.0f);
        
        // Render Nodes
        if (showNodes) {
            glPointSize(5.0f);
            glColor3f(1.0f, 1.0f, 1.0f);
            glBegin(GL_POINTS);
            for (const auto& node : nodes) {
                glVertex3f(node.position.x, node.position.y, node.position.z);
            }
            glEnd();
            glPointSize(1.0f);
        }
    }

    void renderUI() override {
        ImGui::Text("Finite Element Method (Truss)");
        ImGui::Separator();
        
        ImGui::SliderFloat("Load Multiplier", &loadMultiplier, 0.0f, 5.0f);
        ImGui::SliderFloat("Damping", &damping, 0.0f, 1.0f);
        
        if (ImGui::Button("Reset")) {
            resetSimulation();
        }
        ImGui::SameLine();
        if (ImGui::Button("Export VTK")) {
            exportToVTK("fem_output.vtu");
        }
        ImGui::Checkbox("Show Nodes", &showNodes);
    }
    
    void exportToVTK(const std::string& filename) {
        // Prepare node positions
        std::vector<glm::vec3> points;
        for (const auto& node : nodes) {
            points.push_back(node.position);
        }
        
        // Prepare element connectivity
        std::vector<std::vector<int>> cells;
        for (const auto& elem : elements) {
            cells.push_back({elem.nodeIndices[0], elem.nodeIndices[1]});
        }
        
        // Prepare displacement magnitudes as scalar data
        std::vector<float> displacements;
        for (const auto& node : nodes) {
            float disp = glm::length(node.velocity); // Using velocity as displacement indicator
            displacements.push_back(disp);
        }
        
        // Prepare displacement vectors
        std::vector<glm::vec3> dispVectors;
        for (const auto& node : nodes) {
            dispVectors.push_back(node.velocity);
        }
        
        // Export unstructured grid
        bool success = VTKExporter::exportUnstructuredGrid(
            filename,
            points,
            cells,
            displacements,
            "displacement_magnitude",
            dispVectors,
            "displacement"
        );
        
        if (success) {
            std::cout << "VTK export successful: " << filename << std::endl;
        } else {
            std::cerr << "VTK export failed!" << std::endl;
        }
    }    
    
    void renderGraph() override {
        if (ImPlot::BeginPlot("Max Displacement vs Time")) {
            if (!displacementHistory.empty()) {
                std::vector<float> t, d;
                for (const auto& p : displacementHistory) {
                    t.push_back(p.x);
                    d.push_back(p.y);
                }
                ImPlot::PlotLine("Displacement", t.data(), d.data(), t.size());
            }
            ImPlot::EndPlot();
        }
        
        ImVec2 plotMin = ImGui::GetItemRectMin();
        ImVec2 plotMax = ImGui::GetItemRectMax();
        
        if (ImGui::Button("Save Image")) {
            requestPlotCapture("fem_displacement_plot.png", plotMin, plotMax);
        }
    }
    
    json saveState() const override {
        json state;
        state["loadMultiplier"] = loadMultiplier;
        state["damping"] = damping;
        return state;
    }
    
    void loadState(const json& state) override {
        if (state.contains("loadMultiplier")) loadMultiplier = state["loadMultiplier"];
        if (state.contains("damping")) damping = state["damping"];
    }
};
