#pragma once

#include "../ISolver.h"
#include <imgui.h>
#include "implot.h"
#include <GL/glew.h>
#include <vector>
#include <cmath>
#include <algorithm>

// Re-implementation of AdvancedPhysicsEngine logic within the solver for now
// Ideally this should be in a shared physics utility
class SolidMechanicsSolver : public ISolver {
private:
    Material material;
    MaterialState materialState;
    
    float length = 1.0f;     // m
    float area = 0.01f;      // m²
    float appliedForce = 0.0f; // N
    
    // Results
    double stress = 0.0;
    double strain = 0.0;
    double elongation = 0.0;
    
    // Visualization data
    std::vector<glm::vec2> stressStrainData;
    std::vector<glm::vec2> loadingHistory;
    
    // Animation and visualization
    std::vector<glm::vec3> crackLines;
    
public:
    std::string getName() const override { return "Solid Mechanics"; }

    json saveState() const override {
        json state;
        state["appliedForce"] = appliedForce;
        state["length"] = length;
        state["area"] = area;
        state["materialName"] = material.name;
        return state;
    }

    void loadState(const json& state) override {
        if (state.contains("appliedForce")) appliedForce = state["appliedForce"];
        if (state.contains("length")) length = state["length"];
        if (state.contains("area")) area = state["area"];
        // Material loading will be handled by Application setting the material first, 
        // or we need access to DB here. For now, we assume material is set externally 
        // or we just save parameters.
        calculate();
    }

    void setMaterial(const Material& mat) override { 
        material = mat; 
        materialState = MaterialState(); // Reset state when changing materials
        materialState.currentYieldStrength = mat.yieldStrength;
        stressStrainData.clear();
        loadingHistory.clear();
    }
    
    void setDimensions(float l, float a) { length = l; area = a; }
    
    void setForce(float f) { 
        appliedForce = f; 
        calculate(); 
    }
    
    void calculate() {
        // Validate inputs to prevent NaN/inf
        if (area <= 0.0f) area = 0.001f;
        if (length <= 0.0f) length = 0.1f;
        if (material.youngsModulus <= 0.0) {
            std::cerr << "Warning: Invalid Young's Modulus, using default\n";
            return;
        }
        
        // Clamp force to reasonable values
        appliedForce = std::clamp(appliedForce, -1e8f, 1e8f);
        
        strain = calculateStrainFromForce(appliedForce, area, material.youngsModulus);
        
        // Check for numerical issues
        if (std::isnan(strain) || std::isinf(strain)) {
            std::cerr << "Warning: Invalid strain calculated, resetting\n";
            strain = 0.0;
            return;
        }
        
        // Logic from AdvancedPhysicsEngine::calculateStressWithPlasticity
        calculateStressWithPlasticity();
        
        // Validate stress
        if (std::isnan(stress) || std::isinf(stress)) {
            std::cerr << "Warning: Invalid stress calculated, resetting\n";
            stress = 0.0;
            return;
        }
        
        elongation = materialState.totalStrain * length;
        
        // Add to stress-strain curve (with validation)
        if (std::isfinite(materialState.totalStrain) && std::isfinite(stress)) {
            stressStrainData.push_back(glm::vec2(materialState.totalStrain * 100, stress / 1e6));
            loadingHistory.push_back(glm::vec2(ImGui::GetTime(), appliedForce));
        }
        
        if (stressStrainData.size() > 1000) {
            stressStrainData.erase(stressStrainData.begin());
        }
        if (loadingHistory.size() > 1000) {
            loadingHistory.erase(loadingHistory.begin());
        }
        
        // Update crack visualization
        updateCrackVisualization();
    }
    
    void renderUI() override {
        ImGui::Text("Advanced Solid Mechanics Simulation");
        ImGui::Separator();
        
        if (ImGui::SliderFloat("Applied Force (N)", &appliedForce, -100000, 100000)) {
            calculate();
        }
        if (ImGui::SliderFloat("Length (m)", &length, 0.1f, 5.0f)) {
            calculate();
        }
        if (ImGui::SliderFloat("Cross-sectional Area (m²)", &area, 0.001f, 0.1f)) {
            calculate();
        }
        
        // Cyclic loading controls
        ImGui::Separator();
        ImGui::Text("Advanced Loading:");
        
        static float cyclicAmplitude = 10000.0f;
        static float cyclicFrequency = 0.5f;
        static bool enableCyclicLoading = false;
        
        ImGui::Checkbox("Enable Cyclic Loading", &enableCyclicLoading);
        if (enableCyclicLoading) {
            ImGui::SliderFloat("Amplitude (N)", &cyclicAmplitude, 1000, 50000);
            ImGui::SliderFloat("Frequency (Hz)", &cyclicFrequency, 0.1f, 2.0f);
            
            float cyclicForce = cyclicAmplitude * sin(ImGui::GetTime() * 2 * M_PI * cyclicFrequency);
            setForce(cyclicForce);
        }
        
        if (ImGui::Button("Apply Tensile Test")) {
            performTensileTest();
        }
        ImGui::SameLine();
        if (ImGui::Button("Reset Material")) {
            resetMaterial();
        }
        
        ImGui::Separator();
        renderResults();
    }
    
    void renderResults() {
        ImGui::Text("Results:");
        ImGui::Text("Stress: %.2f MPa", stress / 1e6);
        ImGui::Text("Total Strain: %.4f%% (%.4f mm elongation)", 
                   materialState.totalStrain * 100, elongation * 1000);
        
        // Show plastic vs elastic components
        ImGui::Text("Elastic Strain: %.4f%%", materialState.elasticStrain * 100);
        ImGui::Text("Plastic Strain: %.4f%%", materialState.plasticStrain * 100);
        
        // Material state information
        ImGui::Separator();
        ImGui::Text("Material State:");
        
        if (materialState.plasticStrain > 0) {
            ImGui::TextColored(ImVec4(1.0f, 0.6f, 0.0f, 1.0f), "PLASTIC DEFORMATION");
            ImGui::Text("Work Hardening: Current Yield = %.1f MPa", 
                       materialState.currentYieldStrength / 1e6);
        } else {
            ImGui::TextColored(ImVec4(0.0f, 1.0f, 0.0f, 1.0f), "ELASTIC");
        }
        
        if (materialState.isFractured) {
            ImGui::TextColored(ImVec4(1.0f, 0.0f, 0.0f, 1.0f), "FRACTURED");
            ImGui::Text("Damage Parameter: %.2f", materialState.damageParameter);
        }
        
        // Safety factors
        double yieldSafetyFactor = (stress != 0) ? 
            std::abs(materialState.currentYieldStrength / stress) : INFINITY;
        double ultimateSafetyFactor = (stress != 0) ? 
            std::abs(material.ultimateTensileStrength / stress) : INFINITY;
            
        ImGui::Text("Safety Factor (Yield): %.2f", yieldSafetyFactor);
        ImGui::Text("Safety Factor (Ultimate): %.2f", ultimateSafetyFactor);
    }
    
    void renderGraph() override {
        const float plot_height = (ImGui::GetContentRegionAvail().y - ImGui::GetStyle().ItemSpacing.y) * 0.5f;
        const ImVec2 plot_size(-1, plot_height > 1.0f ? plot_height : 1.0f);
        
        // Enable pan/zoom with ImPlot flags
        ImPlotFlags flags = ImPlotFlags_NoTitle;
        
        if (ImPlot::BeginPlot("Advanced Stress-Strain Curve", plot_size, flags)) {
            ImPlot::SetupAxes("Strain (%)", "Stress (MPa)");
            ImPlot::SetupAxisLimits(ImAxis_X1, 0, 10, ImPlotCond_Once);
            ImPlot::SetupAxisLimits(ImAxis_Y1, 0, 500, ImPlotCond_Once);
            
            if (!stressStrainData.empty()) {
                std::vector<float> strainVec, stressVec;
                for (const auto& point : stressStrainData) {
                    strainVec.push_back(point.x);
                    stressVec.push_back(point.y);
                }
                ImPlot::PlotLine("Stress-Strain", strainVec.data(), stressVec.data(), strainVec.size());
                
                // Add material property lines
                // Yield strength line
                ImPlot::PlotInfLines("Yield Point", &material.yieldStrength, 1);
                ImPlot::PlotInfLines("Ultimate Strength", &material.ultimateTensileStrength, 1);
            }
            ImPlot::EndPlot();
        }
        
        ImVec2 plotMin = ImGui::GetItemRectMin();
        ImVec2 plotMax = ImGui::GetItemRectMax();
        
        if (ImGui::Button("Save Image##1")) {
            requestPlotCapture("solid_mechanics_plot.png", plotMin, plotMax);
        }
        ImGui::SameLine();
        if (ImGui::Button("Fit##1")) {
            ImPlot::SetNextAxesToFit();
        }
        ImGui::SameLine();
        if (ImGui::Button("Reset Axes##1")) {
            // This will be handled by ImPlot's context menu or we can force it
        }
        
        // Loading history graph
        if (ImPlot::BeginPlot("Loading History", plot_size, flags)) {
            ImPlot::SetupAxes("Time (s)", "Force (kN)");
            ImPlot::SetupAxisLimits(ImAxis_X1, 0, 10, ImPlotCond_Once);
            ImPlot::SetupAxisLimits(ImAxis_Y1, -100, 100, ImPlotCond_Once);
            
            if (!loadingHistory.empty()) {
                std::vector<float> timeVec, forceVec;
                for (const auto& point : loadingHistory) {
                    timeVec.push_back(point.x);
                    forceVec.push_back(point.y / 1000); // Convert to kN
                }
                ImPlot::PlotLine("Applied Force (kN)", timeVec.data(), forceVec.data(), timeVec.size());
            }
            ImPlot::EndPlot();
        }
        
        if (ImGui::Button("Save Image##2")) {
            ImVec2 plotMin2 = ImGui::GetItemRectMin();
            ImVec2 plotMax2 = ImGui::GetItemRectMax();
            requestPlotCapture("loading_history_plot.png", plotMin2, plotMax2);
        }
        ImGui::SameLine();
        if (ImGui::Button("Fit##2")) {
            ImPlot::SetNextAxesToFit();
        }
    }
    
    void render3D() override {
        glPushMatrix();
        
        // Calculate deformation factors
        double elasticDeformation = 1.0 + materialState.elasticStrain;
        
        // Choose color based on material state
        glm::vec3 currentColor = material.elasticColor;
        
        if (materialState.plasticStrain > 0) {
            // Interpolate between elastic and plastic colors
            float plasticRatio = std::min(materialState.plasticStrain / 0.1, 1.0); // Normalize
            currentColor = material.elasticColor * (1.0f - plasticRatio) + 
                          material.plasticColor * plasticRatio;
        }
        
        if (materialState.isFractured) {
            currentColor = material.fractureColor;
        }
        
        glColor3f(currentColor.x, currentColor.y, currentColor.z);
        
        // Render deformed specimen
        glScalef(elasticDeformation, 1.0f, 1.0f);
        renderDeformedBox(length, 0.1, 0.1);
        
        // Render cracks if fractured
        if (materialState.isFractured) {
            renderCracks();
        }
        
        glPopMatrix();
        
        // Render stress visualization
        renderStressField();
    }
    
private:
    double calculateStrainFromForce(double force, double area, double E) {
        // Prevent division by zero
        if (area <= 1e-10 || E <= 1e-10) {
            return 0.0;
        }
        double stress = force / area;
        double strain = stress / E;
        
        // Clamp to reasonable values
        return std::clamp(strain, -10.0, 10.0); // Max 1000% strain
    }
    
    void calculateStressWithPlasticity() {
        // Prevent division by zero
        if (material.youngsModulus <= 1e-10) {
            stress = 0.0;
            return;
        }
        
        if (strain <= materialState.elasticStrain + materialState.plasticStrain) {
            // Unloading
            materialState.elasticStrain = strain - materialState.plasticStrain;
            stress = materialState.elasticStrain * material.youngsModulus;
        } else {
            // Loading
            materialState.totalStrain = strain;
            double elasticStrainLimit = materialState.currentYieldStrength / material.youngsModulus;

            if (strain <= elasticStrainLimit) {
                // Elastic region
                materialState.elasticStrain = strain;
                materialState.plasticStrain = 0.0;
                stress = strain * material.youngsModulus;
            } else {
                // Plastic region
                materialState.elasticStrain = elasticStrainLimit;
                materialState.plasticStrain = strain - elasticStrainLimit;
                
                // Simplified linear work hardening
                double strainRange = 0.1; // 10% strain range
                double strengthDiff = material.ultimateTensileStrength - material.yieldStrength;
                double hardening = (strainRange > 1e-10) ? (strengthDiff / strainRange) : 0.0;
                stress = materialState.currentYieldStrength + hardening * materialState.plasticStrain;
                
                // Clamp stress to reasonable values
                stress = std::clamp(stress, -1e10, 1e10);
                materialState.currentYieldStrength = stress;

                if (stress > material.ultimateTensileStrength) {
                    materialState.isFractured = true;
                    double denominator = material.ultimateTensileStrength * 0.1;
                    materialState.damageParameter = (denominator > 1e-10) ? 
                        ((stress - material.ultimateTensileStrength) / denominator) : 1.0;
                    materialState.damageParameter = std::clamp(materialState.damageParameter, 0.0, 10.0);
                    if (materialState.crackTips.empty()) {
                        materialState.crackTips.push_back(glm::vec3(0,0,0));
                    }
                }
            }
        }
    }
    
    void performTensileTest() {
        // Automated tensile test - gradually increase force
        static bool testRunning = false;
        static float testForce = 0.0f;
        
        if (!testRunning) {
            testRunning = true;
            testForce = 0.0f;
        }
        
        testForce += 1000.0f; // Increase by 1kN per step
        setForce(testForce);
        
        // Stop test if material fractures
        if (materialState.isFractured) {
            testRunning = false;
            testForce = 0.0f;
        }
    }
    
    void resetMaterial() {
        materialState = MaterialState();
        materialState.currentYieldStrength = material.yieldStrength;
        stressStrainData.clear();
        loadingHistory.clear();
        crackLines.clear();
        setForce(0.0f);
    }
    
    void updateCrackVisualization() {
        if (materialState.isFractured && !materialState.crackTips.empty()) {
            crackLines.clear();
            
            // Simple crack propagation visualization
            for (const auto& crackTip : materialState.crackTips) {
                // Create crack lines
                float crackLength = materialState.damageParameter * length * 0.5f;
                glm::vec3 crackStart = crackTip - glm::vec3(crackLength/2, 0, 0);
                glm::vec3 crackEnd = crackTip + glm::vec3(crackLength/2, 0, 0);
                
                crackLines.push_back(crackStart);
                crackLines.push_back(crackEnd);
            }
        }
    }
    
    void renderDeformedBox(float l, float w, float h) {
        float halfL = l / 2.0f;
        float halfW = w / 2.0f;
        float halfH = h / 2.0f;
        
        // Apply non-uniform deformation for plastic strain
        if (materialState.plasticStrain > 0) {
            // Add necking effect (reduction in cross-section)
            float necking = 1.0f - materialState.plasticStrain * 0.5f;
            glScalef(1.0f, necking, necking);
        }
        
        glBegin(GL_QUADS);
        // Front face
        glVertex3f(-halfL, -halfH, halfW);
        glVertex3f(halfL, -halfH, halfW);
        glVertex3f(halfL, halfH, halfW);
        glVertex3f(-halfL, halfH, halfW);
        // Back face
        glVertex3f(-halfL, -halfH, -halfW);
        glVertex3f(-halfL, halfH, -halfW);
        glVertex3f(halfL, halfH, -halfW);
        glVertex3f(halfL, -halfH, -halfW);
        // Top face
        glVertex3f(-halfL, halfH, -halfW);
        glVertex3f(-halfL, halfH, halfW);
        glVertex3f(halfL, halfH, halfW);
        glVertex3f(halfL, halfH, -halfW);
        // Bottom face
        glVertex3f(-halfL, -halfH, -halfW);
        glVertex3f(halfL, -halfH, -halfW);
        glVertex3f(halfL, -halfH, halfW);
        glVertex3f(-halfL, -halfH, halfW);
        // Right face
        glVertex3f(halfL, -halfH, -halfW);
        glVertex3f(halfL, halfH, -halfW);
        glVertex3f(halfL, halfH, halfW);
        glVertex3f(halfL, -halfH, halfW);
        // Left face
        glVertex3f(-halfL, -halfH, -halfW);
        glVertex3f(-halfL, -halfH, halfW);
        glVertex3f(-halfL, halfH, halfW);
        glVertex3f(-halfL, halfH, -halfW);
        glEnd();
    }
    
    void renderCracks() {
        glLineWidth(3.0f);
        glColor3f(1.0f, 0.0f, 0.0f); // Red cracks
        
        glBegin(GL_LINES);
        for (size_t i = 0; i < crackLines.size(); i += 2) {
            glVertex3f(crackLines[i].x, crackLines[i].y, crackLines[i].z);
            glVertex3f(crackLines[i+1].x, crackLines[i+1].y, crackLines[i+1].z);
        }
        glEnd();
        
        glLineWidth(1.0f);
    }
    
    void renderStressField() {
        // Simple stress field visualization using color-coded particles
        if (stress != 0) {
            float stressLevel = std::abs(stress) / material.ultimateTensileStrength;
            
            glPointSize(3.0f);
            glBegin(GL_POINTS);
            
            for (int i = 0; i < 50; ++i) {
                float x = (i / 49.0f - 0.5f) * length;
                float intensity = stressLevel * (1.0f - std::abs(x) / (length/2));
                
                glColor3f(intensity, 0.2f, 1.0f - intensity);
                glVertex3f(x, 0.2f, 0.0f);
            }
            
            glEnd();
            glPointSize(1.0f);
        }
    }
};
