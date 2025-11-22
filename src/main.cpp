
#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>
#include <GL/glew.h>
#include <imgui.h>
#include "backends/imgui_impl_glfw.h"
#include "backends/imgui_impl_opengl3.h"
#include "implot.h"
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <iostream>
#include <vector>
#include <memory>
#include <cmath>
#include <string>
#include <map>
#include <fstream>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ============================================================================
// UTILITIES
// ============================================================================
void exportToCSV(const std::string& filename, const std::vector<std::string>& headers, const std::vector<std::vector<float>>& data) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file for writing: " << filename << std::endl;
        return;
    }
    
    // Write headers
    for (size_t i = 0; i < headers.size(); ++i) {
        file << headers[i];
        if (i < headers.size() - 1) file << ",";
    }
    file << "\n";
    
    // Write data
    // Assuming data is a vector of columns (each inner vector is a column)
    // We need to iterate rows
    if (data.empty()) return;
    size_t numRows = data[0].size();
    size_t numCols = data.size();
    
    for (size_t i = 0; i < numRows; ++i) {
        for (size_t j = 0; j < numCols; ++j) {
            if (i < data[j].size()) {
                file << data[j][i];
            }
            if (j < numCols - 1) file << ",";
        }
        file << "\n";
    }
    
    std::cout << "Exported data to " << filename << std::endl;
}

// Image Capture State
bool captureNextFrame = false;
std::string captureFilename;
ImVec2 captureMin;
ImVec2 captureMax;

void requestPlotCapture(const std::string& filename) {
    captureNextFrame = true;
    captureFilename = filename;
    // Get the rect of the last item (the plot)
    captureMin = ImGui::GetItemRectMin();
    captureMax = ImGui::GetItemRectMax();
}

// ============================================================================
// CORE PHYSICS CALCULATIONS
// ============================================================================

class PhysicsEngine {
public:
    // Solid Mechanics - Hooke's Law
    static double calculateStress(double force, double area) {
        return force / area; // Pa
    }
    
    static double calculateStrain(double stress, double youngsModulus) {
        return stress / youngsModulus; // dimensionless
    }
    
    static double calculateElongation(double strain, double length) {
        return strain * length; // m
    }
    
    // Fluid Mechanics - Buoyancy
    static double calculateBuoyantForce(double fluidDensity, double volume, double gravity = 9.81) {
        return fluidDensity * volume * gravity; // N
    }
    
    static bool willFloat(double objectDensity, double fluidDensity) {
        return objectDensity < fluidDensity;
    }
    
    // Poiseuille's Law for pipe flow
    static double calculateVelocityProfile(double radius, double r, double pressureGradient, double viscosity) {
        return (pressureGradient / (4.0 * viscosity)) * (radius * radius - r * r);
    }
    
    // Gas Laws - Ideal Gas Law: PV = nRT
    static double idealGasVolume(double pressure, double moles, double temperature, double R = 8.314) {
        return (moles * R * temperature) / pressure; // m³
    }
    
    static double idealGasPressure(double volume, double moles, double temperature, double R = 8.314) {
        return (moles * R * temperature) / volume; // Pa
    }
    
    // Phase Transitions
    static std::string getPhase(double temperature, double meltingPoint = 273.15, double boilingPoint = 373.15) {
        if (temperature < meltingPoint) return "Solid";
        else if (temperature < boilingPoint) return "Liquid";
        else return "Gas";
    }
    
    static double calculateHeatCapacity(double mass, double specificHeatCapacity, double deltaT) {
        return mass * specificHeatCapacity * deltaT; // J
    }
};

// ============================================================================
// MATERIAL DATABASE
// ============================================================================

struct Material {
    std::string name;
    double density;          // kg/m³
    double youngsModulus;    // Pa
    double viscosity;        // Pa·s (for fluids)
    double thermalConductivity; // W/(m·K)
    double specificHeatCapacity; // J/(kg·K)
    double meltingPoint;     // K
    double boilingPoint;     // K
    glm::vec3 color;
    // New fields for EnhancedSolidMechanicsSimulator
    double yieldStrength;
    double ultimateTensileStrength;
    glm::vec3 elasticColor;
    glm::vec3 plasticColor;
    glm::vec3 fractureColor;
};

class MaterialDatabase {
private:
    std::map<std::string, Material> materials;
    
public:
    MaterialDatabase() {
        // Solids
        materials["Steel"] = {"Steel", 7800, 200e9, 0.0, 50.2, 420, 1811, 3273, {0.7f, 0.7f, 0.8f}, 250e6, 400e6, {0.7f, 0.7f, 0.8f}, {0.9f, 0.5f, 0.5f}, {0.2f, 0.2f, 0.2f}};
        materials["Aluminum"] = {"Aluminum", 2700, 70e9, 0.0, 237, 900, 933, 2743, {0.8f, 0.8f, 0.9f}, 95e6, 110e6, {0.8f, 0.8f, 0.9f}, {0.9f, 0.7f, 0.7f}, {0.3f, 0.3f, 0.3f}};
        materials["Rubber"] = {"Rubber", 1500, 0.01e9, 0.0, 0.13, 2000, 200, 500, {0.2f, 0.2f, 0.2f}, 15e6, 20e6, {0.2f, 0.2f, 0.2f}, {0.8f, 0.8f, 0.2f}, {0.1f, 0.1f, 0.1f}};
        materials["Concrete"] = {"Concrete", 2400, 30e9, 0.0, 1.7, 880, 1000, 2000, {0.6f, 0.6f, 0.6f}, 5e6, 5e6, {0.6f, 0.6f, 0.6f}, {0.7f, 0.7f, 0.7f}, {0.4f, 0.4f, 0.4f}};
        
        // Fluids
        materials["Water"] = {"Water", 1000, 0.0, 0.001, 0.606, 4186, 273.15, 373.15, {0.2f, 0.4f, 0.8f}, 0, 0, {}, {}, {}};
        materials["Oil"] = {"Oil", 900, 0.0, 0.1, 0.14, 2000, 250, 573, {0.8f, 0.6f, 0.2f}, 0, 0, {}, {}, {}};
        materials["Honey"] = {"Honey", 1400, 0.0, 10.0, 0.5, 2000, 260, 400, {0.9f, 0.7f, 0.2f}, 0, 0, {}, {}, {}};
        
        // Gases
        materials["Air"] = {"Air", 1.225, 0.0, 1.8e-5, 0.026, 1005, 0, 0, {0.7f, 0.8f, 0.9f}, 0, 0, {}, {}, {}};
        materials["CO2"] = {"CO2", 1.98, 0.0, 1.5e-5, 0.016, 844, 194.7, 194.7, {0.5f, 0.9f, 0.5f}, 0, 0, {}, {}, {}};
        materials["Helium"] = {"Helium", 0.1785, 0.0, 2.0e-5, 0.152, 5193, 1, 4.2, {0.9f, 0.5f, 0.9f}, 0, 0, {}, {}, {}};
    }
    
    Material getMaterial(const std::string& name) {
        if (materials.find(name) != materials.end()) {
            return materials[name];
        }
        return materials["Steel"]; // Default
    }
    
    std::vector<std::string> getMaterialNames() {
        std::vector<std::string> names;
        for (const auto& pair : materials) {
            names.push_back(pair.first);
        }
        return names;
    }
};

// ============================================================================
// SIMULATION SCENARIOS
// ============================================================================

struct MaterialState {
    double totalStrain = 0.0;
    double elasticStrain = 0.0;
    double plasticStrain = 0.0;
    double currentYieldStrength = 0.0;
    bool isFractured = false;
    double damageParameter = 0.0;
    std::vector<glm::vec3> crackTips;
};

class AdvancedPhysicsEngine {
public:
    static std::pair<double, MaterialState> calculateStressWithPlasticity(
        double strain, const Material& material, MaterialState currentState)
    {
        // This is a simplified model for demonstration
        double stress = 0.0;

        if (strain <= currentState.elasticStrain + currentState.plasticStrain) {
            // Unloading
            currentState.elasticStrain = strain - currentState.plasticStrain;
            stress = currentState.elasticStrain * material.youngsModulus;
        } else {
            // Loading
            currentState.totalStrain = strain;
            double elasticStrainLimit = currentState.currentYieldStrength / material.youngsModulus;

            if (strain <= elasticStrainLimit) {
                // Elastic region
                currentState.elasticStrain = strain;
                currentState.plasticStrain = 0.0;
                stress = strain * material.youngsModulus;
            } else {
                // Plastic region
                currentState.elasticStrain = elasticStrainLimit;
                currentState.plasticStrain = strain - elasticStrainLimit;
                
                // Simplified linear work hardening
                double hardening = (material.ultimateTensileStrength - material.yieldStrength) / 0.1; // over 10% strain
                stress = currentState.currentYieldStrength + hardening * currentState.plasticStrain;
                
                currentState.currentYieldStrength = stress;

                if (stress > material.ultimateTensileStrength) {
                    currentState.isFractured = true;
                    currentState.damageParameter = (stress - material.ultimateTensileStrength) / (material.ultimateTensileStrength * 0.1);
                    if (currentState.crackTips.empty()) {
                        currentState.crackTips.push_back(glm::vec3(0,0,0));
                    }
                }
            }
        }
        return {stress, currentState};
    }
};

/* class SolidMechanicsSimulator {
private:
    Material material;
    float length = 1.0f;     // m
    float area = 0.01f;      // m²
    float appliedForce = 0.0f; // N
    
    // Results
    double stress = 0.0;
    double strain = 0.0;
    double elongation = 0.0;
    
    // Visualization data
    std::vector<glm::vec2> stressStrainData;
    
public:
    void setMaterial(const Material& mat) { material = mat; }
    void setDimensions(float l, float a) { length = l; area = a; }
    void setForce(float f) { appliedForce = f; calculate(); }
    
    void calculate() {
        stress = PhysicsEngine::calculateStress(appliedForce, area);
        strain = PhysicsEngine::calculateStrain(stress, material.youngsModulus);
        elongation = PhysicsEngine::calculateElongation(strain, length);
        
        // Add to stress-strain curve
        stressStrainData.push_back(glm::vec2(strain * 100, stress / 1e6)); // % strain, MPa stress
        if (stressStrainData.size() > 100) {
            stressStrainData.erase(stressStrainData.begin());
        }
    }
    
    void renderUI() {
        ImGui::Text("Solid Mechanics Simulation");
        ImGui::Separator();
        
        if (ImGui::SliderFloat("Applied Force (N)", &appliedForce, -50000, 50000)) {
            calculate();
        }
        if (ImGui::SliderFloat("Length (m)", &length, 0.1f, 5.0f)) {
            calculate();
        }
        if (ImGui::SliderFloat("Cross-sectional Area (m²)", &area, 0.001f, 0.1f)) {
            calculate();
        }
        
        if (ImGui::Button("Calculate")) {
            calculate();
        }
        
        ImGui::Separator();
        ImGui::Text("Results:");
        ImGui::Text("Stress: %.2f MPa", stress / 1e6);
        ImGui::Text("Strain: %.4f%%", strain * 100);
        ImGui::Text("Elongation: %.3f mm", elongation * 1000);
        
        // Safety factor
        double yieldStrength = material.youngsModulus * 0.002; // 0.2% offset
        double safetyFactor = (stress != 0) ? yieldStrength / std::abs(stress) : INFINITY;
        ImGui::Text("Safety Factor: %.2f", safetyFactor);
    }
    
    void renderGraph() {
        if (ImPlot::BeginPlot("Stress-Strain Curve")) {
            if (!stressStrainData.empty()) {
                std::vector<float> strainVec, stressVec;
                for (const auto& point : stressStrainData) {
                    strainVec.push_back(point.x);
                    stressVec.push_back(point.y);
                }
                ImPlot::PlotLine("Stress vs Strain", strainVec.data(), stressVec.data(), strainVec.size());
            }
            ImPlot::EndPlot();
        }
    }
    
    void render3D() {
        // Render deformed rod
        glPushMatrix();
        
        // Calculate deformation factor
        double deformationFactor = 1.0 + strain;
        
        // Set color based on stress
        double stressRatio = std::min(std::abs(stress) / (material.youngsModulus * 0.001), 1.0);
        glColor3f(stressRatio, 0.2f, 1.0f - stressRatio);
        
        // Render stretched/compressed box
        glScalef(deformationFactor, 1.0f, 1.0f);
        renderBox(length, 0.1, 0.1);
        
        glPopMatrix();
    }
    
private:
    void renderBox(float length, float width, float height) {
        glBegin(GL_QUADS);
        // Simple box rendering (6 faces)
        // Front face
        glVertex3f(-length/2, -height/2, width/2);
        glVertex3f(length/2, -height/2, width/2);
        glVertex3f(length/2, height/2, width/2);
        glVertex3f(-length/2, height/2, width/2);
        // Back face
        glVertex3f(-length/2, -height/2, -width/2);
        glVertex3f(-length/2, height/2, -width/2);
        glVertex3f(length/2, height/2, -width/2);
        glVertex3f(length/2, -height/2, -width/2);
        // Continue with other faces...
        glEnd();
    }
}; */

class EnhancedSolidMechanicsSimulator {
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
    float damageVisualization = 0.0f;
    
public:
    void setMaterial(const Material& mat) { 
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
        strain = calculateStrainFromForce(appliedForce, area, material.youngsModulus);
        
        // Use enhanced physics engine
        auto result = AdvancedPhysicsEngine::calculateStressWithPlasticity(
            strain, material, materialState);
        
        stress = result.first;
        materialState = result.second;
        elongation = materialState.totalStrain * length;
        
        // Add to stress-strain curve
        stressStrainData.push_back(glm::vec2(materialState.totalStrain * 100, stress / 1e6));
        loadingHistory.push_back(glm::vec2(glfwGetTime(), appliedForce));
        
        if (stressStrainData.size() > 1000) {
            stressStrainData.erase(stressStrainData.begin());
        }
        if (loadingHistory.size() > 1000) {
            loadingHistory.erase(loadingHistory.begin());
        }
        
        // Update crack visualization
        updateCrackVisualization();
    }
    
    void renderUI() {
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
            
            float cyclicForce = cyclicAmplitude * sin(glfwGetTime() * 2 * M_PI * cyclicFrequency);
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
    
    void renderGraph() {
        const float plot_height = (ImGui::GetContentRegionAvail().y - ImGui::GetStyle().ItemSpacing.y) * 0.5f;
        const ImVec2 plot_size(-1, plot_height > 1.0f ? plot_height : 1.0f);
        if (ImPlot::BeginPlot("Advanced Stress-Strain Curve", plot_size)) {
            if (!stressStrainData.empty()) {
                std::vector<float> strainVec, stressVec;
                for (const auto& point : stressStrainData) {
                    strainVec.push_back(point.x);
                    stressVec.push_back(point.y);
                }
                ImPlot::PlotLine("Stress vs Strain", strainVec.data(), stressVec.data(), strainVec.size());
                
                // Add material property lines
                // Yield strength line
                ImPlot::PlotInfLines("Yield Point", &material.yieldStrength, 1);
                ImPlot::PlotInfLines("Ultimate Strength", &material.ultimateTensileStrength, 1);
            }
            ImPlot::EndPlot();
        }
        
        if (ImGui::Button("Export to CSV")) {
            std::vector<std::string> headers = {"Strain (%)", "Stress (MPa)"};
            std::vector<std::vector<float>> data;
            std::vector<float> strainCol, stressCol;
            for (const auto& p : stressStrainData) {
                strainCol.push_back(p.x);
                stressCol.push_back(p.y);
            }
            data.push_back(strainCol);
            data.push_back(stressCol);
            exportToCSV("solid_mechanics_data.csv", headers, data);
        }
        ImGui::SameLine();
        if (ImGui::Button("Save Image")) {
            requestPlotCapture("solid_mechanics_plot.png");
        }
        
        // Loading history graph
        if (ImPlot::BeginPlot("Loading History", plot_size)) {
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
    }
    
    void render3D() {
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
        double stress = force / area;
        return stress / E; // Initial elastic assumption
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

class FluidMechanicsSimulator {
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
    void setObject(const Material& obj) { object = obj; }
    void setFluid(const Material& fl) { fluid = fl; }
    void setObjectProperties(float vol, float rad) { objectVolume = vol; objectRadius = rad; }
    
    void calculate() {
        weight = object.density * objectVolume * 9.81;
        buoyantForce = PhysicsEngine::calculateBuoyantForce(fluid.density, objectVolume);
        floats = PhysicsEngine::willFloat(object.density, fluid.density);
        
        // Calculate terminal velocity (simplified)
        double netForce = buoyantForce - weight;
        terminalVelocity = std::sqrt(std::abs(netForce) / (0.5 * fluid.density * 3.14159 * objectRadius * objectRadius));
        
        // Set target position for animation
        targetPosition = floats ? 0.5f : -0.5f;
    }
    
    void update(float deltaTime) {
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
    
    void renderUI() {
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
    
    void render3D() {
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
        // Render container walls
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
class GasDynamicsSimulator {
private:
    Material gas;
    float initialVolume = 0.001f; // m³
    double initialPressure = 101325; // Pa
    float initialTemperature = 300.0f; // K
    float currentTemperature = 300.0f; // K
    float moles = 0.04f; // mol
    
    // Results
    double currentVolume = 0.001;
    double currentPressure = 101325;
    
    // Graph data
    std::vector<glm::vec2> pvData; // Pressure vs Volume
    std::vector<glm::vec2> tvData; // Temperature vs Volume
    
public:
    void setGas(const Material& g) { gas = g; }
    void setInitialConditions(float vol, double press, float temp, float n) {
        initialVolume = vol;
        initialPressure = press;
        initialTemperature = temp;
        currentTemperature = temp;
        moles = n;
        calculate();
    }
    
    void setTemperature(float temp) {
        currentTemperature = temp;
        calculate();
    }
    
    void calculate() {
        // Using ideal gas law: PV = nRT
        const double R = 8.314; // J/(mol·K)
        
        // Calculate new volume at constant pressure (Charles's Law)
        currentVolume = initialVolume * (currentTemperature / initialTemperature);
        
        // Calculate pressure if volume is constrained
        currentPressure = (moles * R * currentTemperature) / currentVolume;
        
        // Add to graph data
        pvData.push_back(glm::vec2(currentVolume * 1000, currentPressure / 1000)); // L, kPa
        tvData.push_back(glm::vec2(currentTemperature - 273.15, currentVolume * 1000)); // °C, L
        
        if (pvData.size() > 100) {
            pvData.erase(pvData.begin());
            tvData.erase(tvData.begin());
        }
    }
    
    void renderUI() {
        ImGui::Text("Gas Dynamics Simulation");
        ImGui::Separator();
        
        if (ImGui::SliderFloat("Temperature (K)", &currentTemperature, 200.0f, 800.0f, "%.1f")) {
            calculate();
        }
        if (ImGui::SliderFloat("Initial Volume (L)", &initialVolume, 0.1f, 10.0f, "%.3f")) {
            calculate();
        }
        if (ImGui::SliderFloat("Moles", &moles, 0.01f, 0.1f, "%.3f")) {
            calculate();
        }
        
        ImGui::Separator();
        ImGui::Text("Results:");
        ImGui::Text("Current Volume: %.3f L", currentVolume * 1000);
        ImGui::Text("Current Pressure: %.1f kPa", currentPressure / 1000);
        ImGui::Text("Temperature: %.1f K (%.1f °C)", currentTemperature, currentTemperature - 273.15);
        
        // Gas law verification
        double pv_nrt = (currentPressure * currentVolume) / (moles * 8.314 * currentTemperature);
        ImGui::Text("PV/nRT = %.3f (should ≈ 1.0)", pv_nrt);
    }
    
    void renderGraphs() {
        const float plot_height = (ImGui::GetContentRegionAvail().y - ImGui::GetStyle().ItemSpacing.y) * 0.5f;
        const ImVec2 plot_size(-1, plot_height > 1.0f ? plot_height : 1.0f);
        if (ImPlot::BeginPlot("Pressure vs Volume", plot_size)) {
            if (!pvData.empty()) {
                std::vector<float> volumeVec, pressureVec;
                for (const auto& point : pvData) {
                    volumeVec.push_back(point.x);
                    pressureVec.push_back(point.y);
                }
                ImPlot::PlotLine("P vs V", volumeVec.data(), pressureVec.data(), volumeVec.size());
            }
            ImPlot::EndPlot();
        }
        
        if (ImGui::Button("Export PV Data")) {
            std::vector<std::string> headers = {"Volume (L)", "Pressure (kPa)"};
            std::vector<std::vector<float>> data;
            std::vector<float> volCol, pressCol;
            for (const auto& p : pvData) {
                volCol.push_back(p.x);
                pressCol.push_back(p.y);
            }
            data.push_back(volCol);
            data.push_back(pressCol);
            exportToCSV("gas_pv_data.csv", headers, data);
        }
        ImGui::SameLine();
        if (ImGui::Button("Save PV Plot")) {
            requestPlotCapture("gas_pv_plot.png");
        }
        
        if (ImPlot::BeginPlot("Temperature vs Volume", plot_size)) {
            if (!tvData.empty()) {
                std::vector<float> tempVec, volumeVec;
                for (const auto& point : tvData) {
                    tempVec.push_back(point.x);
                    volumeVec.push_back(point.y);
                }
                ImPlot::PlotLine("T vs V", tempVec.data(), volumeVec.data(), tempVec.size());
            }
            ImPlot::EndPlot();
        }
        if (ImGui::Button("Save TV Plot")) {
            requestPlotCapture("gas_tv_plot.png");
        }
    }
    
    void render3D() {
        // Render gas container with particle animation
        glColor3f(gas.color.x, gas.color.y, gas.color.z);
        
        // Container walls
        glColor4f(0.8f, 0.8f, 0.8f, 0.2f);
        renderContainer();
        
        // Animated particles
        int numParticles = 50;
        float containerSize = std::pow(currentVolume, 1.0/3.0) * 2.0f; // Cube root for size
        float speed = std::sqrt(currentTemperature / 300.0f); // Particle speed based on temperature
        
        glColor3f(gas.color.x, gas.color.y, gas.color.z);
        for (int i = 0; i < numParticles; ++i) {
            float time = glfwGetTime();
            float x = cos(time * speed + i) * containerSize * 0.4f;
            float y = sin(time * speed * 1.3f + i) * containerSize * 0.4f;
            float z = cos(time * speed * 0.7f + i) * containerSize * 0.4f;
            
            glPushMatrix();
            glTranslatef(x, y, z);
            renderSphere(0.01f);
            glPopMatrix();
        }
    }
    
private:
    void renderContainer() {
        float s = std::pow(currentVolume, 1.0/3.0) / 2.0f;
        
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        
        // Render container as wireframe box
        glBegin(GL_LINES);
        // Bottom face
        glVertex3f(-s, -s, -s); glVertex3f( s, -s, -s);
        glVertex3f( s, -s, -s); glVertex3f( s, -s,  s);
        glVertex3f( s, -s,  s); glVertex3f(-s, -s,  s);
        glVertex3f(-s, -s,  s); glVertex3f(-s, -s, -s);
        // Top face
        glVertex3f(-s,  s, -s); glVertex3f( s,  s, -s);
        glVertex3f( s,  s, -s); glVertex3f( s,  s,  s);
        glVertex3f( s,  s,  s); glVertex3f(-s,  s,  s);
        glVertex3f(-s,  s,  s); glVertex3f(-s,  s, -s);
        // Vertical edges
        glVertex3f(-s, -s, -s); glVertex3f(-s,  s, -s);
        glVertex3f( s, -s, -s); glVertex3f( s,  s, -s);
        glVertex3f( s, -s,  s); glVertex3f( s,  s,  s);
        glVertex3f(-s, -s,  s); glVertex3f(-s,  s,  s);
        glEnd();
        
        glDisable(GL_BLEND);
    }
    
    void renderSphere(float radius) {
        // Simple particle rendering (low-poly for performance)
        const int stacks = 8;
        const int slices = 8;
        
        for (int i = 0; i < stacks; ++i) {
            float phi1 = M_PI * i / stacks;
            float phi2 = M_PI * (i + 1) / stacks;
            
            glBegin(GL_TRIANGLE_STRIP);
            for (int j = 0; j <= slices; ++j) {
                float theta = 2 * M_PI * j / slices;
                
                glVertex3f(radius * sin(phi1) * cos(theta), radius * cos(phi1), radius * sin(phi1) * sin(theta));
                glVertex3f(radius * sin(phi2) * cos(theta), radius * cos(phi2), radius * sin(phi2) * sin(theta));
            }
            glEnd();
        }
    }
};

class PhaseChangeSimulator {
private:
    Material material;
    float currentTemperature = 273.15f; // K
    float mass = 0.1f; // kg
    double heatInput = 0.0;      // J
    
    std::string currentPhase = "Solid";
    double meltingProgress = 0.0; // 0-1
    double boilingProgress = 0.0; // 0-1
    
    // Graph data
    std::vector<glm::vec2> tempTimeData;
    
public:
    void setMaterial(const Material& mat) { material = mat; }
    void setMass(float m) { mass = m; }
    void setTemperature(float temp) {
        currentTemperature = temp;
        updatePhase();

        // Add to graph to reflect slider changes
        tempTimeData.push_back(glm::vec2(glfwGetTime(), currentTemperature - 273.15));
        if (tempTimeData.size() > 200) {
            tempTimeData.erase(tempTimeData.begin());
        }
    }
    
    void addHeat(double heat) {
        heatInput += heat;
        
        // Calculate temperature change based on heat capacity
        if (currentPhase == "Solid" && currentTemperature < material.meltingPoint) {
            double deltaT = heat / (mass * material.specificHeatCapacity);
            currentTemperature += deltaT;
        } else if (currentPhase == "Liquid" && currentTemperature < material.boilingPoint) {
            double deltaT = heat / (mass * material.specificHeatCapacity);
            currentTemperature += deltaT;
        }
        
        updatePhase();
        
        // Add to graph
        tempTimeData.push_back(glm::vec2(glfwGetTime(), currentTemperature - 273.15));
        if (tempTimeData.size() > 200) {
            tempTimeData.erase(tempTimeData.begin());
        }
    }
    
    void updatePhase() {
        currentPhase = PhysicsEngine::getPhase(currentTemperature, material.meltingPoint, material.boilingPoint);
        
        // Calculate transition progress
        if (currentTemperature >= material.meltingPoint - 10 && currentTemperature <= material.meltingPoint + 10) {
            meltingProgress = (currentTemperature - (material.meltingPoint - 10)) / 20.0;
            meltingProgress = std::max(0.0, std::min(1.0, meltingProgress));
        }
        
        if (currentTemperature >= material.boilingPoint - 10 && currentTemperature <= material.boilingPoint + 10) {
            boilingProgress = (currentTemperature - (material.boilingPoint - 10)) / 20.0;
            boilingProgress = std::max(0.0, std::min(1.0, boilingProgress));
        }
    }
    
    void renderUI() {
        ImGui::Text("Phase Change Simulation");
        ImGui::Separator();
        
        if (ImGui::SliderFloat("Temperature (K)", &currentTemperature, 200.0f, 500.0f, "%.1f")) {
            // We call setTemperature to ensure the phase is updated and data is added to the graph
            setTemperature(currentTemperature);
        }
        if (ImGui::SliderFloat("Mass (kg)", &mass, 0.01f, 1.0f, "%.2f")) {
            // Mass change might affect future heat calculations, no need to recalculate immediately
        }
        
        if (ImGui::Button("Add 1000J Heat")) {
            addHeat(1000.0);
        }
        ImGui::SameLine();
        if (ImGui::Button("Remove 1000J Heat")) {
            addHeat(-1000.0);
        }
        
        ImGui::Separator();
        ImGui::Text("Results:");
        ImGui::Text("Current Phase: %s", currentPhase.c_str());
        ImGui::Text("Temperature: %.1f K (%.1f °C)", currentTemperature, currentTemperature - 273.15);
        ImGui::Text("Melting Point: %.1f K", material.meltingPoint);
        ImGui::Text("Boiling Point: %.1f K", material.boilingPoint);
        ImGui::Text("Total Heat Input: %.0f J", heatInput);
        
        if (meltingProgress > 0 && meltingProgress < 1) {
            ImGui::Text("Melting Progress: %.1f%%", meltingProgress * 100);
        }
        if (boilingProgress > 0 && boilingProgress < 1) {
            ImGui::Text("Boiling Progress: %.1f%%", boilingProgress * 100);
        }
    }
    
    void renderGraph() {
        if (ImPlot::BeginPlot("Temperature vs Time", ImVec2(-1, -FLT_MIN))) {
            if (!tempTimeData.empty()) {
                std::vector<float> timeVec, tempVec;
                for (const auto& point : tempTimeData) {
                    timeVec.push_back(point.x);
                    tempVec.push_back(point.y);
                }
                ImPlot::PlotLine("Temperature (°C)", timeVec.data(), tempVec.data(), timeVec.size());
                
                // Add phase transition lines
                ImPlot::PlotInfLines("Melting Point", &material.meltingPoint, 1);
                ImPlot::PlotInfLines("Boiling Point", &material.boilingPoint, 1);
            }
            ImPlot::EndPlot();
        }
        
        if (ImGui::Button("Export Temp Data")) {
            std::vector<std::string> headers = {"Time (s)", "Temperature (C)"};
            std::vector<std::vector<float>> data;
            std::vector<float> timeCol, tempCol;
            for (const auto& p : tempTimeData) {
                timeCol.push_back(p.x);
                tempCol.push_back(p.y);
            }
            data.push_back(timeCol);
            data.push_back(tempCol);
            exportToCSV("phase_change_data.csv", headers, data);
        }
        ImGui::SameLine();
        if (ImGui::Button("Save Image")) {
            requestPlotCapture("phase_change_plot.png");
        }
    }
    
    void render3D() {
        glPushMatrix();
        
        // Render material based on current phase
        if (currentPhase == "Solid") {
            glColor3f(material.color.x * 0.8f, material.color.y * 0.8f, material.color.z * 0.8f);
            renderCube(0.2f);
        } else if (currentPhase == "Liquid") {
            glColor4f(material.color.x, material.color.y, material.color.z, 0.7f);
            renderLiquid();
        } else if (currentPhase == "Gas") {
            glColor4f(material.color.x, material.color.y, material.color.z, 0.3f);
            renderGas();
        }
        
        glPopMatrix();
    }
    
private:
    void renderCube(float size) {
        float s = size / 2.0f;
        // Render solid cube
        glBegin(GL_QUADS);
        // Front Face
        glVertex3f(-s, -s,  s); glVertex3f( s, -s,  s); glVertex3f( s,  s,  s); glVertex3f(-s,  s,  s);
        // Back Face
        glVertex3f(-s, -s, -s); glVertex3f(-s,  s, -s); glVertex3f( s,  s, -s); glVertex3f( s, -s, -s);
        // Top Face
        glVertex3f(-s,  s, -s); glVertex3f(-s,  s,  s); glVertex3f( s,  s,  s); glVertex3f( s,  s, -s);
        // Bottom Face
        glVertex3f(-s, -s, -s); glVertex3f( s, -s, -s); glVertex3f( s, -s,  s); glVertex3f(-s, -s,  s);
        // Right face
        glVertex3f( s, -s, -s); glVertex3f( s,  s, -s); glVertex3f( s,  s,  s); glVertex3f( s, -s,  s);
        // Left Face
        glVertex3f(-s, -s, -s); glVertex3f(-s, -s,  s); glVertex3f(-s,  s,  s); glVertex3f(-s,  s, -s);
        glEnd();
    }
    
    void renderLiquid() {
        // Render liquid with flowing animation
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        
        // Animated liquid surface
        float time = glfwGetTime();
        const int segments = 20;
        const float size = 0.4f;

        glBegin(GL_TRIANGLE_STRIP);
        for (int i = 0; i <= segments; ++i) {
            float x = (i / (float)segments - 0.5f) * size * 2.0f;
            float wave = sin(time * 2.0f + x * 5.0f) * 0.05f + cos(time * 1.5f + x * 3.0f) * 0.03f;
            
            glVertex3f(x, wave, -size);
            glVertex3f(x, wave,  size);
        }
        glEnd();
        
        glDisable(GL_BLEND);
    }
    
    void renderGas() {
        // Render gas as moving particles
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        
        int numParticles = 30;
        float time = glfwGetTime();
        
        for (int i = 0; i < numParticles; ++i) {
            float x = cos(time * 3 + i) * 0.5f;
            float y = sin(time * 2.5f + i) * 0.5f;
            float z = cos(time * 1.8f + i) * 0.5f;
            
            glPushMatrix();
            glTranslatef(x, y, z);
            renderSphere(0.01f);
            glPopMatrix();
        }
        
        glDisable(GL_BLEND);
    }
    
    void renderSphere(float radius) {
        // Simple sphere rendering (low-poly for performance)
        const int stacks = 8;
        const int slices = 8;
        
        for (int i = 0; i < stacks; ++i) {
            float phi1 = M_PI * i / stacks;
            float phi2 = M_PI * (i + 1) / stacks;
            
            glBegin(GL_TRIANGLE_STRIP);
            for (int j = 0; j <= slices; ++j) {
                float theta = 2 * M_PI * j / slices;
                
                glVertex3f(radius * sin(phi1) * cos(theta), radius * cos(phi1), radius * sin(phi1) * sin(theta));
                glVertex3f(radius * sin(phi2) * cos(theta), radius * cos(phi2), radius * sin(phi2) * sin(theta));
            }
            glEnd();
        }
    }
};

class ThermalConductivitySimulator {
private:
    Material material;
    float length = 1.0f;
    float area = 0.01f;
    float leftTemp = 373.15f; // 100 C
    float rightTemp = 273.15f; // 0 C
    float simulationSpeed = 1.0f;
    
    // Simulation state
    std::vector<float> temperatures;
    int numSegments = 50;
    
    // Graph data
    std::vector<glm::vec2> tempProfileData;
    
public:
    ThermalConductivitySimulator() {
        temperatures.resize(numSegments, 293.15f); // Start at room temp
    }

    void setMaterial(const Material& mat) { material = mat; }
    
    void update(float deltaTime) {
        // 1D Heat Equation: dT/dt = alpha * d^2T/dx^2
        // alpha = k / (rho * cp)
        double alpha = material.thermalConductivity / (material.density * material.specificHeatCapacity);
        
        float dx = length / numSegments;
        float dt = 0.001f * simulationSpeed; // Small time step for stability
        
        // Stability check: dt <= dx^2 / (2*alpha)
        float maxDt = (dx * dx) / (2.0f * alpha);
        if (dt > maxDt) dt = maxDt * 0.9f;
        
        std::vector<float> newTemps = temperatures;
        
        // Boundary conditions
        newTemps[0] = leftTemp;
        newTemps[numSegments-1] = rightTemp;
        
        // Finite difference
        float factor = alpha * dt / (dx * dx);
        
        for (int i = 1; i < numSegments - 1; ++i) {
            newTemps[i] = temperatures[i] + factor * (temperatures[i+1] - 2*temperatures[i] + temperatures[i-1]);
        }
        
        temperatures = newTemps;
        
        // Update graph data
        tempProfileData.clear();
        for (int i = 0; i < numSegments; ++i) {
            float x = (float)i / (numSegments - 1) * length;
            tempProfileData.push_back(glm::vec2(x, temperatures[i] - 273.15f));
        }
    }
    
    void renderUI() {
        ImGui::Text("Thermal Conductivity Simulation");
        ImGui::Separator();
        
        ImGui::SliderFloat("Left Temp (K)", &leftTemp, 200.0f, 500.0f);
        ImGui::SliderFloat("Right Temp (K)", &rightTemp, 200.0f, 500.0f);
        ImGui::SliderFloat("Simulation Speed", &simulationSpeed, 0.0f, 10.0f);
        
        ImGui::Separator();
        ImGui::Text("Material Properties:");
        ImGui::Text("Conductivity: %.2f W/(m·K)", material.thermalConductivity);
        ImGui::Text("Diffusivity: %.2e m²/s", material.thermalConductivity / (material.density * material.specificHeatCapacity));
    }
    
    void renderGraph() {
        if (ImPlot::BeginPlot("Temperature Profile")) {
            if (!tempProfileData.empty()) {
                std::vector<float> xVec, tVec;
                for (const auto& point : tempProfileData) {
                    xVec.push_back(point.x);
                    tVec.push_back(point.y);
                }
                ImPlot::PlotLine("Temp vs Position", xVec.data(), tVec.data(), xVec.size());
            }
            ImPlot::EndPlot();
        }
    }
    
    void render3D() {
        // Render rod with temperature gradient colors
        float halfH = sqrt(area) / 2.0f;
        
        glBegin(GL_QUADS);
        for (int i = 0; i < numSegments - 1; ++i) {
            float x1 = (float)i / numSegments * length - length/2;
            float x2 = (float)(i+1) / numSegments * length - length/2;
            
            glm::vec3 color1 = getTempColor(temperatures[i]);
            glm::vec3 color2 = getTempColor(temperatures[i+1]);
            
            // Front face
            glColor3f(color1.x, color1.y, color1.z);
            glVertex3f(x1, -halfH, halfH);
            glVertex3f(x1, halfH, halfH);
            glColor3f(color2.x, color2.y, color2.z);
            glVertex3f(x2, halfH, halfH);
            glVertex3f(x2, -halfH, halfH);
            
            // Top face
            glColor3f(color1.x, color1.y, color1.z);
            glVertex3f(x1, halfH, halfH);
            glVertex3f(x1, halfH, -halfH);
            glColor3f(color2.x, color2.y, color2.z);
            glVertex3f(x2, halfH, -halfH);
            glVertex3f(x2, halfH, halfH);
            
            // Back face
            glColor3f(color1.x, color1.y, color1.z);
            glVertex3f(x1, halfH, -halfH);
            glVertex3f(x1, -halfH, -halfH);
            glColor3f(color2.x, color2.y, color2.z);
            glVertex3f(x2, -halfH, -halfH);
            glVertex3f(x2, halfH, -halfH);

            // Bottom face
            glColor3f(color1.x, color1.y, color1.z);
            glVertex3f(x1, -halfH, -halfH);
            glVertex3f(x1, -halfH, halfH);
            glColor3f(color2.x, color2.y, color2.z);
            glVertex3f(x2, -halfH, halfH);
            glVertex3f(x2, -halfH, -halfH);
        }
        glEnd();
    }
    
private:
    glm::vec3 getTempColor(float temp) {
        // Map temp to Blue-Red gradient
        float t = (temp - 273.15f) / 100.0f; // 0 to 100 C range
        t = std::max(0.0f, std::min(1.0f, t));
        return glm::vec3(t, 0.0f, 1.0f - t);
    }
};

class ProjectileMotionSimulator {
private:
    Material projectile;
    float mass = 1.0f;
    float radius = 0.1f;
    float launchVelocity = 20.0f;
    float launchAngle = 45.0f;
    float gravity = 9.81f;
    float dragCoefficient = 0.47f; // Sphere
    bool useAirResistance = false;
    
    // Simulation state
    glm::vec3 position;
    glm::vec3 velocity;
    bool isSimulating = false;
    float time = 0.0f;
    
    // Trail
    std::vector<glm::vec3> trail;
    
public:
    ProjectileMotionSimulator() {
        reset();
    }
    
    void setMaterial(const Material& mat) { projectile = mat; }
    
    void reset() {
        position = glm::vec3(-4.0f, 0.0f, 0.0f);
        velocity = glm::vec3(0.0f);
        isSimulating = false;
        time = 0.0f;
        trail.clear();
        trail.push_back(position);
    }
    
    void fire() {
        reset();
        isSimulating = true;
        float angleRad = glm::radians(launchAngle);
        velocity = glm::vec3(launchVelocity * cos(angleRad), launchVelocity * sin(angleRad), 0.0f);
    }
    
    void update(float deltaTime) {
        if (!isSimulating) return;
        
        // Physics update
        glm::vec3 force(0.0f, -mass * gravity, 0.0f);
        
        if (useAirResistance) {
            float v = glm::length(velocity);
            if (v > 0) {
                float airDensity = 1.225f;
                float area = M_PI * radius * radius;
                float dragForce = 0.5f * airDensity * v * v * dragCoefficient * area;
                force -= glm::normalize(velocity) * dragForce;
            }
        }
        
        glm::vec3 acceleration = force / mass;
        velocity += acceleration * deltaTime;
        position += velocity * deltaTime;
        
        // Ground collision
        if (position.y < 0.0f) {
            position.y = 0.0f;
            isSimulating = false;
        }
        
        // Update trail
        if (glm::length(position - trail.back()) > 0.1f) {
            trail.push_back(position);
        }
    }
    
    void renderUI() {
        ImGui::Text("Projectile Motion Simulation");
        ImGui::Separator();
        
        ImGui::SliderFloat("Velocity (m/s)", &launchVelocity, 1.0f, 50.0f);
        ImGui::SliderFloat("Angle (deg)", &launchAngle, 0.0f, 90.0f);
        ImGui::SliderFloat("Gravity (m/s²)", &gravity, 1.0f, 20.0f);
        ImGui::Checkbox("Air Resistance", &useAirResistance);
        
        if (ImGui::Button("FIRE")) {
            fire();
        }
        ImGui::SameLine();
        if (ImGui::Button("Reset")) {
            reset();
        }
        
        ImGui::Separator();
        ImGui::Text("Status:");
        ImGui::Text("Height: %.2f m", position.y);
        ImGui::Text("Distance: %.2f m", position.x - (-4.0f));
        ImGui::Text("Speed: %.2f m/s", glm::length(velocity));
    }
    
    void renderGraph() {
        // Maybe plot height vs time or something?
        // For now, just rely on 3D view
        ImGui::Text("See 3D view for trajectory");
    }
    
    void render3D() {
        // Render ground
        glColor3f(0.2f, 0.5f, 0.2f);
        glBegin(GL_QUADS);
        glVertex3f(-10, -0.01f, -2);
        glVertex3f(10, -0.01f, -2);
        glVertex3f(10, -0.01f, 2);
        glVertex3f(-10, -0.01f, 2);
        glEnd();
        
        // Render projectile
        glPushMatrix();
        glTranslatef(position.x, position.y + radius, position.z);
        glColor3f(projectile.color.x, projectile.color.y, projectile.color.z);
        renderSphere(radius);
        glPopMatrix();
        
        // Render trail
        glLineWidth(2.0f);
        glColor3f(1.0f, 1.0f, 0.0f);
        glBegin(GL_LINE_STRIP);
        for (const auto& pos : trail) {
            glVertex3f(pos.x, pos.y + radius, pos.z);
        }
        // Connect to current pos
        glVertex3f(position.x, position.y + radius, position.z);
        glEnd();
        glLineWidth(1.0f);
        
        // Render cannon/launch point
        glPushMatrix();
        glTranslatef(-4.0f, 0.0f, 0.0f);
        glColor3f(0.4f, 0.4f, 0.4f);
        // Simple cannon base
        renderBox(0.5f, 0.5f, 0.5f);
        // Cannon barrel
        glRotatef(launchAngle, 0.0f, 0.0f, 1.0f);
        glTranslatef(0.5f, 0.0f, 0.0f);
        renderBox(1.0f, 0.2f, 0.2f);
        glPopMatrix();
    }
    
private:
    void renderSphere(float r) {
         const int stacks = 10;
        const int slices = 10;
        for (int i = 0; i < stacks; ++i) {
            float phi1 = M_PI * i / stacks;
            float phi2 = M_PI * (i + 1) / stacks;
            glBegin(GL_TRIANGLE_STRIP);
            for (int j = 0; j <= slices; ++j) {
                float theta = 2 * M_PI * j / slices;
                glVertex3f(r * sin(phi1) * cos(theta), r * cos(phi1), r * sin(phi1) * sin(theta));
                glVertex3f(r * sin(phi2) * cos(theta), r * cos(phi2), r * sin(phi2) * sin(theta));
            }
            glEnd();
        }
    }
    
    void renderBox(float l, float h, float w) {
        float hl = l/2, hh = h/2, hw = w/2;
        glBegin(GL_QUADS);
         glVertex3f(-hl, -hh, hw); glVertex3f(hl, -hh, hw); glVertex3f(hl, hh, hw); glVertex3f(-hl, hh, hw);
        glVertex3f(-hl, -hh, -hw); glVertex3f(-hl, hh, -hw); glVertex3f(hl, hh, -hw); glVertex3f(hl, -hh, -hw);
        glVertex3f(-hl, hh, -hw); glVertex3f(-hl, hh, hw); glVertex3f(hl, hh, hw); glVertex3f(hl, hh, -hw);
        glVertex3f(-hl, -hh, -hw); glVertex3f(hl, -hh, -hw); glVertex3f(hl, -hh, hw); glVertex3f(-hl, -hh, hw);
        glVertex3f(hl, -hh, -hw); glVertex3f(hl, hh, -hw); glVertex3f(hl, hh, hw); glVertex3f(hl, -hh, hw);
        glVertex3f(-hl, -hh, -hw); glVertex3f(-hl, -hh, hw); glVertex3f(-hl, hh, hw); glVertex3f(-hl, hh, -hw);
        glEnd();
    }
};

// ============================================================================
// MAIN APPLICATION CLASS
// ============================================================================

class OrbitalMechanicsSimulator {
private:
    struct Body {
        glm::vec3 position;
        glm::vec3 velocity;
        float mass;
        float radius;
        glm::vec3 color;
        std::vector<glm::vec3> trail;
    };
    
    std::vector<Body> bodies;
    float G = 1.0f; // Gravitational constant (scaled for visual effect)
    float timeScale = 1.0f;
    bool isPaused = false;
    
public:
    OrbitalMechanicsSimulator() {
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
            glm::vec3(1.0f, 1.0f, 0.0f)
        });
        
        // Earth-like
        bodies.push_back({
            glm::vec3(10,0,0),
            glm::vec3(0,0,10), // Velocity for circular orbit: v = sqrt(GM/r) = sqrt(1000/10) = 10
            10.0f,
            0.3f,
            glm::vec3(0.2f, 0.4f, 1.0f)
        });
        
        // Moon-like
        bodies.push_back({
            glm::vec3(11,0,0),
            glm::vec3(0,0,13), // Earth v + Moon v_rel
            1.0f,
            0.1f,
            glm::vec3(0.6f, 0.6f, 0.6f)
        });
    }
    
    void addRandomPlanet() {
        float r = 5.0f + (rand() % 200) / 10.0f;
        float theta = (rand() % 360) * M_PI / 180.0f;
        
        float x = r * cos(theta);
        float z = r * sin(theta);
        
        // Orbital velocity for circular orbit
        float v = sqrt(G * 1000.0f / r);
        float vx = -v * sin(theta);
        float vz = v * cos(theta);
        
        bodies.push_back({
            glm::vec3(x, 0, z),
            glm::vec3(vx, 0, vz),
            1.0f + (rand() % 50) / 10.0f,
            0.2f + (rand() % 30) / 100.0f,
            glm::vec3((rand()%10)/10.0f, (rand()%10)/10.0f, (rand()%10)/10.0f)
        });
    }
    
    void update(float deltaTime) {
        if (isPaused) return;
        
        float dt = deltaTime * timeScale;
        
        // N-body gravity
        for (size_t i = 0; i < bodies.size(); ++i) {
            glm::vec3 force(0.0f);
            
            for (size_t j = 0; j < bodies.size(); ++j) {
                if (i == j) continue;
                
                glm::vec3 diff = bodies[j].position - bodies[i].position;
                float dist = glm::length(diff);
                
                if (dist > 0.1f) { // Avoid singularity
                    float f = G * bodies[i].mass * bodies[j].mass / (dist * dist);
                    force += glm::normalize(diff) * f;
                }
            }
            
            // Symplectic Euler integration
            glm::vec3 acceleration = force / bodies[i].mass;
            bodies[i].velocity += acceleration * dt;
        }
        
        for (auto& body : bodies) {
            body.position += body.velocity * dt;
            
            // Update trail
            if (body.trail.empty() || glm::length(body.position - body.trail.back()) > 0.2f) {
                body.trail.push_back(body.position);
                if (body.trail.size() > 200) {
                    body.trail.erase(body.trail.begin());
                }
            }
        }
    }
    
    void renderUI() {
        ImGui::Text("Orbital Mechanics Simulation");
        ImGui::Separator();
        
        ImGui::SliderFloat("Time Scale", &timeScale, 0.1f, 5.0f);
        ImGui::SliderFloat("Gravitational Constant (G)", &G, 0.1f, 5.0f);
        ImGui::Checkbox("Pause", &isPaused);
        
        if (ImGui::Button("Reset Solar System")) {
            resetSolarSystem();
        }
        if (ImGui::Button("Add Random Planet")) {
            addRandomPlanet();
        }
        
        ImGui::Separator();
        ImGui::Text("Bodies: %zu", bodies.size());
    }
    
    void renderGraph() {
        ImGui::Text("See 3D view for orbits");
    }
    
    void render3D() {
        for (const auto& body : bodies) {
            // Render body
            glPushMatrix();
            glTranslatef(body.position.x, body.position.y, body.position.z);
            glColor3f(body.color.x, body.color.y, body.color.z);
            renderSphere(body.radius);
            glPopMatrix();
            
            // Render trail
            glLineWidth(1.0f);
            glColor3f(body.color.x, body.color.y, body.color.z);
            glBegin(GL_LINE_STRIP);
            for (const auto& pos : body.trail) {
                glVertex3f(pos.x, pos.y, pos.z);
            }
            glVertex3f(body.position.x, body.position.y, body.position.z);
            glEnd();
        }
    }
    
private:
    void renderSphere(float r) {
         const int stacks = 10;
        const int slices = 10;
        for (int i = 0; i < stacks; ++i) {
            float phi1 = M_PI * i / stacks;
            float phi2 = M_PI * (i + 1) / stacks;
            glBegin(GL_TRIANGLE_STRIP);
            for (int j = 0; j <= slices; ++j) {
                float theta = 2 * M_PI * j / slices;
                glVertex3f(r * sin(phi1) * cos(theta), r * cos(phi1), r * sin(phi1) * sin(theta));
                glVertex3f(r * sin(phi2) * cos(theta), r * cos(phi2), r * sin(phi2) * sin(theta));
            }
            glEnd();
        }
    }
};

class PipeFlowSimulator {
private:
    Material fluid;
    float pipeRadius = 0.5f;
    float pipeLength = 4.0f;
    float pressureDiff = 100.0f; // Pa
    
    // Simulation state
    struct Particle {
        glm::vec3 position;
        float speed;
        glm::vec3 color;
    };
    std::vector<Particle> particles;
    float flowRate = 0.0f;
    float maxVelocity = 0.0f;
    
public:
    PipeFlowSimulator() {
        resetParticles();
    }
    
    void setFluid(const Material& f) { fluid = f; }
    
    void resetParticles() {
        particles.clear();
        for (int i = 0; i < 200; ++i) {
            spawnParticle();
        }
    }
    
    void spawnParticle() {
        // Random radial position
        float r = sqrt((float)rand() / RAND_MAX) * pipeRadius;
        float theta = ((float)rand() / RAND_MAX) * 2 * M_PI;
        
        // Random longitudinal position
        float x = ((float)rand() / RAND_MAX) * pipeLength - pipeLength/2;
        
        particles.push_back({
            glm::vec3(x, r * cos(theta), r * sin(theta)),
            0.0f,
            glm::vec3(0.4f, 0.6f, 1.0f)
        });
    }
    
    void update(float deltaTime) {
        // Poiseuille's Law calculations
        // v(r) = (DeltaP / (4 * mu * L)) * (R^2 - r^2)
        
        float viscosity = std::max(fluid.viscosity, 0.0001); // Avoid div by zero
        float factor = pressureDiff / (4.0f * viscosity * pipeLength);
        
        maxVelocity = factor * pipeRadius * pipeRadius;
        flowRate = (M_PI * pipeRadius * pipeRadius * pipeRadius * pipeRadius * pressureDiff) / (8.0f * viscosity * pipeLength);
        
        for (auto& p : particles) {
            float r = sqrt(p.position.y * p.position.y + p.position.z * p.position.z);
            
            // Calculate velocity at this radius
            if (r < pipeRadius) {
                p.speed = factor * (pipeRadius * pipeRadius - r * r);
            } else {
                p.speed = 0.0f;
            }
            
            // Move particle
            p.position.x += p.speed * deltaTime;
            
            // Wrap around
            if (p.position.x > pipeLength/2) {
                p.position.x = -pipeLength/2;
                // Randomize radial position slightly for variety
                float newR = sqrt((float)rand() / RAND_MAX) * pipeRadius;
                float theta = atan2(p.position.z, p.position.y);
                p.position.y = newR * cos(theta);
                p.position.z = newR * sin(theta);
            }
            
            // Color based on speed
            float speedRatio = (maxVelocity > 0) ? p.speed / maxVelocity : 0.0f;
            p.color = glm::vec3(speedRatio, 0.5f, 1.0f - speedRatio);
        }
    }
    
    void renderUI() {
        ImGui::Text("Laminar Pipe Flow Simulation");
        ImGui::Separator();
        
        ImGui::SliderFloat("Pressure Diff (Pa)", &pressureDiff, 0.0f, 1000.0f);
        ImGui::SliderFloat("Pipe Radius (m)", &pipeRadius, 0.1f, 1.0f);
        ImGui::SliderFloat("Pipe Length (m)", &pipeLength, 1.0f, 10.0f);
        
        ImGui::Separator();
        ImGui::Text("Fluid: %s", fluid.name.c_str());
        ImGui::Text("Viscosity: %.4f Pa·s", fluid.viscosity);
        
        ImGui::Separator();
        ImGui::Text("Results:");
        ImGui::Text("Max Velocity: %.2f m/s", maxVelocity);
        ImGui::Text("Flow Rate: %.4f m³/s", flowRate);
        
        // Reynolds Number check
        float density = fluid.density;
        float avgVelocity = maxVelocity / 2.0f;
        float diameter = 2 * pipeRadius;
        float Re = (density * avgVelocity * diameter) / std::max(fluid.viscosity, 0.0001);
        
        ImGui::Text("Reynolds Number: %.0f", Re);
        if (Re > 2300) {
            ImGui::TextColored(ImVec4(1,0,0,1), "WARNING: Flow may be turbulent!");
            ImGui::Text("(Simulation assumes laminar)");
        } else {
            ImGui::TextColored(ImVec4(0,1,0,1), "Flow is Laminar");
        }
    }
    
    void renderGraph() {
        if (ImPlot::BeginPlot("Velocity Profile")) {
            std::vector<float> rVec, vVec;
            for (float r = -pipeRadius; r <= pipeRadius; r += pipeRadius/20.0f) {
                rVec.push_back(r);
                float v = (pressureDiff / (4.0f * std::max(fluid.viscosity, 0.0001) * pipeLength)) * (pipeRadius * pipeRadius - r * r);
                vVec.push_back(v);
            }
            ImPlot::PlotLine("Velocity vs Radius", rVec.data(), vVec.data(), rVec.size());
            ImPlot::EndPlot();
        }
        
        if (ImGui::Button("Export Velocity Data")) {
            std::vector<std::string> headers = {"Radius (m)", "Velocity (m/s)"};
            std::vector<std::vector<float>> data;
            std::vector<float> radCol, velCol;
            for (float r = -pipeRadius; r <= pipeRadius; r += pipeRadius/20.0f) {
                radCol.push_back(r);
                float v = (pressureDiff / (4.0f * std::max(fluid.viscosity, 0.0001) * pipeLength)) * (pipeRadius * pipeRadius - r * r);
                velCol.push_back(v);
            }
            data.push_back(radCol);
            data.push_back(velCol);
            exportToCSV("pipe_flow_data.csv", headers, data);
        }
        ImGui::SameLine();
        if (ImGui::Button("Save Image")) {
            requestPlotCapture("pipe_flow_plot.png");
        }
    }
    
    void render3D() {
        // Render pipe (transparent)
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glColor4f(0.8f, 0.8f, 0.8f, 0.2f);
        
        const int segments = 20;
        glBegin(GL_QUAD_STRIP);
        for (int i = 0; i <= segments; ++i) {
            float theta = 2.0f * M_PI * i / segments;
            float y = pipeRadius * cos(theta);
            float z = pipeRadius * sin(theta);
            
            glVertex3f(-pipeLength/2, y, z);
            glVertex3f(pipeLength/2, y, z);
        }
        glEnd();
        glDisable(GL_BLEND);
        
        // Render particles
        glPointSize(4.0f);
        glBegin(GL_POINTS);
        for (const auto& p : particles) {
            glColor3f(p.color.x, p.color.y, p.color.z);
            glVertex3f(p.position.x, p.position.y, p.position.z);
        }
        glEnd();
        glPointSize(1.0f);
    }
};

class ConvectionSimulator {
private:
    struct Particle {
        glm::vec3 position;
        glm::vec3 velocity;
        float temperature; // 0.0 (Cold) to 1.0 (Hot)
        float size;
    };
    
    std::vector<Particle> particles;
    float boxWidth = 4.0f;
    float boxHeight = 3.0f;
    float boxDepth = 1.0f;
    
    float bottomTemp = 1.0f; // Hot
    float topTemp = 0.0f;    // Cold
    float viscosity = 0.1f;
    float buoyancyStrength = 5.0f;
    
public:
    ConvectionSimulator() {
        resetSimulation();
    }
    
    void resetSimulation() {
        particles.clear();
        for (int i = 0; i < 500; ++i) {
            spawnParticle();
        }
    }
    
    void spawnParticle() {
        float x = ((float)rand() / RAND_MAX) * boxWidth - boxWidth/2;
        float y = ((float)rand() / RAND_MAX) * boxHeight - boxHeight/2;
        float z = ((float)rand() / RAND_MAX) * boxDepth - boxDepth/2;
        
        // Initial temp depends on height
        float normalizedHeight = (y + boxHeight/2) / boxHeight;
        float temp = 1.0f - normalizedHeight; 
        
        particles.push_back({
            glm::vec3(x, y, z),
            glm::vec3(0,0,0),
            temp,
            0.05f + ((float)rand() / RAND_MAX) * 0.05f
        });
    }
    
    void update(float deltaTime) {
        for (auto& p : particles) {
            // 1. Heat Transfer from boundaries
            float distToBottom = p.position.y - (-boxHeight/2);
            float distToTop = (boxHeight/2) - p.position.y;
            
            if (distToBottom < 0.2f) {
                p.temperature += (bottomTemp - p.temperature) * deltaTime * 2.0f;
            }
            if (distToTop < 0.2f) {
                p.temperature += (topTemp - p.temperature) * deltaTime * 2.0f;
            }
            
            // 2. Buoyancy Force (Hot rises, Cold sinks)
            // F_b = (T - T_avg) * g
            float buoyancy = (p.temperature - 0.5f) * buoyancyStrength;
            p.velocity.y += buoyancy * deltaTime;
            
            // 3. Viscosity / Drag
            p.velocity -= p.velocity * viscosity * deltaTime;
            
            // 4. Horizontal movement (simplified convection roll)
            // If hitting top/bottom, move sideways to create cycle
            if (distToTop < 0.5f || distToBottom < 0.5f) {
                // Push away from center to start rolls
                if (p.position.x > 0) p.velocity.x += 1.0f * deltaTime;
                else p.velocity.x -= 1.0f * deltaTime;
            }
            
            // 5. Update Position
            p.position += p.velocity * deltaTime;
            
            // 6. Boundaries
            if (p.position.y > boxHeight/2) {
                p.position.y = boxHeight/2;
                p.velocity.y *= -0.5f;
            }
            if (p.position.y < -boxHeight/2) {
                p.position.y = -boxHeight/2;
                p.velocity.y *= -0.5f;
            }
            if (p.position.x > boxWidth/2) {
                p.position.x = -boxWidth/2; // Wrap around x
            }
            if (p.position.x < -boxWidth/2) {
                p.position.x = boxWidth/2;
            }
            if (p.position.z > boxDepth/2) {
                p.position.z = boxDepth/2;
                p.velocity.z *= -0.5f;
            }
            if (p.position.z < -boxDepth/2) {
                p.position.z = -boxDepth/2;
                p.velocity.z *= -0.5f;
            }
        }
    }
    
    void renderUI() {
        ImGui::Text("Convection Simulation");
        ImGui::Separator();
        
        ImGui::SliderFloat("Bottom Temp (Hot)", &bottomTemp, 0.5f, 2.0f);
        ImGui::SliderFloat("Top Temp (Cold)", &topTemp, -1.0f, 0.5f);
        ImGui::SliderFloat("Viscosity", &viscosity, 0.01f, 2.0f);
        ImGui::SliderFloat("Buoyancy Strength", &buoyancyStrength, 1.0f, 20.0f);
        
        if (ImGui::Button("Reset Simulation")) {
            resetSimulation();
        }
        
        ImGui::Separator();
        ImGui::Text("Particles: %zu", particles.size());
        ImGui::Text("Red = Hot (Rising)");
        ImGui::Text("Blue = Cold (Sinking)");
    }
    
    void renderGraph() {
        ImGui::Text("See 3D view for convection cells");
    }
    
    void render3D() {
        // Debug print
        // std::cout << "Rendering Convection. Particles: " << particles.size() << std::endl;
        if (!particles.empty()) {
            // std::cout << "P0: " << particles[0].position.x << ", " << particles[0].position.y << ", " << particles[0].position.z << std::endl;
        }
        // Render Box
        glColor3f(1.0f, 1.0f, 1.0f);
        glLineWidth(2.0f);
        glBegin(GL_LINES);
        // Bottom rect
        glVertex3f(-boxWidth/2, -boxHeight/2, -boxDepth/2); glVertex3f(boxWidth/2, -boxHeight/2, -boxDepth/2);
        glVertex3f(boxWidth/2, -boxHeight/2, -boxDepth/2); glVertex3f(boxWidth/2, -boxHeight/2, boxDepth/2);
        glVertex3f(boxWidth/2, -boxHeight/2, boxDepth/2); glVertex3f(-boxWidth/2, -boxHeight/2, boxDepth/2);
        glVertex3f(-boxWidth/2, -boxHeight/2, boxDepth/2); glVertex3f(-boxWidth/2, -boxHeight/2, -boxDepth/2);
        // Top rect
        glVertex3f(-boxWidth/2, boxHeight/2, -boxDepth/2); glVertex3f(boxWidth/2, boxHeight/2, -boxDepth/2);
        glVertex3f(boxWidth/2, boxHeight/2, -boxDepth/2); glVertex3f(boxWidth/2, boxHeight/2, boxDepth/2);
        glVertex3f(boxWidth/2, boxHeight/2, boxDepth/2); glVertex3f(-boxWidth/2, boxHeight/2, boxDepth/2);
        glVertex3f(-boxWidth/2, boxHeight/2, boxDepth/2); glVertex3f(-boxWidth/2, boxHeight/2, -boxDepth/2);
        // Connectors
        glVertex3f(-boxWidth/2, -boxHeight/2, -boxDepth/2); glVertex3f(-boxWidth/2, boxHeight/2, -boxDepth/2);
        glVertex3f(boxWidth/2, -boxHeight/2, -boxDepth/2); glVertex3f(boxWidth/2, boxHeight/2, -boxDepth/2);
        glVertex3f(boxWidth/2, -boxHeight/2, boxDepth/2); glVertex3f(boxWidth/2, boxHeight/2, boxDepth/2);
        glVertex3f(-boxWidth/2, -boxHeight/2, boxDepth/2); glVertex3f(-boxWidth/2, boxHeight/2, boxDepth/2);
        glEnd();
        
        // Render Particles
        // Use Quads instead of Points for better visibility/compatibility
        glBegin(GL_QUADS);
        float s = 0.05f; // Size of particle
        for (const auto& p : particles) {
            // Color map: Blue (0.0) -> Red (1.0)
            float t = glm::clamp(p.temperature, 0.0f, 1.0f);
            glColor3f(t, 0.2f, 1.0f - t);
            
            glVertex3f(p.position.x - s, p.position.y - s, p.position.z);
            glVertex3f(p.position.x + s, p.position.y - s, p.position.z);
            glVertex3f(p.position.x + s, p.position.y + s, p.position.z);
            glVertex3f(p.position.x - s, p.position.y + s, p.position.z);
        }
        glEnd();
    }
};

class CoupledHeatSimulator {
private:
    float solidTemp;
    float fluidTemp;
    float solidMass = 1.0f;
    float fluidMass = 5.0f;
    float solidSpecificHeat = 450.0f; // Steel J/kgK
    float fluidSpecificHeat = 4186.0f; // Water J/kgK
    float heatTransferCoeff = 50.0f; // W/m^2K
    float surfaceArea = 0.06f; // 0.1m x 0.1m x 6 faces
    
    std::vector<glm::vec2> solidTempHistory;
    std::vector<glm::vec2> fluidTempHistory;
    float time = 0.0f;
    
public:
    CoupledHeatSimulator() {
        resetSimulation();
    }
    
    void resetSimulation() {
        solidTemp = 100.0f; // Hot
        fluidTemp = 20.0f;  // Cold
        time = 0.0f;
        solidTempHistory.clear();
        fluidTempHistory.clear();
    }
    
    void update(float deltaTime) {
        time += deltaTime;
        
        // Newton's Law of Cooling / Heating
        // Q_dot = h * A * (T_solid - T_fluid)
        // dT_solid/dt = -Q_dot / (m_solid * c_solid)
        // dT_fluid/dt = Q_dot / (m_fluid * c_fluid)
        
        float q_dot = heatTransferCoeff * surfaceArea * (solidTemp - fluidTemp);
        
        float deltaT_solid = -(q_dot / (solidMass * solidSpecificHeat)) * deltaTime;
        float deltaT_fluid = (q_dot / (fluidMass * fluidSpecificHeat)) * deltaTime;
        
        solidTemp += deltaT_solid;
        fluidTemp += deltaT_fluid;
        
        if (time > 0.1f) { // Sample every 0.1s
             solidTempHistory.push_back(glm::vec2(time, solidTemp));
             fluidTempHistory.push_back(glm::vec2(time, fluidTemp));
             if (solidTempHistory.size() > 1000) solidTempHistory.erase(solidTempHistory.begin());
             if (fluidTempHistory.size() > 1000) fluidTempHistory.erase(fluidTempHistory.begin());
        }
    }
    
    void renderUI() {
        ImGui::Text("Coupled Heat Transfer");
        ImGui::Separator();
        
        ImGui::SliderFloat("Solid Temp", &solidTemp, 0.0f, 200.0f);
        ImGui::SliderFloat("Fluid Temp", &fluidTemp, 0.0f, 100.0f);
        ImGui::SliderFloat("Heat Transfer Coeff (h)", &heatTransferCoeff, 1.0f, 500.0f);
        
        if (ImGui::Button("Reset Simulation")) {
            resetSimulation();
        }
        
        ImGui::Separator();
        ImGui::Text("Solid Temp: %.2f C", solidTemp);
        ImGui::Text("Fluid Temp: %.2f C", fluidTemp);
    }
    
    void renderGraph() {
        if (ImPlot::BeginPlot("Temperature vs Time")) {
            ImPlot::SetupAxes("Time (s)", "Temp (C)");
            if (!solidTempHistory.empty()) {
                ImPlot::PlotLine("Solid", &solidTempHistory[0].x, &solidTempHistory[0].y, solidTempHistory.size(), 0, 0, sizeof(glm::vec2));
                ImPlot::PlotLine("Fluid", &fluidTempHistory[0].x, &fluidTempHistory[0].y, fluidTempHistory.size(), 0, 0, sizeof(glm::vec2));
            }
            ImPlot::EndPlot();
        }
        
        if (ImGui::Button("Export Heat Data")) {
            std::vector<std::string> headers = {"Time (s)", "Solid Temp (C)", "Fluid Temp (C)"};
            std::vector<std::vector<float>> data;
            std::vector<float> timeCol, solidCol, fluidCol;
            for (size_t i = 0; i < solidTempHistory.size(); ++i) {
                timeCol.push_back(solidTempHistory[i].x);
                solidCol.push_back(solidTempHistory[i].y);
                if (i < fluidTempHistory.size()) {
                    fluidCol.push_back(fluidTempHistory[i].y);
                } else {
                    fluidCol.push_back(0.0f);
                }
            }
            data.push_back(timeCol);
            data.push_back(solidCol);
            data.push_back(fluidCol);
            exportToCSV("coupled_heat_data.csv", headers, data);
        }
        ImGui::SameLine();
        if (ImGui::Button("Save Image")) {
            requestPlotCapture("coupled_heat_plot.png");
        }
    }
    
    void render3D() {
        // Render Fluid Container (Transparent Box)
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        
        // Fluid Color based on temp (Blue -> Red)
        float t_fluid = glm::clamp((fluidTemp - 20.0f) / 80.0f, 0.0f, 1.0f);
        glColor4f(t_fluid, 0.2f, 1.0f - t_fluid, 0.3f);
        
        float size = 2.0f;
        // Draw Cube (Fluid)
        glBegin(GL_QUADS);
        // Front
        glVertex3f(-size, -size, size); glVertex3f(size, -size, size); glVertex3f(size, size, size); glVertex3f(-size, size, size);
        // Back
        glVertex3f(-size, -size, -size); glVertex3f(-size, size, -size); glVertex3f(size, size, -size); glVertex3f(size, -size, -size);
        // Top
        glVertex3f(-size, size, -size); glVertex3f(-size, size, size); glVertex3f(size, size, size); glVertex3f(size, size, -size);
        // Bottom
        glVertex3f(-size, -size, -size); glVertex3f(size, -size, -size); glVertex3f(size, -size, size); glVertex3f(-size, -size, size);
        // Right
        glVertex3f(size, -size, -size); glVertex3f(size, size, -size); glVertex3f(size, size, size); glVertex3f(size, -size, size);
        // Left
        glVertex3f(-size, -size, -size); glVertex3f(-size, -size, size); glVertex3f(-size, size, size); glVertex3f(-size, size, -size);
        glEnd();
        
        // Render Solid Block (Opaque)
        // Solid Color based on temp (Blue -> Red)
        float t_solid = glm::clamp((solidTemp - 20.0f) / 180.0f, 0.0f, 1.0f);
        glColor4f(t_solid, 0.0f, 1.0f - t_solid, 1.0f);
        
        float solidSize = 0.5f;
        glBegin(GL_QUADS);
        // Front
        glVertex3f(-solidSize, -solidSize, solidSize); glVertex3f(solidSize, -solidSize, solidSize); glVertex3f(solidSize, solidSize, solidSize); glVertex3f(-solidSize, solidSize, solidSize);
        // Back
        glVertex3f(-solidSize, -solidSize, -solidSize); glVertex3f(-solidSize, solidSize, -solidSize); glVertex3f(solidSize, solidSize, -solidSize); glVertex3f(solidSize, -solidSize, -solidSize);
        // Top
        glVertex3f(-solidSize, solidSize, -solidSize); glVertex3f(-solidSize, solidSize, solidSize); glVertex3f(solidSize, solidSize, solidSize); glVertex3f(solidSize, solidSize, -solidSize);
        // Bottom
        glVertex3f(-solidSize, -solidSize, -solidSize); glVertex3f(solidSize, -solidSize, -solidSize); glVertex3f(size, -solidSize, solidSize); glVertex3f(-solidSize, -solidSize, solidSize);
        // Right
        glVertex3f(solidSize, -solidSize, -solidSize); glVertex3f(solidSize, solidSize, -solidSize); glVertex3f(solidSize, solidSize, solidSize); glVertex3f(solidSize, -solidSize, solidSize);
        // Left
        glVertex3f(-solidSize, -solidSize, -solidSize); glVertex3f(-solidSize, -solidSize, solidSize); glVertex3f(-solidSize, solidSize, solidSize); glVertex3f(-solidSize, solidSize, -solidSize);
        glEnd();
        
        glDisable(GL_BLEND);
    }
};

class UniversalSimulator {
private:
    GLFWwindow* window;
    MaterialDatabase materialDB;
    
    // Simulators
    std::unique_ptr<EnhancedSolidMechanicsSimulator> solidSim;
    std::unique_ptr<FluidMechanicsSimulator> fluidSim;
    std::unique_ptr<GasDynamicsSimulator> gasSim;
    std::unique_ptr<PhaseChangeSimulator> phaseSim;
    std::unique_ptr<ThermalConductivitySimulator> thermalSim;
    std::unique_ptr<ProjectileMotionSimulator> projectileSim;
    std::unique_ptr<OrbitalMechanicsSimulator> orbitalSim;
    std::unique_ptr<PipeFlowSimulator> pipeSim;
    std::unique_ptr<ConvectionSimulator> convectionSim;
    std::unique_ptr<CoupledHeatSimulator> coupledSim;
    
    // UI State
    int currentScenario = 0; // 0=Solid, 1=Fluid, 2=Gas, 3=Phase
    bool scientificMode = true; // Toggle between scientific/educational
    std::string selectedMaterial = "Steel";
    std::string selectedFluid = "Water";
    std::string selectedGas = "Air";
    
    // OpenGL state
    glm::mat4 viewMatrix;
    glm::mat4 projectionMatrix;
    
    // Camera State
    float cameraDistance = 10.0f;
    float cameraPitch = 30.0f; // Degrees
    float cameraYaw = 45.0f;   // Degrees
    glm::vec3 cameraTarget = glm::vec3(0.0f, 0.0f, 0.0f);
    
public:
    UniversalSimulator() {
        // Initialize simulators
        solidSim = std::make_unique<EnhancedSolidMechanicsSimulator>();
        fluidSim = std::make_unique<FluidMechanicsSimulator>();
        gasSim = std::make_unique<GasDynamicsSimulator>();
        phaseSim = std::make_unique<PhaseChangeSimulator>();
        thermalSim = std::make_unique<ThermalConductivitySimulator>();
        projectileSim = std::make_unique<ProjectileMotionSimulator>();
        orbitalSim = std::make_unique<OrbitalMechanicsSimulator>();
        pipeSim = std::make_unique<PipeFlowSimulator>();
        convectionSim = std::make_unique<ConvectionSimulator>();
        coupledSim = std::make_unique<CoupledHeatSimulator>();
    }
    
    bool initialize() {
        // Initialize GLFW
        if (!glfwInit()) {
            std::cerr << "Failed to initialize GLFW" << std::endl;
            return false;
        }
        
        // Create window
#ifdef __APPLE__
        // macOS doesn't support OpenGL 3.3 Compatibility Profile.
        // We must use Legacy OpenGL (2.1) for immediate mode rendering (glBegin/glEnd).
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
        // No profile hint needed for 2.1
#else
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_COMPAT_PROFILE);
#endif
        
        window = glfwCreateWindow(1600, 1200, "Universal Material & Fluid Simulator", nullptr, nullptr);
        if (!window) {
            std::cerr << "Failed to create GLFW window" << std::endl;
            glfwTerminate();
            return false;
        }
        
        glfwMakeContextCurrent(window);
        glfwSwapInterval(1); // VSync
        
        // Initialize GLEW
        if (glewInit() != GLEW_OK) {
            std::cerr << "Failed to initialize GLEW" << std::endl;
            return false;
        }
        
        // Setup ImGui
        IMGUI_CHECKVERSION();
        ImGui::CreateContext();
        ImPlot::CreateContext();
        ImGuiIO& io = ImGui::GetIO(); (void)io;
        io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
        
        ImGui::StyleColorsDark();
        ImGui_ImplGlfw_InitForOpenGL(window, true);
        
#ifdef __APPLE__
        ImGui_ImplOpenGL3_Init("#version 120");
#else
        ImGui_ImplOpenGL3_Init("#version 330");
#endif
        
        // Setup OpenGL
        glEnable(GL_DEPTH_TEST);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        
        // Initialize camera
        setupCamera();
        
        // Initialize simulators with default materials
        solidSim->setMaterial(materialDB.getMaterial("Steel"));
        fluidSim->setObject(materialDB.getMaterial("Steel"));
        fluidSim->setFluid(materialDB.getMaterial("Water"));
        gasSim->setGas(materialDB.getMaterial("Air"));
        phaseSim->setMaterial(materialDB.getMaterial("Water"));
        thermalSim->setMaterial(materialDB.getMaterial("Steel"));
        projectileSim->setMaterial(materialDB.getMaterial("Steel"));
        pipeSim->setFluid(materialDB.getMaterial("Water"));
        
        return true;
    }
    
    void setupCamera() {
        projectionMatrix = glm::perspective(glm::radians(45.0f), 1600.0f / 1200.0f, 0.1f, 100.0f);
        viewMatrix = glm::lookAt(
            glm::vec3(3.0f, 2.0f, 3.0f), // Camera position
            glm::vec3(0.0f, 0.0f, 0.0f), // Look at origin
            glm::vec3(0.0f, 1.0f, 0.0f)  // Up vector
        );
    }
    
    void run() {
        while (!glfwWindowShouldClose(window)) {
            glfwPollEvents();
            
            update();
            render();
            
            glfwSwapBuffers(window);
        }
    }
    
    void update() {
        // Handle Camera Input
        ImGuiIO& io = ImGui::GetIO();
        if (!io.WantCaptureMouse) {
            // Zoom
            if (io.MouseWheel != 0.0f) {
                cameraDistance -= io.MouseWheel * 1.0f;
                if (cameraDistance < 1.0f) cameraDistance = 1.0f;
            }
            
            // Orbit (Right Mouse)
            if (ImGui::IsMouseDragging(ImGuiMouseButton_Right)) {
                ImVec2 delta = ImGui::GetMouseDragDelta(ImGuiMouseButton_Right);
                cameraYaw -= delta.x * 0.5f;
                cameraPitch += delta.y * 0.5f;
                
                // Clamp pitch to avoid gimbal lock/flipping
                if (cameraPitch > 89.0f) cameraPitch = 89.0f;
                if (cameraPitch < -89.0f) cameraPitch = -89.0f;
                
                ImGui::ResetMouseDragDelta(ImGuiMouseButton_Right);
            }
            
            // Pan (Middle Mouse)
            if (ImGui::IsMouseDragging(ImGuiMouseButton_Middle)) {
                ImVec2 delta = ImGui::GetMouseDragDelta(ImGuiMouseButton_Middle);
                
                // Calculate camera basis vectors
                float yawRad = glm::radians(cameraYaw);
                float pitchRad = glm::radians(cameraPitch);
                glm::vec3 forward;
                forward.x = cos(pitchRad) * sin(yawRad);
                forward.y = sin(pitchRad);
                forward.z = cos(pitchRad) * cos(yawRad);
                forward = glm::normalize(forward);
                
                glm::vec3 right = glm::normalize(glm::cross(forward, glm::vec3(0, 1, 0))); // Assuming Y-up world
                glm::vec3 up = glm::normalize(glm::cross(right, forward));
                
                float panSpeed = 0.005f * cameraDistance;
                cameraTarget -= right * delta.x * panSpeed;
                cameraTarget += up * delta.y * panSpeed;
                
                ImGui::ResetMouseDragDelta(ImGuiMouseButton_Middle);
            }
        }
        
        // Update Camera Matrix
        float yawRad = glm::radians(cameraYaw);
        float pitchRad = glm::radians(cameraPitch);
        
        glm::vec3 cameraPos;
        cameraPos.x = cameraTarget.x + cameraDistance * cos(pitchRad) * sin(yawRad);
        cameraPos.y = cameraTarget.y + cameraDistance * sin(pitchRad);
        cameraPos.z = cameraTarget.z + cameraDistance * cos(pitchRad) * cos(yawRad);
        
        viewMatrix = glm::lookAt(cameraPos, cameraTarget, glm::vec3(0.0f, 1.0f, 0.0f));
        
        // Update fluid simulation animation
        static float lastTime = 0.0f;
        float currentTime = glfwGetTime();
        float deltaTime = currentTime - lastTime;
        lastTime = currentTime;
        
        
        fluidSim->update(deltaTime);
        thermalSim->update(deltaTime);
        projectileSim->update(deltaTime);
        orbitalSim->update(deltaTime);
        pipeSim->update(deltaTime);
        convectionSim->update(deltaTime);
        coupledSim->update(deltaTime);
    }
    
    void render() {
        // Clear screen
        glClearColor(0.1f, 0.1f, 0.2f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        // Setup matrices
        glMatrixMode(GL_PROJECTION);
        glLoadMatrixf(glm::value_ptr(projectionMatrix));
        glMatrixMode(GL_MODELVIEW);
        glLoadMatrixf(glm::value_ptr(viewMatrix));
        
        // Render 3D scene
        render3DScene();
        
        // Render UI
        renderUI();
    }
    
    void render3DScene() {
        // Render grid
        renderGrid();
        
        // Render current simulation
        switch (currentScenario) {
            case 0: solidSim->render3D(); break;
            case 1: fluidSim->render3D(); break;
            case 2: gasSim->render3D(); break;
            case 3: phaseSim->render3D(); break;
            case 4: thermalSim->render3D(); break;
            case 5: projectileSim->render3D(); break;
            case 6: orbitalSim->render3D(); break;
            case 7: pipeSim->render3D(); break;
            case 8: convectionSim->render3D(); break;
            case 9: coupledSim->render3D(); break;
        }
        
        // Render coordinate axes
        renderAxes();
    }
    
    void renderGrid() {
        glColor3f(0.3f, 0.3f, 0.3f);
        glBegin(GL_LINES);
        for (int i = -5; i <= 5; ++i) {
            // X-axis lines
            glVertex3f(-5.0f, 0.0f, i);
            glVertex3f(5.0f, 0.0f, i);
            // Z-axis lines
            glVertex3f(i, 0.0f, -5.0f);
            glVertex3f(i, 0.0f, 5.0f);
        }
        glEnd();
    }
    
    void renderAxes() {
        glLineWidth(3.0f);
        glBegin(GL_LINES);
        // X-axis (Red)
        glColor3f(1.0f, 0.0f, 0.0f);
        glVertex3f(0.0f, 0.0f, 0.0f);
        glVertex3f(1.0f, 0.0f, 0.0f);
        // Y-axis (Green)
        glColor3f(0.0f, 1.0f, 0.0f);
        glVertex3f(0.0f, 0.0f, 0.0f);
        glVertex3f(0.0f, 1.0f, 0.0f);
        // Z-axis (Blue)
        glColor3f(0.0f, 0.0f, 1.0f);
        glVertex3f(0.0f, 0.0f, 0.0f);
        glVertex3f(0.0f, 0.0f, 1.0f);
        glEnd();
        glLineWidth(1.0f);
    }
    
    void renderUI() {
        // Start ImGui frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        
        // Main control panel
        renderMainPanel();
        
        // Scenario-specific panels
        renderScenarioPanel();
        
        // Graphs panel
        renderGraphsPanel();
        
        // Material selection panel
        renderMaterialPanel();
        
        // Info panel
        renderInfoPanel();
        
        // Render ImGui
        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
    }
    
    void renderMainPanel() {
        ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(300, 200), ImGuiCond_FirstUseEver);
        
        ImGui::Begin("Universal Material & Fluid Simulator");
        
        // Mode toggle
        if (ImGui::Checkbox("Scientific Mode", &scientificMode)) {
            // Toggle between scientific and educational modes
            // In educational mode, show simplified explanations
        }
        
        ImGui::Separator();
        ImGui::Text("Select Simulation Scenario:");
        
        const char* scenarios[] = {"Solid Mechanics", "Fluid Mechanics", "Gas Dynamics", "Phase Change", "Thermal Conductivity", "Projectile Motion", "Orbital Mechanics", "Laminar Pipe Flow", "Fluid Convection", "Coupled Heat Transfer"};
        if (ImGui::Combo("Scenario", &currentScenario, scenarios, 10)) {
            // Switch to selected scenario
            switchScenario(currentScenario);
        }
        
        ImGui::Separator();
        
        // Quick material selection
        std::vector<std::string> materials = materialDB.getMaterialNames();
        static size_t selectedMaterialIndex = 0;
        
        if (ImGui::BeginCombo("Quick Material", selectedMaterial.c_str())) {
            for (size_t i = 0; i < materials.size(); ++i) {
                bool isSelected = (selectedMaterialIndex == i);
                if (ImGui::Selectable(materials[i].c_str(), isSelected)) {
                    selectedMaterialIndex = i;
                    selectedMaterial = materials[i];
                    applyMaterialToCurrentScenario();
                }
                if (isSelected) {
                    ImGui::SetItemDefaultFocus();
                }
            }
            ImGui::EndCombo();
        }
        
        if (ImGui::Button("Reset Simulation")) {
            resetCurrentSimulation();
        }
        
        ImGui::Separator();
        
        // Camera Controls
        if (ImGui::CollapsingHeader("Camera Controls", ImGuiTreeNodeFlags_DefaultOpen)) {
            ImGui::Text("Zoom:");
            ImGui::SameLine();
            if (ImGui::Button("In [+]")) {
                cameraDistance -= 1.0f;
                if (cameraDistance < 1.0f) cameraDistance = 1.0f;
            }
            ImGui::SameLine();
            if (ImGui::Button("Out [-]")) {
                cameraDistance += 1.0f;
            }
            
            ImGui::Text("Orbit:");
            ImGui::SameLine();
            if (ImGui::Button("Left")) cameraYaw -= 5.0f;
            ImGui::SameLine();
            if (ImGui::Button("Right")) cameraYaw += 5.0f;
            
            ImGui::SameLine();
            if (ImGui::Button("Up")) {
                cameraPitch += 5.0f;
                if (cameraPitch > 89.0f) cameraPitch = 89.0f;
            }
            ImGui::SameLine();
            if (ImGui::Button("Down")) {
                cameraPitch -= 5.0f;
                if (cameraPitch < -89.0f) cameraPitch = -89.0f;
            }
            
            if (ImGui::Button("Reset View")) {
                cameraTarget = glm::vec3(0,0,0);
                cameraDistance = 10.0f;
                cameraPitch = 30.0f;
                cameraYaw = 45.0f;
            }
        }
        
        ImGui::End();
    }
    
    void renderScenarioPanel() {
        ImGui::SetNextWindowPos(ImVec2(320, 10), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(400, 500), ImGuiCond_FirstUseEver);
        
        ImGui::Begin("Scenario Controls");
        
        switch (currentScenario) {
            case 0:
                solidSim->renderUI();
                if (!scientificMode) {
                    ImGui::Separator();
                    ImGui::TextWrapped("Educational Mode: Watch how materials stretch or compress when forces are applied. The color shows stress levels!");
                }
                break;
            case 1:
                fluidSim->renderUI();
                if (!scientificMode) {
                    ImGui::Separator();
                    ImGui::TextWrapped("Educational Mode: Drop objects in fluids to see if they float or sink based on density!");
                }
                break;
            case 2:
                gasSim->renderUI();
                if (!scientificMode) {
                    ImGui::Separator();
                    ImGui::TextWrapped("Educational Mode: Heat gases to see particles move faster and containers expand!");
                }
                break;
            case 3:
                phaseSim->renderUI();
                if (!scientificMode) {
                    ImGui::Separator();
                    ImGui::TextWrapped("Educational Mode: Heat materials to watch them melt from solid to liquid to gas!");
                }
                break;
            case 4:
                thermalSim->renderUI();
                if (!scientificMode) {
                    ImGui::Separator();
                    ImGui::TextWrapped("Educational Mode: Watch how heat travels through materials!");
                }
                break;
            case 5:
                projectileSim->renderUI();
                if (!scientificMode) {
                    ImGui::Separator();
                    ImGui::TextWrapped("Educational Mode: Launch things and see how gravity and air resistance affect them!");
                }
                break;
            case 6:
                orbitalSim->renderUI();
                if (!scientificMode) {
                    ImGui::Separator();
                    ImGui::TextWrapped("Educational Mode: Create your own solar system! Gravity pulls planets together.");
                }
                break;
            case 7:
                pipeSim->renderUI();
                if (!scientificMode) {
                    ImGui::Separator();
                    ImGui::TextWrapped("Educational Mode: Watch fluid flow through a pipe. It moves fastest in the middle!");
                }
                break;
            case 8:
                convectionSim->renderUI();
                if (!scientificMode) {
                    ImGui::Separator();
                    ImGui::TextWrapped("Educational Mode: Hot fluid rises (Red), Cold fluid sinks (Blue). This creates a cycle!");
                }
                break;
            case 9:
                coupledSim->renderUI();
                if (!scientificMode) {
                    ImGui::Separator();
                    ImGui::TextWrapped("Educational Mode: The hot block cools down by giving heat to the water, until they are both the same temperature.");
                }
                break;
        }
        
        ImGui::End();
    }
    
    void renderGraphsPanel() {
        ImGui::SetNextWindowPos(ImVec2(730, 10), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(500, 400), ImGuiCond_FirstUseEver);
        
        ImGui::Begin("Data Visualization");
        
        switch (currentScenario) {
            case 0:
                solidSim->renderGraph();
                break;
            case 2:
                gasSim->renderGraphs();
                break;
            case 3:
                phaseSim->renderGraph();
                break;
            case 4:
                thermalSim->renderGraph();
                break;
            case 5:
                projectileSim->renderGraph();
                break;
            case 6:
                orbitalSim->renderGraph();
                break;
            case 7:
                pipeSim->renderGraph();
                break;
            case 8:
                convectionSim->renderGraph();
                break;
            case 9:
                coupledSim->renderGraph();
                break;
            default:
                ImGui::Text("No graphs available for this scenario");
                break;
        }
        
        ImGui::End();
    }
    
    void renderMaterialPanel() {
        ImGui::SetNextWindowPos(ImVec2(10, 220), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(300, 300), ImGuiCond_FirstUseEver);
        
        ImGui::Begin("Material Properties");
        
        Material currentMaterial = materialDB.getMaterial(selectedMaterial);
        
        ImGui::Text("Material: %s", currentMaterial.name.c_str());
        ImGui::Separator();
        
        if (scientificMode) {
            ImGui::Text("Density: %.1f kg/m³", currentMaterial.density);
            ImGui::Text("Young's Modulus: %.1e Pa", currentMaterial.youngsModulus);
            ImGui::Text("Viscosity: %.2e Pa·s", currentMaterial.viscosity);
            ImGui::Text("Thermal Conductivity: %.2f W/(m·K)", currentMaterial.thermalConductivity);
            ImGui::Text("Specific Heat: %.0f J/(kg·K)", currentMaterial.specificHeatCapacity);
            ImGui::Text("Melting Point: %.1f K", currentMaterial.meltingPoint);
            ImGui::Text("Boiling Point: %.1f K", currentMaterial.boilingPoint);
        } else {
            // Educational explanations
            ImGui::TextWrapped("Density: How heavy the material is");
            ImGui::TextWrapped("Elasticity: How much it stretches");
            ImGui::TextWrapped("Viscosity: How thick fluids are");
            ImGui::TextWrapped("Heat Conduction: How fast heat travels through it");
        }
        
        ImGui::Separator();
        
        // Custom material editor
        if (ImGui::CollapsingHeader("Custom Material Editor")) {
            static float customDensity = 1000.0f;
            static float customYoungsModulus = 1e9f;
            static float customViscosity = 0.001f;
            
            ImGui::SliderFloat("Custom Density", &customDensity, 100.0f, 20000.0f, "%.1f kg/m³");
            ImGui::SliderFloat("Custom Young's Modulus", &customYoungsModulus, 1e6f, 1e12f, "%.1e Pa", ImGuiSliderFlags_Logarithmic);
            ImGui::SliderFloat("Custom Viscosity", &customViscosity, 1e-6f, 100.0f, "%.3e Pa·s", ImGuiSliderFlags_Logarithmic);
            
            if (ImGui::Button("Apply Custom Material")) {
                // Create and apply custom material
                Material custom = currentMaterial;
                custom.name = "Custom";
                custom.density = customDensity;
                custom.youngsModulus = customYoungsModulus;
                custom.viscosity = customViscosity;
                
                applyCustomMaterial(custom);
            }
        }
        
        ImGui::End();
    }
    
    void renderInfoPanel() {
        ImGui::SetNextWindowPos(ImVec2(1240, 10), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(350, 600), ImGuiCond_FirstUseEver);
        
        ImGui::Begin("Physics Information");
        
        switch (currentScenario) {
            case 0:
                ImGui::Text("SOLID MECHANICS");
                ImGui::Separator();
                if (scientificMode) {
                    ImGui::TextWrapped("Hooke's Law: σ = E × ε");
                    ImGui::TextWrapped("where σ = stress, E = Young's modulus, ε = strain");
                    ImGui::TextWrapped("Stress = Force / Area");
                    ImGui::TextWrapped("Strain = Change in length / Original length");
                } else {
                    ImGui::TextWrapped("When you pull or push materials, they stretch or compress. Stiffer materials (like steel) don't deform much, while flexible materials (like rubber) deform a lot!");
                }
                break;
                
            case 4:
                ImGui::Text("HEAT TRANSFER");
                ImGui::Separator();
                if (scientificMode) {
                     ImGui::TextWrapped("Fourier's Law: q = -k∇T");
                     ImGui::TextWrapped("Heat flows from hot to cold regions.");
                } else {
                    ImGui::TextWrapped("Heat moves from hot things to cold things. Some materials (like metal) let heat move fast, others (like wood) are slow!");
                }
                break;
                
            case 5:
                ImGui::Text("PROJECTILE MOTION");
                ImGui::Separator();
                if (scientificMode) {
                    ImGui::TextWrapped("y = y0 + vy*t - 0.5*g*t^2");
                    ImGui::TextWrapped("x = x0 + vx*t");
                } else {
                    ImGui::TextWrapped("When you throw something, gravity pulls it down. If there's air, it also slows down!");
                }
                break;
                
            case 6:
                ImGui::Text("ORBITAL MECHANICS");
                ImGui::Separator();
                if (scientificMode) {
                    ImGui::TextWrapped("F = G * m1 * m2 / r^2");
                    ImGui::TextWrapped("Gravity is a force that attracts all objects with mass.");
                } else {
                    ImGui::TextWrapped("Planets orbit stars because of gravity. It's like they are constantly falling but missing the ground!");
                }
                break;
                
            case 7:
                ImGui::Text("LAMINAR PIPE FLOW");
                ImGui::Separator();
                if (scientificMode) {
                    ImGui::TextWrapped("Poiseuille's Law: Q = πR^4ΔP / 8μL");
                    ImGui::TextWrapped("Velocity profile is parabolic.");
                } else {
                    ImGui::TextWrapped("Fluids stick to the walls of the pipe, so they move slower there and faster in the middle!");
                }
                break;
                
            case 8:
                ImGui::Text("FLUID CONVECTION");
                ImGui::Separator();
                if (scientificMode) {
                    ImGui::TextWrapped("Rayleigh-Bénard Convection");
                    ImGui::TextWrapped("Buoyancy driven flow due to temperature gradient.");
                } else {
                    ImGui::TextWrapped("When you heat a fluid from the bottom, it gets lighter and rises. When it cools at the top, it gets heavier and sinks. This creates a loop!");
                }
                break;
                
            case 9:
                ImGui::Text("COUPLED HEAT TRANSFER");
                ImGui::Separator();
                if (scientificMode) {
                    ImGui::TextWrapped("Conjugate Heat Transfer");
                    ImGui::TextWrapped("Newton's Law of Cooling: Q = h*A*(T_obj - T_fluid)");
                } else {
                    ImGui::TextWrapped("Heat flows from the hot object to the cold fluid. The bigger the temperature difference, the faster it flows!");
                }
                break;

            case 1:
                ImGui::Text("FLUID MECHANICS");
                ImGui::Separator();
                if (scientificMode) {
                    ImGui::TextWrapped("Buoyant Force = ρ_fluid × V × g");
                    ImGui::TextWrapped("Object floats if ρ_object < ρ_fluid");
                    ImGui::TextWrapped("Archimedes' Principle governs buoyancy");
                } else {
                    ImGui::TextWrapped("Objects float when they're less dense than the fluid they're in. That's why ice floats on water but rocks sink!");
                }
                break;
                
            case 2:
                ImGui::Text("GAS DYNAMICS");
                ImGui::Separator();
                if (scientificMode) {
                    ImGui::TextWrapped("Ideal Gas Law: PV = nRT");
                    ImGui::TextWrapped("Charles's Law: V/T = constant (at constant P)");
                    ImGui::TextWrapped("Boyle's Law: PV = constant (at constant T)");
                } else {
                    ImGui::TextWrapped("When you heat gases, the tiny particles move faster and need more space, so the gas expands or pressure increases!");
                }
                break;
                
            case 3:
                ImGui::Text("PHASE CHANGES");
                ImGui::Separator();
                if (scientificMode) {
                    ImGui::TextWrapped("Phase transitions occur at specific temperatures");
                    ImGui::TextWrapped("Heat of fusion (solid→liquid)");
                    ImGui::TextWrapped("Heat of vaporization (liquid→gas)");
                } else {
                    ImGui::TextWrapped("Matter changes form when heated or cooled: ice melts to water, water boils to steam. Each material has its own melting and boiling points!");
                }
                break;
        }
        
        ImGui::Separator();
        ImGui::Text("Controls:");
        ImGui::TextWrapped("• Use sliders to adjust parameters");
        ImGui::TextWrapped("• Watch real-time 3D visualization");
        ImGui::TextWrapped("• View graphs and data");
        ImGui::TextWrapped("• Toggle between Scientific/Educational modes");
        
        ImGui::End();
    }
    
    void switchScenario(int /*scenario*/) {
        // Reset camera
        cameraTarget = glm::vec3(0,0,0);
        cameraDistance = 10.0f;
        cameraPitch = 30.0f;
        cameraYaw = 45.0f;
        setupCamera();
        
        // Apply appropriate materials to new scenario
        applyMaterialToCurrentScenario();
    }
    
    void applyMaterialToCurrentScenario() {
        Material mat = materialDB.getMaterial(selectedMaterial);
        
        switch (currentScenario) {
            case 0:
                solidSim->setMaterial(mat);
                break;
            case 1:
                if (mat.viscosity > 0) {
                    fluidSim->setFluid(mat);
                } else {
                    fluidSim->setObject(mat);
                }
                break;
            case 2:
                gasSim->setGas(mat);
                break;
            case 3:
                phaseSim->setMaterial(mat);
                break;
            case 7:
                pipeSim->setFluid(mat);
                break;
        }
    }
    
    void applyCustomMaterial(const Material& custom) {
        switch (currentScenario) {
            case 0:
                solidSim->setMaterial(custom);
                break;
            case 1:
                fluidSim->setObject(custom);
                break;
            case 2:
                gasSim->setGas(custom);
                break;
            case 3:
                phaseSim->setMaterial(custom);
                break;
            case 7:
                pipeSim->setFluid(custom);
                break;
        }
    }
    
    void resetCurrentSimulation() {
        switch (currentScenario) {
            case 0:
                solidSim->setForce(0.0);
                break;
            case 1:
                fluidSim->calculate();
                break;
            case 2:
                gasSim->setTemperature(300.0);
                break;
            case 3:
                phaseSim->setTemperature(273.15);
                break;
            case 7:
                pipeSim->resetParticles();
                break;
            case 8:
                convectionSim->resetSimulation();
                break;
            case 9:
                coupledSim->resetSimulation();
                break;
        }
    }
    
    void cleanup() {
        // Cleanup ImGui
        ImGui_ImplOpenGL3_Shutdown();
        ImGui_ImplGlfw_Shutdown();
        ImPlot::DestroyContext();
        ImGui::DestroyContext();
        
        // Cleanup GLFW
        glfwDestroyWindow(window);
        glfwTerminate();
    }
    
    ~UniversalSimulator() {
        cleanup();
    }
};

// ============================================================================
// MAIN FUNCTION
// ============================================================================

int main() {
    std::cout << "Universal Material & Fluid Simulator" << std::endl;
    std::cout << "====================================" << std::endl;
    
    UniversalSimulator app;
    
    if (!app.initialize()) {
        std::cerr << "Failed to initialize application" << std::endl;
        return -1;
    }
    
    std::cout << "Simulator initialized successfully!" << std::endl;
    std::cout << "Controls:" << std::endl;
    std::cout << "- Select scenario from dropdown" << std::endl;
    std::cout << "- Adjust material properties with sliders" << std::endl;
    std::cout << "- Toggle Scientific/Educational mode" << std::endl;
    std::cout << "- Camera: Right Drag to Orbit, Middle Drag to Pan, Scroll to Zoom" << std::endl;
    
    app.run();
    
    return 0;
}