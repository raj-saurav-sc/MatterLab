#pragma once

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
#include "../stb_image_write.h"
#include "vendor/json.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <string>

#include "../physics/Material.h"
#include "../physics/ISolver.h"
#include "../physics/solvers/SolidMechanicsSolver.h"
#include "../physics/solvers/FluidMechanicsSolver.h"
#include "../physics/solvers/GasDynamicsSimulator.h"
#include "../physics/solvers/PhaseChangeSimulator.h"
#include "../physics/solvers/ThermalConductivitySimulator.h"
#include "../physics/solvers/ProjectileMotionSimulator.h"
#include "../physics/solvers/OrbitalMechanicsSimulator.h"
#include "../physics/solvers/PipeFlowSimulator.h"
#include "../physics/solvers/FluidConvectionSimulator.h"
#include "../physics/solvers/CoupledHeatSimulator.h"
#include "../physics/solvers/FEMSolver.h"
#include "../physics/solvers/CFDSolver.h"

class Application {
private:
    GLFWwindow* window;
    int width = 1280;
    int height = 720;
    
    // Camera
    glm::vec3 cameraTarget{0.0f, 0.0f, 0.0f};
    float cameraDistance = 5.0f;
    float cameraPitch = 30.0f;
    float cameraYaw = 45.0f;
    glm::mat4 viewMatrix;
    glm::mat4 projectionMatrix;
    
    // State
    MaterialDatabase materialDB;
    std::string selectedMaterial = "Steel";
    int currentScenario = 0;
    bool scientificMode = true;
    
    // UI State
    bool showControls = true;
    bool showHelp = false;
    bool showMaterialEditor = false;
    
    // Material Editor
    char materialNameBuffer[64] = "Custom Material";
    Material customMaterial;
    float youngsModulusGPa = 200.0f;
    float yieldStrengthMPa = 250.0f;
    float ultimateStrengthMPa = 400.0f;
    int materialTypePreset = 0; // 0=Metal, 1=Fluid, 2=Gas
    // Temp floats for ImGui sliders (Material uses double)
    float tempDensity = 7800.0f;
    float tempViscosity = 0.0f;
    float tempThermalCond = 50.0f;
    float tempSpecificHeat = 420.0f;
    float tempMeltingPoint = 1811.0f;
    float tempBoilingPoint = 3273.0f;
    
    // Solvers
    std::vector<std::shared_ptr<ISolver>> solvers;
    
public:
    Application() = default;
    
    bool initialize() {
        // Initialize GLFW
        if (!glfwInit()) return false;
        
        // GL 3.0 + GLSL 130
        const char* glsl_version = "#version 130";
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
        
        // Create window
        window = glfwCreateWindow(width, height, "Universal Material & Fluid Simulator (Phase 2)", NULL, NULL);
        if (!window) return false;
        
        glfwMakeContextCurrent(window);
        glfwSwapInterval(1); // Enable vsync
        
        // Initialize GLEW
        if (glewInit() != GLEW_OK) return false;
        
        // Initialize ImGui
        IMGUI_CHECKVERSION();
        ImGui::CreateContext();
        ImPlot::CreateContext();
        ImGuiIO& io = ImGui::GetIO(); (void)io;
        
        ImGui::StyleColorsDark();
        
        ImGui_ImplGlfw_InitForOpenGL(window, true);
        ImGui_ImplOpenGL3_Init(glsl_version);
        
        // Initialize Solvers
        solvers.push_back(std::make_shared<SolidMechanicsSolver>());
        solvers.push_back(std::make_shared<FluidMechanicsSolver>());
        solvers.push_back(std::make_shared<GasDynamicsSimulator>());
        solvers.push_back(std::make_shared<PhaseChangeSimulator>());
        solvers.push_back(std::make_shared<ThermalConductivitySimulator>());
        solvers.push_back(std::make_shared<ProjectileMotionSimulator>());
        solvers.push_back(std::make_shared<OrbitalMechanicsSimulator>());
        solvers.push_back(std::make_shared<PipeFlowSimulator>());
        solvers.push_back(std::make_shared<FluidConvectionSimulator>());
        solvers.push_back(std::make_shared<CoupledHeatSimulator>());
        solvers.push_back(std::make_shared<FEMSolver>());
        solvers.push_back(std::make_shared<CFDSolver>());
        
        for (auto& solver : solvers) {
            solver->initialize();
            solver->setMaterial(materialDB.getMaterial("Steel"));
        }
        
        setupCamera();
        
        return true;
    }
    
    void run() {
        while (!glfwWindowShouldClose(window)) {
            glfwPollEvents();
            
            update();
            render();
            
            // Handle Image Capture
            if (captureNextFrame) {
                captureNextFrame = false;
                
                int x = (int)captureMin.x;
                int y = (int)captureMin.y;
                int w = (int)(captureMax.x - captureMin.x);
                int h = (int)(captureMax.y - captureMin.y);
                
                // Handle high-DPI scaling
                int winW, winH, fbW, fbH;
                glfwGetWindowSize(window, &winW, &winH);
                glfwGetFramebufferSize(window, &fbW, &fbH);
                float scaleX = (float)fbW / winW;
                float scaleY = (float)fbH / winH;
                
                x = (int)(x * scaleX);
                y = (int)(y * scaleY);
                w = (int)(w * scaleX);
                h = (int)(h * scaleY);
                
                int glY = fbH - y - h;
                
                if (w > 0 && h > 0) {
                    glPixelStorei(GL_PACK_ALIGNMENT, 1);
                    std::vector<unsigned char> pixels(w * h * 3);
                    glReadPixels(x, glY, w, h, GL_RGB, GL_UNSIGNED_BYTE, pixels.data());
                    
                    std::vector<unsigned char> flipped(w * h * 3);
                    for (int row = 0; row < h; ++row) {
                        memcpy(&flipped[row * w * 3], &pixels[(h - 1 - row) * w * 3], w * 3);
                    }
                    
                    if (stbi_write_png(captureFilename.c_str(), w, h, 3, flipped.data(), w * 3)) {
                        std::cout << "Saved image to " << captureFilename << std::endl;
                    } else {
                        std::cerr << "Failed to save image to " << captureFilename << std::endl;
                    }
                    glPixelStorei(GL_PACK_ALIGNMENT, 4);
                }
            }
            
            glfwSwapBuffers(window);
        }
    }
    
    ~Application() {
        ImGui_ImplOpenGL3_Shutdown();
        ImGui_ImplGlfw_Shutdown();
        ImPlot::DestroyContext();
        ImGui::DestroyContext();
        
        glfwDestroyWindow(window);
        glfwTerminate();
    }
    
private:
    void update() {
        handleInput();
        
        float deltaTime = ImGui::GetIO().DeltaTime;
        if (currentScenario >= 0 && currentScenario < static_cast<int>(solvers.size())) {
            solvers[currentScenario]->update(deltaTime);
        }
    }
    
    void handleInput() {
        ImGuiIO& io = ImGui::GetIO();
        if (!io.WantCaptureMouse) {
            // Zoom (Mouse Wheel)
            if (io.MouseWheel != 0.0f) {
                cameraDistance -= io.MouseWheel * 1.0f;
                if (cameraDistance < 1.0f) cameraDistance = 1.0f;
                if (cameraDistance > 50.0f) cameraDistance = 50.0f;
            }
            
            // Orbit (Right Mouse)
            if (ImGui::IsMouseDragging(ImGuiMouseButton_Right)) {
                ImVec2 delta = ImGui::GetMouseDragDelta(ImGuiMouseButton_Right);
                cameraYaw -= delta.x * 0.5f;
                cameraPitch += delta.y * 0.5f;
                
                if (cameraPitch > 89.0f) cameraPitch = 89.0f;
                if (cameraPitch < -89.0f) cameraPitch = -89.0f;
                
                ImGui::ResetMouseDragDelta(ImGuiMouseButton_Right);
            }
            
            // Pan (Left Mouse or Middle Mouse)
            if (ImGui::IsMouseDragging(ImGuiMouseButton_Left) || ImGui::IsMouseDragging(ImGuiMouseButton_Middle)) {
                int button = ImGui::IsMouseDragging(ImGuiMouseButton_Left) ? ImGuiMouseButton_Left : ImGuiMouseButton_Middle;
                ImVec2 delta = ImGui::GetMouseDragDelta(button);
                
                // Calculate pan direction based on camera orientation
                float yawRad = glm::radians(cameraYaw);
                float pitchRad = glm::radians(cameraPitch);
                
                glm::vec3 right(cos(yawRad), 0, -sin(yawRad));
                glm::vec3 up = glm::normalize(glm::cross(right, glm::vec3(0, 1, 0)));
                
                float panSpeed = cameraDistance * 0.002f;
                cameraTarget -= right * delta.x * panSpeed;
                cameraTarget += up * delta.y * panSpeed;
                
                ImGui::ResetMouseDragDelta(button);
            }
        }
        
        // Keyboard Controls
        if (!io.WantCaptureKeyboard) {
            float panSpeed = 0.1f;
            float zoomSpeed = 0.5f;
            
            // Pan with Arrow Keys
            if (ImGui::IsKeyPressed(ImGuiKey_LeftArrow)) {
                float yawRad = glm::radians(cameraYaw);
                glm::vec3 right(cos(yawRad), 0, -sin(yawRad));
                cameraTarget -= right * panSpeed;
            }
            if (ImGui::IsKeyPressed(ImGuiKey_RightArrow)) {
                float yawRad = glm::radians(cameraYaw);
                glm::vec3 right(cos(yawRad), 0, -sin(yawRad));
                cameraTarget += right * panSpeed;
            }
            if (ImGui::IsKeyPressed(ImGuiKey_UpArrow)) {
                cameraTarget.y += panSpeed;
            }
            if (ImGui::IsKeyPressed(ImGuiKey_DownArrow)) {
                cameraTarget.y -= panSpeed;
            }
            
            // Zoom with W/S or +/-
            if (ImGui::IsKeyPressed(ImGuiKey_W) || ImGui::IsKeyPressed(ImGuiKey_Equal)) {
                cameraDistance -= zoomSpeed;
                if (cameraDistance < 1.0f) cameraDistance = 1.0f;
            }
            if (ImGui::IsKeyPressed(ImGuiKey_S) || ImGui::IsKeyPressed(ImGuiKey_Minus)) {
                cameraDistance += zoomSpeed;
                if (cameraDistance > 50.0f) cameraDistance = 50.0f;
            }
            
            // Reset with R
            if (ImGui::IsKeyPressed(ImGuiKey_R)) {
                resetCamera();
            }
        }
        
        setupCamera();
    }
    
    void resetCamera() {
        cameraTarget = glm::vec3(0.0f, 0.0f, 0.0f);
        cameraDistance = 5.0f;
        cameraPitch = 30.0f;
        cameraYaw = 45.0f;
    }
    
    void setupCamera() {
        int display_w, display_h;
        glfwGetFramebufferSize(window, &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);
        float aspect = (float)display_w / (float)display_h;
        
        projectionMatrix = glm::perspective(glm::radians(45.0f), aspect, 0.1f, 100.0f);
        
        float yawRad = glm::radians(cameraYaw);
        float pitchRad = glm::radians(cameraPitch);
        
        glm::vec3 offset;
        offset.x = cameraDistance * cos(pitchRad) * sin(yawRad);
        offset.y = cameraDistance * sin(pitchRad);
        offset.z = cameraDistance * cos(pitchRad) * cos(yawRad);
        
        glm::vec3 cameraPos = cameraTarget + offset;
        viewMatrix = glm::lookAt(cameraPos, cameraTarget, glm::vec3(0, 1, 0));
    }
    
    void render() {
        glClearColor(0.1f, 0.1f, 0.2f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glEnable(GL_DEPTH_TEST);
        
        // 3D Rendering
        glMatrixMode(GL_PROJECTION);
        glLoadMatrixf(glm::value_ptr(projectionMatrix));
        glMatrixMode(GL_MODELVIEW);
        glLoadMatrixf(glm::value_ptr(viewMatrix));
        
        renderGrid();
        renderAxes();
        
        if (currentScenario >= 0 && currentScenario < static_cast<int>(solvers.size())) {
            solvers[currentScenario]->render3D();
        }
        
        // UI Rendering
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        
        renderUI();
        
        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
    }
    
    void renderUI() {
        // Main Control Panel
        ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(300, 600), ImGuiCond_FirstUseEver);
        ImGui::Begin("Control Panel");
        
        // Project Management
        if (ImGui::Button("Save Project")) {
            saveProject("project.json");
        }
        ImGui::SameLine();
        if (ImGui::Button("Load Project")) {
            loadProject("project.json");
        }
        ImGui::Separator();
        
        // Scenario Selection
        const char* scenarios[] = {
            "Solid Mechanics", "Fluid Mechanics", "Gas Dynamics", "Phase Changes",
            "Thermal Conductivity", "Projectile Motion", "Orbital Mechanics",
            "Pipe Flow", "Fluid Convection", "Coupled Heat Transfer",
            "Advanced FEM (Truss)", "Advanced CFD (Navier-Stokes)"
        };
        if (ImGui::Combo("Scenario", &currentScenario, scenarios, IM_ARRAYSIZE(scenarios))) {
            // Reset camera or state if needed
        }
        
        // Material Selector
        ImGui::Text("Material:");
        if (ImGui::BeginCombo("##material", selectedMaterial.c_str())) {
            for (const auto& name : materialDB.getMaterialNames()) {
                bool isSelected = (selectedMaterial == name);
                if (ImGui::Selectable(name.c_str(), isSelected)) {
                    selectedMaterial = name;
                    if (currentScenario >= 0 && currentScenario < static_cast<int>(solvers.size())) {
                        solvers[currentScenario]->setMaterial(materialDB.getMaterial(selectedMaterial));
                    }
                }
                if (isSelected) ImGui::SetItemDefaultFocus();
            }
            ImGui::EndCombo();
        }
        
        if (ImGui::Button("Material Editor")) {
            showMaterialEditor = !showMaterialEditor;
        }
        
        ImGui::Separator();
        
        ImGui::Checkbox("Scientific Mode", &scientificMode);
        
        ImGui::Separator();
        
        // Camera Controls
        ImGui::Text("Camera Controls");
        if (ImGui::Button("Zoom In")) {
            cameraDistance -= 1.0f;
            if (cameraDistance < 1.0f) cameraDistance = 1.0f;
        }
        ImGui::SameLine();
        if (ImGui::Button("Zoom Out")) {
            cameraDistance += 1.0f;
            if (cameraDistance > 50.0f) cameraDistance = 50.0f;
        }
        if (ImGui::Button("Reset View")) {
            resetCamera();
        }
        
        ImGui::Text("Zoom: %.1f", cameraDistance);
        
        if (ImGui::CollapsingHeader("Controls Help")) {
            ImGui::TextWrapped("Mouse:");
            ImGui::BulletText("Left Drag: Pan");
            ImGui::BulletText("Right Drag: Rotate");
            ImGui::BulletText("Wheel: Zoom");
            ImGui::TextWrapped("Keyboard:");
            ImGui::BulletText("Arrows: Pan");
            ImGui::BulletText("W/S or +/-: Zoom");
            ImGui::BulletText("R: Reset View");
        }
        
        ImGui::Separator();
        
        // Solver-specific UI
        if (currentScenario >= 0 && currentScenario < static_cast<int>(solvers.size())) {
            solvers[currentScenario]->renderUI();
        }
        
        ImGui::End();
        
        // Graph Panel
        ImGui::SetNextWindowPos(ImVec2(10, 620), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(600, 300), ImGuiCond_FirstUseEver);
        ImGui::Begin("Data Visualization");
        if (currentScenario >= 0 && currentScenario < static_cast<int>(solvers.size())) {
            solvers[currentScenario]->renderGraph();
        }
        ImGui::End();
        
        // Material Editor Window
        renderMaterialEditor();
    }
    
    void renderMaterialEditor() {
        if (!showMaterialEditor) return;
        
        ImGui::SetNextWindowSize(ImVec2(500, 600), ImGuiCond_FirstUseEver);
        ImGui::Begin("Material Editor", &showMaterialEditor);
        
        // Material Name
        ImGui::InputText("Material Name", materialNameBuffer, 64);
        
        ImGui::Separator();
        
        // Preset Templates
        const char* presets[] = {"Metal", "Fluid", "Gas"};
        if (ImGui::Combo("Preset Template", &materialTypePreset, presets, 3)) {
            // Apply preset values
            if (materialTypePreset == 0) { // Metal
                customMaterial.density = 7800;
                youngsModulusGPa = 200.0f;
                customMaterial.viscosity = 0.0;
                customMaterial.thermalConductivity = 50.0;
                customMaterial.specificHeatCapacity = 420;
                yieldStrengthMPa = 250.0f;
                ultimateStrengthMPa = 400.0f;
                customMaterial.color = glm::vec3(0.7f, 0.7f, 0.8f);
            } else if (materialTypePreset == 1) { // Fluid
                customMaterial.density = 1000;
                youngsModulusGPa = 0.0f;
                customMaterial.viscosity = 0.001;
                customMaterial.thermalConductivity = 0.6;
                customMaterial.specificHeatCapacity = 4186;
                yieldStrengthMPa = 0.0f;
                ultimateStrengthMPa = 0.0f;
                customMaterial.color = glm::vec3(0.2f, 0.4f, 0.8f);
            } else { // Gas
                customMaterial.density = 1.225;
                youngsModulusGPa = 0.0f;
                customMaterial.viscosity = 1.8e-5;
                customMaterial.thermalConductivity = 0.026;
                customMaterial.specificHeatCapacity = 1005;
                yieldStrengthMPa = 0.0f;
                ultimateStrengthMPa = 0.0f;
                customMaterial.color = glm::vec3(0.7f, 0.8f, 0.9f);
            }
        }
        
        ImGui::Separator();
        ImGui::Text("Physical Properties");
        
        ImGui::SliderFloat("Density (kg/m³)", &tempDensity, 0.1, 20000, "%.1f", ImGuiSliderFlags_Logarithmic);
        ImGui::SliderFloat("Young's Modulus (GPa)", &youngsModulusGPa, 0.001, 500, "%.3f", ImGuiSliderFlags_Logarithmic);
        ImGui::SliderFloat("Viscosity (Pa·s)", &tempViscosity, 0.0, 10.0, "%.6f", ImGuiSliderFlags_Logarithmic);
        
        ImGui::Separator();
        ImGui::Text("Thermal Properties");
        
        ImGui::SliderFloat("Thermal Conductivity (W/m·K)", &tempThermalCond, 0.01, 500, "%.2f", ImGuiSliderFlags_Logarithmic);
        ImGui::SliderFloat("Specific Heat (J/kg·K)", &tempSpecificHeat, 100, 5000, "%.0f");
        
        ImGui::Separator();
        ImGui::Text("Phase Transition");
        
        ImGui::SliderFloat("Melting Point (K)", &tempMeltingPoint, 0, 4000, "%.0f");
        ImGui::SliderFloat("Boiling Point (K)", &tempBoilingPoint, 0, 6000, "%.0f");
        
        ImGui::Separator();
        ImGui::Text("Mechanical Properties");
        
        ImGui::SliderFloat("Yield Strength (MPa)", &yieldStrengthMPa, 0, 2000, "%.0f");
        ImGui::SliderFloat("Ultimate Strength (MPa)", &ultimateStrengthMPa, 0, 3000, "%.0f");
        
        ImGui::Separator();
        ImGui::Text("Visual");
        
        ImGui::ColorEdit3("Color", &customMaterial.color.x);
        
        ImGui::Separator();
        
        // Save button
        if (ImGui::Button("Save Material", ImVec2(150, 30))) {
            // Copy temp values to Material struct
            customMaterial.density = tempDensity;
            customMaterial.viscosity = tempViscosity;
            customMaterial.thermalConductivity = tempThermalCond;
            customMaterial.specificHeatCapacity = tempSpecificHeat;
            customMaterial.meltingPoint = tempMeltingPoint;
            customMaterial.boilingPoint = tempBoilingPoint;
            
            // Convert units
            customMaterial.name = std::string(materialNameBuffer);
            customMaterial.youngsModulus = youngsModulusGPa * 1e9; // GPa to Pa
            customMaterial.yieldStrength = yieldStrengthMPa * 1e6; // MPa to Pa
            customMaterial.ultimateTensileStrength = ultimateStrengthMPa * 1e6;
            customMaterial.elasticColor = customMaterial.color;
            customMaterial.plasticColor = glm::vec3(0.9f, 0.7f, 0.7f);
            customMaterial.fractureColor = glm::vec3(0.3f, 0.3f, 0.3f);
            customMaterial.vanDerWaalsA = 0.0;
            customMaterial.vanDerWaalsB = 0.0;
            
            materialDB.addMaterial(customMaterial.name, customMaterial);
            selectedMaterial = customMaterial.name;
            
            // Update current solver
            if (currentScenario >= 0 && currentScenario < static_cast<int>(solvers.size())) {
                solvers[currentScenario]->setMaterial(customMaterial);
            }
        }
        
        ImGui::SameLine();
        if (ImGui::Button("Close", ImVec2(100, 30))) {
            showMaterialEditor = false;
        }
        
        ImGui::End();
    }
    
    void saveProject(const std::string& filename) {
        json project;
        project["currentScenario"] = currentScenario;
        project["selectedMaterial"] = selectedMaterial;
        project["scientificMode"] = scientificMode;
        
        // Save all solvers states
        json solversState = json::array();
        for (const auto& solver : solvers) {
            solversState.push_back(solver->saveState());
        }
        project["solvers"] = solversState;
        
        std::ofstream file(filename);
        file << project.dump(4);
        std::cout << "Project saved to " << filename << std::endl;
    }
    
    void loadProject(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Failed to open project file: " << filename << std::endl;
            return;
        }
        
        json project;
        file >> project;
        
        if (project.contains("currentScenario")) currentScenario = project["currentScenario"];
        if (project.contains("selectedMaterial")) selectedMaterial = project["selectedMaterial"];
        if (project.contains("scientificMode")) scientificMode = project["scientificMode"];
        
        if (project.contains("solvers")) {
            json solversState = project["solvers"];
            for (size_t i = 0; i < solvers.size() && i < solversState.size(); ++i) {
                solvers[i]->loadState(solversState[i]);
            }
        }
        
        // Re-apply material to current solver to ensure consistency
        if (currentScenario >= 0 && currentScenario < static_cast<int>(solvers.size())) {
            solvers[currentScenario]->setMaterial(materialDB.getMaterial(selectedMaterial));
        }
        
        std::cout << "Project loaded from " << filename << std::endl;
    }
    
    void renderGrid() {
        glColor3f(0.3f, 0.3f, 0.3f);
        glBegin(GL_LINES);
        for (int i = -5; i <= 5; ++i) {
            glVertex3f(-5.0f, 0.0f, i);
            glVertex3f(5.0f, 0.0f, i);
            glVertex3f(i, 0.0f, -5.0f);
            glVertex3f(i, 0.0f, 5.0f);
        }
        glEnd();
    }
    
    void renderAxes() {
        glLineWidth(3.0f);
        glBegin(GL_LINES);
        glColor3f(1,0,0); glVertex3f(0,0,0); glVertex3f(1,0,0);
        glColor3f(0,1,0); glVertex3f(0,0,0); glVertex3f(0,1,0);
        glColor3f(0,0,1); glVertex3f(0,0,0); glVertex3f(0,0,1);
        glEnd();
        glLineWidth(1.0f);
    }
};
