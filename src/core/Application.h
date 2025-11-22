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

#include <iostream>
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
        
        for (auto& solver : solvers) {
            solver->initialize();
            solver->setMaterial(materialDB.getMaterial(selectedMaterial));
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
        if (currentScenario >= 0 && currentScenario < solvers.size()) {
            solvers[currentScenario]->update(deltaTime);
        }
    }
    
    void handleInput() {
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
                
                if (cameraPitch > 89.0f) cameraPitch = 89.0f;
                if (cameraPitch < -89.0f) cameraPitch = -89.0f;
                
                ImGui::ResetMouseDragDelta(ImGuiMouseButton_Right);
            }
        }
        setupCamera();
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
        
        if (currentScenario >= 0 && currentScenario < solvers.size()) {
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
        
        // Scenario Selection
        const char* scenarios[] = {
            "Solid Mechanics", "Fluid Mechanics", "Gas Dynamics", "Phase Changes",
            "Thermal Conductivity", "Projectile Motion", "Orbital Mechanics",
            "Pipe Flow", "Fluid Convection", "Coupled Heat Transfer"
        };
        if (ImGui::Combo("Scenario", &currentScenario, scenarios, IM_ARRAYSIZE(scenarios))) {
            // Reset camera or state if needed
        }
        
        // Material Selection
        std::vector<std::string> materials = materialDB.getMaterialNames();
        if (ImGui::BeginCombo("Material", selectedMaterial.c_str())) {
            for (const auto& name : materials) {
                bool isSelected = (selectedMaterial == name);
                if (ImGui::Selectable(name.c_str(), isSelected)) {
                    selectedMaterial = name;
                    // Update current solver with new material
                    if (currentScenario >= 0 && currentScenario < solvers.size()) {
                        solvers[currentScenario]->setMaterial(materialDB.getMaterial(selectedMaterial));
                    }
                }
                if (isSelected) ImGui::SetItemDefaultFocus();
            }
            ImGui::EndCombo();
        }
        
        ImGui::Checkbox("Scientific Mode", &scientificMode);
        
        ImGui::Separator();
        
        // Solver-specific UI
        if (currentScenario >= 0 && currentScenario < solvers.size()) {
            solvers[currentScenario]->renderUI();
        }
        
        ImGui::End();
        
        // Graph Panel
        ImGui::SetNextWindowPos(ImVec2(10, 620), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(600, 300), ImGuiCond_FirstUseEver);
        ImGui::Begin("Data Visualization");
        if (currentScenario >= 0 && currentScenario < solvers.size()) {
            solvers[currentScenario]->renderGraph();
        }
        ImGui::End();
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
