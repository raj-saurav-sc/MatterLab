#pragma once

#include <string>
#include "Material.h"
#include "../core/Utils.h"
#include "../core/vendor/json.hpp"

using json = nlohmann::json;

class ISolver {
public:
    virtual ~ISolver() = default;
    
    // Core lifecycle
    virtual void initialize() {}
    virtual void update(float deltaTime) {}
    virtual void reset() {}
    
    // Rendering
    virtual void render3D() {}
    virtual void renderUI() {}
    virtual void renderGraph() {}
    
    // Material handling
    virtual void setMaterial(const Material& material) {}
    
    // State Management
    virtual json saveState() const { return json::object(); }
    virtual void loadState(const json& state) {}
    
    // Identification
    virtual std::string getName() const = 0;
};
