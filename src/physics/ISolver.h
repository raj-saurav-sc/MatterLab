#pragma once

#include <string>
#include "Material.h"
#include "../core/Utils.h"

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
    
    // Identification
    virtual std::string getName() const = 0;
};
