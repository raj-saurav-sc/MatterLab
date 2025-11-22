#pragma once

#include <string>
#include <vector>
#include <map>
#include <glm/glm.hpp>

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

struct MaterialState {
    double totalStrain = 0.0;
    double elasticStrain = 0.0;
    double plasticStrain = 0.0;
    double currentYieldStrength = 0.0;
    bool isFractured = false;
    double damageParameter = 0.0;
    std::vector<glm::vec3> crackTips;
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
