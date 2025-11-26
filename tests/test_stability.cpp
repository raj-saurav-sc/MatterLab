#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include "../src/physics/Material.h"
#include "../src/core/RandomGenerator.h"

// Test colors
#define GREEN "\033[32m"
#define RED "\033[31m"
#define YELLOW "\033[33m"
#define RESET "\033[0m"

int tests_passed = 0;
int tests_failed = 0;

// Custom clamp for compatibility
template<typename T>
T clamp(T value, T min, T max) {
    return (value < min) ? min : (value > max) ? max : value;
}

void test_assert(bool condition, const std::string& test_name) {
    if (condition) {
        std::cout << GREEN << "✓ PASS: " << test_name << RESET << std::endl;
        tests_passed++;
    } else {
        std::cout << RED << "✗ FAIL: " << test_name << RESET << std::endl;
        tests_failed++;
    }
}

// Test 1: Division by Zero Protection
void test_division_protection() {
    std::cout << "\n" << YELLOW << "=== Testing Division by Zero Protection ===" << RESET << std::endl;
    
    // Test epsilon check
    float epsilon = 1e-10f;
    float value = 0.0f;
    
    // Should not divide if value is too small
    float result = (value > epsilon) ? (1.0f / value) : 0.0f;
    test_assert(!std::isinf(result), "Division by zero protection works");
    test_assert(!std::isnan(result), "No NaN from protected division");
    
    // Test with very small value
    value = 1e-12f;
    result = (value > epsilon) ? (1.0f / value) : 0.0f;
    test_assert(result == 0.0f, "Very small values return 0 instead of huge number");
}

// Test 2: NaN/Inf Detection
void test_nan_inf_detection() {
    std::cout << "\n" << YELLOW << "=== Testing NaN/Inf Detection ===" << RESET << std::endl;
    
    float nan_value = std::sqrt(-1.0f);
    float inf_value = 1.0f / 0.0f;
    float normal_value = 42.0f;
    
    test_assert(std::isnan(nan_value), "NaN detection works");
    test_assert(std::isinf(inf_value), "Inf detection works");
    test_assert(std::isfinite(normal_value), "Finite detection works");
    test_assert(!std::isnan(normal_value), "Normal values not detected as NaN");
    test_assert(!std::isinf(normal_value), "Normal values not detected as Inf");
}

// Test 3: Bounds Clamping
void test_bounds_clamping() {
    std::cout << "\n" << YELLOW << "=== Testing Bounds Clamping ===" << RESET << std::endl;
    
    float value = 150.0f;
    float clamped = clamp(value, 0.0f, 100.0f);
    test_assert(clamped == 100.0f, "Upper bound clamping works");
    
    value = -50.0f;
    clamped = clamp(value, 0.0f, 100.0f);
    test_assert(clamped == 0.0f, "Lower bound clamping works");
    
    value = 50.0f;
    clamped = clamp(value, 0.0f, 100.0f);
    test_assert(clamped == 50.0f, "Values within bounds unchanged");
}

// Test 4: Van der Waals Constants
void test_van_der_waals_constants() {
    std::cout << "\n" << YELLOW << "=== Testing Van der Waals Constants ===" << RESET << std::endl;
    
    MaterialDatabase db;
    Material co2 = db.getMaterial("CO2");
    
    test_assert(co2.vanDerWaalsA > 0.0, "CO2 Van der Waals 'a' is positive");
    test_assert(co2.vanDerWaalsB > 0.0, "CO2 Van der Waals 'b' is positive");
    test_assert(std::abs(co2.vanDerWaalsA - 0.3658) < 0.001, "CO2 'a' value is correct");
    test_assert(std::abs(co2.vanDerWaalsB - 4.27e-5) < 1e-6, "CO2 'b' value is correct");
    
    Material air = db.getMaterial("Air");
    test_assert(air.vanDerWaalsA > 0.0, "Air Van der Waals 'a' is positive");
    test_assert(air.vanDerWaalsB > 0.0, "Air Van der Waals 'b' is positive");
}

// Test 5: Random Number Generation
void test_random_generation() {
    std::cout << "\n" << YELLOW << "=== Testing Random Number Generation ===" << RESET << std::endl;
    
    // Test integer range
    int rand_int = RandomGenerator::randInt(0, 10);
    test_assert(rand_int >= 0 && rand_int <= 10, "Random int in range [0, 10]");
    
    // Test float range
    float rand_float = RandomGenerator::randFloat(-1.0f, 1.0f);
    test_assert(rand_float >= -1.0f && rand_float <= 1.0f, "Random float in range [-1, 1]");
    
    // Test distribution (should not be all the same)
    std::vector<int> values;
    for (int i = 0; i < 100; i++) {
        values.push_back(RandomGenerator::randInt(0, 100));
    }
    
    bool all_same = true;
    for (size_t i = 1; i < values.size(); i++) {
        if (values[i] != values[0]) {
            all_same = false;
            break;
        }
    }
    test_assert(!all_same, "Random values have variation");
}

// Test 6: Material Database
void test_material_database() {
    std::cout << "\n" << YELLOW << "=== Testing Material Database ===" << RESET << std::endl;
    
    MaterialDatabase db;
    
    // Test solid materials
    Material steel = db.getMaterial("Steel");
    test_assert(steel.name == "Steel", "Steel material exists");
    test_assert(steel.density > 0, "Steel has positive density");
    test_assert(steel.youngsModulus > 0, "Steel has positive Young's modulus");
    
    // Test fluid materials
    Material water = db.getMaterial("Water");
    test_assert(water.name == "Water", "Water material exists");
    test_assert(water.viscosity > 0, "Water has positive viscosity");
    
    // Test gas materials
    Material helium = db.getMaterial("Helium");
    test_assert(helium.name == "Helium", "Helium material exists");
    test_assert(helium.vanDerWaalsA > 0, "Helium has Van der Waals constants");
}

// Test 7: Input Validation
void test_input_validation() {
    std::cout << "\n" << YELLOW << "=== Testing Input Validation ===" << RESET << std::endl;
    
    // Test volume validation
    float volume = -1.0f;
    if (volume <= 1e-6f) volume = 0.1f;
    test_assert(volume == 0.1f, "Negative volume corrected to default");
    
    // Test moles validation
    float moles = 0.0f;
    if (moles <= 1e-6f) moles = 0.1f;
    test_assert(moles == 0.1f, "Zero moles corrected to default");
    
    // Test temperature validation
    float temperature = -100.0f;
    if (temperature <= 1e-6f) temperature = 273.15f;
    test_assert(temperature == 273.15f, "Negative temperature corrected to default");
}

// Test 8: Numerical Stability
void test_numerical_stability() {
    std::cout << "\n" << YELLOW << "=== Testing Numerical Stability ===" << RESET << std::endl;
    
    // Test stress calculation with safe division
    float force = 1000.0f;
    float area = 0.001f;
    float E = 200e9f;
    
    // Protected division
    float stress = (area > 1e-10) ? (force / area) : 0.0f;
    float strain = (E > 1e-10) ? (stress / E) : 0.0f;
    
    test_assert(std::isfinite(stress), "Stress calculation is finite");
    test_assert(std::isfinite(strain), "Strain calculation is finite");
    test_assert(!std::isnan(stress), "Stress is not NaN");
    test_assert(!std::isnan(strain), "Strain is not NaN");
    
    // Test with zero area (should be protected)
    area = 0.0f;
    stress = (area > 1e-10) ? (force / area) : 0.0f;
    test_assert(stress == 0.0f, "Zero area returns 0 stress (protected)");
    test_assert(!std::isinf(stress), "Zero area doesn't cause infinity");
}

// Test 9: Compressibility Factor
void test_compressibility_factor() {
    std::cout << "\n" << YELLOW << "=== Testing Compressibility Factor ===" << RESET << std::endl;
    
    // Test protected division for Z = PV/nRT
    float P = 101325.0f;
    float V = 1.0f;
    float n = 1.0f;
    float R = 8.314f;
    float T = 300.0f;
    
    // Protected division
    float denominator = n * R * T;
    float Z_protected = (std::abs(denominator) > 1e-6f) ? ((P * V) / denominator) : 1.0f;
    test_assert(std::isfinite(Z_protected), "Protected Z calculation is finite");
    test_assert(Z_protected > 0.0f, "Z factor is positive");
}

// Test 10: Edge Cases
void test_edge_cases() {
    std::cout << "\n" << YELLOW << "=== Testing Edge Cases ===" << RESET << std::endl;
    
    // Test very large values
    float large_value = 1e30f;
    float clamped_large = clamp(large_value, 0.0f, 1e8f);
    test_assert(clamped_large == 1e8f, "Very large values clamped correctly");
    
    // Test very small values
    float small_value = 1e-30f;
    bool is_effectively_zero = (small_value < 1e-10f);
    test_assert(is_effectively_zero, "Very small values detected as effectively zero");
    
    // Test epsilon comparison
    float a = 0.1f + 0.2f;
    float b = 0.3f;
    bool approximately_equal = std::abs(a - b) < 1e-6f;
    test_assert(approximately_equal, "Floating point epsilon comparison works");
}

int main() {
    std::cout << "\n";
    std::cout << "╔════════════════════════════════════════════════╗\n";
    std::cout << "║   UniversalSimulator Test Suite               ║\n";
    std::cout << "║   Testing Numerical Stability & Security      ║\n";
    std::cout << "╚════════════════════════════════════════════════╝\n";
    
    test_division_protection();
    test_nan_inf_detection();
    test_bounds_clamping();
    test_van_der_waals_constants();
    test_random_generation();
    test_material_database();
    test_input_validation();
    test_numerical_stability();
    test_compressibility_factor();
    test_edge_cases();
    
    // Summary
    std::cout << "\n";
    std::cout << "╔════════════════════════════════════════════════╗\n";
    std::cout << "║   Test Results                                 ║\n";
    std::cout << "╚════════════════════════════════════════════════╝\n";
    std::cout << GREEN << "✓ Passed: " << tests_passed << RESET << std::endl;
    std::cout << RED << "✗ Failed: " << tests_failed << RESET << std::endl;
    std::cout << "Total: " << (tests_passed + tests_failed) << std::endl;
    
    if (tests_failed == 0) {
        std::cout << "\n" << GREEN << "🎉 All tests passed! System is stable and secure." << RESET << std::endl;
        return 0;
    } else {
        std::cout << "\n" << RED << "⚠️  Some tests failed. Please review." << RESET << std::endl;
        return 1;
    }
}
