#pragma once

#include <random>

// Modern random number generation utility
class RandomGenerator {
private:
    static std::mt19937& getEngine() {
        static std::random_device rd;
        static std::mt19937 gen(rd());
        return gen;
    }

public:
    // Generate random integer in range [min, max]
    static int randInt(int min, int max) {
        std::uniform_int_distribution<int> dist(min, max);
        return dist(getEngine());
    }
    
    // Generate random float in range [min, max]
    static float randFloat(float min, float max) {
        std::uniform_real_distribution<float> dist(min, max);
        return dist(getEngine());
    }
    
    // Generate random double in range [min, max]
    static double randDouble(double min, double max) {
        std::uniform_real_distribution<double> dist(min, max);
        return dist(getEngine());
    }
    
    // Seed the generator (optional, auto-seeded by default)
    static void seed(unsigned int s) {
        getEngine().seed(s);
    }
};
