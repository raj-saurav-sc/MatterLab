#include "Utils.h"
#include <fstream>
#include <iostream>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../stb_image_write.h"

// Global state for image capture
bool captureNextFrame = false;
std::string captureFilename;
ImVec2 captureMin;
ImVec2 captureMax;

void requestPlotCapture(const std::string& filename, ImVec2 min, ImVec2 max) {
    captureNextFrame = true;
    captureFilename = filename;
    captureMin = min;
    captureMax = max;
}

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
