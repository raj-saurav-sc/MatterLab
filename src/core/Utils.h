#pragma once

#include <string>
#include <vector>
#include <imgui.h>

// Data Export
void exportToCSV(const std::string& filename, const std::vector<std::string>& headers, const std::vector<std::vector<float>>& data);

// Image Capture
extern bool captureNextFrame;
extern std::string captureFilename;
extern ImVec2 captureMin;
extern ImVec2 captureMax;

void requestPlotCapture(const std::string& filename, ImVec2 min, ImVec2 max);
