#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <glm/glm.hpp>

class VTKExporter {
public:
    // Export structured grid (for CFD - regular grid)
    static bool exportStructuredGrid(
        const std::string& filename,
        int nx, int ny, int nz,
        const std::vector<float>& scalarData,
        const std::string& scalarName,
        const std::vector<glm::vec3>& vectorData = {},
        const std::string& vectorName = ""
    ) {
        std::ofstream file(filename);
        if (!file.is_open()) return false;
        
        // VTK XML header
        file << "<?xml version=\"1.0\"?>\n";
        file << "<VTKFile type=\"StructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
        file << "  <StructuredGrid WholeExtent=\"0 " << nx << " 0 " << ny << " 0 " << nz << "\">\n";
        file << "    <Piece Extent=\"0 " << nx << " 0 " << ny << " 0 " << nz << "\">\n";
        
        // Points (grid coordinates)
        file << "      <Points>\n";
        file << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        
        float dx = 10.0f / nx;
        float dy = 10.0f / ny;
        float dz = (nz > 1) ? 10.0f / nz : 1.0f;
        
        for (int k = 0; k <= nz; k++) {
            for (int j = 0; j <= ny; j++) {
                for (int i = 0; i <= nx; i++) {
                    float x = -5.0f + i * dx;
                    float y = -5.0f + j * dy;
                    float z = (nz > 1) ? (-5.0f + k * dz) : 0.0f;
                    file << "          " << x << " " << y << " " << z << "\n";
                }
            }
        }
        
        file << "        </DataArray>\n";
        file << "      </Points>\n";
        
        // Point Data
        std::string scalarsAttr = scalarName.empty() ? "" : " Scalars=\"" + scalarName + "\"";
        std::string vectorsAttr = vectorName.empty() ? "" : " Vectors=\"" + vectorName + "\"";
        file << "      <PointData" << scalarsAttr << vectorsAttr << ">\n";
        
        // Scalar data
        if (!scalarData.empty() && !scalarName.empty()) {
            file << "        <DataArray type=\"Float32\" Name=\"" << scalarName << "\" format=\"ascii\">\n";
            for (size_t i = 0; i < scalarData.size(); i++) {
                file << "          " << scalarData[i];
                if ((i + 1) % 10 == 0) file << "\n";
            }
            file << "\n        </DataArray>\n";
        }
        
        // Vector data
        if (!vectorData.empty() && !vectorName.empty()) {
            file << "        <DataArray type=\"Float32\" Name=\"" << vectorName << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
            for (size_t i = 0; i < vectorData.size(); i++) {
                file << "          " << vectorData[i].x << " " << vectorData[i].y << " " << vectorData[i].z << "\n";
            }
            file << "        </DataArray>\n";
        }
        
        file << "      </PointData>\n";
        file << "    </Piece>\n";
        file << "  </StructuredGrid>\n";
        file << "</VTKFile>\n";
        
        file.close();
        return true;
    }
    
    // Export unstructured grid (for FEM - irregular mesh)
    static bool exportUnstructuredGrid(
        const std::string& filename,
        const std::vector<glm::vec3>& points,
        const std::vector<std::vector<int>>& cells,
        const std::vector<float>& scalarData,
        const std::string& scalarName,
        const std::vector<glm::vec3>& vectorData = {},
        const std::string& vectorName = ""
    ) {
        std::ofstream file(filename);
        if (!file.is_open()) return false;
        
        // VTK XML header
        file << "<?xml version=\"1.0\"?>\n";
        file << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
        file << "  <UnstructuredGrid>\n";
        file << "    <Piece NumberOfPoints=\"" << points.size() << "\" NumberOfCells=\"" << cells.size() << "\">\n";
        
        // Points
        file << "      <Points>\n";
        file << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (const auto& p : points) {
            file << "          " << p.x << " " << p.y << " " << p.z << "\n";
        }
        file << "        </DataArray>\n";
        file << "      </Points>\n";
        
        // Cells
        file << "      <Cells>\n";
        
        // Connectivity
        file << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
        for (const auto& cell : cells) {
            file << "          ";
            for (int idx : cell) {
                file << idx << " ";
            }
            file << "\n";
        }
        file << "        </DataArray>\n";
        
        // Offsets
        file << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
        int offset = 0;
        for (const auto& cell : cells) {
            offset += cell.size();
            file << "          " << offset << "\n";
        }
        file << "        </DataArray>\n";
        
        // Cell types (3 = line, 5 = triangle, 9 = quad, 10 = tetra, 12 = hexahedron)
        file << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
        for (const auto& cell : cells) {
            int type = 3; // Default to line
            if (cell.size() == 2) type = 3;      // Line
            else if (cell.size() == 3) type = 5; // Triangle
            else if (cell.size() == 4) type = 9; // Quad
            file << "          " << type << "\n";
        }
        file << "        </DataArray>\n";
        file << "      </Cells>\n";
        
        // Point Data
        std::string scalarsAttr = scalarName.empty() ? "" : " Scalars=\"" + scalarName + "\"";
        std::string vectorsAttr = vectorName.empty() ? "" : " Vectors=\"" + vectorName + "\"";
        file << "      <PointData" << scalarsAttr << vectorsAttr << ">\n";
        
        // Scalar data
        if (!scalarData.empty() && !scalarName.empty()) {
            file << "        <DataArray type=\"Float32\" Name=\"" << scalarName << "\" format=\"ascii\">\n";
            for (float val : scalarData) {
                file << "          " << val << "\n";
            }
            file << "        </DataArray>\n";
        }
        
        // Vector data
        if (!vectorData.empty() && !vectorName.empty()) {
            file << "        <DataArray type=\"Float32\" Name=\"" << vectorName << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
            for (const auto& v : vectorData) {
                file << "          " << v.x << " " << v.y << " " << v.z << "\n";
            }
            file << "        </DataArray>\n";
        }
        
        file << "      </PointData>\n";
        file << "    </Piece>\n";
        file << "  </UnstructuredGrid>\n";
        file << "</VTKFile>\n";
        
        file.close();
        return true;
    }
};
