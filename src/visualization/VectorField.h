#pragma once

#include <glm/glm.hpp>
#include <GL/glew.h>
#include <vector>
#include <cmath>
#include <algorithm>

class VectorField {
public:
    // Render a single arrow from start point with given direction and magnitude
    static void renderArrow(glm::vec3 start, glm::vec3 direction, float magnitude, 
                           glm::vec3 color, float scale = 1.0f) {
        if (magnitude < 0.001f) return; // Skip tiny vectors
        
        glm::vec3 end = start + direction * magnitude * scale;
        
        // Draw arrow shaft (line)
        glColor3f(color.x, color.y, color.z);
        glBegin(GL_LINES);
        glVertex3f(start.x, start.y, start.z);
        glVertex3f(end.x, end.y, end.z);
        glEnd();
        
        // Draw arrow head (simple triangle)
        float headSize = magnitude * scale * 0.2f;
        if (headSize < 0.01f) headSize = 0.01f;
        if (headSize > 0.3f) headSize = 0.3f;
        
        // Calculate perpendicular vectors for arrow head
        glm::vec3 dir = glm::normalize(direction);
        glm::vec3 perp1, perp2;
        
        if (std::abs(dir.x) < 0.9f) {
            perp1 = glm::normalize(glm::cross(dir, glm::vec3(1, 0, 0)));
        } else {
            perp1 = glm::normalize(glm::cross(dir, glm::vec3(0, 1, 0)));
        }
        perp2 = glm::cross(dir, perp1);
        
        glm::vec3 headBase = end - dir * headSize;
        glm::vec3 p1 = headBase + perp1 * headSize * 0.5f;
        glm::vec3 p2 = headBase + perp2 * headSize * 0.5f;
        glm::vec3 p3 = headBase - perp1 * headSize * 0.5f;
        glm::vec3 p4 = headBase - perp2 * headSize * 0.5f;
        
        glBegin(GL_TRIANGLES);
        // Four triangular faces forming a pyramid
        glVertex3f(end.x, end.y, end.z);
        glVertex3f(p1.x, p1.y, p1.z);
        glVertex3f(p2.x, p2.y, p2.z);
        
        glVertex3f(end.x, end.y, end.z);
        glVertex3f(p2.x, p2.y, p2.z);
        glVertex3f(p3.x, p3.y, p3.z);
        
        glVertex3f(end.x, end.y, end.z);
        glVertex3f(p3.x, p3.y, p3.z);
        glVertex3f(p4.x, p4.y, p4.z);
        
        glVertex3f(end.x, end.y, end.z);
        glVertex3f(p4.x, p4.y, p4.z);
        glVertex3f(p1.x, p1.y, p1.z);
        glEnd();
    }
    
    // Get color based on magnitude (blue -> cyan -> green -> yellow -> red)
    static glm::vec3 getColorFromMagnitude(float magnitude, float maxMagnitude) {
        if (maxMagnitude < 0.001f) return glm::vec3(0.5f, 0.5f, 0.5f);
        
        float normalized = std::min(magnitude / maxMagnitude, 1.0f);
        
        // Color map: blue (0) -> cyan (0.25) -> green (0.5) -> yellow (0.75) -> red (1)
        if (normalized < 0.25f) {
            float t = normalized / 0.25f;
            return glm::vec3(0.0f, t, 1.0f); // Blue to Cyan
        } else if (normalized < 0.5f) {
            float t = (normalized - 0.25f) / 0.25f;
            return glm::vec3(0.0f, 1.0f, 1.0f - t); // Cyan to Green
        } else if (normalized < 0.75f) {
            float t = (normalized - 0.5f) / 0.25f;
            return glm::vec3(t, 1.0f, 0.0f); // Green to Yellow
        } else {
            float t = (normalized - 0.75f) / 0.25f;
            return glm::vec3(1.0f, 1.0f - t, 0.0f); // Yellow to Red
        }
    }
    
    // Render a field of vectors
    static void renderVectorField(const std::vector<glm::vec3>& positions,
                                  const std::vector<glm::vec3>& vectors,
                                  float scale = 1.0f,
                                  bool colorByMagnitude = true,
                                  int stride = 1) {
        if (positions.size() != vectors.size()) return;
        
        // Find max magnitude for color scaling
        float maxMag = 0.0f;
        if (colorByMagnitude) {
            for (const auto& v : vectors) {
                float mag = glm::length(v);
                maxMag = std::max(maxMag, mag);
            }
        }
        
        // Render arrows
        for (size_t i = 0; i < positions.size(); i += stride) {
            glm::vec3 vec = vectors[i];
            float mag = glm::length(vec);
            
            if (mag < 0.001f) continue;
            
            glm::vec3 dir = vec / mag;
            glm::vec3 color = colorByMagnitude ? 
                getColorFromMagnitude(mag, maxMag) : 
                glm::vec3(1.0f, 0.0f, 0.0f);
            
            renderArrow(positions[i], dir, mag, color, scale);
        }
    }
};
