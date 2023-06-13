//
// Created by barth on 12/06/23.
//

#ifndef INC_3DIISPH_PARTICLEINITIALIZATION_HPP
#define INC_3DIISPH_PARTICLEINITIALIZATION_HPP

#include <vector>
#include "Vector.hpp"

void initCube(Vec3f blockPosition, Vec3f blockDimensions, std::vector<Vec3f> &particlePositions) {
    unsigned int initialSize = particlePositions.size();

    const int minX = blockPosition.x;
    const int maxX = blockPosition.x + blockDimensions.x;
    const int minY = blockPosition.y;
    const int maxY = blockPosition.y + blockDimensions.y;
    const int minZ = blockPosition.z;
    const int maxZ = blockPosition.z + blockDimensions.z;

    for (int x = minX; x < maxX; x++) {
        for (int y = minY; y < maxY; y++) {
            for (int z = minZ; z < maxZ; z++) {
                particlePositions.emplace_back(x + 1.25, y + 1.25, z + 1.25);
                particlePositions.emplace_back(x + 1.75, y + 1.25, z + 1.25);
                particlePositions.emplace_back(x + 1.25, y + 1.75, z + 1.25);
                particlePositions.emplace_back(x + 1.75, y + 1.75, z + 1.25);
                particlePositions.emplace_back(x + 1.25, y + 1.25, z + 1.25);
                particlePositions.emplace_back(x + 1.75, y + 1.25, z + 1.25);
                particlePositions.emplace_back(x + 1.25, y + 1.75, z + 1.25);
                particlePositions.emplace_back(x + 1.75, y + 1.75, z + 1.25);
            }
        }
    }
    std::cout << "Simulating " << particlePositions.size() - initialSize << " particles of fluid" << std::endl;
}

#endif //INC_3DIISPH_PARTICLEINITIALIZATION_HPP
