//
// Created by barth on 12/06/23.
//

#ifndef INC_3DIISPH_PARTICLEINITIALIZATION_HPP
#define INC_3DIISPH_PARTICLEINITIALIZATION_HPP

#include <vector>
#include "Vector.hpp"

enum class InitType {
    BLOCK, SPHERE, TORUS
};

void initBlock(Vec3f blockPosition, Vec3f blockDimensions, std::vector<Vec3f> &particlePositions) {
    const int minX = blockPosition.x - blockDimensions.x / 2.0f;
    const int maxX = blockPosition.x + blockDimensions.x / 2.0f;

    const int minY = blockPosition.y - blockDimensions.y / 2.0f;
    const int maxY = blockPosition.y + blockDimensions.y / 2.0f;

    const int minZ = blockPosition.z - blockDimensions.z / 2.0f;
    const int maxZ = blockPosition.z + blockDimensions.z / 2.0f;

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
}

void initSphere(Vec3f spherePosition, float sphereRadius, std::vector<Vec3f> &particlePositions) {
    const int minX = spherePosition.x - sphereRadius;
    const int maxX = spherePosition.x + sphereRadius;

    const int minY = spherePosition.y - sphereRadius;
    const int maxY = spherePosition.y + sphereRadius;

    const int minZ = spherePosition.z - sphereRadius;
    const int maxZ = spherePosition.z + sphereRadius;

    for (int x = minX; x < maxX; x++) {
        for (int y = minY; y < maxY; y++) {
            for (int z = minZ; z < maxZ; z++) {
                if (pow(x - spherePosition.x, 2.0f) + pow(y - spherePosition.y, 2.0f) +
                    pow(z - spherePosition.z, 2.0f) < pow(sphereRadius, 2.0f)) {
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
    }
}

void initTorus(Vec3f torusPosition, float majorRadius, float minorRadius, std::vector<Vec3f> &particlePositions) {
    const int minX = torusPosition.x - majorRadius - minorRadius;
    const int maxX = torusPosition.x + majorRadius + minorRadius;

    const int minY = torusPosition.y - majorRadius - minorRadius;
    const int maxY = torusPosition.y + majorRadius + minorRadius;

    const int minZ = torusPosition.z - minorRadius;
    const int maxZ = torusPosition.z + minorRadius;

    for (int x = minX; x < maxX; x++) {
        for (int y = minY; y < maxY; y++) {
            for (int z = minZ; z < maxZ; z++) {
                float localX = x - torusPosition.x;
                float localY = y - torusPosition.y;
                float localZ = z - torusPosition.z;

                // reference: https://math.stackexchange.com/questions/4380905/how-can-i-determine-if-a-point-x-y-z-is-within-a-torus-r-r
                float u = majorRadius - sqrt(localX * localX + localY * localY);
                float expression = u * u + localZ * localZ - minorRadius * minorRadius;

                if (expression < 0.0f) {
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
    }
}

#endif //INC_3DIISPH_PARTICLEINITIALIZATION_HPP
