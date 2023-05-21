//
// Created by barth on 17/05/23.
//

#ifndef OPENVDBBRIDGE_UTILS_HPP
#define OPENVDBBRIDGE_UTILS_HPP

#include <openvdb/openvdb.h>
#include <string>
#include <fstream>
#include <iostream>

/**
 * Write an openvdb grid to a file in the results folder
 * @param grid the grid to write
 * @param fileName the name of the file
 */
void writeGrid(const openvdb::GridBase::Ptr &grid, const std::string &fileName) {
    // if results folder does not exist, create it
    if (!std::ifstream("./results")) {
        std::cout << "Creating results folder" << std::endl;
        system("mkdir results");
    }

    std::cout << "\nWriting \"" << fileName << "\" to file\n";
    grid->setName(fileName);
    openvdb::GridPtrVec grids;
    grids.push_back(grid);
    openvdb::io::File file("./results/" + fileName + ".vdb");
    file.write(grids);
    file.close();
}

/**
 * Export a list of points and quads to an obj file in the results folder
 * @param points the list of points/vertices of the mesh
 * @param quads the list of quads of the mesh
 * @param triangles the list of triangles of the mesh
 * @param fileName the name of the file
 */
void exportToObj(std::vector<openvdb::Vec3s> &points, std::vector<openvdb::Vec4I> &quads,
                 std::vector<openvdb::Vec3I> &triangles, const std::string &fileName) {
    // if results folder does not exist, create it
    if (!std::ifstream("./results")) {
        std::cout << "Creating results folder" << std::endl;
        system("mkdir results");
    }

    std::ofstream objFile;
    std::cout << "\nWriting \"" << fileName << "\" to file\n";
    objFile.open("./results/" + fileName + ".obj");
    for (auto &p: points) {
        objFile << "v " << p.x() << " " << p.y() << " " << p.z() << "\n";
    }
    for (auto &q: quads) {
        objFile << "f " << q.x() << " " << q.y() << " " << q.z() << " " << q.w() << "\n";
    }
    for (auto &t: triangles) {
        objFile << "f " << t.x() << " " << t.y() << " " << t.z() << "\n";
    }
    objFile.close();
}

#endif //OPENVDBBRIDGE_UTILS_HPP
