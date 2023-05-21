//
// Created by barth on 17/05/23.
//

#ifndef OPENVDBBRIDGE_UTILS_HPP
#define OPENVDBBRIDGE_UTILS_HPP

#include <openvdb/openvdb.h>
#include <string>
#include <fstream>

void writeGrid(const openvdb::GridBase::Ptr &grid, const std::string &fileName) {
    std::cout << "\nWriting \"" << fileName << "\" to file\n";
    grid->setName(fileName);
    openvdb::GridPtrVec grids;
    grids.push_back(grid);
    openvdb::io::File file("./results/" + fileName + ".vdb");
    file.write(grids);
    file.close();
}

void exportToObj(std::vector<openvdb::Vec3s> &points, std::vector<openvdb::Vec4I> &quads,
                 std::vector<openvdb::Vec3I> &triangles, const std::string &fileName) {
    std::ofstream objFile;
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
