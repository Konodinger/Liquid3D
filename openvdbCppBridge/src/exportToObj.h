//
// Created by barth on 21/05/23.
//

#ifndef OPENVDBBRIDGE_EXPORTTOOBJ_H
#define OPENVDBBRIDGE_EXPORTTOOBJ_H

#include <openvdb/openvdb.h>
#include <fstream>

void exportToObj(std::vector<openvdb::Vec3s> &points, std::vector<openvdb::Vec4I> &quads,
                 std::vector<openvdb::Vec3I> &triangles) {
    std::ofstream objFile;
    objFile.open("output.obj");
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

#endif //OPENVDBBRIDGE_EXPORTTOOBJ_H
