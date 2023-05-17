//
// Created by barth on 17/05/23.
//

#ifndef OPENVDBBRIDGE_UTILS_HPP
#define OPENVDBBRIDGE_UTILS_HPP

#include <openvdb/openvdb.h>
#include <string>

void writeGrid(const openvdb::GridBase::Ptr &grid, const std::string &fileName) {
    std::cout << "\nWriting \"" << fileName << "\" to file\n";
    grid->setName(fileName);
    openvdb::GridPtrVec grids;
    grids.push_back(grid);
    openvdb::io::File file(fileName + ".vdb");
    file.write(grids);
    file.close();
}

#endif //OPENVDBBRIDGE_UTILS_HPP
