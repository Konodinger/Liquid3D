#include <openvdb/openvdb.h>
#include <openvdb/io/Stream.h>

#include <iostream>
#include <fstream>

using namespace openvdb;

int main() {
    // Initialize the OpenVDB library.  This must be called at least
    // once per program and may safely be called multiple times.
    initialize();

    std::vector<openvdb::Vec3R> particles;

    // read particles from file particles.txt
    // each line contains 3 float numbers separated by space
    std::ifstream infile("../particles.txt");
    if (!infile) {
        std::cout << "Error opening particles file" << std::endl;
        return 1;
    }

    float x, y, z;
    while (infile >> x >> y >> z) {
        particles.emplace_back(x, y, z);
    }

    std::cout << "Number of particles: " << particles.size() << std::endl;
    std::cout << "First particle: " << particles[0] << std::endl;

    // Create an empty floating-point grid with background value 0.
    FloatGrid::Ptr grid = FloatGrid::create(0.0);
    std::cout << "Testing random access:" << std::endl;
    // Get an accessor for coordinate-based access to voxels.
    FloatGrid::Accessor accessor = grid->getAccessor();
    // Define a coordinate with large signed indices.
    Coord xyz(1000, -200000000, 30000000);
    // Set the voxel value at (1000, -200000000, 30000000) to 1.
    accessor.setValue(xyz, 1.0);
    // Verify that the voxel value at (1000, -200000000, 30000000) is 1.
    std::cout << "Grid" << xyz << " = " << accessor.getValue(xyz) << std::endl;
    // Reset the coordinates to those of a different voxel.
    xyz.reset(1000, 200000000, -30000000);
    // Verify that the voxel value at (1000, 200000000, -30000000) is
    // the background value, 0.
    std::cout << "Grid" << xyz << " = " << accessor.getValue(xyz) << std::endl;
    // Set the voxel value at (1000, 200000000, -30000000) to 2.
    accessor.setValue(xyz, 2.0);
    // Set the voxels at the two extremes of the available coordinate space.
    // For 32-bit signed coordinates these are (-2147483648, -2147483648, -2147483648)
    // and (2147483647, 2147483647, 2147483647).
    accessor.setValue(Coord::min(), 3.0f);
    accessor.setValue(Coord::max(), 4.0f);
    std::cout << "Testing sequential access:" << std::endl;
    // Print all active ("on") voxels by means of an iterator.

    /*for (FloatGrid::ValueOnCIter iter = grid->cbeginValueOn(); iter; ++iter) {
        std::cout << "Grid" << iter.getCoord() << " = " << *iter << std::endl;
    }*/

    // multiple grids structure
    GridPtrVecPtr grids(new GridPtrVec);
    grids->push_back(grid);

    // Stream the grids to a file.
    std::ofstream ofile("mygrids.vdb", std::ios_base::binary);
    io::Stream(ofile).write(*grids);
    ofile.close();
}