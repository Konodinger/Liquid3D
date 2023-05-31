# Liquid3D

IGR205 project

The goal of this project is to simulate fluids using the IISPH technique in order to use the OpenVDB library to create a mesh out of the fluid particles. The result can then be imported inside Blender.

## Build

### Everything at once

If OpenVDB is installed on your computer in your home directory you can build both the simulation and the openvdb bridge by entering the following commands:

```bash
./build.sh
```

If OpenVDB is installed somewhere else you can specify the path to the OpenVDB installation in the file `openvdbCppBridge/CMakeLists.txt` and then run the build script.

### Simulation

The 3D IISPH Simulation can be built using CMake or using Make whichever you prefer. The 2D IISPH simulation can only be built using Make.

### OpenVDB C++ bridge

The C++ bridge can be built using CMake or using Docker. The build instructions are detailled in the `README.md` file in the folder. Should it fail, you can use our python fallback.

## Run

### Simulation

You can run the 3D IISPH simulation by running the following:

```bash
make run

# or

./3DIISPH
```

The 3D IISPH outputs a `.txt` file containing the positions of the particles at each frame. This file is then to be used using the OpenVDB bridges as explained below.

### OpenVDB bridges

To run the C++ bridge you can use the following command:

```bash
./openvdbCppBridge/openvdbCppBridge -f <path_to_particle_file>
```

If this should fail, you can fallback on the python implementation of the bridge. It is located in the `openvdbPyBridge` folder. Note that the python implmentation works around limitations of the `pyopenvdb` module and that the resulting `.vdb` files are not as good as the ones generated using the C++ bridge. Long simulations will also take ages to convert using the python bridge.

### Blender

You can use the `vdb_blender.blend` file to import the `.vdb` files generated using the bridge. The file contains a material to visualize fluids.


## Useful links

Example vdb files: https://www.openvdb.org/download/
