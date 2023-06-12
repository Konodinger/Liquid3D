# Liquid3D

[![3D IISPH Makefile CI](https://github.com/Konodinger/Liquid3D/actions/workflows/makefile.yml/badge.svg)](https://github.com/Konodinger/Liquid3D/actions/workflows/makefile.yml)
[![3D IISPH CMake](https://github.com/Konodinger/Liquid3D/actions/workflows/cmake.yml/badge.svg)](https://github.com/Konodinger/Liquid3D/actions/workflows/cmake.yml)

This is the main repository for our IGR205 project at Telecom Paris.

The goal of this project is to simulate fluids using the IISPH technique in order to use the OpenVDB library to create a mesh out of the fluid particles. The result can then be imported inside Blender.

![Demo Image](./cover.gif)

## Build

### Simulation

The 3D IISPH Simulation can be built using CMake or using Make whichever you prefer. The 2D IISPH simulation can only be built using Make.

### OpenVDB C++ bridge

The C++ bridge can be built using CMake or using Docker. The build instructions are detailled in the `README.md` file in the `openVdbCppBridge` folder. Should it fail, you can use our python fallback.

### Everything at once

If OpenVDB is installed on your computer in your home directory you can build both the simulation and the openvdb bridge by entering the following commands:

```bash
./build.sh
```

If OpenVDB is installed somewhere else than your home directory you must specify the path to the OpenVDB installation in the file `openvdbCppBridge/CMakeLists.txt` and then run the build script.

## Run

### Simulation

You can run the 3D IISPH simulation by running the following command inside the `simulation` folder:

```bash
make run

# or

./3DIISPH
```

The 3D IISPH outputs a `.txt` file containing the positions of the particles at each frame. This file is then to be used using the OpenVDB bridges as explained below.

### OpenVDB bridge

To run the C++ bridge you can use the following command:

```bash
./openvdbCppBridge/openvdbCppBridge -f <path_to_particle_file>
```

More options are available, you can see them by running the command with the `-h` flag.

If this should fail, you can fallback on the python implementation of the bridge. It is located in the `pythonFallback` folder. You will obtain the same result as with the C++ bridge but you will notice a significant performance drop (at least 100x slower).

### Blender

You can use the `vdb_blender.blend` file to import the `.vdb` files generated using the bridge. The file contains a material to visualize fluids.
