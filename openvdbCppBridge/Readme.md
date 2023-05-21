# OpenVDB C++ Bridge

As the python module is not capable of handling a set of particles, we developped a new version using the C++ API of
OpenVDB

## Installation

Go to the git repository of OpenVDB and follow the instructions to install it.

When it is done, go to `CmakeLists.txt` and change the path to the OpenVDB library to match your installation.

## Build

You can use cmake to build the project. The executable will be in the `build` folder.

## Usage

The executable takes 0 or 1 argument.

When no arguments are given, the program uses a debug file to lookup the particles.

You can however give a path to a file containing the particle data. The file must be formatted as follow:

```
x1 y1 z1
x2 y2 z2
...
```

The program will then create 2 `.vdb` files in the `results` folder containing the point data grid and its SDF
representation.
Moreover, a `.obj` file will be created in the `results` folder containing the mesh generated from the SDF.

Point data grids are not compatible with Blender but the SDF and the mesh can be directly imported.