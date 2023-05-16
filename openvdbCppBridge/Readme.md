# OpenVDB C++ Bridge

As the python module is not capable of handling a set of particles, we developped a new version using the C++ API of
OpenVDB

## Installation

Go to the git repository of OpenVDB and follow the instructions to install it.

When it is done, go to `CmakeLists.txt` and change the path to the OpenVDB library to match your installation.

## Usage

You can use cmake to build the project. The executable will be in the `build` folder.

## Findings

Point grids can be exported to `.vdb` files. However, they are not compatible with Blender which only supports voxel
grids.