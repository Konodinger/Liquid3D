# Liquid3D

IGR205 project

The goal of this project is to simulate fluids using the IISPH technique in order to use the OpenVDB library to create a mesh out of the fluid particles. The result can then be imported inside Blender.

## OpenVDB bridges

The project provides 2 OpenVDB bridges. The first written in python is not working yet because PyOpenVDB does not expose low-level SDF functions of the OpenVDB API. We will try to overcome this challenge using other libraries. The instructions to run it are available in its folder with some examples of what can be done with it.

The second bridge is written in C++. Contrary to the python bridge, the C++ bridge is working as expected. It takes a particle file containing the positions of the particles, converts it to a SDF and then makes a mesh out of it. The bridge outputs a `.vdb` file containing the SDF and a `.obj` file containing the mesh. Both can be imported in Blender easily.

To compile the C++ bridge on your machine, make sure to change the path to your own OpenVDB installation.


## Useful links

Example vdb files: https://www.openvdb.org/download/
