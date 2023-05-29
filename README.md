# Liquid3D

IGR205 project

The goal of this project is to simulate fluids using the IISPH technique in order to use the OpenVDB library to create a mesh out of the fluid particles. The result can then be imported inside Blender.

## Simulation

The project includes 2D and 3D IISPH simulations that can be compiled and run by entering the following commands:

```bash
make
make run
```

The 3D IISPH outputs a `.txt` file containing the positions of the particles at each frame. This file is then to be used using the OpenVDB bridges as explained below.

## OpenVDB bridges

You can build and run the C++ bridge in the `OpenvdbCppBridge` folder using cmake or docker. The build instructions are detailled in the `README.md` file in the folder.

You can fallback on the python implementation of the bridge if you can't build the C++ bridge. It is located in the `openvdbPyBridge` folder. Note that the python implmentation works around limitations of the `pyopenvdb` module and that the resulting `.vdb` files are not as good as the ones generated using the C++ bridge. Long simulations will also take ages to convert using the python bridge.

## Blender

You can use the `vdb_blender.blend` file to import the `.vdb` files generated using the bridge. The file contains a material to visualize fluids.


## Useful links

Example vdb files: https://www.openvdb.org/download/
