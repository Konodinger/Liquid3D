# OpenVDB C++ Bridge

This is the C++ bridge for OpenVDB. It takes as input a file containing the positions of particles for different time
steps and outputs `.vdb` SDF and point data grids files as well as `.obj` mesh files in the `results` folder.

## Installation

This program needs OpenVDB to be installed to work.
Visit [the official repository of OpenVDB](https://github.com/AcademySoftwareFoundation/openvdb) and follow the
installation guide.

**The bridge assumes OpenVDB is installed in your home folder. If this is not the case, edit `CMakeLists.txt` and change
the path to the OpenVDB library to match your installation.**

If you have any issue with this step, please refer to the official documentation of
OpenVDB: https://www.openvdb.org/documentation/doxygen/build.html (Building with OpenVDB)

### Docker

If you can't install OpenVDB on your machine, you can use Docker to run the program.

You can pull the docker image from Docker Hub with the following command:

```
docker pull barthpaleologue/openvdb_bridge
```

You can also build the image yourself using the given `Dockerfile`.

The docker image contains this program and OpenVDB ready to be used.

## Build

You can use CMake to build the project.

```
mkdir build
cd build
cmake ..
make
```

## Usage

The program can be used with the following command:

```
./OpenVdbBridge -f <fileName>
```

The options are:

```
Options:
	-f <fileName>	Specify the file to read
	-o <fileName>	Specify the output file name
	-h		Print this help message
	-p	    Generate point grids as well
	--obj	Generate obj files as well
	--particleRadius <radius>	Specify the radius of the particles
	--voxelSize <size>	Specify the size of the voxels
	--halfWidth <width>	Specify the half width of the SDF
```

The file must be formatted as follows:

```
gridSizeX gridSizeY gridSizeZ
dt nbTimeSteps
nbParticles
x1 y1 z1
x2 y2 z2
...
```

You can see an example in the `particles_example.txt` file.

## Output

The program will generate one `.vdb` file per timestep containing a level set representation of the fluid.

You can also generate a `.vdb` containing a point data grid by adding the `-p` flag, or generate `.obj` files by adding
the `--obj` flag.

Point data grids are not compatible with Blender but the SDF and the mesh can be directly imported.