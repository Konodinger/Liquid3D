# OpenVDB C++ Bridge

This is the C++ bridge for OpenVDB. It takes as input a file containing the positions of particles for different time
steps and outputs `.vdb` SDF and point data grids files as well as `.obj` mesh files in the `results` folder.

## Installation

This program needs OpenVDB to be installed to work.
Visit [the official repository of OpenVDB](https://github.com/AcademySoftwareFoundation/openvdb) and follow the
installation guide.

The bridge assumes OpenVDB is installed in your home folder. If this is not the case, edit `CMakeLists.txt` and change
the path to the OpenVDB library to match your installation.
nbBoundaryParticles
xb1 yb1 zb1
xb2 yb2 zb2
If you have any issue with this step, please refer to the official documentation of
OpenVDB: https://www.openvdb.org/documentation/doxygen/build.html (Building with OpenVDB)

### Docker

If you can't install OpenVDB on your machine, you can use Docker to run the program.

You can pull the docker image from Docker Hub with the following command:

```
docker pull barthpaleologue/openvdb_bridge
```

You can also build the image with the script `buildDocker.sh` since the Dockerfile is included.
The docker image is based on Ubuntu 22.04 and contains OpenVDB and all the dependencies needed to run the program.

You can then run it with `runDocker.sh`.
You will need to build the bridge inside the container. The C++ project is mounted at
`/openvdb_bridge`. See the next section for more information.

## Build

You can use CMake to build the project.

```
mkdir build
cmake -B build -S .
cmake --build build
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
	-h		Print this help message
	-p	    Generate point grids as well
	--obj	Generate obj files as well
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

The program will generate one `.vdb` file and one `.obj` file per time step. The `.vdb` file contains the SDF and the
`.obj` file contains the mesh.

You can also generate a `.vdb` containing a point data grid by adding the `-p` flag.

Point data grids are not compatible with Blender but the SDF and the mesh can be directly imported.

https://docs.blender.org/manual/en/latest/modeling/volumes/introduction.html

![image](https://github.com/Konodinger/Liquid3D/assets/31370477/70bf1565-6bd1-4577-8520-5b444dda864d)
