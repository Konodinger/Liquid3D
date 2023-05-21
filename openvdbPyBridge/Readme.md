# OpenVDB Python Bridge

We use the `pyopenvdb` package to read and write OpenVDB files. This package has very restrictive installation requirements such as a very specific version of libc therefore we use a docker container to run the code.

To build the docker image run:

```sh
docker build -t fluid_simu .

#or 

./buildDocker.sh
```

> - `-t`: we tag the builded image with `fluid_simu`

Then to run the docker container:

```sh
docker run -it --rm -v $(pwd):/scripts fluid_simu bash

# or

./runDockerUnix.sh
```

> - `-it`: interactive with terminal
> - `-rm`: delete at the end
> - `-v`: to bind the volume on `$(pwd):/scripts`
> - `bash`: run a bash after starting the container

## MacOs M1

M1 is a different architecture

```sh
docker build . -t fluid_simu --platform=linux/amd64
docker run --rm  -it -v $(pwd):/scripts --platform linux/amd64 fluid_simu bash
```

> - we need to add `--platform=linux/amd64` !

## Test the installed `pyopenvdb` lib

```sh
python -c 'import pyopenvdb as vdb'
```

## Usage

The python scripts will be located at `/scripts` inside the docker container (the volume on the computer is binded).
You can run them with:

```sh
cd /scripts
python main.py
```

The pyopenvdb is online again: https://academysoftwarefoundation.github.io/openvdb/python/index.html

Some examples can be found here: https://www.openvdb.org/documentation/doxygen/python.html

More documentation: https://artifacts.aswf.io/io/aswf/openvdb/openvdb_toolset_2013/1.0.0/openvdb_toolset_2013-1.0.0.pdf

As far as we understand it, the python API does not exposes low-level features required to create level sets from points, but it can create meshes out of level sets.

Either we find a way to create the level set from the points ourselves, or we use the C++ API.
The C++ version of this project is available in the folder openvdbCppBridge.

It is possible to perform the conversion using python according to this stackoverflow page: https://stackoverflow.com/questions/56965268/how-do-i-convert-a-3d-point-cloud-ply-into-a-mesh-with-faces-and-vertices
