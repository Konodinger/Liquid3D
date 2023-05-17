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

The pyopenvdb doc is no longer online, but you can view its backed up version here: https://web.archive.org/web/20190102171432/http://www.openvdb.org:80/documentation/doxygen/python/index.html

Some examples can be found here: https://www.openvdb.org/documentation/doxygen/python.html

More documentation: https://artifacts.aswf.io/io/aswf/openvdb/openvdb_toolset_2013/1.0.0/openvdb_toolset_2013-1.0.0.pdf
