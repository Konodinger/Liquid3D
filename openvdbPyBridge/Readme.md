# OpenVDB Python Bridge

We use the `pyopenvdb` package to read and write OpenVDB files. This package has very restrictive installation requirements such as a very specific version of libc therefore we use a docker container to run the code.

The docker is based on this repository: https://github.com/surf-visualization/pyopenvdb-docker-image

To build the docker image run:

```bash
./buildDocker.sh

# or

docker build -t pyopenvdb_img .
```

Then to run the docker container:

```bash
./runDocker.sh

# or

docker run -it --rm -v `pwd`:/scripts pyopenvdb_img bash
```

The python scripts will be located at `/scripts` inside the docker container. You can run them with:

```bash
cd /scripts
python3.7 main.py
```

The pyopenvdb doc is no longer online, but you can view its backed up version here: https://web.archive.org/web/20190102171432/http://www.openvdb.org:80/documentation/doxygen/python/index.html

Some examples can be found here: https://www.openvdb.org/documentation/doxygen/python.html