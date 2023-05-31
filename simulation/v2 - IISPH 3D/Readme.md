# 3D IISPH

This is an implementation of the Implicit Incompressible SPH method for 3D fluids, as described in the
paper [Implicit Incompressible SPH](https://cg.informatik.uni-freiburg.de/publications/2013_TVCG_IISPH.pdf).

## Build

### Using CMake

The project can be built using CMake. The following commands can be used to build the project:

```bash
mkdir build
cmake -B build -S .
cmake --build build
```

### Using Make

The project can also be built using Make. The following commands can be used to build the project:

```bash
make
```

## Run

The project can be run using the following command:

```bash
./3DIISPH
```

or

```bash
make run
```