import pyopenvdb as vdb
import os
import numpy as np

# this file code uses code from https://www.openvdb.org/documentation/doxygen/python.html to mesh a point cloud

# we must make a triangular mesh out of the points

points = np.loadtxt('points_float.txt')

grid = vdb.FloatGrid()
grid = grid.createLevelSetFromPolygons(points)

print(grid.evalActiveVoxelBoundingBox())
grid.name = 'grid'

# those are the same
points, triangles, quads = grid.convertToPolygons(adaptivity=0.5)
mesh = grid.convertToPolygons(adaptivity=0.8)

# if not results folder exists, create it
if not os.path.exists('./results'):
    os.makedirs('./results')

# Write to VDB file
vdb.write('./results/polygon_levelset.vdb', grids=[grid])