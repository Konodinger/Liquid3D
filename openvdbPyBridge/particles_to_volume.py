import pyopenvdb as vdb
import numpy as np
import os

# read file points.txt that contains one 3d point on each line plus a float value
points = np.loadtxt('points_float.txt')

max_x = np.max(points[:, 0])
max_y = np.max(points[:, 1])
max_z = np.max(points[:, 2])

min_x = np.min(points[:, 0])
min_y = np.min(points[:, 1])
min_z = np.min(points[:, 2])

resolution_x = (max_x - min_x) / 10
resolution_y = (max_y - min_y) / 10
resolution_z = (max_z - min_z) / 10

# Create a FloatGrid with background value 0.0
grid = vdb.FloatGrid()

# for each voxel of the grid, check if there is a point inside it
# if there is, set the value to 1
for x in np.arange(min_x, max_x, resolution_x):
    for y in np.arange(min_y, max_y, resolution_y):
        for z in np.arange(min_z, max_z, resolution_z):
            for point in points:
                if x <= point[0] < x + resolution_x and y <= point[1] < y + resolution_y and z <= point[2] < z + resolution_z:
                    grid.getAccessor().setValueOn((int(x), int(y), int(z)), 1.0)


print(grid.evalActiveVoxelBoundingBox())
grid.name = 'grid'

# Write to VDB file
vdb.write('./results/volume.vdb', grids=[grid])