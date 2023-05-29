# this file code uses code from https://www.openvdb.org/documentation/doxygen/python.html to mesh a point cloud

import pyopenvdb as vdb
import os
import numpy as np
from time import time

DEBUG = True
print("\nStart\n")
t0 = time()

def from_points_to_voxel(array, resolution=10):
    """
    resolution is an arbitrary value
    """
    # Define dimension of the volume
    max_x = np.max(array[:, 0])
    max_y = np.max(array[:, 1])
    max_z = np.max(array[:, 2])
    min_x = np.min(array[:, 0])
    min_y = np.min(array[:, 1])
    min_z = np.min(array[:, 2])

    voxel_size = min(max_x-min_x,max_y-min_y,max_z-min_z)/resolution

    # create empty grid
    resolution_x = int(np.ceil((max_x-min_x)/voxel_size))
    resolution_y = int(np.ceil((max_y-min_y)/voxel_size))
    resolution_z = int(np.ceil((max_z-min_z)/voxel_size))
    array_grid = np.zeros((resolution_x, resolution_y, resolution_z))

    if DEBUG:
        print("resolution_x:", resolution_x)
        print("resolution_y:", resolution_y)
        print("resolution_z:", resolution_z)

    # fill the grid with the number of particles in each voxel
    particle_number = len(array[:,0])
    print("There are", particle_number, "particles")

    particle_count = 0
    for x in np.arange(min_x, max_x + voxel_size, voxel_size):
        #if DEBUG: print("x:", x)
        for y in np.arange(min_y, max_y + voxel_size, voxel_size):
            #if DEBUG: print("y:", y)
            for z in np.arange(min_z, max_z + voxel_size, voxel_size):
                #if DEBUG: print("z:", z)
                for point in array:
                    if ((x <= point[0] < x + voxel_size) and (y <= point[1] < y + voxel_size) and (z <= point[2] < z + voxel_size)):
                        x_ind, y_ind, z_ind = int(x/voxel_size), int(y/voxel_size), int(z/voxel_size)
                        array_grid[x_ind][y_ind][z_ind] += 1
                        particle_count += 1
                        delta_t = time() - t0
                        if DEBUG: print("Progress:", int(particle_count*10000/particle_number)/100, "% in", int(delta_t*100)/100, "s")
    if DEBUG: 
        if (particle_count==particle_number):
            print("All particles have been added in point_grid")
        else:
            print("There are missing particles in point_grid...")
    return array_grid


# Load the point cloud
points = np.loadtxt('points_float.txt')
if DEBUG : print("points.shape:", points.shape)

# Create the vdb grid
grid = vdb.FloatGrid()
grid.name = "density" # to use it in blender shaders

# Create voxel grid
point_grid = from_points_to_voxel(points, resolution=50)

if DEBUG:
    for i in range(point_grid.shape[0]):
        for j in range(point_grid.shape[1]):
            for k in range(point_grid.shape[2]):
                if (point_grid[i][j][k] > 0):
                    print("There are", point_grid[i][j][k], "particles in voxel", i, j, k)

# Convert the numpy array into vdb grid

# Convert the numpy voxel grid into vdb grid
grid.copyFromArray(point_grid)
print("There are",grid.activeVoxelCount(), "active voxels in the grid")

# Convert the vdb grid into a polygon mesh
mesh = grid.convertToPolygons(adaptivity=0.8)

# if not results folder exists, create it
if not os.path.exists('./results'):
    os.makedirs('./results')

# Write to VDB file
vdb.write('./results/polygon_levelset.vdb', grids=[grid])

print("\nDone in", int((time()-t0)*100)/100, "s\n")

