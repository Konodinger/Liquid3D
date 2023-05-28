import pyopenvdb as vdb
import os
import numpy as np

# this file code uses code from https://www.openvdb.org/documentation/doxygen/python.html to mesh a point cloud

DEBUG = True

# Write to OBJ file
def writeObjFile(filename, points, triangles=[], quads=[]):
    """Write out a triangulated or quad mesh OBJ file.
    points is an array of (x,y,z) coordinates.
    triangles is an array of (i,j,k) tuples, where i,j,k are the
    indices into points of the three vertices that make up a triangle.
    quads is an array of (i,j,k,l) tuples, where i,j,k,l are the
    indices into points of the four vertices that make up a quad.

    src: https://www.openvdb.org/documentation/doxygen/python.html
    """
    f = open(filename, 'w')
    # Output vertices.
    for xyz in points:
        f.write('v %g %g %g\n' % tuple(xyz))
    f.write('\n')
    # Output faces.
    for ijk in triangles:
        f.write('f %d %d %d\n' %
            (ijk[0]+1, ijk[1]+1, ijk[2]+1)) # offset vertex indices by one
    for ijkl in quads:
        f.write('f %d %d %d %d\n' %
            (ijkl[0]+1, ijkl[1]+1, ijkl[2]+1, ijkl[3]+1))
    f.close()


print("Start\n")

points = np.loadtxt('points_float.txt')
if DEBUG : print("points.shape:", points.shape)
grid = vdb.FloatGrid()

# Define dimension of the volume
max_x = np.max(points[:, 0])
max_y = np.max(points[:, 1])
max_z = np.max(points[:, 2])

min_x = np.min(points[:, 0])
min_y = np.min(points[:, 1])
min_z = np.min(points[:, 2])

if DEBUG:
    print("max_x:", max_x)
    print("max_y:", max_y)
    print("max_z:", max_z)
    print("min_x:", min_x)
    print("min_y:", min_y)
    print("min_z:", min_z)

resolution = 10 #arbitrary value
voxel_size = min(max_x-min_x,max_y-min_y,max_z-min_z)/resolution

# create empty grid

resolution_x = int(np.ceil((max_x-min_x)/voxel_size))
resolution_y = int(np.ceil((max_y-min_y)/voxel_size))
resolution_z = int(np.ceil((max_z-min_z)/voxel_size))
point_grid = np.zeros((resolution_x, resolution_y, resolution_z))

if DEBUG:
    print("resolution_x:", resolution_x)
    print("resolution_y:", resolution_y)
    print("resolution_z:", resolution_z)

# fill the grid with the number of particles in each voxel

particle_number = len(points[:,0])
print("There are", particle_number, "particles")

particle_count = 0
for x in np.arange(min_x, max_x + voxel_size, voxel_size):
    if DEBUG: print("x:", x)
    for y in np.arange(min_y, max_y + voxel_size, voxel_size):
        for z in np.arange(min_z, max_z + voxel_size, voxel_size):
            for point in points:
                if ((x <= point[0] < x + voxel_size) and (y <= point[1] < y + voxel_size) and (z <= point[2] < z + voxel_size)):
                    x_ind, y_ind, z_ind = int(x/voxel_size), int(y/voxel_size), int(z/voxel_size)
                    point_grid[x_ind][y_ind][z_ind] += 1
                    particle_count += 1

    

if (particle_count==particle_number):
    print("All particles have been added in point_grid")
else:
    print("There are missing particles in point_grid...")

for i in range(point_grid.shape[0]):
    for j in range(point_grid.shape[1]):
        for k in range(point_grid.shape[2]):
            if (point_grid[i][j][k] > 0):
                print("There are", point_grid[i][j][k], "particles in voxel", i, j, k)

# Convert the numpy array into vdb grid

grid.copyFromArray(point_grid)
print(grid.activeVoxelCount())



# those are the same
#points, triangles, quads = grid.convertToPolygons(adaptivity=0.5)
mesh = grid.convertToPolygons(adaptivity=0.8)
print(mesh)

# if not results folder exists, create it
if not os.path.exists('./results'):
    os.makedirs('./results')

writeObjFile('./results/polygon_levelset.obj', *mesh)

# Write to VDB file
vdb.write('./results/polygon_levelset.vdb', grids=[grid])

print("\nDone")

