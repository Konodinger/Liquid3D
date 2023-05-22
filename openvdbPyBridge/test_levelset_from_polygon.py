import pyopenvdb as vdb
import os
import numpy as np

# this file code uses code from https://www.openvdb.org/documentation/doxygen/python.html to mesh a point cloud

points = np.loadtxt('points_float.txt')
grid = vdb.FloatGrid()

# Define dimension of the volume
max_x = np.max(points[:, 0])
max_y = np.max(points[:, 1])
max_z = np.max(points[:, 2])

min_x = np.min(points[:, 0])
min_y = np.min(points[:, 1])
min_z = np.min(points[:, 2])

resolution = 500 #arbitrary value
voxel_size = min(max_x-min_x,max_y-min_y,max_z-min_z)/resolution
print(voxel_size)

# create empty grid

resolution_x = int(np.ceil((max_x-min_x)/voxel_size))
resolution_y = int(np.ceil((max_y-min_y)/voxel_size))
resolution_z = int(np.ceil((max_z-min_z)/voxel_size))
point_grid = np.zeros((resolution_x, resolution_y, resolution_z))

# fill the grid with the number of particles in each voxel

particle_number = len(points[:,0])
print("there is ", particle_number, "particles")

particle_count = 0
for x in np.arange(min_x, max_x, resolution_x):
    for y in np.arange(min_y, max_y, resolution_y):
        for z in np.arange(min_z, max_z, resolution_z):
            for point in points:
                if (x <= point[0] < x + resolution_x and y <= point[1] < y + resolution_y and z <= point[2] < z + resolution_z):
                    x_ind, y_ind, z_ind = int(min_x/resolution_x), int(min_y/resolution_y), int(min_z/resolution_z)
                    point_grid[x_ind][y_ind][z_ind] += 1
                    particle_count += 1

if (particle_count==particle_number):
    print("All particles have been added")
else:
    print("There are missing particles...")

print(point_grid)
print(np.sum(point_grid))



 



grid.copyFromArray(point_grid)

grid.name = 'grid'

# those are the same
points, triangles, quads = grid.convertToPolygons(adaptivity=0.5)
mesh = grid.convertToPolygons(adaptivity=0.8)











# if not results folder exists, create it
if not os.path.exists('./results'):
    os.makedirs('./results')

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

writeObjFile('./results/polygon_levelset.obj', *mesh)

# Write to VDB file
vdb.write('./results/polygon_levelset.vdb', grids=[grid])

