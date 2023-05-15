import pyopenvdb as vdb
import numpy as np
import os

# this file code uses code from https://www.openvdb.org/documentation/doxygen/python.html to mesh a point cloud

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

# read file points.txt that contains one 3d point on each line plus a float value
# points = np.loadtxt('points_float.txt')

points = np.random.rand(50, 3)
print(points)

# we must make a triangular mesh out of the points

# Create a FloatGrid with background value 0.0
grid = vdb.FloatGrid()

# set the value to 1 at each point
#for point in points:
    #help(grid.getAccessor())
    #grid.getAccessor().setValueOn(point, 1)
    #grid.getAccessor().setValueOn((float(point[0]), float(point[1]), float(point[2])), 1.0)

#grid.copyFromArray(values)

print(grid.evalActiveVoxelBoundingBox())
grid.name = 'grid'

# those are the same
points, triangles, quads = grid.convertToPolygons(adaptivity=0.5)
mesh = grid.convertToPolygons(adaptivity=0.8)

# if not results folder exists, create it
if not os.path.exists('./results'):
    os.makedirs('./results')

writeObjFile('./results/result.obj', *mesh)

# Write to VDB file
vdb.write('./results/grid.vdb', grids=[grid])

# write the mesh to a VDB file
# according to the documentation https://web.archive.org/web/20190102171432/http://www.openvdb.org:80/documentation/doxygen/python/index.html
# we can only write grids to a VDB file
# vdb.write('mesh.vdb', points, triangles, quads)
