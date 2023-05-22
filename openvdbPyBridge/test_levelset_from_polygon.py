import pyopenvdb as vdb
import os
import numpy as np

# this file code uses code from https://www.openvdb.org/documentation/doxygen/python.html to mesh a point cloud

# we must make a triangular mesh out of the points

points = np.loadtxt('points_float.txt')
print(points)
print(points.shape)
grid = vdb.FloatGrid()


#grid = grid.createLevelSetFromPolygons(points)
points = np.random.rand(200, 200, 200)
print(points.shape)
print(points)
grid.copyFromArray(points)
"""At this point we have a problem of dimensions, 
we have to figure out what shape the array has to have"""

print(grid.evalActiveVoxelBoundingBox())
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

