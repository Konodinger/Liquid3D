from parse_file import parse_file
from particles_to_levelset import particles_to_levelset

import os
import pyopenvdb as vdb

def convert_to_vdb(filename):
    # check if file exists
    if not os.path.isfile(filename):
        print("File does not exist:", filename)
        return
    
    (grid_size, dt, nb_time_steps, particles_per_time_step) = parse_file(filename)

    for i, particles in enumerate(particles_per_time_step):
        print("Processing time step", i)
        levelset = particles_to_levelset(particles, particle_radius=0.5, voxel_size=0.3)
        vdb.write(f"./results/levelset{i}.vdb", levelset)

convert_to_vdb("./liquidPointCloud.txt")