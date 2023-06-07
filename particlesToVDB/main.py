from parse_file import parse_file
from particles_to_levelset import particles_to_levelset

import pyopenvdb as vdb

(grid_size, dt, nb_time_steps, particles_per_time_step) = parse_file("./liquidPointCloud.txt")

for i, particles in enumerate(particles_per_time_step):
    print("Processing time step", i)
    levelset = particles_to_levelset(particles, particle_radius=1.0, voxel_size=0.5)
    vdb.write(f"./results/levelset{i}.vdb", levelset)