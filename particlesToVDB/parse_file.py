# file format:

# grid_size_x grid_size_y grid_size_z
# dt nb_time_steps
# nb_particles
# x y z
# ...

import numpy as np

def parse_file(filename: str) -> tuple[np.ndarray, float, int, list[np.ndarray]]:
    with open(filename, 'r') as f:
        print(f"Reading file: {filename}")
        lines = f.readlines()
        
        grid_size = np.array([int(x) for x in lines[0].split()])
        print(f"grid_size: {grid_size}")

        dt, nb_time_steps = [float(x) for x in lines[1].split()]
        print(f"dt: {dt}, nb_time_steps: {nb_time_steps}")

        nb_particles = int(lines[2])
        print(f"nb_particles: {nb_particles}")

        all_particles = np.array([[float(x) for x in line.split()] for line in lines[3:]])
        print(f"particles: {all_particles}")

        particles_per_time_step = np.array_split(all_particles, nb_time_steps)

        return grid_size, dt, nb_time_steps, particles_per_time_step