import pyopenvdb as vdb
from scipy.spatial import cKDTree
import numpy as np

def kernel(sqr_d, sqr_h):
        return np.maximum(0, (1 - sqr_d / sqr_h)**3)

def particles_to_levelset(particles: np.ndarray, particle_radius: float, voxel_size: float) -> vdb.FloatGrid:
    """
    creates a density grid from a set of particles
    based on https://github.com/isl-org/DeepLagrangianFluids/
    """
    if not particle_radius >= voxel_size:
        raise ValueError(f"Error: particle_radius ({particle_radius}) must exceed voxel_size ({voxel_size})")
    if not voxel_size > 0:
        raise ValueError(f"Error: voxel_size must be > 0 but got {voxel_size}")

    point_radius = particle_radius / voxel_size
    points = particles / voxel_size
    tree = cKDTree(points)

    grid = vdb.FloatGrid()
    accessor = grid.getAccessor()

    visited_indices = set()
    unvisited_indices = set()
    visited_points = np.zeros(shape=points.shape[0:1], dtype=np.uint8)

    def compute_value(ijk):
        nn = tree.query_ball_point(ijk, point_radius)
        if nn:
            sqr_dist = np.sum((points[nn] - np.asarray(ijk))**2, axis=-1)
            # use negative values to get normals pointing outwards
            value = np.sum(kernel(sqr_dist, point_radius**2))
            visited_points[nn] = 1
        else:
            value = 0
        accessor.setValueOn(ijk, value=value)
        return value

    # 6-connected 3D grid
    neighbors = np.array([
        [-1, 0, 0],
        [1, 0, 0],
        [0, 1, 0],
        [0, -1, 0],
        [0, 0, 1],
        [0, 0, -1],
    ], dtype=np.int32)

    while np.count_nonzero(visited_points) != len(visited_points):
        unvisited_points_mask = visited_points == 0
        idx = tuple(np.round(points[unvisited_points_mask][0]).astype(np.int32))
        visited_points[np.argwhere(unvisited_points_mask)[0]] = 1
        unvisited_indices.add(idx)

        while len(unvisited_indices):
            idx = unvisited_indices.pop()
            visited_indices.add(idx)
            value = compute_value(idx)
            if value:
                new_indices = neighbors + idx
                for x in new_indices:
                    i = tuple(x)
                    if not i in visited_indices:
                        unvisited_indices.add(i)

    grid.signedFloodFill()

    return grid
