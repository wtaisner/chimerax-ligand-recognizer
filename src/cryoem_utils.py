"""Adapted from https://github.com/dabrze/cryo-em-ligand-cutter/tree/main"""
from math import sqrt

import numpy as np
import scipy as sp  # type: ignore
from scipy import signal
from scipy.stats import mode, norm  # type: ignore


def create_binary_kernel(radius):
    """
    Creates a binary kernel of a given radius.

    Args:
        radius (int): The radius of the kernel.

    Returns:
        numpy.ndarray: A binary kernel of shape (2*radius+1, 2*radius+1, 2*radius+1).
    """
    boxsize = 2 * radius + 1
    kern_sphere = np.zeros(shape=(boxsize, boxsize, boxsize), dtype="float")
    kx = ky = kz = boxsize
    center = boxsize // 2

    r1 = center
    for i in range(kx):
        for j in range(ky):
            for k in range(kz):
                dist = sqrt((i - center) ** 2 + (j - center) ** 2 + (k - center) ** 2)
                if dist < r1:
                    kern_sphere[i, j, k] = 1

    return kern_sphere


def get_ligand_mask(atom_radius, unit_cell, map_array, origin, ligand_coords):
    """
    Creates a binary mask of the ligand in the given map_array.

    Args:
        atom_radius (float): The radius of the atoms in the ligand.
        unit_cell (tuple): The dimensions of the unit cell.
        map_array (numpy.ndarray): The 3D density map.
        origin (tuple): The origin of the map.
        ligand_coords (list): The coordinates of the atoms in the ligand.

    Returns:
        numpy.ndarray: A binary mask of the ligand in the given map_array.
    """
    x_ligand = np.array(ligand_coords)[:, 0] - origin[0]
    y_ligand = np.array(ligand_coords)[:, 1] - origin[1]
    z_ligand = np.array(ligand_coords)[:, 2] - origin[2]

    grid_3d = np.zeros((map_array.shape), dtype="float")
    x = x_ligand * map_array.shape[0] / unit_cell[0]
    y = y_ligand * map_array.shape[1] / unit_cell[1]
    z = z_ligand * map_array.shape[2] / unit_cell[2]

    for ix, iy, iz in zip(x, y, z):
        grid_3d[int(round(iz)), int(round(iy)), int(round(ix))] = 1.0

    pixsize = unit_cell[0] / map_array.shape[0]
    kern_rad = round(atom_radius / pixsize)
    if kern_rad > 0:
        grid2 = signal.fftconvolve(grid_3d, create_binary_kernel(kern_rad), "same")
        grid2_binary = grid2 > 1e-5
        dilate = sp.ndimage.binary_dilation(grid2_binary, iterations=1)
        mask = dilate
    else:
        mask = grid_3d

    mask = mask * (mask >= 1.0e-5)
    mask = np.where(grid2_binary, 1.0, mask)
    shift_z = origin[0]
    shift_y = origin[1]
    shift_x = origin[2]
    mask = np.roll(
        np.roll(np.roll(mask, -shift_z, axis=0), -shift_y, axis=1),
        -shift_x,
        axis=2,
    )

    return mask


def get_mask_bounding_box(masked_array):
    """
    Returns the bounding box of a masked array, defined as the minimum and maximum indices of nonzero elements
    in each dimension.

    Args:
        masked_array (numpy.ndarray): A 3D numpy array with boolean values indicating the masked voxels.

    Returns:
        tuple: A tuple containing the minimum and maximum indices of nonzero elements in each dimension, in the
        following order: (min_x, max_x, min_y, max_y, min_z, max_z).
    """
    nonzero_indices = np.nonzero(masked_array)

    min_x, max_x = np.min(nonzero_indices[0]), np.max(nonzero_indices[0])
    min_y, max_y = np.min(nonzero_indices[1]), np.max(nonzero_indices[1])
    min_z, max_z = np.min(nonzero_indices[2]), np.max(nonzero_indices[2])

    return min_x, max_x, min_y, max_y, min_z, max_z


def resample_blob(blob, target_voxel_size, unit_cell, map_array):
    """
    Resamples a given blob to a target voxel size using the provided unit cell and map array.

    Args:
        blob (numpy.ndarray): The blob to be resampled.
        target_voxel_size (float): The target voxel size (in Angstroms).
        unit_cell (tuple): The unit cell dimensions (in Angstroms).
        map_array (numpy.ndarray): The map array.

    Returns:
        numpy.ndarray: The resampled blob.
    """
    blob = sp.ndimage.zoom(
        blob,
        [
            unit_cell[0] / target_voxel_size / map_array.shape[0],
            unit_cell[1] / target_voxel_size / map_array.shape[1],
            unit_cell[2] / target_voxel_size / map_array.shape[2],
        ],
        prefilter=False,
    )

    return blob


def get_sphere_volume(radius):
    """
    Calculates the volume of a sphere given its radius.

    Args:
        radius (float): The radius of the sphere.

    Returns:
        float: The volume of the sphere.
    """
    return 4.0 / 3.0 * 3.14 * (radius**3)


def get_blob_volume(voxel_count, voxel_size):
    """
    Calculates the volume of a blob given the number of voxels and the size of each voxel.

    Args:
        voxel_count (int): The number of voxels in the blob.
        voxel_size (float): The size of each voxel in angstroms.

    Returns:
        float: The volume of the blob in cubic angstroms.
    """
    return voxel_count * (voxel_size**3)


MAP_VALUE_MAPPER = {
    1.0: 0.66,
    1.1: 0.63,
    1.2: 0.57,
    1.3: 0.57,
    1.4: 0.54,
    1.5: 0.50,
    1.6: 0.48,
    1.7: 0.44,
    1.8: 0.42,
    1.9: 0.39,
    2.0: 0.36,
    2.1: 0.33,
    2.2: 0.31,
    2.3: 0.30,
    2.4: 0.28,
    2.5: 0.25,
    2.6: 0.25,
    2.7: 0.23,
    2.8: 0.21,
    2.9: 0.21,
    3.0: 0.20,
    3.1: 0.18,
    3.2: 0.18,
    3.3: 0.17,
    3.4: 0.15,
    3.5: 0.16,
    3.6: 0.14,
    3.7: 0.12,
    3.8: 0.14,
    3.9: 0.15,
    4.0: 0.17,
}


def extract_ligand(
    density_threshold,
    min_blob_radius,
    atom_radius,
    target_voxel_size,
    resolution,
    res_cov_threshold,
    blob_cov_threshold,
    padding,
    unit_cell,
    map_array,
    origin,
    ligand_coords,
):
    """
    Extracts a ligand blob from a given map array and saves it as a compressed numpy file.

    Args:
        density_threshold (float): Density threshold for the ligand blob.
        min_blob_radius (float): Minimum radius of the ligand blob.
        atom_radius (float): Radius of the atoms in the ligand.
        target_voxel_size (float): Target voxel size for the resampled blob.
        resolution (float): Resolution of the map.
        res_cov_threshold (float): Minimum coverage threshold for the model.
        blob_cov_threshold (float): Minimum coverage threshold for the blob.
        padding (int): Padding size for the blob.
        unit_cell (np.ndarray): Unit cell dimensions.
        map_array (np.ndarray): Map array.
        origin (np.ndarray): Origin of the map.
        ligand_coords (np.ndarray): Coordinates of the ligand.

    Returns:
        np.ndarray: cut ligand or None
    """
    mask = get_ligand_mask(atom_radius, unit_cell, map_array, origin, ligand_coords)
    min_x, max_x, min_y, max_y, min_z, max_z = get_mask_bounding_box(mask)
    blob = map_array * mask
    blob = blob[
        min_x - padding : max_x + 1 + padding,
        min_y - padding : max_y + 1 + padding,
        min_z - padding : max_z + 1 + padding,
    ]

    # Resampling the map to target voxel size
    blob = resample_blob(blob, target_voxel_size, unit_cell, map_array)
    blob[blob < density_threshold] = 0
    blob_volume = get_blob_volume(np.sum(blob != 0), target_voxel_size)

    if blob_volume >= get_sphere_volume(min_blob_radius):
        fragment_mask = mask[
            min_x - padding : max_x + 1 + padding,
            min_y - padding : max_y + 1 + padding,
            min_z - padding : max_z + 1 + padding,
        ]
        fragment_mask = resample_blob(
            fragment_mask, target_voxel_size, unit_cell, map_array
        )
        res_voxels = fragment_mask > 0
        blob_voxels = blob > 0
        res_cov_frac = np.sum(res_voxels & blob_voxels) / np.sum(res_voxels)
        blob_cov_frac = np.sum(res_voxels & blob_voxels) / np.sum(blob_voxels)

        if res_cov_frac >= res_cov_threshold and blob_cov_frac >= blob_cov_threshold:
            # rescaling density values
            blob = blob * (MAP_VALUE_MAPPER[resolution] / blob[blob > 0].min())

            # print(
            #     f"Dimensions: {blob.shape}, Blob min value: {blob[blob > 0].min():.3f}, "
            #     + f"Blob max value: {blob.max():.3f}, Non-zero: {np.sum(blob != 0):,}, "
            #     + f"Zero: {np.sum(blob == 0):,}, NA count: {np.sum(np.isnan(blob)):,}, "
            #     + f"Blob volume: {blob_volume:.3f}, Model coverage: {res_cov_frac:.2f}, "
            #     + f"Blob coverage: {blob_cov_frac:.2f}"
            # )
            return blob
        else:
            # print(
            #     f"Model coverage: {res_cov_frac:.2f}, "
            #     + f"Blob coverage: {blob_cov_frac:.2f}. Not enough coverage."
            # )
            return None
    else:
        print(f"Not enough density.")
        return None


def read_map(map_model):
    file_header = map_model.data.file_header

    map_array = np.asarray(map_model.full_matrix(), dtype="float")
    # chimerax reads the ccp4 file in 'right' order, so we don't need to change it

    unit_cell = np.zeros(6, dtype="float")
    cell = np.array(
        [(file_header["xlen"], file_header["ylen"], file_header["zlen"])],
        dtype=[("x", "<f4"), ("y", "<f4"), ("z", "<f4")],
    )
    unit_cell[:3] = cell.astype([("x", "<f4"), ("y", "<f4"), ("z", "<f4")]).view(
        ("<f4", 3)
    )

    # swapping a and c to compatible with ZYX convension
    unit_cell[0], unit_cell[2] = unit_cell[2], unit_cell[0]
    unit_cell[3:] = float(90)
    origin = [
        1 * file_header["ncstart"],  # nxstart
        1 * file_header["nrstart"],  # nystart
        1 * file_header["nsstart"],  # nzstart
    ]

    return unit_cell, map_array, origin


def em_stats(cif_model):
    try:
        em_3d_keys = cif_model.metadata["em_3d_reconstruction"][1:]
        em_3d_values = cif_model.metadata["em_3d_reconstruction data"]
        em_3d_dict = dict(zip(em_3d_keys, em_3d_values))
    except Exception:
        em_3d_dict = {}
    try:
        resolution = round(float(em_3d_dict["resolution"]), 1)
    except Exception:
        resolution = None
    try:
        num_particles = int(em_3d_dict["num_particles"])
    except Exception:
        num_particles = None

    if resolution is not None:
        if resolution > 4.0:
            resolution = 4.0
        elif resolution < 1.0:
            resolution = 1.0

    return resolution, num_particles


def extract_ligand_coords(cif_model, residue):
    resolution, num_particles = em_stats(cif_model)
    ligand_coords = np.array([atom.coord for atom in residue.atoms], dtype=np.float32)

    return ligand_coords, resolution, num_particles


def cut_ligand(
    map_model,
    cif_model,
    residue,
    density_std_threshold=2.8,
    min_blob_radius=0.8,
    atom_radius=1.5,
    target_voxel_size=0.2,
    res_cov_threshold=0.02,
    blob_cov_threshold=0.01,
    padding=2,
):
    ligand_coords, resolution, num_particles = extract_ligand_coords(cif_model, residue)

    if resolution is None or num_particles is None:
        resolution = 3.0

    unit_cell, map_array, origin = read_map(map_model)

    map_median = np.median(map_array)
    map_std = np.std(map_array)
    value_mask = (map_array < map_median - 0.5 * map_std) | (
        map_array > map_median + 0.5 * map_std
    )

    quantile_threshold = norm.cdf(density_std_threshold)
    density_threshold = np.quantile(map_array[value_mask], quantile_threshold)

    blob = extract_ligand(
        density_threshold,
        min_blob_radius,
        atom_radius,
        target_voxel_size,
        resolution,
        res_cov_threshold,
        blob_cov_threshold,
        padding,
        unit_cell,
        map_array,
        origin,
        ligand_coords,
    )

    return blob
