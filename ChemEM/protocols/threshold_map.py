#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 12:42:06 2024

@author: aaron.sweeney
"""

import numpy as np 
import os
from scipy.ndimage import binary_dilation, gaussian_filter
from ChemEM.tools.rdtools import RDTools
from scipy.stats import norm

class ConfidenceMap:
    def __init__(self, system):
        self.system = system 
        self.radius_in_angstrom = 6.0
        self.use_model = False
        
    def get_position_input(self):
        
        
        self.atoms = [i for i in self.protein_openff_structure.atoms if i.element > 1]
        self.positions = np.array([[i.xx, i.xy, i.xz] for i in self.atoms ])
        
        self.atom_radii = np.array( [RDTools.get_van_der_waals_radius(i.element_name) for i in self.atoms] )
        #pdb.set_trace()
    
    
    def get_protein_atom_radius(self):
        #improve this to include the resolution and atomm element
        self.voxel_radius = int(np.ceil(self.radius_in_angstrom / self.densmap.apix[0]))
        
    def get_grid_coords(self):
        self.grid_coords = coords_to_grid(self.positions,
                         self.densmap.origin,
                         self.densmap.apix) #x,y,z
    
    def get_protein_structure(self):
        #TODO! system ligands!!
        self.protein_openff_structure =  self.system.proteins[0].openff_structure
        self.densmap = self.system.maps[0]    
    
    def get_model_mask(self):
        
        
        model_mask = create_mask_from_model(self.densmap.density_map.shape, #z,y,x
                                                 self.grid_coords,
                                                 self.voxel_radius)
        
        self.model_mask = refine_mask(model_mask)
    
    def get_background_by_mask(self):
        self.bg_mean, self.bg_std = estimate_background_stats(self.densmap.density_map, self.model_mask)
    
    def estimate_background_from_cubes(self):
        self.bg_mean, self.bg_std  = estimate_background_distribution(self.densmap.density_map, cube_size=10, num_cubes=4)
    
    def compute_confidence_map(self, cube_size=10, num_cubes=4, desired_ppv=0.99):
        """
        Compute a confidence map from the given density map following the described procedure.
        """
       
        
        # 2. Compute p-values for each voxel
        pvals = compute_p_values(self.densmap.density_map, self.bg_mean, self.bg_std)
        
        # 3. Apply Benjamini-Yekutieli to get q-values
        qvals = benjamini_yekutieli(pvals)
        
        # 4. Convert q-values to PPVs: PPV = 1 - q
        ppvs = 1 - qvals
        
        # 5. Threshold the PPV map based on desired PPV (e.g., 0.99)
        # This thresholded version can be visualized. For general visualization,
        # you might store the PPV map to a file format supported by your visualization software.
        thresholded_map = (ppvs >= desired_ppv).astype(np.float32)
        
        
        outpath = os.path.join(self.system.output, 'fdr_maps')
        try:
            os.mkdir(outpath)
        except FileExistsError:
            pass
            
        confidence_map = self.densmap.copy() 
        confidence_map.density_map = ppvs 
        out_file = os.path.join(outpath, 'confidence_map.mrc')
        confidence_map.write_mrc(out_file)
        
        masked_map = self.densmap.copy() 
        masked_map.density_map = masked_map.density_map * thresholded_map
        out_file = os.path.join(outpath, 'masked_map.mrc')
        masked_map.write_mrc(out_file)
        self.system.confidence_map = masked_map
        #return ppvs, thresholded_map
    
    def run(self):
        self.get_protein_structure() 
        
        
        if self.use_model:
            self.get_position_input()
            self.get_grid_coords()
            self.get_protein_atom_radius()
            self.get_model_mask() #replace this with a blurred map!!
            self.get_background_by_mask()
        else:
            self.estimate_background_from_cubes()
        
        self.compute_confidence_map()
        
        


def compute_p_values(density_map, bg_mean, bg_std):
    """
    Compute one-sided p-values for each voxel assuming Gaussian noise distribution.
    p-value = P(X >= value) under N(bg_mean, bg_std^2).
    """
    # Z-score each voxel
    z_map = (density_map - bg_mean) / (bg_std + 1e-12)
    
    # p-values = 1 - CDF(z)
    p_values = 1 - norm.cdf(z_map)
    return p_values


def benjamini_yekutieli(pvals, alpha=0.05):
    """
    Benjamini-Yekutieli procedure for FDR control.
    This method handles arbitrary dependencies.
    
    Returns array of q-values.
    """
    # Flatten p-values and sort
    pvals_flat = pvals.flatten()
    n = len(pvals_flat)
    sort_indices = np.argsort(pvals_flat)
    sorted_pvals = pvals_flat[sort_indices]
    
    # BY correction factor for arbitrary dependency
    # harmonic_n = sum_{i=1}^n (1/i)
    harmonic_n = np.sum(1.0 / np.arange(1, n+1))
    
    # Compute thresholded q-values
    qvals = np.zeros_like(sorted_pvals)
    min_ratio = 1.0
    for i in reversed(range(n)):
        current = sorted_pvals[i] * n * harmonic_n / (i+1)
        if current < min_ratio:
            min_ratio = current
        qvals[i] = min_ratio
    
    # Revert to original ordering
    qvals_original = np.empty_like(qvals)
    qvals_original[sort_indices] = qvals
    return qvals_original.reshape(pvals.shape)

def estimate_background_distribution(density_map, cube_size=10, num_cubes=4):
    """
    Estimate the background noise distribution parameters (mean, std).
    For simplicity, we pick `num_cubes` regions along the principal axes far from the particle center.
    Assumes the particle is roughly centered.
    """
    nx, ny, nz = density_map.shape
    half_x, half_y, half_z = nx//2, ny//2, nz//2
    
    # Define cube corners far outside the center. For instance:
    # 1) along +x direction, 2) -x direction, 3) +y direction, 4) -y direction.
    # Adjust these based on known map orientation and empty regions.
    coords = [
        (nx - cube_size, nx, half_y, half_y+cube_size, half_z, half_z+cube_size),
        (0, cube_size, half_y, half_y+cube_size, half_z, half_z+cube_size),
        (half_x, half_x+cube_size, ny - cube_size, ny, half_z, half_z+cube_size),
        (half_x, half_x+cube_size, 0, cube_size, half_z, half_z+cube_size),
    ]
    
    # Extract background voxels
    background_voxels = []
    for (xstart, xend, ystart, yend, zstart, zend) in coords:
        background_voxels.append(density_map[xstart:xend, ystart:yend, zstart:zend].ravel())
    background_voxels = np.concatenate(background_voxels)
    
    # Estimate mean and std from these voxels
    bg_mean = np.mean(background_voxels)
    bg_std = np.std(background_voxels)
    return bg_mean, bg_std

def coords_to_grid(coords, origin, voxel_size):
    """
    Convert atomic coordinates in Å to voxel indices in the map.
    
    Arguments:
    - coords: Nx3 array of atomic coordinates in Å.
    - origin: The origin of the map in Å (x0, y0, z0) corresponding to map[0,0,0].
    - voxel_size: The size of each voxel in Å/pixel (scalar or tuple).
    
    Returns:
    - grid_coords: Nx3 array of integer voxel indices.
    """
    # If voxel_size is scalar, make it a tuple (vx, vy, vz)
    if np.isscalar(voxel_size):
        voxel_size = (voxel_size, voxel_size, voxel_size)
        
    grid_coords = (coords - origin) / np.array(voxel_size)
    return np.round(grid_coords).astype(int)

def create_mask_from_model(density_shape, atom_indices, radius_pixels):
    """
    Create a binary mask for the protein by marking voxels within a certain 
    radius (in pixels) of any atom position.
    
    Arguments:
    - density_shape: shape of the 3D map (nz, ny, nx)
    - atom_indices: Nx3 array of voxel indices of atoms (v,y,z)
    - radius_pixels: the radius around each atom (in voxel units) to consider inside
    
    Returns:
    - mask: 3D boolean array, True = inside protein envelope
    """
    mask = np.zeros(density_shape, dtype=bool)
    
    # A simple approach: for each atom, mark a sphere of radius_pixels as True.
    # For large maps, this might be slow. Consider more efficient methods or spatial data structures.
    
    nz, ny, nx = density_shape
    rr = radius_pixels**2
    for x_a, y_a, z_a in atom_indices:
        # Determine bounding box for the sphere
        x_min = max(x_a - radius_pixels, 0)
        x_max = min(x_a + radius_pixels + 1, nx)
        y_min = max(y_a - radius_pixels, 0)
        y_max = min(y_a + radius_pixels + 1, ny)
        z_min = max(z_a - radius_pixels, 0)
        z_max = min(z_a + radius_pixels + 1, nz)
        
        # Compute distances within this sub-box
        xs, ys, zs = np.ogrid[z_min:z_max, y_min:y_max, x_min:x_max]
        dist_sq = (xs - x_a)**2 + (ys - y_a)**2 + (zs - z_a)**2
        mask[z_min:z_max, y_min:y_max, x_min:x_max] |= (dist_sq <= rr)
    
    return mask


def refine_mask(mask, dilation_iterations=2, smooth_sigma=1.0):
    """
    Refine the mask by dilating and smoothing, which can help ensure the entire molecular
    envelope is captured. 
    
    - binary_dilation expands the mask slightly
    - gaussian_filter on a float version and thresholding can create a smoother boundary
    
    Returns a refined boolean mask.
    """
    refined = binary_dilation(mask, iterations=dilation_iterations)
    # Smooth the refined mask by converting to float, applying Gaussian, and thresholding
    smoothed = gaussian_filter(refined.astype(float), sigma=smooth_sigma)
    # Threshold at 0.5 to get a cleaner envelope
    refined_mask = (smoothed > 0.5)
    return refined_mask

def estimate_background_stats(density_map, protein_mask):
    """
    Estimate background noise statistics from voxels outside the protein mask.
    """
    outside_voxels = density_map[~protein_mask].ravel()
    bg_mean = np.mean(outside_voxels)
    bg_std = np.std(outside_voxels)
    return bg_mean, bg_std
    
    
