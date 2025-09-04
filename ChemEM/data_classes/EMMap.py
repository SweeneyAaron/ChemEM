# This file is part of the ChemEM software.
#
# Copyright (c) 2023 - Topf Group & Leibniz Institute for Virology (LIV),
# Hamburg, Germany.
#
# This module was developed by:
#   Aaron Sweeney    <aaron.sweeney AT cssb-hamburg.de>

import numpy as np
import math
import mrcfile 
from scipy.fftpack import fftn, ifftn
from scipy.ndimage import fourier_gaussian
import copy


class EMMap:
    def __init__(self,
                origin,
                apix, 
                density_map,
                resolution):
       
        self.origin = origin
        self.apix = apix
        self.density_map = density_map 
        self.resolution = resolution
        self.map_contour = 0.0
    
   
    @property 
    def box_size(self):
        # box_size = z,y,z
        return self.density_map.shape
    @property 
    def x_size(self):
        return self.box_size[2]
    @property 
    def y_size(self):
        return self.box_size[1] 
    @property
    def z_size(self):
        return self.box_size[0]
    @property 
    def x_origin(self):
        return self.origin[0]
    @property 
    def y_origin(self):
        return self.origin[1]
    @property 
    def z_origin(self):
        return self.origin[2]
    
    @property 
    def std(self):
        return np.std(self.density_map)
    
    @property 
    def mean(self):
        return np.mean(self.density_map)
    
    def voxel_to_point(self,x,y,z):
        
        # Calculate the real-world coordinates
        real_world_x = self.origin[0] + x * self.apix[0]
        real_world_y = self.origin[1] + y * self.apix[1]
        real_world_z = self.origin[2] + z * self.apix[2]
        
        return (real_world_x, real_world_y, real_world_z)
    
    def normalise(self):
        
        
        if self.std == 0:
            pass
        else:
            self.density_map = (self.density_map - self.mean) / self.std
    
    def set_map_contour(self):
        self.map_contour = 0.0
        
    
    def write_mrc(self, outfile):
        
        if not outfile.endswith('.mrc'):
            outfile += '.mrc'
        
        data = self.density_map.astype('float32')
        with mrcfile.new(outfile, overwrite=True) as mrc:
            mrc.set_data(data)
            
            mrc.header.nxstart = 0
            mrc.header.nystart = 0
            mrc.header.nzstart = 0
            
            mrc.header.mx = self.x_size
            mrc.header.my = self.y_size
            mrc.header.mz = self.z_size
            
            mrc.header.mapc = 1
            mrc.header.mapr = 2
            mrc.header.maps = 3
            
            mrc.header.cellb.alpha = 90
            mrc.header.cellb.beta = 90
            mrc.header.cellb.gamma = 90
            
            mrc.header.origin.x = self.origin[0]
            mrc.header.origin.y = self.origin[1]
            mrc.header.origin.z = self.origin[2]
            
            mrc.voxel_size = tuple(self.apix)
            
            if hasattr(self, 'ispg'):
                mrc.header.ispg = self.ispg
            if hasattr(self, 'extra1 '):
                mrc.header.extra1 = self.extra1
            if hasattr(self, 'extra2 '):
                mrc.header.extra2 = self.extra2 
            if hasattr(self, 'exttyp'):
                mrc.header.exttyp = self.exttyp
            if hasattr(self, 'extended_header'):
                mrc.set_extended_header(self.extended_header)
    
    def center_of_mass(self):
        """
        Calculate the center of mass (COM) of non-zero voxels in a 3D numpy array.
        
        Parameters:
            array (np.ndarray): 3D array with voxel values.
            origin (tuple): The origin coordinates (x, y, z) in real space.
            apix (float): The spacing between grid points (angstroms per pixel).
        
        Returns:
            tuple: Center of mass coordinates in real-space (x, y, z).
        """
        # Ensure the array is a 3D NumPy array
        if not isinstance(self.density_map, np.ndarray) or self.density_map.ndim != 3:
            raise ValueError("Input array must be a 3D NumPy array.")
        
        # Get the indices of the non-zero voxels
        indices = np.array(np.nonzero(self.density_map)).T  # shape: (N, 3) where N is the number of non-zero voxels
        
        if indices.size == 0:
            raise ValueError("The input array contains no non-zero voxels.")
        
        # Calculate the weights (values of the non-zero voxels)
        #weights = self.density_map[self.density_map > 0]
        weights = self.density_map[indices[:, 0], indices[:, 1], indices[:, 2]]
        # Compute the center of mass in grid coordinates
        com_grid = np.average(indices, axis=0, weights=weights)  # Weighted average of indices
        com_grid = com_grid[::-1] # convert to z,y,x
        # Convert grid coordinates to real-space coordinates
        com_real = np.array(self.origin) + self.apix * com_grid
        
        return tuple(com_real)
            
    def copy(self):
        return copy.deepcopy(self)
            


    
    
    
    
    
    
