# This file is part of the ChemEM software.
#
# Copyright (c) 2023 - Topf Group & Leibniz Institute for Virology (LIV),
# Hamburg, Germany.
#
# This module was developed by:
#   Aaron Sweeney    <aaron.sweeney AT cssb-hamburg.de>

import numpy as np 
from scipy.ndimage import measurements
from scipy.ndimage.morphology import generate_binary_structure
from scipy.spatial import distance_matrix 
from scipy.ndimage import grey_opening, grey_closing
from scipy.ndimage import gaussian_filter
from ChemEM.tools.map_tools import MapTools
from ChemEM.tools.rdtools import RDTools
from ChemEM.messages import Messages
import os
from skimage import filters
from scipy.ndimage import label
from scipy.ndimage.interpolation import map_coordinates



class DifferenceMapOriginal:
    
    def __init__(self, system):
        self.system = system 
        self.t = 1.5 
        self.scale  = True
        self.shell_width = 0.02
        self.filt = True
        self.contour_mask = True
        self.dust = True
        self.softdust = True 
        self.randsize=0.2
        self.apix=None
        
    def get_data(self):
        self.map = self.system.segment.maps[0]
        self.blur_model()
       
    
    def blur_model(self):
        
        mol_map_input = []
        for protein in self.system.segment._protein_mols:
            mol_map_input += MapTools._prepare_molmap_input(protein)
        
        self.mol_map = MapTools.gaussian_blur(mol_map_input,
                                              self.map.resolution,
                                              self.map)
        
        
    def difference_map_input(self):
        c1 = MapTools.map_contour(self.map, t = self.t)
        
        if self.map.resolution > 20.0 :
            mt = 2.0
        elif self.map.resolution > 10.0:
            mt = 1.0
        elif self.map.resolution > 6.0:
            mt = 0.5
        else:
            mt = 0.1
         
        c2 = MapTools.model_contour(self.mol_map, t = mt)
         
        refsc = True
        self.c1 = c1 
        self.c2 = c2 
        self.refsc = refsc
        
    
    def _alignment_box(self, self_map, map2, s):
        (ox, oy, oz) = (
            self_map.origin[0],
            self_map.origin[1],
            self_map.origin[2]
        )
        (o1x, o1y, o1z) = (
            map2.origin[0],
            map2.origin[1],
            map2.origin[2]
        )
        offset = (
            o1x - ox,
            o1y - oy,
            o1z - oz
        )
        (m1x, m1y, m1z) = (
            ox + self_map.density_map.shape[2] * self_map.apix[0],
            oy + self_map.density_map.shape[1] * self_map.apix[1],
            oz + self_map.density_map.shape[0] * self_map.apix[2]
        )
        (m2x, m2y, m2z) = (
            o1x + map2.density_map.shape[2] * map2.apix[0],
            o1y + map2.density_map.shape[1] * map2.apix[1],
            o1z + map2.density_map.shape[0] * map2.apix[2]
        )
        (nx, ny, nz) = (o1x, o1y, o1z)
        if offset[0] > 0:
            nx = ox
        if offset[1] > 0:
            ny = oy
        if offset[2] > 0:
            nz = oz

        (lz, ly, lx) = (
            (m2z-nz) / float(s),
            (m2y-ny) / float(s),
            (m2x-nx) / float(s)
        )
        
        if m2x < m1x:
            lx = (m1x - nx) / float(s)
        if m2y < m1y:
            ly = (m1y - ny) / float(s)
        if m2z < m1z:
            lz = (m1z - nz) / float(s)
        
        gridshape = (int(lz), int(ly), int(lx))
        new_origin = (nx, ny, nz)
        return gridshape, new_origin


    def generate_map(self):
        
        
        c1 = (self.c1 - self.map.density_map.min())
        c2 = (self.c2 - self.mol_map.density_map.min())
        
        
        dens_map = self.map.density_map.copy() - self.map.density_map.min()
        mol_map = self.mol_map.density_map.copy() - self.mol_map.density_map.min()
        
        dens_map_copy = dens_map.copy()
        mol_map_copy = mol_map.copy()
        
        #-----new-----
        if sum( self.map.apix) > sum(self.mol_map.apix):
            spacing =  self.dens_map.apix
        else:
            spacing = self.mol_map.apix
        
        grid_shape, new_ori =self._alignment_box(self.map, self.mol_map, spacing[0])
        #------end new
        
        
        if self.scale:
            dens_map_copy, mol_map_copy = self.amplitude_match(dens_map, mol_map,0,0,self.shell_width,0,0, lpfiltb=self.filt, ref=self.refsc)
        
        
        
        #-----new
        #resample scaled maps to common grid    
        if self.apix == None:
            spacing = [self.map.resolution*0.33, self.map.resolution*0.33, self.map.resolution*0.33]
        else:
            spacing = self.apix
        
        apix_ratio = self.map.apix/spacing
        
        diff1 = self._interpolate_to_grid(self.map, dens_map_copy, grid_shape, spacing, new_ori, 1)
        diff2 = self._interpolate_to_grid(self.mol_map, mol_map_copy, grid_shape, spacing, new_ori, 1)
        
        
        #diff1 = dens_map_copy.copy()
        #diff2 = mol_map_copy.copy()
        
        dens_map_copy = (dens_map > c1) *1
        mol_map_copy = (mol_map > c2)*1

        mask1 = self._interpolate_to_grid(self.map, dens_map_copy , grid_shape, spacing, new_ori, 1, 'zero')
        mask2 = self._interpolate_to_grid(self.map, mol_map_copy, grid_shape, spacing, new_ori, 1, 'zero')
        
        
        mask1.density_map = mask1.density_map > 0.1 
        mask2.density_map = mask2.density_map > 0.1 
        
        
        
        min1 = diff1.density_map.min()
        min2 = diff2.density_map.min()
        min_scaled_maps = min(min1,min2)
        
        diff1.density_map = diff1.density_map - min_scaled_maps
        diff2.density_map = diff2.density_map - min_scaled_maps
        
        
        
        min1 = np.amin(diff1.density_map[mask1.density_map])
        diffc1 = min1+0.10*(np.amax(diff1.density_map)-min1)
        
        min2 = np.amin(diff2.density_map[mask2.density_map])
        diffc2 = min2+0.10*(np.amax(diff2.density_map)-min2)
        
        
        #diffmap = diff1.copy()
        diff1.density_map = diff1.density_map - diff2.density_map
        
        
        #apply mask to difference maps
        if self.contour_mask:
            diff1.density_map = diff1.density_map * mask1.density_map
            
            if self.dust:
                if self.softdust:
                    #a = 1
                    diffc1 = min1 + 0.5 * (diffc1 - min1)
                    
                diff1.density_map =self.label_patches(diff1.density_map ,diffc1, prob=self.randsize)[0]
                
        
        
        #spacing = np.array([self.map.resolution * 0.33, self.map.resolution * 0.33, self.map.resolution * 0.33])
        
        #final_mask = self._interpolation(diff1, spacing, self.map.apix, self.map.origin)
        
        #print('D1 diff c1:',diffc1)
        #new_map =  self._interpolate_to_grid(self.map, diff1.density_map , self.map.density_map.shape, self.map.apix, self.map.origin, 1, 'zero')
        new_map =  self._interpolate_to_grid(diff1, diff1.density_map , self.map.density_map.shape, self.map.apix, self.map.origin, 1, 'zero')
        

        #new_map = self.map.copy()
        #new_map.density_map = diff1
        
        self.system.difference_maps.append(new_map)
    
    
    def _interpolate_to_grid(  # noqa:F811
            self,
            self_map,
            arr,
            gridshape,
            s,
            ori,
            order_spline=3,
            fill='min'
    ):
        """
        Spline interpolation to a grid.

        Arguments:
            *gridshape*
                shape of new grid array (z,y,x)
            *s*
                new grid spacing
            *ori*
                origin of the new grid
            *order_spine*
                order of the spline used for interpolation

        Return:
            Interpolated map
        """
        (ox, oy, oz) = (
            self_map.origin[0],
            self_map.origin[1],
            self_map.origin[2],
        )
        (o1x, o1y, o1z) = (
            float(ori[0]),
            float(ori[1]),
            float(ori[2])
        )
        scale = float(s[0]) / self_map.apix[0]
        offset = (o1x - ox, o1y - oy, o1z - oz)

        new_map_origin = (o1x, o1y, o1z)
        grid_indices = np.indices(gridshape)
        z_ind = grid_indices[0]
        z_ind.ravel()
        y_ind = grid_indices[1]
        y_ind.ravel()
        x_ind = grid_indices[2]
        x_ind.ravel()
        z_ind = ((offset[2]) / self_map.apix[0]) + scale * z_ind
        y_ind = ((offset[1]) / self_map.apix[1]) + scale * y_ind
        x_ind = ((offset[0]) / self_map.apix[2]) + scale * x_ind
        
        
        filtered_array = arr
        
        
        if fill == 'zero':
            fillval = 0.0
        else:
            fillval = arr.min()
            
        new_array = map_coordinates(
            filtered_array,
            [z_ind, y_ind, x_ind],
            cval=fillval,
            order=order_spline,
            prefilter=False,
        )
        
        new_map = self_map.copy()
        new_map.origin = new_map_origin 
        new_map.apix = s 
        new_map.density_map = new_array.reshape(gridshape)
        
        return new_map
    


    def _interpolation(self, arr, spacing, apix, origin):
        
        (ox, oy, oz) = (
            origin[0],
            origin[1],
            origin[2],
        )
        
        (o1x, o1y, o1z) = (
            float(origin[0]),
            float(origin[1]),
            float(origin[2])
        )
        
        scale = apix / spacing
        offset = (o1x - ox, o1y - oy, o1z - oz)
        
        
        gridshape = arr.shape
        grid_indices = np.indices(gridshape)
        z_ind = grid_indices[0]
        z_ind.ravel()
        y_ind = grid_indices[1]
        y_ind.ravel()
        x_ind = grid_indices[2]
        x_ind.ravel()
        z_ind = ((offset[2]) / spacing[2]) + scale[2] * z_ind
        y_ind = ((offset[1]) / spacing[1]) + scale[1] * y_ind
        x_ind = ((offset[0]) / spacing[0]) + scale[0] * x_ind
        
        new_array = map_coordinates(
            arr,
            [z_ind, y_ind, x_ind],
            cval=0.0,
            order=1,
            prefilter=False,
        )
        
        return new_array
    
    def grid_footprint(self):
        a = np.zeros((3, 3, 3))
        a[1, 1, 1] = 1
        a[0, 1, 1] = 1
        a[1, 0, 1] = 1
        a[1, 1, 0] = 1
        a[2, 1, 1] = 1
        a[1, 2, 1] = 1
        a[1, 1, 2] = 1

        return a
    
    def label_patches(self, dens_map, c, prob=0.2):
        
        fp = self.grid_footprint()
        bin_map = dens_map > float(c)
        label_array, labels = measurements.label(
            dens_map * bin_map,
            structure=fp
        )
        
        sizes = measurements.sum(bin_map, label_array, range(labels + 1))
        if labels < 10:
            m_array = sizes < 0.05 * sizes.max()
            ct_remove = np.sum(m_array)
            remove_points = m_array[label_array]
            label_array[remove_points] = 0
            
            return (
                (label_array > 0) * (dens_map * bin_map),
                labels - ct_remove + 1
            )
        
        
        freq, bins = np.histogram(sizes[1:], 20)
        m_array = np.zeros(len(sizes))
        ct_remove = 0
        for i in range(len(freq)):
            fr = freq[i]
            s2 = bins[i + 1]
            s1 = bins[i]
            p_size = float(fr) / float(np.sum(freq))
            if s2 < 10 or p_size > prob:
                m_array = m_array + ((sizes >= s1) & (sizes < s2))
                ct_remove += 1
        m_array = m_array > 0
        remove_points = m_array[label_array]
        label_array[remove_points] = 0
        return (
            (label_array > 0) * (dens_map * bin_map),
            labels - ct_remove
        )
    
    
    def amplitude_match(self, arr1, arr2,
                        shellmin,
                        shellmax,
                        step=0.005,
                        c1=0,
                        c2=0,
                        lpfiltb=False,
                        lpfilta=False,
                        ref=False):
        
        
        
        ft1 = self.fft_array(arr1) #density map
        ft2 = self.fft_array(arr2) #mol map
        
        
       
        if self.map.resolution is not None:
            cutoff1 = self.map.apix[0] / self.map.resolution 
            cutoff2 = self.mol_map.apix[0] / self.map.resolution 
        
            if lpfiltb and not lpfilta:
                self.tanh_lowpass(ft1, cutoff1, fall=0.2)
                self.tanh_lowpass(ft2, cutoff2, fall=0.2)
            
        
        dist1 = self.make_fourier_shell(arr1) / self.map.apix[0]
        dist2 = self.make_fourier_shell(arr2) / self.mol_map.apix[0]
        
        ft1_avg = []
        ft2_avg = []
        ft1_avg_new = []
        lfreq = []
        # select max spatial frequency to iterate to. low resolution map
        maxlevel = 0.5 / np.max((self.map.apix[0],  self.mol_map.apix[0]), axis=0) #change second density map to molmap
        nc = 0
        x = 0.0
        highlevel = x + step
        while (x < maxlevel):
            
            # print x,highlevel, maxlevel
            # indices between upper and lower shell bound
            fshells1 = ((dist1 < min(maxlevel, highlevel)) & (dist1 >= x))
            
            # radial average
            shellvec1 = ft1[fshells1]
            
            # indices between upper and lower shell bound
            fshells2 = ((dist2 < min(maxlevel, highlevel)) & (dist2 >= x))
            
            # radial average
            shellvec2 = ft2[fshells2]
            abs1 = abs(shellvec1)
            abs2 = abs(shellvec2)
            ns1 = len(np.nonzero(abs1)[0])
            ns2 = len(np.nonzero(abs2)[0])
            
            
            
            if ns1 < 10 or ns2 < 10:
                nc += 1
                highlevel = min(maxlevel, x + (nc + 1) * step)
                x = max(0.0, x - nc * step)
                
                continue
            else:
                nc = 0
            
            mft1 = np.mean(abs1)
            mft2 = np.mean(abs2)
            if mft1 == 0.0 and mft2 == 0.0:
                continue
            # sq of radial avg amplitude
            ft1_avg.append(np.log10(np.mean(np.square(abs1))))
            ft2_avg.append(np.log10(np.mean(np.square(abs2))))

            # scale to amplitudes of the ref map
            if ref:
                if mft1 == 0.0:
                    continue
                ft1[fshells1] = shellvec1 * (mft2 / mft1)
            else:
                # replace with avg amplitudes for the two maps
                ft1[fshells1] = shellvec1 * (mft2 + mft1) / (2 * mft1)
                ft2[fshells2] = shellvec2 * (mft2 + mft1) / (2 * mft2)

            # new radial average (to check)
            mft1 = np.mean(abs(ft1[fshells1]))
            ft1_avg_new.append(
                np.log10(
                    np.mean(
                        np.square(abs(ft1[fshells1]))
                    )
                )
            )
            lfreq.append(highlevel)

            sampling_frq = highlevel

            cutoff_freq = min((1.0 / self.map.resolution) + 0.25, maxlevel)

            # scale the rest and break after relevant frequencies
            if sampling_frq > cutoff_freq:
                fshells1 = (dist1 >= highlevel)
                shellvec1 = ft1[fshells1]
                mft1 = np.mean(abs(shellvec1))
                fshells2 = (dist2 >= highlevel)
                shellvec2 = ft2[fshells2]
                mft2 = np.mean(abs(shellvec2))
                if mft1 == 0.0 and mft2 == 0.0:
                    break
                ft1_avg.append(np.log10(np.mean(np.square(abs(shellvec1)))))
                ft2_avg.append(np.log10(np.mean(np.square(abs(shellvec2)))))

                if ref:
                    if mft1 == 0.0:
                        break
                    ft1[fshells1] = shellvec1*(mft2/mft1)
                else:
                    ft1[fshells1] = shellvec1*(mft2+mft1)/(2*mft1)
                    ft2[fshells2] = shellvec2*(mft2+mft1)/(2*mft2)

                mft1 = np.mean(abs(ft1[fshells1]))
                ft1_avg_new.append(
                    np.log10(
                        np.mean(
                            np.square(abs(ft1[fshells1]))
                        )
                    )
                )
                lfreq.append((highlevel + step / 2))
                break
            x = highlevel
            highlevel = x + step
        
        if self.map.resolution is not None:
            if lpfilta and not lpfiltb:
                self.tanh_lowpass(ft1, cutoff1, fall=0.2)
                self.tanh_lowpass(ft2, cutoff2,fall=0.2)
        
        density_map_filtered = self.ifft_array(ft1)
        mol_map_filtered = self.ifft_array(ft2)
        return density_map_filtered, mol_map_filtered
               
           
    def fft_array(self, density_map):
        density_copy = density_map.copy()
        ft1 = np.fft.fftshift(np.fft.fftn(density_copy))
        return ft1
    
    def ifft_array(self, density_map):
        density_map = np.real(np.fft.ifftn(np.fft.ifftshift(density_map)))
        return density_map
    
    def tanh_lowpass(self, array, cutoff, fall=0.2):
        
        cutoff = min(float(cutoff), 0.5)
        drop = np.pi / (2 * cutoff * fall)
        dist = self.make_fourier_shell(array)
        array *= 0.5 * (np.tanh(drop * (dist + cutoff)) - np.tanh(drop * (dist - cutoff)))

    def make_fourier_shell(self, array, fill=1.0):
        # Computing floor and ceil values once
        z_floor, z_ceil = np.floor(array.shape[0] / 2.0), np.ceil(array.shape[0] / 2.0)
        y_floor, y_ceil = np.floor(array.shape[1] / 2.0), np.ceil(array.shape[1] / 2.0)
        x_floor, x_ceil = np.floor(array.shape[2] / 2.0), np.ceil(array.shape[2] / 2.0)
    
        # Creating radial grids for each dimension
        rad_z = np.arange(-z_floor, z_ceil) / array.shape[0]
        rad_y = np.arange(-y_floor, y_ceil) / array.shape[1]
        rad_x = np.arange(-x_floor, x_ceil) / array.shape[2]
    
        # Squaring each component and adding them using broadcasting
        dist = np.sqrt(rad_z[:, None, None]**2 + rad_y[None, :, None]**2 + rad_x[None, None, :]**2)
        
        return dist
    
    def write_map(self):
        for num, dens_map in enumerate(self.system.difference_maps):
            
            outfile = os.path.join(self.system.preprocessing_output, f'Difference_map_{num}.mrc')
            dens_map.write_mrc(outfile)
            
    def run(self):
        
        print(Messages.difference_map())
        self.get_data()
        self.difference_map_input()
        self.generate_map()
        self.write_map()
    
class DifferenceMap:
    #better difference map !!
    def __init__(self, system):
        self.system = system 
        self.percent_above_otsu_threshold = 0.0 # values between 0. and 1.
        self.non_zero_voxel_count_threshold = 50 #basically size !!! 
    
    
    def get_data(self):
        self.map = self.system.segment.maps[0].copy()
        
        vectorised_voxels, map_indices = self.vectorise_points(self.map.density_map.shape, 
                                                               np.array(self.map.origin), 
                                                               np.array(self.map.apix))
        
        
       
        positions = []
        atom_radii = []
        
        for residue_id in self.system.segment.proteins[0]:
            residue = self.system.proteins[0].get_residue(residue_id)
            for atom in residue:
                positions.append(atom.coords)
                atom_radii.append(RDTools.get_van_der_waals_radius(atom.element))
                
           
        self.positions = np.array(positions)
        self.atom_radii = np.array(atom_radii)
        points_outside_mask = self.find_points_outside_atoms(self.positions, self.atom_radii, 
                                                         vectorised_voxels)
        map_indices = map_indices[points_outside_mask]
        
        copy_map = np.zeros(self.map.density_map.shape) 
        for index in map_indices:
            copy_map[index[0],index[1],index[2]] = self.map.density_map[index[0],index[1],index[2]]
        
        self.map.density_map = copy_map
        
        
    
    def get_disconected_densities(self, image, threshold, struct = None):
        if struct is None:
            struct = generate_binary_structure(3, 1)
            
        bool_image = (image > threshold)
        labels, num_features = label(bool_image, structure = struct)
        return labels
    
    def get_sigma(self, resolution, sigma_coeff=0.356):
        return resolution * sigma_coeff
    
    def smooth_image(self,image, sigma):
        smoothed_image = gaussian_filter(image, sigma=sigma)
        return smoothed_image
    
    def  get_percent_over_threshold(self, masked_map, threshold):
        
        count_above_zero = np.sum(masked_map > 0)
        
        count_above_otsu = np.sum(masked_map > threshold)
        percent = round(count_above_otsu / count_above_zero , 2)
        
        return count_above_zero, count_above_otsu, percent
    
    
    def vectorise_points(self, shape, origin, apix):
        z_indices, y_indices, x_indices = np.indices(shape)
        
        real_world_coords_x = origin[0] + x_indices * apix[0]
        real_world_coords_y = origin[1] + y_indices * apix[1]
        real_world_coords_z = origin[2] + z_indices * apix[2]
        combined_coords = np.stack([real_world_coords_x, real_world_coords_y, real_world_coords_z], axis=-1)
        vectorized_voxels = combined_coords.reshape(-1, 3)
        indices = np.stack([z_indices, y_indices, x_indices], axis=-1).reshape(-1, 3)
        
        return vectorized_voxels, indices
    
    def find_points_outside_atoms(self, atomic_coords_array, atomic_radii_array,
                                  vectorised_voxels):
        
        dist_matrix = distance_matrix(vectorised_voxels, atomic_coords_array)
        is_outside_matrix = dist_matrix > atomic_radii_array
        points_outside_mask = np.all(is_outside_matrix, axis=1)
        return points_outside_mask
    
    def significant_features(self):
        feature_maps = []
        full_map = self.system.maps[0].density_map * (self.system.maps[0].density_map > 0) 
        full_map_threshold = filters.threshold_otsu(full_map.ravel())
       
        #self.map.density_map = self.map.density_map * (self.map.density_map > 0)
        map_copy = self.map.density_map.copy()
        structuring_element = generate_binary_structure(3, 1)
       
        image_3d = map_copy * (map_copy  > 0)
        eroded_image = grey_opening(image_3d, footprint=structuring_element)
        flattened = eroded_image.ravel()
        global_thresh = filters.threshold_otsu(flattened)
        segmented_density = self.get_disconected_densities(eroded_image, global_thresh)
        sigma = self.get_sigma(self.map.resolution)
        for num in np.unique(segmented_density)[1:]:
            mask_1 = (segmented_density == num)
            map_1 = eroded_image * mask_1
            map_1_closed = grey_closing(map_1, footprint=structuring_element)
            smooth_mask = self.smooth_image(map_1_closed , sigma)
            flattened = smooth_mask.ravel()
            smooth_thresh = filters.threshold_otsu(flattened)
            threshold_image = smooth_mask * (smooth_mask > smooth_thresh)
            masked_region = map_copy * (threshold_image > 0)
            n_voxels_above_zero, n_voxels_above_otsu, percent_above_otsu = self.get_percent_over_threshold( masked_region, full_map_threshold)
            feature_maps.append([ ( n_voxels_above_zero, n_voxels_above_otsu, percent_above_otsu),
                               masked_region])
       
        self.feature_maps = feature_maps
        
    #def temp_save(self):
        
    #    for i in
    
    def filter_maps(self):
        
        current_filtered_maps = []
        new_map  = np.zeros(self.map.density_map.shape)
        
        for num, sig_feat in enumerate(self.feature_maps):
            #filter for voxels % above threshold
            if sig_feat[0][2] > self.percent_above_otsu_threshold:
                #filter for voxels % above threshold 
                #TODO! make a maximum_threshold and slider.
                if sig_feat[0][0] > self.non_zero_voxel_count_threshold:
                    new_map += sig_feat[1]
        
        
        new_map_object = self.map.copy() 
        new_map_object.density_map = new_map
        new_map_object.normalise()
        self.system.difference_maps.append(new_map_object)
    
    def write_maps(self):
        for num, dens_map in enumerate(self.system.difference_maps):
            
            outfile = os.path.join(self.system.preprocessing_output, f'Difference_map_{num}.mrc')
            dens_map.write_mrc(outfile)
            
    
    def run(self):
        print(Messages.difference_map())
        self.get_data()
        self.significant_features()
        self.filter_maps()
        self.write_maps()
        
        
 
        