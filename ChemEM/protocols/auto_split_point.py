# This file is part of the ChemEM software.
#
# Copyright (c) 2023 - Topf Group & Leibniz Institute for Virology (LIV),
# Hamburg, Germany.
#
# This module was developed by:
#   Aaron Sweeney    <aaron.sweeney AT cssb-hamburg.de>

import numpy as np 
from ChemEM.tools.math import MathTools
from ChemEM.messages import Messages
import os

class AutoSplitPoint:
    
    def __init__(self, system):
        self.system = system 
    
    
    def get_closest_labels(self):
        
        
        
        masked_difference_maps = []
        labeled_map = self.system.labeled_maps[0] #At the moment only 1 map
        features = np.unique(labeled_map.density_map)[1:]

        for point in self.system.centroid:
            min_dist = 9999999999.00
            min_label = -1
            
            for label in features:
                dist = self.centroid_centroid_dist(label, point, labeled_map)
                if dist < min_dist:
                    min_dist = dist
                    min_label = label
        
            
            
            
            masked_map = self.apply_mask(labeled_map, [min_label])
            masked_difference_maps.append(masked_map)
        
        self.system.difference_maps = masked_difference_maps
        
        
    
    def centroid_centroid_dist(self, label, point, labeled_map):
        points = []
        for z in range(labeled_map.density_map.shape[0]):
            for y in range(labeled_map.density_map.shape[1]):
                for x in range(labeled_map.density_map.shape[2]):
                    if labeled_map.density_map[z,y,x] == label:
                        points.append(list(labeled_map.voxel_to_point(x,y,z)))

        map_point = MathTools.calc_centroid(points)
        dist = MathTools.euclidean_distance(point, map_point)
        return dist
    
    def apply_mask(self, labeled_map, label):
        new_lab_map = labeled_map.copy() 
        new_lab_map.density_map = new_lab_map.density_map * np.isin(new_lab_map.density_map, label)
        new_lab_map.density_map = new_lab_map.density_map > 0
        
        new_diff_map = self.system.difference_maps[0].copy()
        new_diff_map.density_map = new_diff_map.density_map * new_lab_map.density_map 
       
        return new_diff_map
        
   
    def write_maps(self):
        for num, densmap in enumerate(self.system.difference_maps):
            fn = f'System_{self.system.system_id}_autosplit_point_map_{num}.mrc'
            out_file = os.path.join(self.system.preprocessing_output, fn)
            densmap.write_mrc(out_file)
        
        
        
    def run(self):
        print(Messages.auto_split_point())
        self.get_closest_labels()
        self.write_maps()
        
    