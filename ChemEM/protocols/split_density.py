# This file is part of the ChemEM software.
#
# Copyright (c) 2023 - Topf Group & Leibniz Institute for Virology (LIV),
# Hamburg, Germany.
#
# This module was developed by:
#   Aaron Sweeney    <aaron.sweeney AT cssb-hamburg.de>


from ChemEM.tools.map_tools import MapTools
from ChemEM.messages import Messages
import os

class SplitDensity:
    
    
    def __init__(self, system):
        self.system = system 
    
    def split_density(self):
        labeled_maps = []
        for diff_map in self.system.difference_maps:
            masked_map, labeled_map, num_features =  MapTools.split_density(diff_map,
                                                                            label_threshold_sigma = self.system.label_threshold_sigma,
                                                                            label_threshold = self.system.label_threshold,
                                                                            struct = None)
            
            
            labeled_maps.append(labeled_map)
            
        self.system.labeled_maps = labeled_maps
    
    def write_maps(self):
        for num, densmap in enumerate(self.system.labeled_maps):
            fn = f'System_{self.system.system_id}_labeled_map_{num}.mrc'
            out_file = os.path.join(self.system.preprocessing_output, fn)
            densmap.write_mrc(out_file)
        
    
    def run(self):
        print(Messages.split_density())
        self.split_density()
        self.write_maps()
        
