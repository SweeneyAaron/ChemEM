# This file is part of the ChemEM software.
#
# Copyright (c) 2023 - Topf Group & Leibniz Institute for Virology (LIV),
# Hamburg, Germany.
#
# This module was developed by:
#   Aaron Sweeney    <aaron.sweeney AT cssb-hamburg.de>


from ChemEM.messages import Messages
from ChemEM.scoring_functions import ScoringFunctions
from ChemEM.tools.map_tools import MapTools 
from ChemEM.tools.system_tools import SystemTools
from pathlib import Path
import os

class Rescore:
    
    def __init__(self, system):
        self.system = system 
    
    def mk_output_dir(self):
        self.system.rescore_output = os.path.join(self.system.output, 'rescore')
        try:
            os.mkdir(self.system.rescore_output)
        except FileExistsError:
            Messages.overwrite(self.system.rescore_output, self.system.overwrite)
    
    
    def get_aco_object(self):
        if not self.system.dock_only:
            if self.system.mi_weight is None:
                self.system.mi_weight = MapTools.get_mi_weight(self.system.maps[0].resolution)
 
        self.docking_run = SystemTools.dock_prep('D0',
                              [self.system.ligands[0]], 
                              self.system.proteins[0],
                              self.system.docking_radius,
                              self.system.centroid,
                              self.system.segment.maps[0],
                              self.system.system_data,
                              platform = self.system.platform,
                              mi_weight =self.system.mi_weight,
                              global_k = self.system.global_k,
                              cutoff = self.system.cutoff,
                              flexible_side_chains = self.system.flexible_side_chains,
                              solvent = self.system.solvent,
                              n_cpu = self.system.n_cpu)
    
    def score(self):
        all_chemem_score = []
        for lig in self.system.ligands:
            chemem_score = ScoringFunctions.ChemEM_score([lig],
                                                         self.docking_run,
                                                         difference_map = [self.system.difference_maps[0]]  )  
            
            full_path = Path(lig.lig_id)
            lig_nam = full_path.name 
            line = f'{lig_nam} : {chemem_score}\n'
            all_chemem_score.append((line, chemem_score))
            
        all_chemem_score = sorted(all_chemem_score, key = lambda x: x[1], reverse = True)
        self.all_chemem_score = [i[0] for i in all_chemem_score]
        
    def write_file(self):
        out_file = os.path.join(self.system.rescore_output, 'rescore.txt')
        with open(out_file, 'w') as f:
            f.writelines(self.all_chemem_score)

    def run(self):
        print(Messages.rescore())
        self.mk_output_dir()
        self.get_aco_object()
        self.score()
        self.write_file()
        