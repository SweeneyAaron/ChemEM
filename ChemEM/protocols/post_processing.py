# This file is part of the ChemEM software.
#
# Copyright (c) 2023 - Topf Group & Leibniz Institute for Virology (LIV),
# Hamburg, Germany.
#
# This module was developed by:
#   Aaron Sweeney    <aaron.sweeney AT cssb-hamburg.de>

from ChemEM.messages import Messages
from ChemEM.tools.map_tools import MapTools 
from ChemEM.tools.openff_tools import  OpenFFTools
from ChemEM.tools.dock_tools import Solution
from rdkit import Chem
import os

class PostProcessing:
    def __init__(self, system):
        self.system = system 
    
    def mk_output_dir(self):
        self.system.post_process_output = os.path.join(self.system.output, 'post_processing')
        try:
            os.mkdir(self.system.post_process_output)
        except FileExistsError:
            Messages.overwrite(self.system.post_process_output, self.system.overwrite)
    
    def get_ligands(self):
        if self.system.docking_jobs[0].max_iterations:
            self.post_docking = 1
        else:
            self.post_docking = 0
            
    def get_map_bounds(self):
        
        if self.system.segment.maps:
            self.map_bounds =  MapTools.get_map_bounds(self.system.segment.maps[0])
            self.densmap = self.system.segment.maps[0]
            
        else:
            
            self.map_bounds = None 
            self.densmap = None

    def get_global_k(self):
        
        self.system.global_k = OpenFFTools.get_global_k(self.densmap.resolution)
   
    def get_post_process_complex_system(self):

        if self.system.hold_fragment:
            if self.post_docking:
                print(Messages.no_hold_fragment())
                exclude = None
            else:
                exclude = OpenFFTools.get_excluded_ligand_atoms(self.system.docking_jobs[0], self.system.hold_fragment)
        else:
            exclude = None
        
        
        self.post_process_complex_system = OpenFFTools.openmm_complex_system(self.system.proteins[0],
                                                                self.system.ligands,
                                                                self.system.centroid,
                                                                dens_map = self.densmap,
                                                                map_bounds = self.map_bounds,
                                                                global_k = self.system.global_k,
                                                                cutoff= self.system.post_process_radius,
                                                                flexible_side_chains = self.system.refine_side_chains,
                                                                solvent = self.system.solvent,
                                                                platform = self.system.platform,
                                                                exclude = exclude)
        
    
    
    def get_solutions_from_system(self):
        
        mols = [Chem.Mol(i) for i in self.system.ligands]
        self.system.docking_jobs[0].best_solutions.append(
            Solution(iteration  = 0,
                    chemem_solution = mols,
                    chemem_score = -1,
                    refined_solution = mols,
                    refined_score = -1,
                    protein_positions =  self.system.docking_jobs[0].complex_system.positions) )
        
    def post_process_solutions(self):
        import pdb 
        pdb.set_trace()
        self.system.post_processed_solutions = []
        for num, solution in enumerate(self.system.docking_jobs[0].best_solutions):
            if num + 1 > self.system.post_process_num_solutions:
                break
            
            out_file = os.path.join(self.system.post_process_output, f'Ligand_{solution.iteration}')
            
            post_processed_solutions = OpenFFTools.simulated_anneling(
                self.post_process_complex_system,
                solution,
                missing_res = self.system.docking_jobs[0].protein._missing_residues,
                out_file = out_file,
                cycles = self.system.cycles,
                start_temp = self.system.start_temp,
                norm_temp = self.system.norm_temp,
                temperature_step = self.system.temperature_step,
                pressure = self.system.pressure,
                barostatInterval = self.system.barostatInterval,
                initial_heating_interval = self.system.initial_heating_interval,
                heating_interval = self.system.heating_interval,
                steps = self.system.steps,
                aco_object = self.system.docking_jobs[0],
                difference_map = self.system.difference_maps)
            self.system.post_processed_solutions.append(post_processed_solutions)
    
    def write_results(self):
        results_out = os.path.join(self.system.post_process_output, 'results.txt')
        results_file = []
        
        for sol_group in self.system.post_processed_solutions:
            results_file.append(f'{sol_group[0].ligand_id}\n')
            for sol in sol_group:
                results_file.append(f'\tCycle {sol.cycle} : {sol.chemem_score}\n')
        
        
        with open(results_out, 'w') as f:
            f.writelines(results_file)

    def run(self):
        
        print(Messages.post_processing())
        self.mk_output_dir()
        self.get_ligands()
        self.get_map_bounds()
        if self.densmap is not None:
            self.get_global_k()
        self.get_post_process_complex_system()
       
        if not self.post_docking:
            self.get_solutions_from_system()
        
        self.post_process_solutions()
        self.write_results()
        
        



        
        
        
        
        
        