# This file is part of the ChemEM software.
#
# Copyright (c) 2023 - Topf Group & Leibniz Institute for Virology (LIV),
# Hamburg, Germany.
#
# This module was developed by:
#   Aaron Sweeney    <aaron.sweeney AT cssb-hamburg.de>


from ChemEM.messages import Messages
from ChemEM.tools.rdtools import RDTools
from ChemEM.tools.map_tools import MapTools 
from ChemEM.tools.dock_tools import DockTools, Solution
from ChemEM.tools.system_tools import SystemTools
from ChemEM.tools.openff_tools import  OpenFFTools
import os 

class Fitting:
    def __init__(self, system):
        self.system = system 
    
    
    def mk_output_dir(self):
        self.system.fitting_output = os.path.join(self.system.output, 'fitting')
        try:
            os.mkdir(self.system.fitting_output)
        except FileExistsError:
            Messages.overwrite(self.system.fitting_output, self.system.overwrite)
    
    def verify_ligands_and_centroids(self):
    
        self.system.centroid = self._extend_or_retract(self.system.ligands,  self.system.centroid)
        if not self.system.dock_only:
            self.system.difference_maps = self._extend_or_retract(self.system.ligands,  self.system.difference_maps)
            
    def _extend_or_retract(self, list1, list2):
         
        if len(list2) < len(list1):
            length_difference = len(list1) - len(list2)
            list2.extend([list2[0]] * length_difference)
            
        elif len(list2) > len(list1):
            list2[:len(list1)]
        
        return list2 
    
    def dock_ligand(self):
        
        for centroid, lig in zip(self.system.centroid, self.system.ligands):

            RDTools.translate_mol_to_centroid(lig, centroid)
        
        if not self.system.dock_only:
            if self.system.mi_weight is None:
                self.system.mi_weight = MapTools.get_mi_weight(self.system.maps[0].resolution)
        
       
        if not self.system.segment.maps:
            densmap = None
        else:
            densmap = self.system.segment.maps[0]
       
        #change back to docking radius!!
        if self.system.docking_radius is None:
            self.system.docking_radius = self.system.segment_dimensions
            
        elif type(self.system.docking_radius) in [float, int]:
            
            self.system.docking_radius = [self.system.docking_radius] * 3
        
        docking_run = SystemTools.dock_prep('D0',
                              self.system.ligands, 
                              self.system.proteins[0],
                              self.system.docking_radius, #cahnged fromself.system.docking_radius
                              self.system.centroid,
                              densmap,
                              self.system.system_data,
                              platform = self.system.platform,
                              mi_weight =self.system.mi_weight,
                              global_k = self.system.global_k,
                              cutoff = self.system.cutoff,
                              flexible_side_chains = self.system.flexible_side_chains,
                              solvent = self.system.solvent,
                              n_cpu = self.system.n_cpu,
                              sasa_cutoff = self.system.sasa_cutoff)
        
       
        
        if self.system.max_iterations is not None:
            docking_run.max_iterations = self.system.max_iterations
        
        if self.system.n_ants is not None:
            docking_run.n_ants = self.system.n_ants
        
        if self.system.theta is not None:
            docking_run.theta = self.system.theta 
        
        if self.system.rho is not None:
            docking_run.rho = self.system.rho 
        
        if self.system.sigma is not None:
            docking_run.sigma = self.system.sigma
        
        docking_run.generate_diverse_solutions = self.system.generate_diverse_solutions
        DockTools.ACO(docking_run, self.system.difference_maps)
        
        out_file = []
        protein_path = os.path.join(self.system.fitting_output, 'PDB')
        try:
            os.mkdir(protein_path)
        except:
            pass
        
        copy_of_coords = docking_run.complex_system.positions.copy()
        
        if self.system.generate_diverse_solutions:
            diverse, mapping, removed = self.get_diverse_solutions(docking_run.best_solutions, cutoff = self.system.generate_diverse_solutions )
            docking_run.best_solutions = diverse + removed
        
        for result in docking_run.best_solutions:
            result.write_refined_solution(self.system.fitting_output)
            #result.write_chemem_solution(self.system.fitting_output)
            out_file.append(f'{result.iteration} :  {result.refined_score} \n')
            protein_out = os.path.join(protein_path, f'Ligand_{result.iteration}.pdb')
            docking_run.complex_system.simulation.context.setPositions(result.protein_positions)
            OpenFFTools.write_complex_system(docking_run.complex_system.simulation.context.getState(getPositions=True).getPositions(),
                                             docking_run.complex_system.simulation.topology,
                                             protein_out, 
                                             missing_res = docking_run.protein._missing_residues)
            
        docking_run.complex_system.simulation.context.setPositions(copy_of_coords)
        of = os.path.join(self.system.fitting_output , 'results.txt')
        with open(of, 'w') as f:
            f.writelines(out_file)
        
        
        self.system.docking_jobs.append(docking_run)
        
        
    def get_diverse_solutions(self,solutions, cutoff = 2.0):
        diverse = []
        removed = []
        mapping = {}
        for sol in solutions:
            for sol2 in diverse:
                
                rmsd = RDTools.diverse_rmsd(sol,sol2)

                if rmsd <= cutoff:
                    mapping[sol2.iteration].append(sol.iteration)
                    removed.append(sol)
                    break
            else:
                diverse.append(sol)
                mapping[sol.iteration] = []
        
        return diverse, mapping, removed
    

    def run(self):
        if not self.system.max_iterations == 0:
            print(Messages.fitting())
        self.mk_output_dir()
        self.verify_ligands_and_centroids()
        self.dock_ligand()
    