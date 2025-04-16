# This file is part of the ChemEM software.
#
# Copyright (c) 2023 - Topf Group & Leibniz Institute for Virology (LIV),
# Hamburg, Germany.
#
# This module was developed by:
#   Aaron Sweeney    <aaron.sweeney AT cssb-hamburg.de>

from ChemEM.messages import Messages
import sys
import os 
from ChemEM.tools.rdtools import RDTools
from ChemEM.tools.system_tools import SystemTools

class PreProcess:
    def __init__(self, system):
        self.system = system 
    
    
    def centroid_check(self):
        if not self.system.centroid:
            Messages.no_class_attribute_found(self.system.__class__.__name__, 'centroid')
            sys.exit()
    
    def mk_output_dir(self):
        self.system.preprocessing_output = os.path.join(self.system.output, 'preprocessing')
        try:
            os.mkdir(self.system.preprocessing_output)
        except FileExistsError:
            Messages.overwrite(self.system.preprocessing_output, self.system.overwrite)
        
    def set_dock_only(self):
        
        if hasattr(self.system, 'dock_only'):
            pass
        elif not self.system.maps:
            self.system.dock_only  = True
        else:
            self.system.dock_only = False
    
    def set_up_rd_fragments(self):
        try:
            for mol in self.system.ligands:
                name = mol.residue_name
                RDTools.BRICS_decomposition(mol)
                fragment_file_name = os.path.join(self.system.preprocessing_output, f'{name}_fragment_file.txt')
                RDTools.fragment_file_writter(mol, fragment_file_name )
                atom_indices_file = os.path.join(self.system.preprocessing_output, f'{name}_mol_indeces.svg')
                RDTools.draw_molecule_with_atom_indices(mol, atom_indices_file)
                frag_image_file = os.path.join(self.system.preprocessing_output, f'{name}_fragments.svg')
                RDTools.save_smiles_with_index(mol, frag_image_file)
        
        except Exception as e:
            print(Messages.chemem_warning(self.__class__, 'set_up_rd_fragments', e))
            
            
    def segment_system(self):
        #add exclude !!!
         
        segment_id, segment  = SystemTools.segment_system(self.system)
        self.system.segment = segment
        if not self.system.dock_only:
            for num, mp in enumerate(self.system.segment.maps):
                file = os.path.join(self.system.preprocessing_output, f"Segment_{num}_full_map.mrc")
                mp.write_mrc(file)
                                    
                
        
    def run(self):
        print(Messages.preprocess())
        self.centroid_check()
        self.mk_output_dir()
        self.set_dock_only()
        self.set_up_rd_fragments()
        self.segment_system()
        
        
        