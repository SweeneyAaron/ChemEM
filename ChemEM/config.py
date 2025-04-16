# This file is part of the ChemEM software.
#
# Copyright (c) 2023 - Topf Group & Leibniz Institute for Virology (LIV),
# Hamburg, Germany.
#
# This module was developed by:
#   Aaron Sweeney    <aaron.sweeney AT cssb-hamburg.de>

#get args from config!! 
from ChemEM.system import System
from ChemEM.data_classes.data import BasicData, SystemAttrs
from ChemEM.parsers import Parsers
from ChemEM.protocol import protocol_selector
from dataclasses import dataclass, field 
from typing import List
import ast
import multiprocessing 
import os

@dataclass
class Config:
    pre_process : bool = field(default=True)
    fitting : bool = field(default=False)
    post_process : bool = field(default=False)
    simulation : bool = field(default=False)
    rescore : bool = field(default=False)
    protein : List = field(default_factory=list)
    ligand : List = field(default_factory=list)
    system_ligand_file : List =  field(default_factory=list)
    densmap : List = field(default_factory=list)
    resolution : List = field(default_factory=list)
    centroid : List = field(default_factory=list)
    output : str = field(default="./")
    load = None
    segment_centroid = None
    map_contour : List = field(default_factory=list)
    local_contour : List = field(default_factory=list)
    full_map_dimensions : List[int] = field(default=(30, 30, 30))
    segment_dimesions : List[int] = field(default=(30, 30, 30))
    exclude : List= field(default_factory=list)
    ligands_from_dir : str = field(default=None)
    system_ligands_from_dir : str = field(default=None)
    difference_map : List = field(default_factory=list)
    local_resolution : List = field(default_factory=list)
    full_map_id : str = field(default=None)
    platform : str = 'OpenCL' #get best platform
    cutoff : float = 20.0
    flexible_side_chains : bool = field(default=False)
    solvent : bool = field(default=False)
    n_cpu : int =  None
    post_process_solution : List = field(default_factory=list)
    hold_fragment : List = field(default_factory=list)
    protonation : bool = field(default=True)
    chirality : bool = field(default=True)
    rings : bool = field(default=True)
    pH : List[float] = field(default=(6.4, 8.4))
    pKa_prec : float = field(default=1.0) #change this to pKa_prec
    forcefield: List[str] = field(default_factory=list)
    sasa_cutoff : float = field(default=None)
    difference_map_protocol: str =  field(default="version_2")
    __list_ids__ = ['protein', 'ligand', 'densmap', 'resolution', 'centroid', 'local_contour', 'local_resolution',
                    'map_contour', 'exclude', 'difference_map', 'post_process_solution' ,'hold_fragment', 'system_ligand_file',
                    'forcefield']
    
    
    def _process_line(self, line):
        line = line.split('=', maxsplit=1)
        if len(line) < 2:
            return

        attr_id, value = line[0].strip(), line[1].strip()
        try:
            value = ast.literal_eval(value)
        except:
            pass
        
        if attr_id in self.__list_ids__:
            getattr(self, attr_id).append(value)
            
        else:
            setattr(self, attr_id, value)
        
        
        print(f"\t{attr_id}, {value}, {type(value)}")

    def load_config(self, config_file):
        with open(config_file, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    self._process_line(line.strip())
        
        if self.n_cpu is None:
            cpus = round(multiprocessing.cpu_count() / 2)
            if cpus > 20:
                cpus = 20
            
            self.n_cpus = cpus
        
        return self.create_system()
    
    def create_system(self):
        system = System() 
        system_data = BasicData()
        system.system_data = system_data 
        system.__system_data_ref__ = id(system_data)
        
        #----------------add protein 
        for num, protein_file in enumerate(self.protein):
            
            system.proteins.append(self.add_protein(protein_file , num + 1, system.__system_data_ref__, self.forcefield))
        
        
        #----------------add ligands
        ligand_residue_number = max([i.num_residues for i in system.proteins]) + 1
        
        if self.ligands_from_dir is not None:
            self.ligand += self.get_ligands_from_dir(self.ligands_from_dir)
            
        for num, ligand in enumerate(self.ligand):
            system.ligands.append( self.add_ligand(ligand, 
                                                   num, 
                                                   system.__system_data_ref__,
                                                   residue_number = num + ligand_residue_number) )
        
        ligand_residue_number += num
        #----------------add  system ligands
        
        if self.system_ligands_from_dir is not None:
            self.system_ligand_file += self.get_ligands_from_dir(self.system_ligands_from_dir)
        
        for num, ligand in enumerate(self.system_ligand_file):
            system_ligand = self.add_ligand(ligand, 
                                            num, 
                                            system.__system_data_ref__,
                                            residue_number = num + ligand_residue_number) 
            
            #TODO!! expand to add to system rather than protein 
            system.proteins[0].add_system_ligand(system_ligand)
        
        #----------------add maps

        for (density_map, map_contour), resolution in zip(self.match_files_with_contours(self.densmap, self.map_contour), self.resolution):
            
            system.maps.append( self.add_map(density_map, resolution,  map_contour))
        
        #----------------add difference_maps 
        if len(self.difference_map) > 0:
            combined_diffmap_resolution = self.match_files_with_contours(self.difference_map, self.local_resolution, value = self.resolution[0])
            combined_diffmap_contour = self.match_files_with_contours(self.difference_map, self.local_contour)
            
            for (density_map, resolution), (_, map_countour) in zip(combined_diffmap_resolution, combined_diffmap_contour):
                system.difference_maps.append(self.add_map(density_map, resolution, map_contour))

        #----------------set system attrs
        sys_attrs = SystemAttrs()
        for attr in sys_attrs.system_attrs:
            value = getattr(self, attr, None)
            if value is not None:
                setattr(system, attr, value)
        
        #---------------Set Protocols
        
        system.add_protocol(protocol_selector('pre_process'))
        
        if system.maps:
            #add difference mapping protocols for now only differnece map!!
            if not system.difference_maps:
                
                if self.difference_map_protocol == 'original':
                    system.add_protocol(protocol_selector('difference_map_original'))
                    
                elif self.difference_map_protocol == 'version_2':
                    system.add_protocol(protocol_selector('difference_map'))

                else:
                    system.add_protocol(protocol_selector('difference_map'))
                #---split densirty
                if self.pre_process_split_density:
                    system.add_protocol(protocol_selector('pre_process_split_density'))
                
                #----autosplit point
                if self.auto_split_point:
                    if not self.pre_process_split_density:
                        system.add_protocol(protocol_selector('pre_process_split_density'))
                    system.add_protocol(protocol_selector('auto_split_point'))
                    
                #----autosplit point
                #enforce one or the other
                elif self.auto_split_zone:
                    if not self.pre_process_split_density:
                        system.add_protocol(protocol_selector('pre_process_split_density'))
                    system.add_protocol(protocol_selector('auto_split_zone'))
                
        #the set up for this will change for version 2 
        if self.fitting:
            
            system.add_protocol(protocol_selector('fitting'))

        if self.post_process:
            if not self.fitting:
                system.add_protocol(protocol_selector('fitting'))
                system.max_iterations = 0
            system.add_protocol(protocol_selector('post_process'))
    
        
        if self.rescore:
            
            system.add_protocol(protocol_selector('rescore'))

        return system
    
    
    def add_protein(self, protein_file, protein_id, system_data_ref , forcefield = []):
        
        
        protein_id = f'P{protein_id}'
        
        protein = Parsers.protein_parser(protein_file, protein_id, system_data_ref , forcefield = forcefield)
        return protein
        
    
    def add_ligand(self, ligand_input, 
                   ligand_id , 
                   system_data_ref,
                   chain = '',
                   residue_number = None,
                   ):
        
        
        ligand_id = f'L{ligand_id}'
        
        if residue_number is None:
            residue_number = 0
            for prot in self.proteins.values():
                residue_number += prot.num_residues
            
            residue_number += len(self.ligands)
            
        ligand = Parsers.ligand_parser(ligand_input,
                                       system_data_ref,
                                       ligand_id = ligand_id,
                                       chain = chain,
                                       residue_number = residue_number,
                                       protonation = self.protonation,
                                       chirality= self.chirality,
                                       rings = self.rings,
                                       pH = self.pH,
                                       n = self.pKa_prec
                                      )
        
        return ligand
    
    def get_ligands_from_dir(self, path):
        return [os.path.join(path , i) for i in os.listdir(path) if i.endswith('.sdf') or i.endswith('.mol2')]
    
    def match_files_with_contours(self, list1, list2 , value = None):
    
        # Determine the length difference
        length_difference = len(list1) - len(list2)
        #length_difference = len(self.densmap) - len(self.map_contour)
        
        # Adjust the size of contour_list based on the length difference
        if length_difference > 0:
            # If file_list is longer, extend contour_list with None
            list2.extend([value] * length_difference)
        else:
            # If contour_list is longer, reduce its size to match file_list
            list2 = list2[:len(list1)]
        
        # Combine the file paths and map contours into a list of tuples
        combined_list = list(zip(list1, list2))
        
        return combined_list
    
    def add_map(self, file,resolution,  contour = None):
        
        density_map = Parsers.map_parser(file, resolution, contour)
        return density_map

#build included 
