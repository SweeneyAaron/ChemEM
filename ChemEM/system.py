# This file is part of the ChemEM software.
#
# Copyright (c) 2023 - Topf Group & Leibniz Institute for Virology (LIV),
# Hamburg, Germany.
#
# This module was developed by:
#   Aaron Sweeney    <aaron.sweeney AT cssb-hamburg.de>

import numpy as np
import multiprocessing
import math
from ChemEM.messages import Messages
import sys 

class System():
    def __init__(self):
        
        self.system_id = 0
        self.system_data = None
        self.overwrite = True
        
        self.segment = None
        self.docking_jobs = []
        
        #map properties
        self.maps = []
        self.difference_maps = []
        #self.map_contours = [] #move map contours to map objects!!
        
        #protein props
        self.proteins = []
        self.ligands = []
        self.system_ligands = []
        self.num_system_ligands = 0
        
        #centroid data
        self.centroid = []
        self.segment_centroid = None
        self.output = './'
        
        #preprocessing data
        self.auto_split_radius = 5.0
        self.label_threshold = None 
        self.label_threshold_sigma = None
        self.pH = [6.4, 8.4]
        self.pKa_prec = 1.0
        #exclude data
        self.exclude = []
        #self.include_ligands_in_diff_map = 0 will include system ligands now!!
        self.segment_dimensions = (20, 20, 20)
        
        #fitting options
        self.mi_weight = None 
        self.global_k = None 
        self.docking_radius = None
        
        self.platform = 'OpenCL'
        self.cutoff = 15.0       
        self.flexible_side_chains = 0
        self.solvent = False
        self.n_cpu = min(math.ceil(multiprocessing.cpu_count() / 2), 20)
        self.n_ants = None
        self.theta  = None
        self.rho = None
        self.sigma = None
        self.max_iterations = None
        self.generate_diverse_solutions = None
        self.sasa_cutoff = 5.0
        
        #post processing options
        self.post_process_num_solutions = 15
        self.post_process_solution = None
        self.refine_side_chains = 1
        self.cycles = 4
        self.start_temp = 0
        self.norm_temp = 300
        self.top_temp = 315
        self.temperature_step = 1
        self.pressure = 1
        self.barostatInterval = 10
        self.initial_heating_interval = 10
        self.heating_interval = 100
        self.steps = 1000
        self.hold_fragment = None
        self.protocols = []
        self.post_process_radius = 15.0
        
   
    def __repr__(self):
        return Messages.system_repr(self)
    
    def add_protocol(self, protocol):
        
        self.protocols.append(protocol(self))
        
    def run(self):
        for protocol in self.protocols:
            try:
                protocol.run()
            except Exception as e:
                print(Messages.fatal_exception(protocol.__class__, e))
                sys.exit()
                
    def run_protocol(self, protocol):
         protocol(self).run()
         
    


        
       
        
       