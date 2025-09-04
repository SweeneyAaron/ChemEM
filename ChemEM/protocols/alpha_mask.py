# This file is part of the ChemEM software.
#
# Copyright (c) 2023 - Topf Group & Leibniz Institute for Virology (LIV),
# Hamburg, Germany.
#
# This module was developed by:
#   Aaron Sweeney    <aaron.sweeney AT cssb-hamburg.de>
"""
ChemEM alpha-sphere SASA masking (backport)

Summary
-------
Back-ported from ChemEM v2.0 to extend usefulness for ChemEM >=0.6.6 and <2.0.0.
Masks a cryo-EM map to identify candidate ligand densities within SASA-exposed
pocket volumes while excluding bulk solvent.

Method
------
1) Convert the protein to alpha-spheres:
   - small spheres → protein (core)
   - intermediate spheres → SASA / pocket boundary
   - large spheres → bulk solvent
2) Restrict the map to the SASA region and seed voxels whose mean density
   exceeds the local baseline.
3) Grow regions on a Sobel-filtered map until the rate of change in intensity
   falls below 0.2.

Outputs
-------
- Masked map(s) of unoccupied/candidate ligand density
- Segmented regions and their centroids
"""



import os 
import sys 
import json
import numpy as np
from scipy.spatial import Delaunay
from ChemEM.tools.rdtools import RDTools
from ChemEM.tools.math import MathTools
from ChemEM.tools.map_tools import MapTools
from ChemEM.protocols.threshold_map import ConfidenceMap
from ChemEM.data_classes.EMMap import EMMap
from ChemEM.messages import Messages

class AlphaMask:
    
    def __init__(self, system):
        self.system = system 
        self.PROBE_SPHERE_MIN = 3.0  # Minimum radius for binding site detection
        self.PROBE_SPHERE_MAX = 6.0
        self.pad = self.system.alpha_pad# TODO needs a new option in main
        self.am_sigma_thr = self.system.am_sigma_thr
    
    def mk_output(self):
        self.output = os.path.join(self.system.output, 'alpha_mask')
        try:
            os.mkdir(self.output)
        except FileExistsError:
            pass
        
    def set_grid_spacing(self):
        
        if self.system.maps[0] is not None:
            self.grid_spacing = self.system.maps[0].apix[0]
            
        else:
            print('ChemEM Fatal Error: No density map found')
            sys.exit()    
    
    def get_protein_structure(self):
        #TODO! system ligands!!
        self.protein_openff_structure =  self.system.proteins[0].openff_structure
    
    def get_position_input(self):
        
        
        self.atoms = [i for i in self.protein_openff_structure.atoms if i.element > 1]
        self.positions = np.array([[i.xx, i.xy, i.xz] for i in self.atoms ])
        self.atom_radii = np.array( [RDTools.get_van_der_waals_radius(i.element_name) for i in self.atoms] )
        
    def get_delaunay(self):
        self.delaunay = Delaunay(self.positions)
        self.tetrahedra = self.delaunay.simplices 
    
    
    def calculate_circumsphere(self):
        self.circumspheres = []
        for tetra in self.tetrahedra:
            vertices = self.positions[tetra]
            
            center, radius = MathTools.compute_circumsphere(vertices)
            
            self.circumspheres.append((center, radius))
    
    def filter_circumspheres(self):
        self.candidate_tetrahedra = []
        self.candidate_centers = []
        self.candidate_radii = []

        
        for i, (center, radius) in enumerate(self.circumspheres):
            if self.PROBE_SPHERE_MIN <= radius <= self.PROBE_SPHERE_MAX:
                self.candidate_tetrahedra.append(self.tetrahedra[i])
                self.candidate_centers.append(center)
                self.candidate_radii.append(radius)
        
    
    def make_confidence_map(self):
        thr_map = ConfidenceMap(self.system)
        thr_map.run()
        
    def get_mask(self):
        
        copy_map = self.system.confidence_map.copy()
        copy_map = MapTools.trim_map_to_atoms(copy_map, self.positions, self.pad)
        
        shape = copy_map.density_map.shape
        
        
        site_mask = MapTools.create_binding_site_mask(shape, 
                                             self.candidate_centers, 
                                             self.candidate_radii, 
                                             copy_map.apix,
                                             copy_map.origin)
        
        
        
        #copy_map.density_map = site_mask
        #copy_map.write_mrc(out_fn)
        
        masked_map = copy_map.copy()
        masked_map.density_map = copy_map.density_map * site_mask
    
        self.mask_map = copy_map 
        self.masked_density = masked_map
        
        boundry_mask = MapTools.find_binding_site_boundary(site_mask)
        distance_map = MapTools.compute_distance_map(boundry_mask, copy_map.apix)
        distance_map *= site_mask
        
        distance_copy = copy_map.copy() 
        distance_copy.density_map = distance_map 
        self.distance_map = distance_copy
        
    
    def get_features(self):
        
        
        
        ligand_dens, ligand_feat = MapTools.extract_ligand_density(self.masked_density.density_map, 
                                                          self.distance_map.density_map, 
                                                          self.masked_density.apix,
                                                          high_threshold_sigma=self.am_sigma_thr)
        
        
        
        combined_map = np.zeros(self.masked_density.density_map.shape)
        apix = np.array(self.masked_density.apix, dtype=float)    # [ax, ay, az]
        base_origin = np.array(self.masked_density.origin, dtype=float)  # [ox, oy, oz]
        
        
        centroids = {}
        num_sites = 0
       
        
        system_sites = {}
       
        
        
        
        for dens, feat in zip(ligand_dens, ligand_feat):
            combined_map += dens
            
            bbox = MapTools.extract_min_bbox(dens, thr=0.0, pad=1.0)
            if bbox is None:
                
                continue

            subvol, (z0, y0, x0), (z1, y1, x1) = bbox

            

            # Array axes: [Z, Y, X] -> indices (z0,y0,x0)
            # Real-space origin/apix order: [X, Y, Z]
            new_origin = base_origin + np.array([
                x0 * apix[0],  # X
                y0 * apix[1],  # Y
                z0 * apix[2],  # Z
            ], dtype=float)
            
            sub_map = EMMap(new_origin, 
                           apix,
                           subvol,
                           self.masked_density.resolution)
            
           
            
            sub_centroid = sub_map.center_of_mass()
            
            
            
            fid = feat.get('feature_id', None)
            
            
            out_name = f"feature_{fid:04d}_bbox.mrc" if fid is not None else f"feature_{x0}_{y0}_{z0}_bbox.mrc"
            
            centroids[out_name] = sub_centroid
            out_fn = os.path.join(self.output, out_name)
            sub_map.write_mrc(out_fn)
            
            
            
        
        copy_map = self.masked_density.copy()
        copy_map.density_map = combined_map 
        out_fn = os.path.join(self.output, 'combined_map.mrc')
        copy_map.write_mrc(out_fn)
        
        
        with open(os.path.join(self.output, 'centroids.json'), 'w') as f:
            json.dump(centroids, f, ensure_ascii=False, indent=2)
       

    def run(self):
        print(Messages.alpha_mask())
        self.mk_output() 
        self.get_protein_structure()
        self.set_grid_spacing()
        self.get_position_input()
        self.get_delaunay()
        self.calculate_circumsphere()
        self.filter_circumspheres()
        self.make_confidence_map()
        self.get_mask()
        self.get_features()
        
       
        
        