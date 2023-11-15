#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 31 11:13:51 2023

@author: aaron.sweeney
"""

#script to run ChemEM

import sys
from ChemEM.system import System, SystemFactory, Config
import os

def test():
    
    if len(sys.argv) == 2:
        
        platform = sys.argv[1]
        if platform not in ['CPU', 'OpenCL', 'CUDA']:
            print("Usage: chemem.test <platform>")
            print('Currently included platforms: CPU, OpenCL, CUDA')
            sys.exit(1)
   
    elif len(sys.argv) == 1:
        platform = 'OpenCL'
    
    else:
        print("Usage: chemem.test <platform>")
        print('Currently included platforms: CPU, OpenCL, CUDA')
        sys.exit(1)
        
        
    pdb_id = '7jjo'
   
        
    data_path = f'{os.path.dirname(os.path.abspath(__file__))}/test_data/'
    conf_path = os.path.join(data_path, f'{pdb_id.lower()}_conf.txt')
    new_output_dir = f'./ChemEM_test_{pdb_id}' 
    os.mkdir(new_output_dir)
    new_conf_path = get__test_conf_file(conf_path, new_output_dir, data_path, pdb_id, platform)
    
    try:
        # read conf file
        config = Config()
        config.load_config(new_conf_path)
    
    except FileNotFoundError:
        print(f'FileNotFoundError no file named: { new_conf_path}')
        sys.exit(1)
    # set up system
    new_system = SystemFactory()
    new_system.create_system(config)
    # run commands
    new_system.make_run()
    
    sys.exit(1)

def get__test_conf_file(file, output, data_path, pdb_id, platform):
    with open(file, 'r') as f:
        f = f.read().splitlines()
    new_conf = []
    for i in f:
        if i.startswith('protein'):
            prot = pdb_id + '_protein.pdb'
            new_conf.append(f'protein = {os.path.join(data_path, prot )}\n')
        elif i.startswith('densmap'):
            densmap = pdb_id + '.mrc'
            new_conf.append(f'densmap = {os.path.join(data_path, densmap )}\n')
        elif i.startswith('output'):
            new_conf.append(f'output = { output }\n')

        else:
            line = i + '\n'
            new_conf.append(line)
    
    new_conf.append(f'platform = {platform}\n')
    
    fn = pdb_id + '_conf.txt'
    new_conf_path = os.path.join(output,fn)
    with open(new_conf_path, 'w') as f:
        f.writelines(new_conf)
    return new_conf_path
        
        
    
    

def main():
    if len(sys.argv) != 2:
        print("Usage: chemem <config_file>")
        sys.exit(1)
    
    conf_file = sys.argv[1]
    if conf_file in ['-h', '--help']:
        print('Usage chemem <config file>')
        sys.exit(1)
    
    if conf_file == '7jjo':
        print('Usage chemem.test <pdb_id>')
        sys.exit(1)
    
    if conf_file == '6tti':
        print('Usage chemem.test <pdb_id>')
        sys.exit(1)
        pass
    try:
        # read conf file
        config = Config()
        config.load_config(conf_file)
    
    except FileNotFoundError:
        print(f'FileNotFoundError no file named: { conf_file }')
        sys.exit(1)
    # set up system
    new_system = SystemFactory()
    new_system.create_system(config)
    # run commands
    new_system.make_run()

if __name__ == "__main__":
    main()

