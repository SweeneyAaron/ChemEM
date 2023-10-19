#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 31 11:13:51 2023

@author: aaron.sweeney
"""

#script to run ChemEM

import sys
from ChemEM.system import System, SystemFactory, Config

def main():
    if len(sys.argv) != 2:
        print("Usage: chemem.py <config_file>")
        sys.exit(1)
    
    conf_file = sys.argv[1]

    # read conf file
    config = Config()
    config.load_config(conf_file)
    # set up system
    sys = SystemFactory()
    sys.create_system(config)
    # run commands
    sys.make_run()

if __name__ == "__main__":
    main()

