# This file is part of the ChemEM software.
#
# Copyright (c) 2023 - Topf Group & Leibniz Institute for Virology (LIV),
# Hamburg, Germany.
#
# This module was developed by:
#   Aaron Sweeney    <aaron.sweeney AT cssb-hamburg.de>


import sys
import os
import ChemEM
from ChemEM.config import Config
from ChemEM.messages import Messages
from ChemEM.tools.rdtools import RDTools


def test():
    
    print(Messages.test())
    
    sys.exit(0)


def protonate():
    
    if len(sys.argv) != 2:
        Messages.protonation_usage()
        sys.exit(1)
    
    conf_file = sys.argv[1]
    if os.path.exists(conf_file):
        #is conf file 
        try:
            print(Messages.load_config())
            # read conf file
            config = Config()
            new_system = config.load_config(conf_file)
            
            ligands = config.ligand
            pH = new_system.pH 
            n = new_system.pKa_prec
            output = new_system.output
            protonation_states = {}
            
        
        except Exception as e:
            Messages.fatal_exception('chemem.protonation', e)
        
        #if the file was read!
        for lig in ligands:
            try:
                found_protonation_states = RDTools.protonate(lig, pH=pH, n=n)
                protonation_states[lig] = found_protonation_states 
                
            except Exception as e:
                Messages.chemem_warning('chemem.protonation', lig, e)
                sys.exit()
        
        
        if protonation_states:
            out_file = ''
            for key, value in protonation_states.items():
                out_file += f'Input SMILES: {key}\nFound states at pH range {pH[0]} - {pH[1]}:\n'
                for smi in value:
                    out_file += f'\t{smi}\n'
            
            out_fn = os.path.join(output, 'protonation.txt')
            with open( out_fn, 'w') as f:
                f.write(out_file)
            
            print(f'Protonation states written to file: {out_fn}')
            print("To indicate protonation states, indicate your prefered state in the configuration file along with 'protonation = False'")
                
        
    else:
        #is smiles
        try:
            protonation_states = RDTools.protonate(conf_file)
            print('Protonation states identified:')
            for state in protonation_states:
                print(state)
            
                
        except Exception as e:
            Messages.fatal_exception('chemem.protonation', e)
            

def main():
    print(Messages.intro(ChemEM.__version__))
    if len(sys.argv) != 2:
        print("Usage: chemem <config_file>")
        sys.exit(1)
    
    conf_file = sys.argv[1]
    if conf_file in ['-h', '--help']:
        print('Usage chemem <config file>')
        sys.exit(1)
    try:
        print(Messages.load_config())
        # read conf file
        config = Config()
        new_system = config.load_config(conf_file)
        
        
    except Exception as e:
        print(Messages.fatal_exception('Main', e))
        sys.exit(1)
    #run ChemEM
    new_system.run()

if __name__ == "__main__":
    #export_simulation()
    main()
    
    