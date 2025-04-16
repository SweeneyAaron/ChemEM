# This file is part of the ChemEM software.
#
# Copyright (c) 2023 - Topf Group & Leibniz Institute for Virology (LIV),
# Hamburg, Germany.
#
# This module was developed by:
#   Aaron Sweeney    <aaron.sweeney AT cssb-hamburg.de>

class Messages:
    
    @staticmethod
    def intro(version):
        return f"""
*******************************************************************************
*                              ChemEM {version}                                   *
*                                                                             *
* Copyright (c) 2023 - Topf Group & Leibniz Institute for Virology (LIV),     *
* Hamburg, Germany.                                                           *
* Developed by Aaron Sweeney <aaron.sweeney AT cssb-hamburg.de>               *
*                                                                             *
* If you use ChemEM in your work, please cite:                                *
* ChemEM: Flexible Docking of Small Molecules in Cryo-EM Structures           *
* Aaron Sweeney, Thomas Mulvaney, Mauro Maiorca, and Maya Topf                *
* Journal of Medicinal Chemistry 2024, 67 (1), 199-212                        *
* DOI: 10.1021/acs.jmedchem.3c01134                                           *
*******************************************************************************
\n"""

    @staticmethod 
    def test():
        return '''
              Deprecation Warning for ChemEM v0.0.2

    Please note that the test function has been deprecated as of ChemEM version 0.0.2. \n
    We now require users to download the latest test data and example configuration files.

    For your convenience, these files can be accessed and downloaded from the following URL:\n
    https://gitlab.com/topf-lab/chemem_test_data/-/archive/main/chemem_test_data-main.zip \n
    Alternatively, if you are familiar with Git, you can clone the repository using the following command in your terminal or command prompt:\n
     git clone https://gitlab.com/topf-lab/chemem_test_data.git'''

    @staticmethod 
    def load_config():
        return """
╔═════════════════════════════════════════════════════════════════════════════╗
║                             Loading Config Data                             ║
╚═════════════════════════════════════════════════════════════════════════════╝
"""

    @staticmethod 
    def preprocess():
        return """
╔═════════════════════════════════════════════════════════════════════════════╗
║                                Pre-Processing                               ║
╚═════════════════════════════════════════════════════════════════════════════╝
"""
    @staticmethod 
    def difference_map():
        return """
╔═════════════════════════════════════════════════════════════════════════════╗
║                              Difference mapping                             ║
╚═════════════════════════════════════════════════════════════════════════════╝
"""
    @staticmethod 
    def split_density():
        return """
╔═════════════════════════════════════════════════════════════════════════════╗
║                                  Split Density                              ║
╚═════════════════════════════════════════════════════════════════════════════╝
"""

    @staticmethod 
    def auto_split_point():
        return """
╔═════════════════════════════════════════════════════════════════════════════╗
║                                Auto-split Point                             ║
╚═════════════════════════════════════════════════════════════════════════════╝
"""

    @staticmethod 
    def auto_split_zone():
        return """
╔═════════════════════════════════════════════════════════════════════════════╗
║                                Auto-split Zone                              ║
╚═════════════════════════════════════════════════════════════════════════════╝
"""

    @staticmethod 
    def fitting():
        return """
╔═════════════════════════════════════════════════════════════════════════════╗
║                                    Fitting                                  ║
╚═════════════════════════════════════════════════════════════════════════════╝
"""
    @staticmethod 
    def post_processing():
        return """
╔═════════════════════════════════════════════════════════════════════════════╗
║                               Post-Processing                               ║
╚═════════════════════════════════════════════════════════════════════════════╝
"""
    @staticmethod 
    def rescore():
        return """
╔═════════════════════════════════════════════════════════════════════════════╗
║                                   Rescore                                   ║
╚═════════════════════════════════════════════════════════════════════════════╝
"""
    @staticmethod
    def system_repr(system):
        return f'System(num_proteins={len(system.proteins)}, num_ligands={len(system.ligands)}, num_system_ligands={len(system.system_ligands)} num_maps={len(system.maps)})'
    
    @staticmethod 
    def no_hold_fragment():
        return 'Hold Fragments option is not currently enabled with fitting solutions.'
    
    @staticmethod 
    def chemem_warning(class_name,func_name, error):
        return f"ChemEM-Warning: Non-fatal error in {class_name}, skipping running {func_name}. \nFull Error: {error}"
    @staticmethod 
    def fatal_exception(class_name, error):
        return f"ChemEM-Error: Fatal error in {class_name}. \n Full Error: {error}"
    
    @staticmethod 
    def no_class_attribute_found(class_name, attr):
        return f"ChemEM-AttrError: '{class_name}' object has no attribute '{attr}'"
    
    @staticmethod 
    def overwrite(path, overwrite):
        if overwrite:
            print( f"""ChemEM-Warning: Overwriting data in {path}. 
               To stop automatic overwriting, set overwrite = 0 in the configuration file.""")
        else:
            print(f""""ChemEM-OverWriteError: {path} exists.
                  To avoid accidentally overwriting files please choose another directory or set overwrite to 1 in your configuration file""")
    @staticmethod
    def protonation_usage():
        
        print("""Usage: 
              \tchemem.protonate <config_file>
              \tchemem.protonate <ligand smiles>""")
        
        