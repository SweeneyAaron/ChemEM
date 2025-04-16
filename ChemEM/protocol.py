# This file is part of the ChemEM software.
#
# Copyright (c) 2023 - Topf Group & Leibniz Institute for Virology (LIV),
# Hamburg, Germany.
#
# This module was developed by:
#   Aaron Sweeney    <aaron.sweeney AT cssb-hamburg.de>

from ChemEM.protocols.preprocess import PreProcess
from ChemEM.protocols.difference_map import DifferenceMap
from ChemEM.protocols.difference_map import DifferenceMapOriginal
from ChemEM.protocols.split_density import SplitDensity
from ChemEM.protocols.auto_split_point import AutoSplitPoint
from ChemEM.protocols.auto_split_zone import AutoSplitZone
from ChemEM.protocols.fitting import Fitting 
from ChemEM.protocols.post_processing import PostProcessing
from ChemEM.protocols.rescore import Rescore

def protocol_selector(protocol_name):
    if protocol_name == 'pre_process':
        return PreProcess
    elif protocol_name == 'difference_map':
        return DifferenceMap
    
    elif protocol_name == 'difference_map_original':
        return DifferenceMapOriginal
    
    elif protocol_name == 'pre_process_split_density':
        return SplitDensity
    
    elif protocol_name == 'auto_split_point':
        return AutoSplitPoint
    
    elif protocol_name == 'auto_split_zone':
        return AutoSplitZone
    
    elif protocol_name == 'fitting':
        return Fitting
    
    elif protocol_name == 'post_process':
        return PostProcessing
    
    elif protocol_name == 'rescore':
        return Rescore
    
    
        
    
    