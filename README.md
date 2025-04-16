# ChemEM

## Description

**ChemEM** is a software package designed for fitting small molecules into Cryo-Electron Microscopy (Cryo-EM) data. It helps place small molecule structures accurately into electron density maps, supporting research in structural biology and drug discovery.

## Installation

### Recommended Installation (via Conda)

The recommended method is to install ChemEM using conda, which ensures all dependencies are met. For best results, please use a Python 3.11 version prior to 3.11.12.

1. **Create a new conda environment with Python 3.11.11:**

   
  <code> conda create -n chemem_env python=3.11.11
   conda activate chemem_env </code> 


2. **Install ChemEM from conda-forge:**

   <code> conda install -c conda-forge chemem </code>

2.1 **Alternatively, if you prefer to install from the ChemEM Anaconda channel:**

   <code> conda install -c conde-forge -c chemem chemem

## Additional Documentation and Citation 

For more detailed documentation, including downloads, tutorials, and usage examples, please visit the [ChemEM Documentation Page](https://chemem.topf-group.com/download.html).

If you use ChemEM in your research, please cite our work: 

Sweeney A., Mulvaney T., Maiorca M., & Topf M. (2024). ChemEM: flexible docking of small molecules in Cryo-EM structures using difference maps. Journal of Medicinal Chemistry, 67(1), 199-212 
