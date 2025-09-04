# ChemEM

ChemEM fits small molecules into Cryo-EM maps using difference-map–guided docking and optional MD-based refinement.

- **Python:** ≥ 3.11 (recommended: **3.11.9**)
- **OS:** Linux/macOS; on Windows use **WSL/WSL2**
- **Docs:** full guide, tutorials, and FAQs at **https://chemem.topf-group.com**

## Installation

### Fast (mamba, recommended)
```bash
# install Mambaforge if needed, then:
mamba create -n chemem python=3.11.11
mamba activate chemem
mamba install -c conda-forge chemem
```

### Conda
```bash
conda create -n chemem python=3.11.11
conda activate chemem
conda install -c conda-forge chemem
```

### Windows
Use **WSL/WSL2 (Ubuntu)**, then follow the Linux instructions above.

## Quick start

ChemEM runs from a **configuration file**. Minimum example:

```ini
# ChemEM config file
protein   = ~/test_data/7jjo_protein.pdb
ligand    = CC(C)NC[C@@H](c1ccc(c(c1)O)O)O
centroid  = (134.087, 133.507, 174.180)
output    = ~/some_directory
densmap   = ~/test_data/7jjo.mrc
resolution= 2.6

# enabled stages (1 = enable, 0 = disable)
pre_process               = 1
pre_process_split_density = 1
auto_split_point          = 1
auto_split_zone           = 0
fitting                   = 1
dock_only                 = 0
post_process              = 0
rescore                   = 0
```

Then run the appropriate ChemEM command for your workflow (see **Usage** below).

## Usage overview

### Stages (toggle via config)
- **pre_process** – Prepares inputs (always enable unless loading a saved ChemEM object).
  - **pre_process_split_density**: split the difference map into disconnected regions and label them.
  - **auto_split_point**: pick the disconnected region closest to the **centroid**.
  - **auto_split_zone** (+ `auto_split_radius`): keep regions whose centers fall within the radius of **centroid**.
- **fitting** – ChemEM/ChemDock Stage 1 (difference-map–guided docking).
- **dock_only** – Dock without Cryo-EM data (ChemDock score only).
- **post_process** – ChemEM/ChemDock Stage 2 (MD simulated annealing of top N; set `post_process_num_solutions`).
- **rescore** – Rescore one or more solutions; writes `rescore.txt` at the output root.

### Extra inputs (optional)
- **Multiple ligands**: add multiple `ligand = ...` lines **or** use `ligands_from_dir = ~/test_data/`
- **System ligands (fixed)**: `system_ligand_file = sys_lig.sdf` **or** `system_ligands_from_dir = ./sys_lig_dir`
- **User-provided difference map**: `difference_map = /path/to/map.mrc`
  - *If set, ChemEM ignores auto difference-map generation and auto-split options.*

### Protonation, rings, chirality
- **Protonation (Dimorphite-DL)**
  - CLI helper: `chemem.protonate <config_file | 'SMILES'>`
  - Example options: `pH = [6.4, 7.4]`, `pKa_prec = 1.0`
  - To skip auto-protonation once finalized: `protonation = 0`
- **Rings from SMILES (RDKit)**: disable if needed with `rings = 0`
- **Chirality**: ensure all centers are defined; to disable auto-assignment: `chirality = 0`

### Force fields (defaults shown)
```ini
forcefield = amber14        # required when specifying force fields
Protein    = amber14/protein.ff14SB.xml
Ions       = amber14/tip3pfb.xml
Implicit   = amber14/tip3pfb.xml
# Supported ions: amber14/tip3p.xml, amber14/tip3pfb.xml, amber14/tip4pew.xml, amber14/tip4pfb.xml, amber14/spce.xml, amber14/opc.xml, amber14/opc3.xml
# Supported implicit: implicit/hct.xml, implicit/obc1.xml, implicit/obc2.xml, implicit/gbn.xml, implicit/gbn2.xml
```

## Alpha-mask (v0.0.6+)

Generate a map of **unaccounted densities** (α-sphere SASA masking):
```bash
chemem.alpha_mask <config_file>
```
Outputs include masked maps of candidate ligand density and segmented centroids.

## Documentation & support

- Full docs, parameter reference, and examples: **https://chemem.topf-group.com**
- Please open issues and feature requests in this repository.

## Citation

Sweeney A., Mulvaney T., Maiorca M., & Topf M. (2024). *ChemEM: flexible docking of small molecules in Cryo-EM structures using difference maps.* **J. Med. Chem.** 67(1), 199–212.

---

*Tested best with Python 3.11.11. ChemEM requires Python ≥ 3.11.*
