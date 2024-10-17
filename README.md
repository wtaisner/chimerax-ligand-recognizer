[![Generic badge](https://img.shields.io/badge/Jingle-Awesome-<COLOR>.svg)](https://www.udio.com/songs/sDBuoBLtEKxWree8wvSSjh)
[![Streamlit - Demo](https://img.shields.io/badge/Streamlit-Demo-green)](https://ligands.cs.put.poznan.pl)
[![bioRxiv - Preprint](https://img.shields.io/badge/bioRxiv-Preprint-red)](https://www.biorxiv.org/content/10.1101/2024.08.27.610022v1)

# ChimeraX bundle for density map classification.

This bundle provides tools for classifying density maps in ChimeraX. Refer to `bundle_info.xml` for more information. 
This work is based on [this repository](https://github.com/jkarolczak/ligand-classification) and the corresponding publication.

## Usage
This tool assumes, that both PDB structure and density map are loaded into ChimeraX.
i.e. `open 8SOR/8sor.cif; open 8SOR/map_model_difference_1.ccp4;`
Then, it can be run from the command line by pointing to a specific ligand in the PDB structure.

The tool implements new commands:
- blob validate CHAIN: validates class of a ligand -- (e.g. `blob validate /A:1401`)
- blobus validatus: validates class of a ligand but with more magic involved (functionality is the same though) -- (e.g. `blobus validatus /A:1401`)

## Installation

TODO: add from toolshed once it's available.

You can install the bundle manually by downloading the repository and running the following commands in the repository directory:
```bash
# set ChimeraX alias -- change this to your ChimeraX path
alias chimerax=/usr/bin/chimerax
# build the package
chimerax --nogui --cmd 'devel build . exit true'
# install the package
chimerax --nogui --cmd 'devel install . exit true'
```

## Citation
```bibtex
@article {Karolczak2024.08.27.610022,
	author = {Karolczak, Jacek and Przyby{\l}owska, Anna and Szewczyk, Konrad and Taisner, Witold and Heumann, John M. and Stowell, Michael H.B. and Nowicki, Micha{\l} and Brzezinski, Dariusz},
	title = {Ligand Identification using Deep Learning},
	elocation-id = {2024.08.27.610022},
	year = {2024},
	doi = {10.1101/2024.08.27.610022},
	publisher = {Cold Spring Harbor Laboratory},
	abstract = {Motivation Accurately identifying ligands plays a crucial role in the process of structure-guided drug design. Based on density maps from X-ray diffraction or cryogenic-sample electron microscopy (cryoEM), scientists verify whether small-molecule ligands bind to active sites of interest. However, the interpretation of density maps is challenging, and cognitive bias can sometimes mislead investigators into modeling fictitious compounds. Ligand identification can be aided by automatic methods, but existing approaches are available only for X-ray diffraction and are based on iterative fitting or feature-engineered machine learning rather than end-to-end deep learning.Results Here, we propose to identify ligands using a deep learning approach that treats density maps as 3D point clouds. We show that the proposed model is on par with existing machine learning methods for X-ray crystallography while also being applicable to cryoEM density maps. Our study demonstrates that electron density map fragments can be used to train models that can be applied to cryoEM structures, but also highlights challenges associated with the standardization of electron microscopy maps and the quality assessment of cryoEM ligands.Availability Code and model weights are available on GitHub at https://github.com/jkarolczak/ligands-classification. Datasets used for training and testing are hosted at Zenodo: 10.5281/zenodo.10908325.Contact dariusz.brzezinski{at}cs.put.poznan.plCompeting Interest StatementThe authors have declared no competing interest.},
	URL = {https://www.biorxiv.org/content/early/2024/08/28/2024.08.27.610022},
	eprint = {https://www.biorxiv.org/content/early/2024/08/28/2024.08.27.610022.full.pdf},
	journal = {bioRxiv}
}
```