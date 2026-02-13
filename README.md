[![GitHub Org](https://img.shields.io/badge/GitHub-HESCOR-blue?logo=github&logoColor=white)](https://github.com/HESCOR)

# HEP-Model – Human Existence Potential Model

## Overview

The **Human Existence Potential Model (HEP)** quantifies the capacity of a society, within a specific cultural and environmental context, to sustain human life. Inspired by **species distribution modelling**, the model integrates archaeological presence/absence data with climate variables in a **logistic regression framework**.

The HEP-Model uses:
- Archaeological site locations as presence and absence points  
- BioClim variables to describe climate in a biologically meaningful way  

The main output of the model is a spatially explicit **HEP map**, representing the relative potential for human existence across a study region.


## Scientific Context: HESCOR

The HEP-Model was developed within the **HESCOR project** at the University of Cologne.

HESCOR is an interdisciplinary research initiative investigating how environmental and earth system processes interact and co-evolve with human societies. The project aims to develop integrated frameworks for:

- Modeling human–environment feedbacks  
- Assessing data limitations and uncertainties  
- Translating insights between natural and social sciences  

More information: https://www.hescor-project.com/



## Repository Contents

- `README.md` – Project overview and usage instructions  
- `pyvenv_list.txt` – Python package dependencies  
- `setup.sh` – Script for creating a local virtual environment  
- `configure.py` – Configuration settings for model execution  
- `ehep_run.py` – Main execution script  
- `ehep_methods.py` – Model calculation functions  
- `ehep_inout.py` – Input/output handling functions  
- `ehep_util.py` – Utility functions  
- `plot_hep.py` – Script for plotting the HEP map  
- `routine_overview.md` – Detailed overview of the model workflow  

## Prerequisites

- Python **3.8** or newer  
- Unix-like operating system (Linux or macOS recommended)  
- Sufficient disk space and memory for raster climate data  
- External datasets (e.g. BioClim variables, archaeological site data)

Further details on data preparation and assumptions are described in `routine_overview.md`.

## Quick Setup

A setup script is provided to create an isolated Python virtual environment and install all required dependencies.

### 1. Create the virtual environment

Run the provided setup script to create a local virtual environment and install dependencies:

```bash
./setup.sh
```

The script will:

* Create (or reuse) a `venv/` directory in the project root
* Upgrade `pip` inside the environment
* Install all packages listed in `pyvenv_list.txt`

### 2. Activate the environment

After the setup completes, activate the environment manually:

```bash
source venv/bin/activate
```

To deactivate later:

```bash
deactivate
```


## Usage

After activating the virtual environment, the model can be executed via:

```bash
python3 ehep_run.py
```

Model configuration parameters (e.g. input data paths, output locations) are defined in `configure.py`.

The resulting HEP maps can be visualized using:

```bash
python3 plot_hep.py
```



## Outputs

Typical outputs include:

* Raster maps representing Human Existence Potential
* Intermediate data products generated during model execution
* Visualizations of model results

Output formats and locations are configurable via `configure.py`.



## Related Publications 
- Shao Y, Klein K, Wegener C, Weniger G-C (2025) Pathways at the Iberian crossroads: Dynamic modeling of the middle–upper paleolithic transition. PLoS One 20(12): e0339184. https://doi.org/10.1371/journal.pone.0339184
- Shao, Y., Wegener, C., Klein, K. et al. Reconstruction of human dispersal during Aurignacian on pan-European scale. Nat Commun 15, 7406 (2024). https://doi.org/10.1038/s41467-024-51349-y
- Yaping Shao, Heiko Limberg, Konstantin Klein, Christian Wegener, Isabell Schmidt, Gerd-Christian Weniger, Andreas Hense, Masoud Rostami, Human-existence probability of the Aurignacian techno-complex under extreme climate conditions, Quaternary Science Reviews, Volume 263, 2021, 106995, ISSN 0277-3791, https://doi.org/10.1016/j.quascirev.2021.106995.
- Rostami, M., Klein, K., Wegener, C., Shao, Y., and Weniger, G.-C.: Impacts of Heinrich events upon Human existence potential in Europe, EGU General Assembly 2020, Online, 4–8 May 2020, EGU2020-7778, https://doi.org/10.5194/egusphere-egu2020-7778, 2020.

## Citation

If you use this code in academic work, please cite the associated publication and reference the HESCOR project.


## License

This project is distributed under the terms of the MIT License. See [LICENSE](LICENSE) for details.
