# DynaBench
The Dynabench package consists of 2 parts: analysis and visualization. Quality control, Residue-based, and Interaction-based analysis are performed. RMSD, RG, and RMSF analysis under the Quality Control; SASA, biophysical type, and energy analysis under the Residue-Based; hydrogen, hydrophobic, and ionic bond analysis are performed under the Interaction-Based analysis. Outputs are in the CSV form under the tables folder. The package accepts .pdb and .dcd inputs. Outputs from Quality Control, Residue Based and Interaction Based analysis are used for visualization.

## Usage
### System Dependencies
* python3 (3.8 or higher) or anaconda3
### Python Dependencies
* ipywidgets
* ipykernel
* pdbfixer
* seaborn
* mdanalysis
* jsonpickle
* moviepy
* matplotlib
* pdb-tools
* freesasa

### Installation
#### Note for Windows users:
  **DynaBench uses EvoEF for energy calculations. Before installing the package, please install the EvoEF (https://github.com/tommyhuangthu/EvoEF) by yourself and configure to your laptop.After installing the EvoeF, you can follow the installation steps.**
  
**Windows users may need to download the WSL and, C++ compiler. The programs can be downloaded by the following links:**

  C++ Compiler:https://visualstudio.microsoft.com/visual-cpp-build-tools/
    
  WSL: https://learn.microsoft.com/en-us/windows/wsl/install

**We recommend you use git clone instead of download with .zip while installing the package. Downloading with .zip may cause problems in test pdbs. If you don't have one, please download the git from here:https://git-scm.com/downloads**

### Clone the repository
```
git clone https://github.com/Atakanzsn/DynaBench.git
```
```
cd DynaBench
```
### Create DynaBench environment
```
conda env create -f requirements.yml
```
```
conda activate DynaBench
```
### Create IpyKernel for DynaBench
```
python -m ipykernel install --user --name dynabench
```
### Build the DynaBench
```
python setup.py build
python setup.py install
```
### Build EvoEF for energy analysis
#### Linuc or MacOS users
```
python build_evoef.py
```
#### Windows users
```
python build_evoef.py --evoef_path=absolute_path_of_EvoEF_folder
```
## Tests
Please run the jupyter notebooks in the **tests** folder to test the package. They will create folders with *tests* at the end. Check outputs with the folder without the *tests* at the end. Outputs should be identical.

## Run
You can run from a script or shell. To run from the script, you can check the test jupyter notebooks. To run from the terminal, you can either run with JSON file or direct commands. If running from the terminal, please do not run inside the DynaBench folder. It may cause import errors.

### Run from terminal
Complete run (Quality control, Residue-Based, and Interaction-Based analysis with all visualizations) with input file from the terminal:
#### With .pdb input
```
dynabench -input_file=input_file.pdb --commands=all_analysis,all_plots
```
#### With .dcd input
```
dynabench -input_file=input_file.dcd --commands=all_analysis,all_plots --dcd_pdb=input_pdb.pdb
```

#### To run with JSON file from the terminal:
```
dynabench --table_json=table_params.json --plot_json=plot_params.json
```

#### For more commands to run from terminal: 
```
dynabench -h
```
