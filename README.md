# DynaBench
The Dynabench package consists of 2 parts: analysis and visualization. Quality control, Residue-based, and Interaction-based analysis are performed. RMSD, RG, and RMSF analysis under the Quality Control; SASA, biophysical type, and energy analysis under the Residue-Based; hydrogen, hydrophobic, and ionic bond analysis are performed under the Interaction-Based analysis. Outputs are in the CSV form under the tables folder. The package accepts .pdb and .dcd inputs. For visualization, outputs from Quality Control, Residue Based and Interaction Based analysis are used.
## Installation
**Windows users may need to download C++ compiler. You can download from here: https://visualstudio.microsoft.com/visual-cpp-build-tools/

Also windows users should install WSL, ones can install by the following link:https://learn.microsoft.com/en-us/windows/wsl/install**

**We recommend you to use git clone instead of download with .zip while installing the package. Downloading with .zip may cause problems in test pdbs.**

**If you don't have, plase download the git from here:https://git-scm.com/downloads**

To install the package, inside the DynaBench folder, run the following commands in the shell.
```
git clone https://github.com/Atakanzsn/DynaBench.git
cd DynaBench
conda env create -f requirements.yml
conda activate DynaBench
python -m ipykernel install --user --name dynabench
python setup.py build
python setup.py install
python build_evoef.py
```
## Tests
Please run the jupyter notebooks in the **tests** folder to test the package. They will create folders with *tests* at the end. Check outputs with the folder without the *tests* at the end. Outputs should be identical.

## Run
You can run from a script or shell. To run from the script, you can check the test jupyter notebooks. To run from the terminal, you can either run with JSON file or direct commands. If running from the terminal, please do not run inside the DynaBench folder. It may cause import errors.

### Run from terminal
Complete run (Quality control, Residue-Based, and Interaction-Based analysis with all visualizations) with input file from the terminal:
```
# for .pdb input
dynabench -input_file=input_file.pdb --commands=all_analysis,all_plots

#for .dcd input
dynabench -input_file=input_file.dcd --commands=all_analysis,all_plots --dcd_pdb=input_pdb.pdb
```
To run with JSON file from the terminal:
```
dynabench --table_json=table_params.json --plot_json=plot_params.json
```

For help in running from the terminal with direct commands: 
```
dynabench -h
```
