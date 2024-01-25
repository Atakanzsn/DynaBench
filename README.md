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
  >**DynaBench uses EvoEF for energy calculations. Before installing the package, please install the [EvoEF](https://github.com/tommyhuangthu/EvoEF) by yourself and configure to your laptop. After installing the EvoeF, you can follow the installation steps.**
>**Windows users also may need to download the [WSL](https://learn.microsoft.com/en-us/windows/wsl/install) and, [C++ compiler](https://visualstudio.microsoft.com/visual-cpp-build-tools/).**

**We recommend you use [git](https://git-scm.com/downloads) clone instead of download with .zip while installing the package. Downloading with .zip may cause problems in test pdbs.**

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
You can run from a script or shell. To run from the terminal, you can either run with JSON file or you can choose options from the shell. If running from the terminal, please do not run inside the DynaBench folder. It may cause import errors.
#### Example Script
```
import DynaBench

#load trajectory .dcd and topology .pdb file
mol = DynaBench.dynabench(inp_file='your_trajectory_path', dcd_pdb='your_topology_path',
 split_models=False, show_time_as='Frame')

#define plotter class
draw = DynaBench.Plotter(job_name=mol.job_name)

#Quality control analysis and visualization
mol.run_quality_control(rmsd_data={'ref_struc':None, 'ref_frame':0})

draw.plot_rmsd(path=None)
draw.plot_rg(path=None)
#RMSF plot can be plotted after residue-based analysis due to marking of interface residues on the plot.

#Residue Based analysis and visualization
mol.run_res_based()

draw.plot_rmsf(rmsf_path=None, intf_path=None)
draw.plot_int_energy(thereshold=50.0, res_path=None, intf_path=None)
draw.plot_biophys(path=None)

#Interaction Based analysis and visualization
mol.run_inter_based()

draw.plot_bond_freq_barplot(path=None)

#Get .json files for analysis and plots.
mol._get_params_()
draw._get_params_()
```

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
