# DynaBench
Dynabench package consists of 2 parts. Analysis and Vizualization. Quality control, Residue based and Interation based analysis are performed. RMSD, RG, and RMSF analysis under the Quality Control; SASA, biophysical type, and ebergy analysis under the Residue Based; hydrogen, hydrophobic and ionic bond analysis are performed under the Interaction Based analysis. Outputs are in the csv form under the *tables* folder. Package accepts .pdb and .dcd inputs.
For visualization, outputs from Quality Control, Residue Based and Interaction Based analysis are used.
## Installation
To install to package, inside the DynaBench folder, run the following commands in the shell.
```
conda env create -f requirements.yml
conda activate DynaBench
python -m ipykernel install --user --name dynabench
python setup.py build
python setup.py install
pytjon build_evoef.py
```
## Tests
To test the package, please run the jupyter notebooks in the **tests** folder. They will create folders with *tests* at the end. Check outputs with the folder without the *tests* at the end. Outputs should be identical.

## Run
You can run from script or shell. To run from script, you can check the test jupyter notebooks. To run from terminal, you can either run json file or direct commands. If you are running from the terminal, please do not run inside the DynaBench folder. It may cause import errors.

### Run from terminal
Complete run (Quality control, Residue Based and Interaction based analysis with all visualizations) with input file from terminal:
```
# for .pdb input
dynabench -input_file=input_file.pdb --commands=all_tables,all_plots

#for .dcd input
dynabench -input_file=input_file.dcd --commands=all_tables,all_plots --dcd_pdb=input_pdb.pdb
```
To run with json file from terminal:
```
dynabench --table_json=table_params.json --plot_json=plot_params.json
```

For help in running from terminal with direct commands: 
```
dynabench -h
```
