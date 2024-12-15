import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from DynaBench.handling import tables_errors
import json
from matplotlib import cm




class Plotter:
    def __init__(self, job_name=None):
        """A class to visualize several analyses from dynabench class. Reads csv files and performs visualizations. Saves figures into the folder named 'figures', with 300 dpi and .png extension.
        
        Return: return_description
        """
        self.handler = tables_errors()

        #params
        self._job_name = job_name
        self._rmsd = False
        self._rmsd_path = None
        self._rg = False
        self._rg_path = None
        self._rmsf = False
        self._rmsf_path = None
        self._rmsf_intf_path = None
        self._biophys = False
        self._biophys_path = None
        self._biophys_palette = None
        self._pie = False
        self._pie_hbond = None
        self._pie_hph = None
        self._pie_ionic = None
        self._pie_path = None
        self._pie_palette = None
        self._bar = False
        self._bar_path = None
        self._bar_palette = None
        self._ene = False
        self._ene_thr = None
        self._ene_intf = None
        self._ene_path = None
        self._plot_SASA = False
        self._sasa_path = None
        self._plot_irmsd = False
        self._irmsd_path = None
        self._plot_fnonnat = False
        self._fnonnat_path = None
        self._plot_lrmsd = False
        self._lrmsd_path = None
        self._plot_dockq = False
        self._dockq_path = None
        self._plot_fnat = False
        self._fnat_path = None
        self._plot_dssp = False
        self._dssp_path = None
        self._dssp_intf_path = None
        self._dssp_thereshold = 50.0

        sns.set_style('whitegrid')

        if job_name is None:
            job_name = input('Plase provide job name, or existing job file where analysis folder is located:\n')
            self.job_path = job_name

        else:
            self.job_path = os.path.join(os.getcwd(), job_name)

        self.target_path = os.path.join(self.job_path, "figures")
        self.table_path = os.path.join(self.job_path, "tables")
        if not os.path.exists(self.target_path):
            os.mkdir(self.target_path)

        self.chain_colors = ['#009E73', '#E69F00', '#F0E442', '#CC79A7']
        self.complex_color = '#c26a77'

    def plot_rmsd(self, path=None):
        """A function to perform line plot RMSD visualization for each chain and the overall complex. 'Reads QualityControl-overtime.csv' file.
        
        Return: None
        """
        
        
        if path:
            self.handler.test_inp_path(path)
            df = pd.read_csv(path)
        else:
            df = pd.read_csv(os.path.join(self.table_path, "QualityControl-overtime.csv"))

        self._rmsd = True
        self._rmsd_path = path

        time_name = df.columns[0]
        time_var = df.iloc[:, 0].values


        cols_use = [x for x in df.columns.tolist() if "RMSD" in x]
        chains = cols_use[:-1]
        complex = cols_use[-1]

        fig, ax = plt.subplots(figsize=(5, 2.7), layout='constrained')
        for i, col in enumerate(chains):
            mol_name = col.split("RMSD")
            ax.plot(time_var, df[col], label=mol_name[0], color=self.chain_colors[i])
        
        ax.plot(time_var, df[complex], label='Complex', color=self.complex_color)
        ax.set_xlabel(time_name)
        ax.set_ylabel(f'RMSD (Å)')

        ax.set_title(f"RMSD Analysis")
        ax.legend()
        fig.savefig(os.path.join(self.target_path, f'RMSD_Analysis.png'), dpi=300)

    def plot_rg(self, path=None):
        """A function to perform line plot RG visualization for each chain and the overall complex. 'Reads QualityControl-overtime.csv' file.
        
        Return: None
        """
        if path:
            self.handler.test_inp_path(path)
            df = pd.read_csv(path)
        else:
            df = pd.read_csv(os.path.join(self.table_path, "QualityControl-overtime.csv"))

        self._rg = True
        self._rg_path = path

        cols_use = [x for x in df.columns.tolist() if "RDG" in x]
        chains = cols_use[:-1]
        complex = cols_use[-1]

        time_name = df.columns[0]
        time_var = df.iloc[:, 0].values

        fig, ax = plt.subplots(figsize=(5, 2.7), layout='constrained')
        for i, col in enumerate(chains):
            mol_name = col.split("RDG")
            ax.plot(time_var, df[col], label=mol_name[0], color=self.chain_colors[i])

        ax.plot(time_var, df[complex], label='Complex', color=self.complex_color)
        ax.set_xlabel(time_name)
        ax.set_ylabel(f'RDG (Å)')

        ax.set_title(f"RG Analysis")
        ax.legend()
        fig.savefig(os.path.join(self.target_path, f'RG_Analysis.png'), dpi=300)

    def plot_rmsf(self, rmsf_path=None, intf_path=None):
        """A function to perform line plot RMSF visualization for each chain, with interface residues marked. Reads 'interface_label_perc.csv' and 'QualityControl-overres' files.
        
        Return: None
        """

        if rmsf_path:
            self.handler.test_inp_path(rmsf_path)
            df1 = pd.read_csv(rmsf_path)

        if intf_path:
            self.handler.test_inp_path(intf_path)
            df = pd.read_csv(intf_path)

        else:
            df1 = pd.read_csv(os.path.join(self.table_path, "QualityControl-overres.csv"))
            if os.path.exists(os.path.join(self.table_path, "interface_label_perc.csv")):
                intf = True
                df = pd.read_csv(os.path.join(self.table_path, "interface_label_perc.csv"))
            else:
                intf=False

        self._rmsf = True
        self._rmsf_path = rmsf_path
        self._rmsf_intf_path = intf_path

        if intf:

            int_df = df[(df["Interface Label"] == 4) | (df["Interface Label"] == 2) | (df["Interface Label"] == 3)]

            int_df = int_df[int_df["Percentage"] >= 50]
            
            g2 = int_df.groupby(["Chain"])

        max_y = np.max(df1['RMSF'])

        groups = df1.groupby(["Molecule"])
        
        fig, axes = plt.subplots(1, len(groups), figsize=(10, 6), sharey=True)
        
        def get_jumpseq(data):
            sequence = list()
            res_list = data["Residue Number"].tolist()
            for index,row in data.iterrows():
                resnum = row['Residue Number']
                
                res_num_idx = res_list.index(resnum)
                latter = res_list[res_num_idx-1]
                diff = resnum - latter
                if diff != 1:
                    for i in range(diff-1):
                        sequence.append('.')
                sequence.append(resnum)
            return sequence
        
        def get_jumped_markers(data, seq_list, jumped_data):
            markers = list() #interface data alıyor
            for x in data['Residue']:
                res_num = int(x[3:])
                try:
                    ix = jumped_data[jumped_data['Residue Number'] == res_num].index[0]
                except: continue
                markers.append(ix)

            return markers


        for ind, (ax, x) in enumerate(zip(axes, np.unique(df1["Molecule"]))):
            data = groups.get_group(x)
            first_resnum = data.iloc[0]['Residue Number']
            last_resnum = data.iloc[-1]['Residue Number']
            resnum_len = last_resnum - first_resnum

            if resnum_len+1 != len(data): #check the jump

                sequence = get_jumpseq(data)

                ax.set_xlim([1,len(sequence)])

                jumped_data = pd.DataFrame(columns=['Molecule', 'Residue Number', 'RMSF'])
                colss = []
                idxs = []

                overall_list = []

                for idx,res in enumerate(sequence):
                    if res != ".":

                        if len(idxs) > 1:
                            colss.append(idxs)
                            idxs = []

                        jumped_data.loc[len(jumped_data)] = data[data['Residue Number']== res].values[0]
                        complete_end=False
                        part_end = False

                        try:
                            if sequence[idx+1] == '.':
                                part_end = True
                        except IndexError:
                            complete_end=True

                        if complete_end or part_end:
                            if intf:
                                int_data = g2.get_group(x)
                                markers = get_jumped_markers(int_data, sequence, jumped_data)
                                overall_list.append([jumped_data, markers])

                            else:
                                overall_list.append([jumped_data])

                            jumped_data = pd.DataFrame(columns=['Molecule', 'Residue Number', 'RMSF']) 

                    else:
                        idxs.append(idx + 1)

                for el in overall_list:
                    if len(el) > 1:                        
                        ax.plot(el[0]["Residue Number"], el[0]["RMSF"], '-o', markevery=el[1], markersize=3.5,
                            color=self.chain_colors[ind], linewidth=3)
                        ax.plot(el[0]["Residue Number"], el[0]["RMSF"], 'o', markevery=el[1],
                            label='Interface Residues', markersize=3.5, color="red", linewidth=3)
                    else:
                        data = el[0]
                        ax.plot(data["Residue Number"], data["RMSF"],
                                color=self.chain_colors[ind], linewidth=3)

                for c in colss:
                    ax.axvspan(np.min(c), np.max(c), facecolor="black", hatch="///", edgecolor="white", linewidth=0.0, alpha=0.45)


            else:
                if intf: #if no jmp, markers, intf yes

                    int_data = g2.get_group(x)
                    
                    res_list = int_data["Residue"].tolist()
                    markerss = [int(x[3:]) -1 for x in res_list]
                    l = data["Residue Number"].tolist()
                    markers = [l.index(i) for i in markerss]

                    ax.plot(data["Residue Number"], data["RMSF"],markevery=markers, markersize=3.5,
                            color=self.chain_colors[ind])
                    ax.plot(data["Residue Number"], data["RMSF"], 'o', markevery=markers,
                            label='Interface Residues', markersize=3.5, color="red")  #jump yok intf var
                    
                else:
                    
                    ax.plot(data["Residue Number"], data["RMSF"],
                            color=self.chain_colors[ind]) #jump yok intf yok
                
            ax.set_xlabel("Residue Number")
            ax.set_ylabel(f'RMSF (Å)')
            ax.set_title(f"Chain {x}", fontsize=16)
            ax.xaxis.label.set_size(12)
            ax.yaxis.label.set_size(12)
        plt.suptitle("RMSF Analysis", fontsize=18)
        if intf:
            plt.legend()
        fig.savefig(os.path.join(self.target_path, f'RMSF-proteinCA.png'), dpi=300)

    def plot_irmsd(self, path=None):
        """A function to perform lineplot visualization of interface RMSD through the simulation. Reads 'dockq_results.csv' file.
        
        Return: None
        """

        self._plot_irmsd = True

        if path:
            self.handler.test_inp_path(path)
            df = pd.read_csv(path)
        else:
            df = pd.read_csv(os.path.join(self.table_path, "dockq_results.csv"))

        
        self._irmsd_path = path
        
        d = df.groupby("mapping")

        fig,ax  = plt.subplots(figsize=(5,2.7), layout='constrained')

        for i in d.groups:
            data = d.get_group(i)
            ax.plot(data.columns[0], "iRMSD",data=data, label=i)
        
        ax.set_xlabel(df.columns[0])
        ax.set_ylabel('iRMSD (Å)')
        ax.set_ylim(bottom=0)

        ax.set_title('iRMSD Analysis')
        ax.legend(title="Mapping")
        fig.savefig(os.path.join(self.target_path, f'irmsd_analysis.png'), dpi=300)

    def plot_lrmsd(self, path=None):
        """A function to perform lineplot visualization of ligand RMSD through the simulation. Reads 'dockq_results.csv' file.
        
        Return: None
        """

        self._plot_lrmsd = True

        if path:
            self.handler.test_inp_path(path)
            df = pd.read_csv(path)
        else:
            df = pd.read_csv(os.path.join(self.table_path, "dockq_results.csv"))

        
        self._lrmsd_path = path

        d = df.groupby("mapping")

        fig,ax  = plt.subplots(figsize=(5,2.7), layout='constrained')

        for i in d.groups:
            data = d.get_group(i)
            ax.plot(data.columns[0], "lRMSD",data=data, label=i)

        ax.set_xlabel(df.columns[0])
        ax.set_ylabel('lRMSD (Å)')
        ax.set_ylim(bottom=0)


        ax.set_title('lRMSD Analysis')
        ax.legend(title="Mapping")
        fig.savefig(os.path.join(self.target_path, f'lrmsd_analysis.png'), dpi=300)

    def plot_dockq(self, path=None):
        """A function to perform lineplot visualization of ligand RMSD through the simulation. Reads 'dockq_results.csv' file.
        
        Return: None
        """

        self._plot_dockq = True

        if path:
            self.handler.test_inp_path(path)
            df = pd.read_csv(path)
        else:
            df = pd.read_csv(os.path.join(self.table_path, "dockq_results.csv"))

        
        self._dockq_path = path
        d = df.groupby("mapping")

        fig,ax  = plt.subplots(figsize=(5,2.7), layout='constrained')

        for i in d.groups:
            data = d.get_group(i)
            ax.plot(data.columns[0], "Total", data=data, label=i)

        ax.set_xlabel(df.columns[0])
        ax.set_ylabel('DockQ Score')
        ax.set_ylim(bottom=0)

        ax.legend(title="Mapping")
        ax.set_title('DockQ Score Analysis')
        fig.savefig(os.path.join(self.target_path, f'dockq_score.png'), dpi=300)

    def plot_fnonnat(self, path=None):
        """A function to perform lineplot visualization of fraction of native contacts through the simulation. Reads 'dockq_results.csv' file.
        
        Return: None
        """

        self._plot_fnonnat = True

        if path:
            self.handler.test_inp_path(path)
            df = pd.read_csv(path)
        else:
            df = pd.read_csv(os.path.join(self.table_path, "dockq_results.csv"))

        self._fnonnat_path = path
        d = df.groupby("mapping")

        fig,ax  = plt.subplots(figsize=(5,2.7), layout='constrained')

        for i in d.groups:
            data = d.get_group(i)
            ax.plot(data.columns[0], "fnonnat",data=data, label=i)

        ax.set_xlabel(df.columns[0])
        ax.set_ylabel('Fraction of Non-Native Contacts')
        ax.set_ylim(bottom=0)

        ax.legend(title="Mapping")
        ax.set_title('Fraction of Non-Native Contacts')
        fig.savefig(os.path.join(self.target_path, f'fnonnat_analysis.png'), dpi=300)

    def plot_fnat(self, path=None):
        """A function to perform lineplot visualization of fraction of native contacts through the simulation. Reads 'dockq_results.csv' file.
        
        Return: None
        """

        self._plot_fnat = True

        if path:
            self.handler.test_inp_path(path)
            df = pd.read_csv(path)
        else:
            df = pd.read_csv(os.path.join(self.table_path, "dockq_results.csv"))

        self._fnat_path = path
        d = df.groupby("mapping")

        fig,ax  = plt.subplots(figsize=(5,2.7), layout='constrained')

        for i in d.groups:
            data = d.get_group(i)
            ax.plot(data.columns[0], "fnat",data=data, label=i)

        ax.set_xlabel(df.columns[0])
        ax.set_ylabel('Fraction of Native Contacts')
        ax.set_ylim(bottom=0)

        ax.legend(title="Mapping")
        ax.set_title('Fraction of Native Contacts')
        fig.savefig(os.path.join(self.target_path, f'fnat_analysis.png'), dpi=300)

    def plot_biophys(self, path=None):
        """A function to perform barplot core-rim and biophysical classification visualization. Reads 'inetrface_label.csv' file.
        
        Return: None
        """
        if not self._biophys_palette:
            self._biophys_palette = sns.color_palette("Dark2", 4)
        

        if path:
            self.handler.test_inp_path(path)
            intf_df = pd.read_csv(path)
        else:    
            intf_df = pd.read_csv(os.path.join(self.table_path, "interface_label_perc.csv"))

        self._biophys = True
        self._biophys_path = path

        def _calc_biophys_type(res_type):
            res_dict = {"GLY": 'Hydrophobic', "ALA": 'Hydrophobic', "PRO": 'Hydrophobic', "VAL": 'Hydrophobic',
                        "LEU": 'Hydrophobic', "ILE": 'Hydrophobic', "MET": 'Hydrophobic', "TRP": 'Hydrophobic',
                        "PHE": 'Hydrophobic',
                        "SER": 'polar',
                        "THR": 'polar', "TYR": 'polar', "ASN": 'polar', "GLN": 'polar', "CYS": 'polar',
                        "LYS": '+ly charged', "ARG": '+ly charged', "HIS": '+ly charged', "ASP": '-ly charged',
                        "GLU": '-ly charged'}
            return res_dict[res_type]

        biophy = list()
        for index, row in intf_df.iterrows():
            biophy.append(_calc_biophys_type(row["Residue"][:3]))

        intf_df["Biophysical Type"] = biophy

        plot_df = intf_df.groupby(["Biophysical Type", "Interface Label"]).count().drop(["Residue", "Percentage"],
                                                                                        axis=1).reset_index()
        ty = list()
        for ind, r in plot_df.iterrows():
            if r["Interface Label"] == 2:
                ty.append("Support")
            elif r["Interface Label"] == 3:
                ty.append("Rim")
            elif r["Interface Label"] == 4:
                ty.append("Core")
            elif r["Interface Label"] == 1:
                ty.append("Surface")
            elif r["Interface Label"] == 0:
                ty.append("Interior")

        plot_df["Interface Text"] = ty

        pl = list()
        total = sum(plot_df["Chain"])
        for index, row in plot_df.iterrows():
            pl.append(row["Chain"] / total * 100)

        plot_df["Count"] = pl
        plot_df.drop("Chain", axis=1, inplace=True)

        plot_df.sort_values("Interface Label", inplace=True)

        fig, ax = plt.subplots()
        sns.barplot(x="Interface Text", data=plot_df, hue="Biophysical Type", y="Count", ax=ax, palette=sns.color_palette(self._biophys_palette, 4))
        l = ["Support", "Rim", "Core"]
        for ind, i in enumerate(ax.get_xticklabels()):
            if ax.get_xticklabels()[ind].get_text() in l:
                plt.setp(ax.get_xticklabels()[ind], color='red', style="italic", weight='bold')

        plt.ylabel("Percentage of Total Residue Types")
        plt.xlabel("Interface Label")

        plt.title("Biophysical Classification Counts of Residues")
        fig.savefig(os.path.join(self.target_path, f'Biophys_count.png'), dpi=300)

    def plot_DSSP(self, thereshold=50.0, intf_path=None, path=None):
        self._plot_DSSP = True 
        self._dssp_thereshold=thereshold
        

        if intf_path:
            self.handler.test_inp_path(intf_path)
            df = pd.read_csv(intf_path)

        if path:
            self.handler.test_inp_path(path)
            df = pd.read_csv(path)
        else:
            df = pd.read_csv(os.path.join(self.table_path, "residue_based_tbl.csv"))
            intf_df = pd.read_csv(os.path.join(self.table_path, "interface_label_perc.csv"))

        self._dssp_path = path  
        self._dssp_intf_path = intf_path

        time_name = df.columns[0]

        int_df = intf_df[(intf_df["Interface Label"] == 4) | (intf_df["Interface Label"] == 2) | (intf_df["Interface Label"] == 3) & (intf_df["Percentage"] >= thereshold)]

        df = df.loc[:, [time_name,"Chain",'Residue', 'Residue Number', 'Secondary Structure']]

        df["Residue"] = [a + str(b) for a, b in zip(df["Residue"], df["Residue Number"])]
        mylabels = {
            'H': 'Alpha Helix',
            'B': 'Beta Bridge',
            'E': 'Strand',
            'G': 'Helix-3',
            'I': 'Helix-5',
            'T': 'Turn',
            'S': 'Bend',
            'Null': 'Null',
            'P': 'Coil',
            'C': 'Coil',
            '!': 'Chain Break'
        }

        # Load data
        interface_groupped = int_df.groupby('Chain')

        groups = df.groupby(['Chain'])
        fig, axes = plt.subplots(len(groups), 1, figsize=(14, 10))
        
        for ax, chain, interface in zip(axes, groups, interface_groupped):
        
            chain_df = chain[1]

            intf_ch = interface[1]

            for index,row in chain_df.iterrows():
                if row['Residue'] not in intf_ch['Residue'].tolist():
                    chain_df.drop(index, inplace=True)

            sns.scatterplot(data=chain_df, x=chain_df.columns[0], y="Residue", hue="Secondary Structure", palette="deep", marker="s", ax=ax,s=50)

            ax.set_xlabel(f'{time_name}')
            ax.set_ylabel('Residue')
            ax.set_title(f"Chain {chain[0][0]}", fontsize=14)
            ax.tick_params(labelsize=9)
            ax.get_legend().remove()
            ax.xaxis.label.set_size(12)
            ax.yaxis.label.set_size(12)
            handles, labels = ax.get_legend_handles_labels()
    
        # Add labels and legend
        plt.suptitle("DSSP Analysis for Interface Residues", fontsize=18)
        fig.legend(handles, [mylabels[i] for i in labels ],ncol=3, loc='upper left') ## ortak legend
        plt.tight_layout()
        fig.subplots_adjust(hspace=0.20)
        fig.savefig(os.path.join(self.target_path, f'dssp_analysis.png'), dpi=300)

    def plot_SASA(self, path=None):
        """A function to perform line plot Interface Area in A^2 visualization for each chain and the overall complex. 'Reads residue_based_tbl.csv' file.
        
        Return: None
        """

        self.plot_SASA = True

        if path:
            self.handler.test_inp_path(path)
            df = pd.read_csv(path)
        else:
            df = pd.read_csv(os.path.join(self.table_path, "residue_based_tbl.csv"))

        self._sasa_path = path
        df = df.iloc[:, [0,1,2,3,7,8]]

        time_name = df.columns[0]

        plot_df = pd.DataFrame(columns=[time_name, "Chain", "SASA"])


        int_df = df[(df["Interface Label"] == 4) | (df["Interface Label"] == 2) | (df["Interface Label"] == 3)]
            
        g2 = int_df.groupby(["Chain"])

        fig, ax = plt.subplots(figsize=(5, 2.7), layout='constrained')
        for chain in g2.groups:
            chain_df = g2.get_group(chain)
            frames = chain_df.groupby([time_name])
            
            for frame in frames.groups:
                frame_df = frames.get_group(frame)
                frame_sasa = frame_df['SASA'].sum()
                plot_df.loc[len(plot_df)] = [frame, chain, frame_sasa]

        for i, chain_ in enumerate(plot_df.groupby('Chain').groups):
            ax.plot(plot_df.groupby('Chain').get_group(chain_)[time_name], plot_df.groupby('Chain').get_group(chain_)['SASA'], label=chain_, color=self.chain_colors[i])

        sums = [plot_df.groupby(time_name).get_group(group)['SASA'].sum() for group in plot_df.groupby(time_name).groups]

        ax.plot(plot_df[time_name].unique() ,sums, color=self.complex_color, label="Complex")
         
        ax.set_xlabel(time_name)
        ax.set_ylabel(f'Interface Area (Å^2)')

        ax.set_title(f"Interaction Surface Area throughout Simulation")
        ax.legend(title="Chains")
        fig.savefig(os.path.join(self.target_path, f'interface_areas.png'), dpi=300)

    def plot_pairwise_freq(self, path=None):
        """A function to perform barplot visualization of interaction frequency in simulation of residue pairs with locations (side chain or backbone) of atoms that participated in interaction. Reads 'int_based_table.csv' file.
        
        Return: None
        """
        if path:
            self.handler.test_inp_path(path)
            df = pd.read_csv(path)
            time_name = df.columns[0]
            time_var = df.iloc[:,0]
            sim_time = len(np.unique(time_var.tolist()))
            df = df[[df.columns[0],"itype", "pairwise", "atom_a", "atom_b"]].groupby('itype')

        else:
            df = pd.read_csv(os.path.join(self.table_path, "int_based_table.csv"))
            time_name = df.columns[0]
            time_var = df.iloc[:,0]
            sim_time = len(np.unique(time_var.tolist()))
            df = df[[df.columns[0],"itype", "pairwise", "atom_a", "atom_b"]].groupby('itype')
            
        self._bar = True
        self._bar_path = path

        if not self._bar_palette:
            self._bar_palette = [['#fee8c9','#fc8e5a','#b30000'], ['#c3dceb', '#68a9cf','#133345' ], ['#efedf5', '#807dba', '#3a1a63']]

        
        backbones = ["HN", "N", "CA", "HA", "C", "O"]

        #sim_time = len(np.unique(df.iloc[:,0].tolist()))

        #df = df.groupby('itype')
        #first, find percentage of pairwises
        for x,palet in zip(df.groups, self._bar_palette):
            ll = dict()
            z = df.get_group(x).copy()
            z_gr = z.groupby('pairwise')
            for pair in z_gr.groups:
                df2 = z_gr.get_group(pair)
                l = len(np.unique(df2.iloc[:,0].tolist()))
                percentage = l/sim_time*100
                ll[pair] = percentage

            locs = list()
            for i in z.iterrows():
                lo = list()
                if i[1]["atom_a"] in backbones:
                    lo.append("bb")
                else:
                    lo.append("sc")

                if i[1]["atom_b"] in backbones:
                    lo.append("bb")
                else:
                    lo.append("sc")

                locs.append("-".join(lo))

            z["locs"] = locs

            #find counting locations
            plot_dict = dict()
            for index, row in z.iterrows():
                if row["pairwise"] not in plot_dict.keys():
                    plot_dict[row["pairwise"]] = [0, 0, 0]
                if row["locs"] == "bb-bb":
                    plot_dict[row["pairwise"]][0] += 1
                if row["locs"] == "bb-sc" or row["locs"] == "sc-bb":
                    plot_dict[row["pairwise"]][1] += 1
                if row["locs"] == "sc-sc":
                    plot_dict[row["pairwise"]][2] += 1

            pairs = list(plot_dict.keys())

            #portion of locations for each pairwise
            for pair, val in plot_dict.items():
                total = sum(val)
                for i in [0,1,2]:
                    val[i] = val[i]*ll[pair]/total

            plott = dict(sorted(plot_dict.items(), key=lambda a: sum(a[1]), reverse=True))

            #plot
            fig, ax = plt.subplots(figsize=(20, 12))

            c = list(plott.keys())
            v = np.array(list(plott.values()))
            
            p = ax.bar(pairs, v[:, 0], color=palet[0])
            bottom = v[:, 0]

            for i in [1,2]:
                    p = ax.bar(pairs, v[:, i], color=palet[i], bottom=bottom)
                    bottom += v[:, i]

            ax.legend(labels=("bb-bb", "bb-sc or sc-bb", "sc-sc"), fontsize=20)


            plt.xlabel("Pairwise", fontweight="bold", size=20)
            plt.ylabel(f"Percentage to Total Simulation Time", fontweight="bold", size=20)

            plt.xticks(rotation='vertical')

            plt.xticks(fontsize=12, fontweight='bold', rotation=85)
            plt.yticks([x for x in range(0,101,10)], fontsize=13, fontweight='bold')

            plt.title(f"General {x.strip('bond').capitalize()}-bond Percentage to Simulation Time", fontweight='bold', size=20)

            plt.savefig(os.path.join(self.target_path, f"Pairwise_{x}-frequency.png"), dpi=300, bbox_inches='tight', format="png")                  

    def plot_int_energy(self, thereshold=50.0, intf_path=None, res_path=None):
        """A function to perform boxplot visualization of interaction energy variation of residues that has constant interface label at least 50%-by default- of all simulation. Reads 'interface_label_perc' and 'residue_based_tbl.csv' files.
        
        Keyword arguments:
        thereshold -- int: The desired percentage value of the residual interface label to be constant throughout the simulation
        Return: None
        """

        if intf_path:
            self.handler.test_inp_path(intf_path)
            df = pd.read_csv(intf_path)

        if res_path:
            self.handler.test_inp_path(res_path)
            df = pd.read_csv(res_path)

        else:
            df = pd.read_csv(os.path.join(self.table_path, "interface_label_perc.csv"))
            energy_df = pd.read_csv(os.path.join(self.table_path, "residue_based_tbl.csv"), usecols=["Chain", "Residue", "Residue Number", "Total Residue Energy"])

        self._ene = True
        self._ene_thr = thereshold
        self._ene_intf = intf_path
        self._ene_path = res_path

        int_df = df[(df["Interface Label"] == 4) | (df["Interface Label"] == 2) | (df["Interface Label"] == 3) & (
                df["Percentage"] >= thereshold)]

        energy_df["Residue"] = [a + str(b) for a, b in zip(energy_df["Residue"], energy_df["Residue Number"])]
        g = energy_df.groupby(["Chain", "Residue"])
        n_df = pd.DataFrame(columns=["Chain", "Residue", "Total Residue Energy"])
        my_palette = dict()
        ch_num = 1
        for index, row in int_df.iterrows():
            mg = (row["Chain"], row["Residue"])
            if mg[0] not in my_palette.keys():
                my_palette[mg[0]] = self.chain_colors[ch_num - 1]
                ch_num += 1
            data = g.get_group(mg)
            n_df = pd.concat([n_df, data])

        fig, ax = plt.subplots(figsize=(15, 8))
        sns.boxplot(data=n_df.sort_values(["Chain", "Residue Number"]), x="Residue", y="Total Residue Energy", hue="Chain",
                    palette=my_palette)
        plt.xticks(rotation=90)
        plt.title("Interface Residue Based Energy Evulation for Complete Simulation")
        plt.savefig(os.path.join(self.target_path, "interface_energy_eval.png"), dpi=300)

    def _get_params_(self):
        params = {
            'job_name':self._job_name,
            'output_path': self.target_path,
            'table_path': self.table_path,
            'chain_color_palette': self.chain_colors,
            'complex_color': self.complex_color,
            'PlotRMSD': {
                'Run': self._rmsd,
                'rmsd_table_path':self._rmsd_path
                },
            'PlotRG': {
                'Run': self._rg,
                'rg_table_path':self._rg_path
                },
            'PlotRMSF': {
                'Run': self._rmsf,
                'rmsf_table_path':self._rmsf_path,
                'rmsf_intf_res_table': self._rmsf_intf_path
                },
            'PlotBiophys': {
                'Run': self._biophys,
                'biophys_table_path':self._biophys_path,
                'biophys_palette': self._biophys_palette
                },
            'PlotSASA' : {
                'Run': self._plot_SASA,
                'sasa_path': self._sasa_path,
            },
            'PlotiRMSD': {
                'Run': self._plot_irmsd,
                'irmsd_path': self._irmsd_path,

            },
            'PlotlRMSD': {
                'Run': self._plot_lrmsd,
                'irmsd_path': self._lrmsd_path,

            },
            'PlotDockQ': {
                'Run': self._plot_dockq,
                'irmsd_path': self._dockq_path,

            },
            'PlotFnonnat': {
                'Run': self._plot_fnonnat,
                'fnonnat_path': self._fnonnat_path
            },
            'PlotFnat': {
                'Run': self._plot_fnat,
                'fnat_path': self._fnat_path
            },
            'PlotPairwiseFreq': {
                'Run': self._bar,
                'bar_table_path':self._bar_path,
                'bar_palette': self._bar_palette
                },
            'PlotResEne': {
                'Run': self._ene,
                'interface_th':self._ene_thr,
                'interface_table_path': self._ene_intf,
                'residue_based_table': self._ene_path
                },
            'PlotDSSP': {
                'Run': self._plot_dssp,
                'dssp_path': self._dssp_path,
                'dssp_thereshold':self._dssp_thereshold,
                'dssp_intf_path':self._dssp_intf_path
            },

            }
        json_path = os.path.join(self.job_path, 'plot_params.json')
        with open(json_path, 'w+') as ofh:
            json.dump(params, ofh)