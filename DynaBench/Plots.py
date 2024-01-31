import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from DynaBench.handling import tables_errors
import json



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

        sns.set_style('whitegrid')

        if job_name is None:
            #job_name = input('Plase provide job name:\n')
            for file in os.listdir(os.getcwd()):
                if '_db' in file:
                    self.job_path = os.path.join(os.getcwd(), file)
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

        #self.handler.test_intf_path(os.path.join(self.table_path, "interface_label_perc.csv"))

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

            int_df = df[(df["Interface Label"] == 4) | (df["Interface Label"] == 2) | (df["Interface Label"] == 3) & (
                    df["Percentage"] >= 50)]
            
            g2 = int_df.groupby(["Chain"])

        groups = df1.groupby(["Molecule"])
        
        try:
            fig, axes = plt.subplots(1, len(groups), sharex="all", figsize=(10, 6), sharey=True)
            for ind, (ax, x) in enumerate(zip(axes, np.unique(df1["Molecule"]))):
                data = groups.get_group(x)
                if intf:
                    int_data = g2.get_group(x)
                    markers = [int(x[3:]) for x in int_data["Residue"]]
                    ax.plot(data["Residue Number"], data["RMSF"], '-o', markevery=markers, markersize=3.5,
                            color=self.chain_colors[ind])
                    ax.plot(data["Residue Number"], data["RMSF"], 'o', markevery=markers,
                            label='Interface Residues', markersize=3.5, color="red")
                else:
                    ax.plot(data["Residue Number"], data["RMSF"],
                            color=self.chain_colors[ind])

                ax.set_xlabel("Residue Number")
                ax.set_ylabel(f'RMSF (Å)')
                ax.set_title(f"Chain {x}", fontsize=16)
                ax.xaxis.label.set_size(12)
                ax.yaxis.label.set_size(12)
            plt.suptitle("RMSF Analysis", fontsize=18)
            if intf:
                plt.legend()
            fig.savefig(os.path.join(self.target_path, f'RMSF-proteinCA.png'), dpi=300)

        except:
            fig, axes = plt.subplots(1, len(groups), sharex="all", figsize=(10, 6), sharey=True)
            for ind, (ax, x) in enumerate(zip(axes, np.unique(df1["Molecule"]))):
                data = groups.get_group(x)
                ax.plot(data["Residue Number"], data["RMSF"],
                        color=self.chain_colors[ind])

                ax.set_xlabel("Residue Number")
                ax.set_ylabel(f'RMSF (Å)')
                ax.set_title(f"Chain {x}", fontsize=16)
                ax.xaxis.label.set_size(12)
                ax.yaxis.label.set_size(12)
            plt.suptitle("RMSF Analysis", fontsize=18)
            if intf:
                plt.legend()
            fig.savefig(os.path.join(self.target_path, f'RMSF-proteinCA.png'), dpi=300)


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
                ty.append("Interior")
            elif r["Interface Label"] == 0:
                ty.append("Surface")

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
        # ax.get_xticklabels[3:].set_color("red")
        l = ["Support", "Rim", "Core"]
        for ind, i in enumerate(ax.get_xticklabels()):
            if ax.get_xticklabels()[ind].get_text() in l:
                plt.setp(ax.get_xticklabels()[ind], color='red', style="italic", weight='bold')

        plt.ylabel("Percentage of Total Residue Types")
        plt.xlabel("Interface Label")

        plt.title("Biophysical Classification Counts of Residues")
        fig.savefig(os.path.join(self.target_path, f'Biophys_count.png'), dpi=300)

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

            plt.xticks(fontsize=10, fontweight='bold', rotation=85)
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

            }
        json_path = os.path.join(self.job_path, 'plot_params.json')
        with open(json_path, 'w+') as ofh:
            json.dump(params, ofh)