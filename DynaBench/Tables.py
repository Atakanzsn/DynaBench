import os

import MDAnalysis
import freesasa
import MDAnalysis as mda
import pandas as pd
from MDAnalysis.analysis.rms import RMSD, RMSF
from MDAnalysis.analysis import align
import interfacea as ia
from openmm import app
import numpy as np
import DynaBench.handling as hp
import DynaBench.pdb_tool_modified as ptm
import random
import json
import subprocess
import shutil

freesasa.setVerbosity(freesasa.silent)
handler = hp.tables_errors()

class dynabench:
    def __init__(self, inp_file, stride=1, split_models=False, chains=None, job_name=None, dcd_pdb=None, show_time_as="Frame", timestep=None, time_unit=None):
        """A class to perform Quality Control, Residue Based and Interaction Based analyses on MD simulations. Results of the analyses are printed as .csv files under a folder named 'tables'. MD outputs with any exception are transformed to .pdb file and the analyses are run through this .pdb file. Number of frames that will be fetched from initial input file to be analysed can be set by stride value.
        
        Keyword arguments:
        inp_file -- Path to the input file
        stride -- Number of frames to be skipped. By default, 1 which means no frames to be skipped.
        split_models -- Boolean. By default, False. If True, all the models will be splitted into a folder named 'models'.
        Return: None
        """
        handler.test_inp_path(inp_file)
        handler.test_stride(stride)
        handler.test_split_models(split_models)

        if timestep is None:
            self.timestep = 1.0
        else:
            self.timestep = timestep 

        #params
        self.time_Type = show_time_as
        self.time_unit = time_unit
        self.inp_file_ = inp_file
        self.split_models_ = split_models
        self.chains_ = chains
        self.qc_flag = False
        self.rb_flag = False
        self.ib_flag = False
        self.dcd_pdb = dcd_pdb
        self.stride = stride
        self.rmsd_data= None
        self.get_all_hph = None
        self.foldx_path = None

        if job_name is None:
            file_name = inp_file.split('.')[0]
            job_name = file_name + '_db' + str(random.randrange(100,999))

        
        self.job_name = job_name #

        self.job_path = os.path.join(os.getcwd(), self.job_name) #

        if not os.path.exists(self.job_path): #create job folder
            os.mkdir(self.job_path)

        self.target_path = os.path.join(self.job_path, 'tables') #define tables folder path

        if not  os.path.exists(self.target_path): #create tables folder
            os.mkdir(self.target_path)

        #check for preprocess ,define pdb_file
        file_name = inp_file.split(".")[0]
        file_ext = inp_file.split(".")[-1]

        if file_ext == 'dcd':
            name = file_name.split("\\")[-1].split(".")[0]
            out_file = os.path.join(self.job_path, f"{name}.pdb")
            if not self.dcd_pdb:
                self.dcd_pdb = input('Please provide input pdb file:\n')

            self.pdb_file = os.path.join(self.job_path, self._preprocess_dcd(self.dcd_pdb, out_file, stride, inp_file))

        elif file_ext == 'pdb':
            if self.stride != 1:
                name = file_name.split("\\")[-1]
                out_file = os.path.join(self.job_path, f"{name}_stride{stride}.pdb")
                u = mda.Universe(inp_file)

                with mda.Writer(out_file, u.atoms.n_atoms) as W:
                    for ts in u.trajectory[::stride]:
                        W.write(u.atoms)

                self.pdb_file = out_file
            else:
                self.pdb_file = inp_file

        #check for split models
        if split_models:
            self._split_models(self.pdb_file, self.job_path)

        #check for split chains
        if chains:
            handler.test_chain_sel(chains, self.pdb_file)
            self.pdb_file = self._sel_chain(self.pdb_file, chains, self.job_path)


    def _sel_chain(self, pdb, chains, job_path):
        u = mda.Universe(pdb)
        new_name = os.path.join(job_path, f"{pdb.split('.')[0]}_{','.join(chains)}.pdb")

        sel = []
        for el in chains:
            sel.append(f"chainID {el.upper()}")

        t = u.select_atoms(' or '.join(sel))
        t.write(new_name, frames='all')

        return new_name

    def _get_params_(self):
        params = {
            'job_name':self.job_name,
            'input_file': self.inp_file_,
            'dcd_pdb': self.dcd_pdb,
            'stride': self.stride,
            'split_models':self.split_models_,
            'chains': self.chains_,
            'pdb_file': self.pdb_file,
            'show_time_as': self.time_Type,
            'timestep': self.timestep,
            'timeunit': self.time_unit,

            'QualityControl': {
                'Run': self.qc_flag,
                'rmsd_data':self.rmsd_data
                },
            'ResidueBased':{
                'Run': self.rb_flag,
                'FoldX_path': self.foldx_path
            },
            'InteractionBased': { 
                'Run': self.ib_flag,
                'get_all_hph': self.get_all_hph
            }

        }

        json_path = os.path.join(self.job_path, 'table_params.json')
        with open(json_path, 'w+') as ofh:
            json.dump(params, ofh)


    @staticmethod
    def _split_models(file, job_path):

        """ Creates folder 'models' and splits the trajectory into the models folder.
        
        Keyword arguments:
        file -- trajectory file
        Return: None
        """
        models_path = os.path.join(job_path, 'models')
        if not os.path.exists(models_path):
            os.mkdir(models_path)
        os.chdir(models_path)
        os.system(f"pdb_splitmodel ..\\{file}")
        os.chdir(job_path)

    @staticmethod
    def _preprocess_dcd(inp_pdb, output_file, stride, ab_file):
        """ Transforms input trajectory file into the .pdb file for given stride value.
        
        Keyword arguments:
        inp_file -- Input trajectory file
        output_file -- Output .pdb file
        stride -- Number of frames to be skipped.
        Return: None
        """
        u = mda.Universe(inp_pdb, ab_file)
        name = output_file.split("\\")[-1].split(".")[0]

        with mda.Writer(output_file, u.atoms.n_atoms) as W:
            for ts in u.trajectory[::stride]:
                W.write(u.atoms)

        ptm.chains(output_file, ['V', 'S'])
        #os.system(f'python pdb_tool_modified.py -V,S {output_file} > {name}_chained.pdb')

        return name + "_chained.pdb"
        

    def run_quality_control(self, rmsd_data={'ref_struc':None, 'ref_frame':0}):
        """ The function to run Quality Control analyses in one line. Defines Quality Control class and runs the functions in it.
        
        Return: None
        """
        handler.test_rmsd_dict(rmsd_data)

        self.rmsd_data = rmsd_data
        self.qc_flag = True

        c = self.QualityControl(self.pdb_file, self.target_path, rmsd_data=rmsd_data,
                                timestep= self.timestep, time_unit = self.time_unit, time_type=self.time_Type, stride=self.stride)
        c.quality_tbls_overtime()
        c.quality_tbls_overres()

    def run_res_based(self, foldx_path):
        """ The function to run Residue Based analyses in one line. Defines Residue Based class and runs the functions in it.
        
        Return: None
        """
        self.foldx_path = foldx_path
        self.rb_flag = True
        c = self.ResidueBased(self.pdb_file, self.target_path, timestep= self.timestep, timeunit = self.time_unit, time_type=self.time_Type, stride=self.stride, job_path=self.job_path, foldx_path=foldx_path)
        c.res_based_tbl()
        c.interface_table()

    def run_inter_based(self, get_all_hph=False):
        """The function to run Interaction Based analyses in one line. Defines Interaction Based class and runs the functions in it.
        
        Keyword arguments:
        get_all_hph -- True or False. Default is False. If True, the function calculates all the possible hydrophobic interactions.
        Return: None
        """
        self.ib_flag = True
        self.get_all_hph=get_all_hph
        c = self.InteractionBased(self.pdb_file, self.target_path, get_all_hph=get_all_hph, timestep=self.timestep, timeunit=self.time_unit, time_type=self.time_Type, stride=self.stride)
        c.int_based_tbl()

    
    class QualityControl:
        def __init__(self, pdb_file, target_path, rmsd_data, timestep, time_unit, time_type, stride):
            """ A class to perform Quality Control analyses (RMSD, RG, and RMSF) of the trajectory. Initially runs private methods and defines RMSD, RG, and RMSF results of the trajectory.

            Return: None
            
            """
            self.stride = stride
            self.timetype = time_type
            self.timestep = timestep
            if time_type == 'Time':
                self.u = mda.Universe(pdb_file, dt=timestep)
                if time_unit.lower() == 'ns' or time_unit.lower() == 'nanosecond':
                    self.u.trajectory.units = {'time':'nanosecond', 'length': 'Angstrom'}
            else:
                self.u = mda.Universe(pdb_file, in_memory=True)

            handler.test_rmsd_refstruc(rmsd_data['ref_struc'])
            handler.test_rmsd_refframe(rmsd_data['ref_frame'], len(self.u.trajectory))
        

            self.target_path = target_path
            self.segs = self.u.segments
            self.protein = self.u.select_atoms('protein')

            self.rmsd_ref_struc = rmsd_data['ref_struc']
            self.rmsd_ref_frame = rmsd_data['ref_frame']
            self.rmsd_dict = self._calc_rmsd()
            self.rdg_dict = self._calc_rdg()
            self.rmsf_dict = self._calc_rmsf()

            if time_type == 'Time':
                self.header_overtime = ["Time (ns)"]
            else:
                self.header_overtime = ["Frame"]
            self.header_overres = ["Molecule", "Residue Number", "RMSF"]

        #
        # Private Methods
        #

        def _calc_rmsd(self):
            """ Calculates RMSD of each chain and the complex during the simulation.
            

            Return: Dict: A dictionary that contains the chains and complexes as keys and their corresponding RMSD values as values.
            """
            
            analysis = {}
            for chain in self.segs:
                chainid = chain.segid
                backbone = self.protein.select_atoms(f'chainid {chainid} and backbone')
                R1 = RMSD(backbone, reference=self.rmsd_ref_struc, ref_frame=self.rmsd_ref_frame)
                R1.run()
                analysis[f"Chain {chainid}"] = [f"{x:.03f}" for x in R1.results.rmsd.T[2]]

            backbone = self.protein.select_atoms('backbone')
            R1 = RMSD(backbone, reference=self.rmsd_ref_struc, ref_frame=self.rmsd_ref_frame)
            R1.run()
            analysis["Complex"] = [f"{x:.03f}" for x in R1.results.rmsd.T[2]]

            return analysis

        def _calc_rdg(self):
            """Calculates RG of each chain and the complex during the simulation.
            
            Return: Dict: A dictionary that contains the chains and complexes as keys and their corresponding RMSD values as values.
            """
            
            analysis = {}
            proteins = {}
            for chain in self.segs:
                chainid = chain.segid
                #if chainid != "SOLV" and chainid != "IONS":

                proteins[chainid] = self.u.select_atoms(f"chainid {chainid} and protein")
                analysis[f"Chain {chainid}"] = []

            proteins["Complex"] = self.protein
            analysis["Complex"] = []

            for _ in self.u.trajectory:
                for key, protein in proteins.items():
                    try:
                        analysis[f"Chain {key}"].append(f"{protein.radius_of_gyration():.03f}")
                    except KeyError:
                        analysis[f"Complex"].append(f"{protein.radius_of_gyration():.03f}")

            return analysis

        def _calc_rmsf(self):
            """ Calculates the RMSF of each residue for all simulation time.

            Return: Dict: A dictionary that contains the chain ids as keys and a list that has residue number in the index 0 and the corresponding RMSF values in the index 1 as values.
            """
            
            analysis = {}
            for chain in self.segs:
                chainid = chain.segid
                #if chainid != "SOLV" and chainid != "IONS":
                analysis[chainid] = []

                calphas_p = self.protein.select_atoms(f'chainid {chainid} and name CA')
                rmsf = RMSF(calphas_p)
                rmsf.run()

                analysis[chainid].append(calphas_p.resnums)
                analysis[chainid].append([f"{x:.03f}" for x in rmsf.results.rmsf])

            return analysis

        #
        # Public Methods
        #

        def quality_tbls_overtime(self):
            """ Reads RMSd and RG results and writes them to the csv file with respect to frame number.
            
            Return: None
            """
            

            file = open(os.path.join(self.target_path, "QualityControl-overtime.csv"), "w")

            for key in list(self.rmsd_dict.keys()):
                self.header_overtime.append(f"{key} RMSD")
                self.header_overtime.append(f"{key} RDG")

            print(",".join(self.header_overtime), end="\n", file=file)

            for i in range(0, len(list(self.rmsd_dict.values())[0])):  # stride'ı burada kullan
                row = [str(i * self.stride * float(self.timestep))]

                rmsd_values = [self.rmsd_dict[key][i] for key in self.rmsd_dict.keys()]
                rdg_values = [self.rdg_dict[key][i] for key in self.rdg_dict.keys()]

                for r1, r2 in zip(rmsd_values, rdg_values):
                    row.append(str(r1))
                    row.append(str(r2))

                print(",".join(row), end="\n", file=file)
            file.close()

        def quality_tbls_overres(self):
            """Reads RMSF results and writes them to the csv file with respect to chain id and residue number.

            Return: None
            """
            

            file = open(os.path.join(self.target_path, "QualityControl-overres.csv"), "w")

            print(",".join(self.header_overres), end="\n", file=file)

            for key, value in self.rmsf_dict.items():
                for index, val in enumerate(value[0]):
                    rms = value[1][index]

                    row = f"{key},{val},{rms}\n"

                    file.write(row)

            file.close()

    class ResidueBased:
        def __init__(self, pdb_path, target_path, time_type, timestep, timeunit, stride, job_path, foldx_path):
            """ A class to perform Residue Based analyses such as core-rim, and  biophysical type classifications; Van der Waals, electrostatic, desolvation, and hydrogen bond energies either between same-chain or different chain residues. Also prints the interface class that residues had with the highest percentage for all simulation time. Initially calculates rASA values and residue energies.
            
            Keyword arguments:
            pdb_path -- Input .pdb file
            target_path -- Path to 'tables' folder
            Return: None
            """
            self.time_type = time_type
            self.timestep = timestep
            self.timeunit = timeunit
            self.stride = stride
            self.foldx_path=foldx_path

            if self.time_type and 'Time' in self.time_type:
                if self.timeunit.lower() == 'ns' or self.timeunit.lower() == 'nanosecond':
                    t = 'Time (ns)'
            else:
                t = 'Frame'

            self.target_path = target_path
            self.header = [f"{t}", "Chain", "Residue", "Residue Number", "rASAc", "rASAm", "delta rASA",
                           "Interface Label", "Residue Biophysical Type", "Backbone Hbond Energy","Sidechain Hbond Energy", "Van der Waals Energy", "Electrostatic Energy",
                           "Total Residue Energy","\n"]

            self.rasam_array = freesasa.structureArray(pdb_path,
                                                       {'separate-chains': True,
                                                        'separate-models': True})  # rasc, len=402
            self.rasac_array = freesasa.structureArray(pdb_path,
                                                       {'separate-chains': False,
                                                        'separate-models': True})  # rasm, len=201

            self.energies = self._res_en(pdb_path, job_path, foldx_path)

        #
        # Private Methods
        #

        @staticmethod
        def _calc_interface_class(delta_ras, rasc, rasm):
            """ Residues are classified according to the their rASA (relative accessible solvent area) values for both monomer (rASAm) and complex (rASAc) conformations by Levy et al. For more information about core-rim classification, please visit https://doi.org/10.1016/j.jmb.2010.09.028.

            This function classifies the residue according to the given rASA values. The classification labels are:
                0=Surface, 1=Interior, 2=Support, 3=Rim, 4=Core

            Keyword arguments:
            delta_ras -- rASAc - rASAm
            rasc -- rASA value for complex conformation
            rasm -- rASA value for monomer conformation
            Return: int: Interface Class label.
            """
            label = 0  # label will be used to extract static/dynamic interfaces
            if delta_ras == 0:
                if rasc < 0.25:
                    label = 0
                else:
                    label = 1

            if delta_ras > 0:
                if rasm < 0.25:
                    label = 2

                elif rasc > 0.25:
                    label = 3

                elif rasm > 0.25 > rasc:
                    label = 4

            return label

        @staticmethod
        def _calc_biophys_type(res_type):
            """ Classifies residues as +/-ly charged, hydrophobic or polar. Output labels are:
                hydrophobic: 0
                +ly charged: 1
                -ly charged: 2
                polar: 3

            Keyword arguments:
            res_type -- 3 letter residue code.
            Return: int: Biophysical Type label.
            """
            res_dict = {"GLY": 0, "ALA": 0, "PRO": 0, "VAL": 0, "LEU": 0, "ILE": 0, "MET": 0, "TRP": 0, "PHE": 0,
                        "SER": 3,
                        "THR": 3, "TYR": 3, "ASN": 3, "GLN": 3, "CYS": 3, "LYS": 1, "ARG": 1, "HIS": 1, "ASP": 2,
                        "GLU": 2}

            return res_dict[res_type]

        @staticmethod
        def _res_en(inp_file, job_path, foldx_path):
            """ Calculates the residue energies (Van der Waals, electrostatic, desolvation, and hydrogen bond for both same chain and different chain interactions) by running EvoEF1 for each frame.
            
            Keyword arguments:
            inp_file -- Input .pdb file
            Return: dict: Dictionary-in-dictionary with the order frame_number-residue-energy_values.
            """
            

            class ResEnergy:
                def __init__(self, name, resnum, chain, total, bb_hbond, sc_hbond, vdw, elec):
                    """A class to collect energy values of a residue. This class is defined for each residue. The energy values that ends with the letter 's' describes same-chain interaction, while letter 'd' stands for different chain interations.
                    
                    Keyword arguments:
                    vdwatt -- Van der Waals attraction energy
                    vdwrep -- Van der Waals repulsion energy
                    elec -- Electrostatic energy
                    HB -- description
                    desH -- Hydrophobic desolvation energy
                    desP -- Polar desolvation energy

                    Return: None
                    """
                    
                    self.res_name = name
                    self.resnum= resnum
                    self.chain = chain
                    self.total = total
                    self.bb_hbond = bb_hbond
                    self.sc_hbond = sc_hbond
                    self.vdw = vdw
                    self.elec = elec

                    self.hbond = self.bb_hbond + self.sc_hbond

            result = dict()

            #run foldx

            inp_path = os.path.abspath(inp_file)

            current = os.getcwd()

            for file in os.listdir(foldx_path):
                if 'foldx' in file:
                    foldx_exe_path = os.path.abspath(os.path.join(foldx_path, file))

            if not os.path.exists(os.path.join(job_path, 'models')):
                os.mkdir(os.path.join(job_path, 'models'))
            
            os.mkdir(os.path.join(job_path, 'foldx_outputs'))


            models_path = os.path.join(job_path, 'models')

            output_path = os.path.join(job_path, 'foldx_outputs')
            
            os.chdir(models_path)
            
            os.system(f"pdb_splitmodel {inp_path}")

            os.chdir(current)

            with open(os.path.join(job_path, 'pdb_list.out'), 'w+') as fh:
                for file in os.listdir(models_path):
                    fh.write(f"{file}\n")

            subprocess.run(f"{foldx_exe_path} --command=SequenceDetail --output-dir={output_path} --pdb-dir={models_path} --pdb-list={os.path.join(job_path,'pdb_list.out')}", stdout=subprocess.DEVNULL)

            os.remove(os.path.join(job_path, 'pdb_list.out'))

            #read outputs

            for file in os.listdir(output_path):
                frame = int(file.split('.')[0].split('_')[-1]) - 1
                if frame not in result.keys():
                    result[frame] = dict()
                with open(os.path.join(output_path, file), 'r+') as fh:
                    for r in fh:
                        r = r.rstrip("\n")
                        splitted = r.split("\t")

                        try:
                            resname = splitted[1]
                            chain = splitted[2]
                            res_num = splitted[3]

                            total= float(splitted[8])
                            bb_hbond = float(splitted[9])
                            sc_hbond=float(splitted[10])
                            vdw = float(splitted[11])
                            elec = float(splitted[12])

                            res= ResEnergy(resname, res_num, chain, total, bb_hbond, sc_hbond, vdw, elec)

                            result[frame][f"{res.chain}{res.resnum}{res.res_name}"] = res
                        
                        except IndexError:
                            continue

            return result

        #
        # Public Methods
        #

        def res_based_tbl(self):
            """Merges the residue based analyses results and writes to .csv file.
            
            
            Return: None
            """
            
            file = open(os.path.join(self.target_path, "residue_based_tbl.csv"), "w")  # open the output file

            file.write(",".join(self.header))

            i = 0
            for frame, frame_obj in enumerate(self.rasac_array):  # stride'ı frame'i düzenlemek için kullan
                rasc_res_frame = freesasa.calc(frame_obj).residueAreas()  # first, get the first frame from rasc

                chain_len = len(rasc_res_frame.keys())  # get the total chain num

                monomer_list = []
                for x in range(chain_len):
                    rasm_res_frame = freesasa.calc(
                        self.rasam_array[i]).residueAreas()  # get the chains of this frame from rasm.
                    monomer_list.append(rasm_res_frame)
                    i += 1

                for key, value in rasc_res_frame.items():  # iterate through chains by using rasc dict., key=chain
                    rasm_chain_dict = [x for x in monomer_list if key in x.keys()][0]  # get the rasm dict of this chain

                    for resnum, val in value.items():
                        rasc = val.relativeTotal
                        rasm = rasm_chain_dict[key][resnum].relativeTotal
                        delta_ras = rasm - rasc

                        label = self._calc_interface_class(delta_ras, rasc, rasm)
                        try:
                            biophy_class = self._calc_biophys_type(val.residueType)
                        except:
                            biophy_class = None

                        name = f"{key.upper()}{val.residueNumber}{val.residueType}"
                        try:
                            obj = self.energies[frame][name]
                        except:
                            continue
                        
                        if self.time_type == 'Time':
                            t = frame * self.stride * float(self.timestep)
                        else:
                            t = frame
                        row = f"{t},{key},{val.residueType},{resnum},{f'{rasc:.03f}'},{f'{rasm:.03f}'}, {f'{delta_ras:.03f}'},{label},{biophy_class},{f'{obj.bb_hbond:.03f}'},{f'{obj.sc_hbond:.03f}'},{f'{obj.vdw:.03f}'},{f'{obj.elec:.03f}'},{f'{obj.total:.03f}'}\n"

                        file.write(row)
            file.close()



        def interface_table(self):
            """Reads residue based csv file and writes the highest precentage of interface label of whole simulation for each residue. 
            
            Return: None
            """
            
            df = pd.read_csv(os.path.join(self.target_path, "residue_based_tbl.csv"),
                             usecols=["Chain", "Residue", "Residue Number", "Interface Label"])
            df["Residue Name"] = [a + str(b) for a, b in zip(df["Residue"], df["Residue Number"])]
            df.drop(["Residue", "Residue Number"], axis=1, inplace=True)
            df2 = pd.DataFrame(columns=["Chain", "Residue", "Interface Label", "Percentage"])
            a = df.groupby(["Residue Name", "Chain"])
            for g in a.groups:
                data = a.get_group(g)
                length_all = len(data["Interface Label"])
                a2 = data.groupby("Interface Label")
                lbl_list = list()
                for g2 in a2.groups:
                    data2 = a2.get_group(g2)
                    len_label = len(data2)
                    perc = len_label / length_all * 100
                    lbl_list.append((g2, perc))
                label = sorted(lbl_list, key=lambda x: x[1], reverse=True)
                df2.loc[len(df2.index)] = [np.unique(data["Chain"])[0], g[0], label[0][0], format(label[0][1], ".2f")]
                df2.sort_values(by="Chain", inplace=True)

            df2.to_csv(os.path.join(self.target_path, "interface_label_perc.csv"), index=False, mode="w+")

    class InteractionBased:
        def __init__(self, pdb_path, target_path, timestep, timeunit, time_type, stride, get_all_hph=False):
            """ A class to calculate Hydrogen, Electrostatic, and Hydrophobic interactions with atom informations that contribute to the interaction, and writes interaction information for all frames to a .csv file.
            
            Keyword arguments:
            pdb_path -- Input .pdb file.
            target_path -- Path to the 'tables' fodler.
            get_all_hph -- Boolean. If True, prints all possible hydrophobic interactions. By default, prints only the closest interaction. 
            Return: None
            """
            self.timeunit = timeunit
            self.time_type = time_type
            self.stride = stride
            self.timestep = timestep

            self.hph_status = get_all_hph
            self.pdb_path = pdb_path
            self._struct = app.PDBFile(pdb_path)
            self.target_path = target_path

        def _calc_interactions(self):
            """Calculates the interactions by Interfacea package.
            
            Return: list: A list of lists such as [interaction_table, frame_number].
            """
            
            return_list = []
            for i in range(0, len(self._struct._positions)):
                self._struct.positions = self._struct._positions[i]
                self._struct.topology.createDisulfideBonds(self._struct.positions)
                mol = ia.Structure("pp.pdb", self._struct)

                analyzer = ia.InteractionAnalyzer(mol)
                analyzer.get_hbonds(strict=True)
                analyzer.get_hydrophobic(get_all_atoms=self.hph_status)
                analyzer.get_ionic()
                table = analyzer.itable._table

                table.drop_duplicates(inplace=True)

                cha = table["chain_a"].tolist()
                chb = table["chain_b"].tolist()
                resa = table["resname_a"].tolist()
                resb = table["resname_b"].tolist()
                ida = table["resid_a"].tolist()
                idb = table["resid_b"].tolist()

                pairwises = [f"{z}{x}.{c}-{v}{b}.{n}" for z, x, c, v, b, n in zip(resa, ida, cha, resb, idb, chb)]
                resnamesa = [f"{a}{b}" for a, b in zip(resa, ida)]
                resnamesb = [f"{a}{b}" for a, b in zip(resb, idb)]

                table["pairwise"] = pairwises
                table.drop(["resname_a", "resname_b", "resid_a", "resid_b"], axis=1, inplace=True)

                table.insert(3, "residue_a", resnamesa)
                table.insert(4, "residue_b", resnamesb)

                return_list.append((analyzer.itable._table, i))  # stride'ı burada kullan
            return return_list

        def int_based_tbl(self):
            """Reads list from _calc_interactions function and writes to .csv file with respect to the frame number.
            
            Return: None
            """
            
            handle = open(os.path.join(self.target_path, "int_based_table.csv"), "w+", newline="")
            return_list = self._calc_interactions()
            if self.time_type == 'Time':
                if self.timeunit.lower() == 'ns' or self.timeunit.lower() == 'nanosecond':
                    t = 'Time (ns),'
            else:
                t = 'Frame,'
            handle.write(t + ",".join(return_list[0][0].columns.tolist()) + "\n")

            for tup in return_list:
                df = tup[0]
                frame = tup[1]

                if t != 'Frame,':
                    frame = float(frame) * self.timestep * self.stride

                df.insert(0, t.rstrip(','), frame)
                df.to_csv(handle, index=False, mode="w+", header=False)
            handle.close()

