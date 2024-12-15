import os
import MDAnalysis as mda
import subprocess

class tables_errors():
    def __init__(self) -> None:
        pass

    def test_inp_path(self, v):
        if not os.path.exists(v):
            raise Exception('PathNotFound: Please provide an existing path.')

    def test_stride(self, stride):
        if not isinstance(stride, int):
            raise TypeError('Please provide an integer value.')

    def test_split_models(self, v):
        if not isinstance(v, bool):
            raise TypeError('Please provide boolean split_model value.')
        
    def test_rmsd_refstruc(self, v):
        if isinstance(v, str):
            self.test_inp_path(v)

    def test_rmsd_refframe(self, v, l):
        self.test_stride(v)
        if l < v:
            raise Exception('Reference frame not within the simulation. Please provide a reference frame in the simulation')
        
    def test_rmsd_dict(self, v):
        if not isinstance(v, dict):
            raise Exception('Please provide an dictionary for RMSD Reference Data.')
        
    def test_table_float(self, v):
        if not isinstance(v, float) or not isinstance(v, int):
            raise TypeError('Please provide an integer or float value for stride.')
        if v<0 or v>=100:
            raise Exception('Thereshold values must be between 0 and 100.')
        
    def test_intf_path(self, v):
        if not os.path.exists(v):
            raise Exception('FileNotFound: Please provide an existing path for interface residues or run residue_based analysis.')
        
    def test_chain_sel(self, v, pdb):
        u = mda.Universe(pdb)
        c = list(u.segments.segids)
        if not isinstance(v, list):
            raise TypeError('Please provide chains in a list.')
        for el in v:
            if el not in c:
                raise Exception(f"Chain {el} could not found in pdb. Plase provide a valid chain ID.")

    def check_dssp(self, myos):
        commands = ['which', 'mkdssp']
        if myos == 'win32':
            commands.insert(0,'wsl')
        result = subprocess.run(
                commands,
                stdout = subprocess.PIPE,
            )
        if b"not found" in result.stdout:
            raise Exception("mkdssp not found. Please install mkdssp by 'sudo apt-get mkdssp'")
        
    def check_wsl(self):
        result = subprocess.run(
                ['where', 'wsl'],
                stdout = subprocess.PIPE,
            )
        if b"not found" in result.stdout:
            raise Exception("WSL not found. Please install WSL by 'wsl --install'")
        


