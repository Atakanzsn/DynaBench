# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 11:42:38 2024

@author: inant
"""

import numpy as np
import MDAnalysis as mda
import os

cutoff = 4.5  # cutoff distance
no_contact = 15000  # assumption for the number of contacts

# Load the trajectory using MDAnalysis
def run_ContactMatrix(pdb_path, dcd_path, stride, output_path): ## Dynabench inputs
    u= mda.Universe(pdb_path, dcd_path, all_coordinates = False, continuous = False)
    print(f"MDAnalysis Universe created with PDB: {pdb_path} and DCD: {dcd_path}")

# Iterate over each frame in the trajectory
    for ts in u.trajectory [::stride]:
        print(f"Processing frame: {ts.frame}")
    
    # Initialize dictionaries to store coordinates
        residue_coordinatesX = {}
        residue_coordinatesY = {}
        residue_coordinatesZ = {}
        residue_atom_counts = {}
        residue_chainIDs = {}
    
    # Extract coordinates of atoms belonging to each residue
        for residue in u.residues:
            residue_number = (residue.segid, residue.resid)
            chain_id = residue.segid
            if residue_number not in residue_coordinatesX:
                residue_coordinatesX[residue_number] = []
                residue_coordinatesY[residue_number] = []
                residue_coordinatesZ[residue_number] = []
                residue_atom_counts[residue_number] = 0
                residue_chainIDs[residue_number] = chain_id
        
            for atom in residue.atoms:
                residue_coordinatesX[residue_number].append(atom.position[0])
                residue_coordinatesY[residue_number].append(atom.position[1])
                residue_coordinatesZ[residue_number].append(atom.position[2])
                residue_atom_counts[residue_number] += 1
    
    # Convert the dictionary values to numpy arrays
        for residue_number, coords in residue_coordinatesX.items():
            residue_coordinatesX[residue_number] = np.array(coords)
        for residue_number, coords in residue_coordinatesY.items():
            residue_coordinatesY[residue_number] = np.array(coords)
        for residue_number, coords in residue_coordinatesZ.items():
            residue_coordinatesZ[residue_number] = np.array(coords)
    
    # Create arrays to store coordinates
        max_atomsX = max(len(coords) for coords in residue_coordinatesX.values())
        max_atomsY = max(len(coords) for coords in residue_coordinatesY.values())
        max_atomsZ = max(len(coords) for coords in residue_coordinatesZ.values())
    
        residueX_array = np.full((len(residue_coordinatesX), max_atomsX), np.nan)
        residueY_array = np.full((len(residue_coordinatesY), max_atomsY), np.nan)
        residueZ_array = np.full((len(residue_coordinatesZ), max_atomsZ), np.nan)
    
        #residueX_array[np.isnan(residueX_array)] = 0
        #residueY_array[np.isnan(residueY_array)] = 0
        #residueZ_array[np.isnan(residueZ_array)] = 0
    
    # Fill the arrays with coordinates
        for i, coords in enumerate(residue_coordinatesX.values()):
            residueX_array[i, :len(coords)] = coords
        for i, coords in enumerate(residue_coordinatesY.values()):
            residueY_array[i, :len(coords)] = coords
        for i, coords in enumerate(residue_coordinatesZ.values()):
            residueZ_array[i, :len(coords)] = coords
    
        residue_array = np.array(list(residue_atom_counts.values()))
    
        total_residues = len(residue_array)
    
    # Calculation of Contact Matrix
        contact = np.zeros((no_contact, 100), dtype=float)  # zeros matrix for contacts calculation
        dx = []
        dy = []
        dz = []
        weight = []
        matrix_final = np.zeros((no_contact, 5), dtype=object)
        # Updated dtype to object to accommodate chain IDs
    
        for i in range(0, total_residues):
            for j in range(i + 1, total_residues):
                contact = 0
                for k in range(residue_array[i]):
                    for l in range(residue_array[j]):
                        dx.append(residueX_array[i, k] - residueX_array[j, l])
                        dy.append(residueY_array[i, k] - residueY_array[j, l])
                        dz.append(residueZ_array[i, k] - residueZ_array[j, l])
                        r2 = dx[-1] * dx[-1] + dy[-1] * dy[-1] + dz[-1] * dz[-1]
                        if i != j and r2 <= cutoff * cutoff:
                            contact += 1
            
                if contact != 0:
                    weight.append(contact / (residue_array[i] ** 0.5) / (residue_array[j] ** 0.5))
                    matrix_final[len(weight) - 1] = [residue_chainIDs[list(residue_coordinatesX.keys())[i]], list(residue_coordinatesX.keys())[i], residue_chainIDs[list(residue_coordinatesX.keys())[j]], list(residue_coordinatesX.keys())[j], weight[-1]]
    
        filtered_matrix = np.array([row for row in matrix_final if row[1] != 0 and row[3] != 0 and row[4] != 0], dtype=object)
        output_files = os.path.join(target_path, "matrix_final_{ts.frame}.out")
        np.savetxt(output_files, filtered_matrix, fmt=['%s', '%d', '%s', '%d', '%.6f'], delimiter='\t')
    return output_files
