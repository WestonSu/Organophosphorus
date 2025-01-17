# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 13:33:51 2024

@author: Wenyuan Su
"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem import MolStandardize, Descriptors
from pathlib import Path
# Define the path to the input CSV file and the output CSV file
input_csv_path = Path("OPC.csv")
output_csv_path = Path("filtered_smiles_OPC.csv")
# Read the input CSV file. Assuming there's no header and the SMILES are in the second column (index 1)
# Adjust the usecols and names parameters to read the first column as well
df = pd.read_csv(input_csv_path, header=None, usecols=[0, 1], names=['id', 'smiles'])
# Initialize the LargestFragmentChooser
lfc = MolStandardize.fragment.LargestFragmentChooser()
# Define the substructures we are looking for as SMILES strings
substructures_smiles = [
    "COP(O)(OC)=O", "COP(OC)(OC)=S", "COP(SC)(OC)=O", "COP(OC)(C)=O", "CP(C)(C)=O", "C[P+](C)(C)C"  
] # OPEs, OTPEs, OTPEs, OPNs, POs, QPCs
# Convert the SMILES strings to RDKit Mol objects for substructure searching
substructures_mols = [Chem.MolFromSmiles(smiles) for smiles in substructures_smiles]
def process_smiles(smiles):
    """Keep only the largest fragment of a given SMILES, standardize, handle errors, and check molecular weight."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mol2 = lfc.choose(mol)
            standardized_smiles = Chem.MolToSmiles(mol2, isomericSmiles=True)
            mol_weight = Descriptors.MolWt(mol2)
            if 100 <= mol_weight <= 800:
                return standardized_smiles
    except Exception as e:
        print(f"Could not process SMILES {smiles}: {str(e)}")
    return None
def contains_any_substructure(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return mol is not None and any(mol.HasSubstructMatch(substructure) for substructure in substructures_mols)
# Apply the processing function
df['processed_smiles'] = df['smiles'].apply(process_smiles)
df.dropna(subset=['processed_smiles'], inplace=True)
# Deduplicate based on processed_smiles
df_unique = df.drop_duplicates(subset=['processed_smiles'])
# Filter for substructures and check if the original SMILES (not just the processed ones) contain any substructures
df_filtered = df_unique[df_unique['processed_smiles'].apply(contains_any_substructure)]
# Select the original ID and SMILES columns for the output, instead of just the processed SMILES
output_df = df_filtered[['id', 'smiles']]
# Save the filtered dataframe to a new CSV file
output_df.to_csv(output_csv_path, index=False, header=False)
