# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 16:31:58 2023

@author: Wenyuan Su
"""

# Importing required modules from RDKit and other libraries
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import numpy as np

# Reading semi-quantitative compound IDs and SMILES from the file
# The file contains tab-separated values (IDs and SMILES strings)
f = open(r"C:\Users\Zz\Desktop\OPC_tanimoto_Semi.txt")
line = f.readline()
Semi_list_ID = []
while line:
      num = list(map(str,line.split('\t')))
      Semi_list_ID.append(num)
      line = f.readline()
f.close()
Semi_data = np.array(Semi_list_ID) # Convert the list to a NumPy array for easier handling

# Reading target compound IDs and SMILES from the file
# The file contains tab-separated values (IDs and SMILES strings)
f = open(r"C:\Users\Zz\Desktop\OPC_tanimoto_Target.txt")
line = f.readline()
Target_list_ID = []
while line:
      num = list(map(str,line.split('\t')))
      Target_list_ID.append(num)
      line = f.readline()
f.close()
Target_data = np.array(Target_list_ID) # Convert the list to a NumPy array for easier handling

# Open a file for writing Tanimoto similarity results
with open('TanimotoSimilarity1.txt', 'w') as f:
    # Write the header (Target compound IDs) in the first row
    for i in range(1,Target_data.shape[0]):
        f.write("\t" + Target_data[i,0])
    f.write("\n")
    # Iterate over all semi-quantitative compounds
    for j in range(1,Semi_data.shape[0]):
        f.write( Semi_data[j,0]) #写入行名
        print(Semi_data[j,1])
        # Convert SMILES to molecular structure for the semi-quantitative compound
        m1 = Chem.MolFromSmiles(Semi_data[j,1])
        # Generate a Morgan fingerprint (circular fingerprint) for the semi-quantitative compound
        # Radius is set to 2, nBits=64 for the size of the fingerprint vector
        fp1 = AllChem.GetMorganFingerprintAsBitVect(m1, 2, nBits=64, useFeatures=True)
        # Iterate over all target compounds
        for k in range(1,Target_data.shape[0]):
            print(Target_data[k,1])
            m2 = Chem.MolFromSmiles(Target_data[k,1])
            fp2 = AllChem.GetMorganFingerprintAsBitVect(m2, 2, nBits=64, useFeatures=True)
            # Calculate the Tanimoto similarity between the two fingerprints (fp1 and fp2)
            Similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
            print(Similarity)
            # Write the similarity value to the file, separating by tab
            f.write("\t" + str(Similarity))
        f.write("\n")
