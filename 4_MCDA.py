# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 13:32:43 2024

@author: Zz
"""

import pandas as pd
import numpy as np
# Load the dataset
df = pd.read_csv('ELECTRE III_PBTMQD.csv', sep=",") 
col_names = df.columns
# Get the indices of the special rows (indifference, preference, and veto)
indif_index = df[(df.value == 'indifference')].index.values[0]
prefer_index = df[(df.value == 'preference')].index.values[0] 
veto_index = df[(df.value == 'veto')].index.values[0]

# Get the column indices for specific evaluation criteria
T_start = df.columns.get_loc("Mutagenicity")  # T start index
Q_start = df.columns.get_loc("conc")          # Q start index
# Negate certain criteria values for correct evaluation (higher is safer, lower is more harmful)
df.iloc[0:indif_index, T_start:Q_start] *= -1

# Concordance matrix calculation
df_dict = dict()
for t in range(1, df.shape[1]):  # Iterate over columns
    a = [[] for i in range(0, indif_index)]
    df_dict[df.columns[t]] = []
    for i in range(0, indif_index):  # Iterate over rows
        for j in range(0, indif_index):  # Iterate over columns
            if df.iloc[i, t] - df.iloc[j, t] > df.iloc[prefer_index, t]:
                a[i].append(0)
            elif df.iloc[i, t] - df.iloc[j, t] < df.iloc[indif_index, t]:
                a[i].append(1)
            else:
                a[i].append((df.iloc[prefer_index, t] - df.iloc[i, t] + df.iloc[j, t]) / (df.iloc[prefer_index, t] - df.iloc[indif_index, t]))
    df_dict[df.columns[t]].append(a)

# Equally weighted criteria calculation
b = [[] for i in range(0, indif_index)]
P_end = df.columns.get_loc("Biodegradebility")  # P criteria end index
T_end = df.columns.get_loc('Algae (log mg/L)')  # T  criteria end index
for i in range(0, indif_index):  # Iterate over rows
    for j in range(0, indif_index):  # Iterate over columns
        element_P = 0
        element_B = 0
        element_T = 0
        element_M = 0
        element_Q = 0
        element_D = 0
        # Calculate weights for Persistence (P) criteria
        for t in range(0, P_end):
            element_P += df_dict[list(df_dict.keys())[t]][0][i][j]
        # Calculate weight for Biodegradability (B) criteria
        element_B = df_dict[list(df_dict.keys())[P_end]][0][i][j]
        # Calculate weights for Toxicity (T) criteria
        for t in range(T_start-1, T_end):
            element_T += df_dict[list(df_dict.keys())[t]][0][i][j]
        # Calculate weight for Mutagenicity (M) criteria
        element_M = df_dict[list(df_dict.keys())[-3]][0][i][j]
        # Calculate weight for Concentration (Q) criteria
        element_Q = df_dict[list(df_dict.keys())[-2]][0][i][j]
        # Calculate weight for Degradation (D) criteria
        element_D = df_dict[list(df_dict.keys())[-1]][0][i][j]
        # Combine weights for the final evaluation
        b[i].append((element_P/4 + element_B  + element_M + element_T/13 + element_Q + element_D)/6)  #PBMTQD

# Discordance matrix calculation
discordance_dict = dict()
for t in range(1, df.shape[1]):  # Iterate over columns
    a = [[] for i in range(0, indif_index)]
    discordance_dict[df.columns[t]] = []
    for i in range(0, indif_index):  # Iterate over rows
        for j in range(0, indif_index):  # Iterate over columns
            if df.iloc[i, t] - df.iloc[j, t] < df.iloc[prefer_index, t]:
                a[i].append(0)
            elif df.iloc[i, t] - df.iloc[j, t] > df.iloc[veto_index, t]:
                a[i].append(1)
            else:
                a[i].append((df.iloc[prefer_index, t] - df.iloc[i, t] + df.iloc[j, t]) / (df.iloc[prefer_index, t] - df.iloc[veto_index, t]))
    discordance_dict[df.columns[t]].append(a)

# Tj(a,b) calculation
Tj_dict = dict()
for t in range(1, df.shape[1]):  # Iterate over columns
    a = [[] for i in range(0, indif_index)]
    Tj_dict[df.columns[t]] = []
    for i in range(0, indif_index):  # Iterate over rows
        for j in range(0, indif_index):  # Iterate over columns
            m = discordance_dict[df.columns[t]][0][i][j]
            n = b[i][j]
            if m > n:
                a[i].append((1-m) / (1-n))
            else:
                a[i].append(1)
    Tj_dict[df.columns[t]].append(a)

# S(a,b) calculation
S = [[] for i in range(0, indif_index)]
for i in range(0, indif_index):  # Iterate over rows
    for j in range(0, indif_index):  # Iterate over columns
        if i == j:
            S[i].append(0)
            continue
        element = b[i][j]
        for key in Tj_dict:
            element *= Tj_dict[key][0][i][j]
        S[i].append(element)

# Convert S to numpy array and calculate row and column sums
S_array = np.array(S)
S_sumrows = np.sum(S_array, axis=1)
S_sumcols = np.sum(S_array, axis=0)

# Lambda calculation
lambda_S = max(map(max, S))
# g(lambda) calculation
g_lambda = -0.15 * lambda_S + 0.3

# T(a,b) calculation
T = [[] for i in range(0, indif_index)]
for i in range(0, indif_index):  # Iterate over rows
    for j in range(0, indif_index):  # Iterate over columns
        if i == j:
            T[i].append(0)
            continue
        elif S[i][j] >= (lambda_S - g_lambda):
            T[i].append(1)
        else:
            T[i].append(0)
# Convert T to numpy array and calculate row and column sums
T_array = np.array(T)
T_sumrows = np.sum(T_array, axis=1)
T_sumcols = np.sum(T_array, axis=0)

# Rank calculation using raw scores
from scipy.stats import rankdata
raw_score = S_sumrows - S_sumcols
raw_rank = len(raw_score) - rankdata(raw_score, method='max') + 1

# Rank calculation using defuzzified scores
defuzzified_score = T_sumrows - T_sumcols
defuzzified_rank = len(defuzzified_score) - rankdata(defuzzified_score, method='max') + 1

# Write the results to a text file
with open('ELECTRE_III_Rank_PBMTQD.txt', 'w') as f:
    f.write('Name\t' + 'raw_score\t' + 'raw_rank\t' + 'defuzzified_score\t' + 'defuzzified_rank\n')
    for i in range(0, indif_index):
        f.write(df.iloc[i, 0] + '\t' + str(raw_score[i]) + '\t' + str(raw_rank[i]) + '\t' + str(defuzzified_score[i]) + '\t' + str(defuzzified_rank[i]) + '\n')

