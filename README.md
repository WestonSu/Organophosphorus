# Identification and Prioritization of Emerging Organophosphorus Compounds Beyond Organophosphate Esters in Chinese Estuarine Waters (Unpublished paper)
<p align="left">
<img src="https://img.shields.io/badge/Python-3776AB.svg?style&logo=Python&logoColor=white" alt="Python" />
<img src="https://img.shields.io/badge/RStudio-75AADB.svg?style&logo=RStudio&logoColor=white" alt="RStudio" />
<img src="https://img.shields.io/badge/license-MIT-blue.svg" alt="MIT License" />
</p>
This repository provides a rational workflow for the identification and prioritization of pollutants in environmental samples by the combination of nontarget screening and MCDA prioritization approach.

---

## Overview

### 1. **Industrial Chemical Database** ([1_Industrial_chemicals.csv](https://github.com/WestonSu/Organophosphorus/blob/main/1_Industrial_chemicals.csv))
This database consists of **92,955** registered or pre-registered substances from:
- The Inventory of Existing Chemical Substances in China (IECSC)
- The Toxic Substances Control Act Inventory (TSCA) of the United States
- The Domestic Substances List (DSL) of Canada
- The Regulation for Registration, Evaluation, Authorisation, and Restriction of Chemicals (REACH) of the European Union

### 2. **OPC Extraction** ([2_Extract_OPCs.py](https://github.com/WestonSu/Organophosphorus/blob/main/2_Extract_OPCs.py))
The `2_Extract_OPCs.py` script extracts OPCs by:
- Filtering out compounds containing counterions (e.g., Na+/K+/Cl−/Br−) and metal/metalloid-containing compounds.
- Standardizing the resulting list into an "MS-ready" format for HRMS analysis.

### 3. **HRMS Data Processing** ([3_HRMS_data_processing.R](https://github.com/WestonSu/Organophosphorus/blob/main/3_HRMS_data_processing.R))
HRMS raw files are batch processed using the `patRoon` R package.

### 4. **Multi-Criteria Decision Analysis (MCDA)** ([4_MCDA.py](https://github.com/WestonSu/Organophosphorus/blob/main/4_MCDA.py))
This Python script implements the **ELECTRE III** strategy to prioritize chemicals based on their hazards. ELECTRE III is a multi-criteria decision analysis (MCDA) method used to rank options under uncertainty and complexity.

### 5. **Structural Similarity Calculation** ([5_Tanimoto.py](https://github.com/WestonSu/Organophosphorus/blob/main/5_Tanimoto.py))
This script calculates the structural similarity between compounds using the **Tanimoto coefficient**, which is commonly used in cheminformatics for molecular similarity assessments.


## Repository Contents
- `1_Industrial_chemicals.csv`: The combined industrial chemical database.
- `2_Extract_OPCs.py`: Script to extract OPCs and generate an MS-ready suspect list.
- `3_HRMS_data_processing.R`: R script for processing HRMS data.
- `4_MCDA.py`: Python implementation for chemical hazard prioritization using ELECTRE III.
- `5_Tanimoto.py`: Python script for calculating structural similarity between chemicals.

### Python Packages:
- `pandas`
- `numpy`
- `scikit-learn`

### R Packages:
- `patRoon`

## Citations 
We ask users to directly cite the following paper:

Su, W. et al. Identification and Prioritization of Emerging Organophosphorus Compounds Beyond Organophosphate Esters in Chinese Estuarine Waters. 

This project also builds on a number of other projects, algorithms and ideas. Please consider citing the following full list of papers when relevant: 

1. Su, W.; Li, P.; Zhong, L.; Liang, W.; Li, T.; Liu, J.; Ruan, T.; Jiang, G. Occurrence and distribution of antibacterial quaternary ammonium compounds in Chinese estuaries revealed by machine learning-assisted mass spectrometric analysis. Environ. Sci. Technol. 2024. 58 (26), 11707−11717. DOI: 10.1021/acs.est.4c02380
2. Helmus, R.; ter Laak, T. L.; van Wezel, A. P.; de Voogt, P.; Schymanski, E. L. patRoon: open source software platform for environmental mass spectrometry based non-target screening. J. Cheminform. 2021, 13 (1), 1. DOI: 10.1186/s13321-020-00477-w
3. Ruttkies, C.; Schymanski, E. L.; Wolf, S.; Hollender, J.; Neumann, S. MetFrag relaunched: incorporating strategies beyond in silico fragmentation. J. Cheminform 2016, 8 (1), 3. DOI: 10.1186/S13321-016-0115-9
4. Dührkop, K.; Fleischauer, M.; Ludwig, M.; Aksenov, A. A.; Melnik, A. V.; Meusel, M.; Dorrestein, P. C.; Rousu, J.; Böcker, S. SIRIUS 4: a rapid tool for turning tandem mass spectra into metabolite structure information. Nat. Methods 2019, 16 (4), 299−302. DOI: 10.1038/s41592-019-0344-8

