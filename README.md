# summer-urop
## Computer-aided discovery of synthetic routes to high-value molecules from renewable feedstocks

Files posted here contain relevant programs to perfrom the first product expansion from the initial list of starting compounds.
Necessary data files will be shared here (https://www.dropbox.com/sh/ojx2hrvknry7vov/AAAaA-_2HGcGANMiBpwr8s1Ia?dl=0)
Target compound SMILES were downloaded from https://zinc.docking.org. Compounds categorized as FDA approved were used.
Enzymatic reaction rules used in the first expansion were from a NorthWestern study that automatically derived rules from KEGG and MetaCyc. The group created new SMARTS that reduced the redundancy and faults of the named databases (paper: https://www.sciencedirect.com/science/article/abs/pii/S1096717621000288).
Note: I was not able to include my main notebook due to the upload size constraints.
___

### Files:

Post Expansion NB.ipynb: notebook containing code I used to get the similarity scores

one_step_reaction.py: expands the starting material from ~50 compounds to ~40,000

start-target-tscores.py: compares the SMILES of the starting ~50 compounds to ~1000 targets (FDA-approved compounds) using Tanimoto Similarity Scoring

similar-prodcuts-target.py: compares the SMILES of the products of small-scale expansion from 15 starting compounds (Tscores > 0.75) to ~1000 targets (FDA-approved compounds) using Tanimoto Similarity Scoring
