from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
import pandas as pd

#accessing list of smiles string of both feedstock and fda approved drugs
start = pd.read_excel("Molecule Database.xlsx", sheet_name = "Start")
fda = pd.read_csv("fda.csv")

fda_smiles = fda["smiles"]
start_smiles = start["SMILES"]

#creating list of mols
fda_mol = [Chem.MolFromSmiles(f) for f in fda_smiles]
start_mol = [Chem.MolFromSmiles(s) for s in start_smiles]

#creating list of fingerprints
fda_fps = [FingerprintMols.FingerprintMol(f) for f in fda_mol]
start_fps = [FingerprintMols.FingerprintMol(s) for s in start_mol]

#intiliazes lists for dataframes
st, tg, sim = [], [], []

for i, m in enumerate(start_fps):
    #bulk comparing a start mol fp value with the entire list of fda drugs
    values = DataStructs.BulkTanimotoSimilarity(m, fda_fps)
    
    #collect the smiles and sim values
    for j, v in enumerate(values):
        sim.append(v)
        st.append(start_smiles[i])
        tg.append(fda_smiles[j])

#creating dataframe and sorting by sim values
d = {'Start':st, 'Target':tg, 'Similarity':sim}
df_final = pd.DataFrame(data=d)
df_final = df_final.sort_values('Similarity', ascending=False)

### Uncomement the follwoing to create file:
#df_final.to_excel('similarity.xlsx', index = False)
