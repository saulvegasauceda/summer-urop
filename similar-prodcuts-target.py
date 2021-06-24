from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
import pandas as pd

#accessing list of smiles string of both feedstock and fda approved drugs
start = pd.read_excel("Molecule Database.xlsx", sheet_name = "Start")
top_sim = pd.read_excel("similarity.xlsx", sheet_name = ">0.75")

prom_start = set(top_sim["Start"])
name = start["Name"]
st_smiles = start["SMILES"]
#removing duplicates
prom_start = list(set(prom_start))

#creating dictionary between start smiles and name
#Note: This was not used to make the file but will be later on
st = dict(zip(st_smiles, name))


rules_df = pd.read_excel("enz_rules.xlsx", sheet_name = "Method A")
start_df = pd.read_excel("Molecule Database.xlsx", sheet_name = "Start")
cf = pd.read_excel("et_rules.xlsx", sheet_name = "co_factors")

mol_smiles = start_df["SMILES"]
#list of cof smiles per reaction
raw_cof_smiles = rules_df["SMILES"]
#lists of cof smiles
cof = cf["SMILES"]

# splits the cofactors SMILES strings into lists
cof_smiles = []
for c in raw_cof_smiles:
    cof = c.split(';')
    cof_smiles.append(cof)
    
rxn_smarts = rules_df["SMARTS"]
num_react = rules_df["Reactant Number"]

#creates a list of tuples w/ the following format:
# (smiles string of cofactors, number of reactants, smarts of reaction)
enz_rxns = list(zip(cof_smiles, num_react, rxn_smarts))

def one_step_reaction(molecule_smiles, enzyme_reaction_tuples, cutoff):
    """
    parameters:
    molecule_smiles: list of string
        SMILES strings of the starting molecules
    enzyme_reaction_tuples: list of tuples
        tuples contain the following information for each reaction: 
        (SMILES string of necessary cofactors, # of reactants, SMARTS patterns for the reaction)
    cutoff: molecular weight cutoff in g/mol
    
    return:
    nodes: set of string
    SMILES of all of the new products and initial molecules
    edges: set of tuples
    (SMILES of reactants, SMILES of products, SMARTS of rxn)
    """
    nodes = set()
    edges = set()
    #run expansion loop
    
    for mol in molecule_smiles:
        for rxn_tup in (enzyme_reaction_tuples):
            #initiliazes rxn object
            rxn = AllChem.ReactionFromSmarts(rxn_tup[2])
            #get number of reactants necessary
            rct_num = rxn_tup[1]
            #initializing reactants list
            rct = []
            if rxn_tup[0] != 'none':
                #adds cofactors when necessary
                rct = rxn_tup[0].copy()
            for i in range(rct_num):
                #adds starting molecule
                rct.append(mol)
            #removing any empty cells
            rct_list = [r for r in rct if r != "none"]
            
            
            reactants = [Chem.MolFromSmiles(r) for r in rct_list]
            for outcome in rxn.RunReactants(reactants):
                #can not test for mol validity at this point!!!
                for prod_mol in outcome:
                    #unmap atom
                    for atom in prod_mol.GetAtoms():
                        atom.SetAtomMapNum(0)
                    prod_smiles = Chem.MolToSmiles(prod_mol)
                    prod = prod_smiles.split(".")
                    for p in prod:
                        if p not in cof:
                            p_mol = Chem.MolFromSmiles(p)
                            #checking for valid mol and MWC
                            if p_mol is not None and Chem.rdMolDescriptors.CalcExactMolWt(p_mol) < cutoff:
                                nodes.add(p)
                                edges.add((mol, p, rxn_tup[2]))

    #include starting materials in products
    return nodes.union(molecule_smiles), edges

prod, rxn = one_step_reaction(prom_start, enz_rxns, 1000)

#accessing smiles
fda = pd.read_csv("fda.csv")
fda_smiles = fda["smiles"]
prod_smiles = list(prod)

#creating list of mols
fda_mol = [Chem.MolFromSmiles(f) for f in fda_smiles]
prod_mol = [Chem.MolFromSmiles(p) for p in prod]

#creating list of fingerprints
fda_fps = [FingerprintMols.FingerprintMol(f) for f in fda_mol]
prod_fps = [FingerprintMols.FingerprintMol(p) for p in prod_mol]

#intiliazes lists for dataframes
pr, tg, sim = [], [], []

for i, m in enumerate(prod_fps):
    #bulk comparing a start mol fp value with the entire list of fda drugs
    values = DataStructs.BulkTanimotoSimilarity(m, fda_fps)
    
    #collect the smiles and sim values (only recording values surpassing similarity threshold)
    for j, v in enumerate(values):
        if v > 0.85:
            sim.append(v)
            pr.append(prod_smiles[i])
            tg.append(fda_smiles[j])

#creating dataframe and sorting by sim values
df = pd.DataFrame({'Product':pr, 'Target': tg, 'Similarity':sim})
df = df.sort_values('Similarity', ascending=False

### Uncomement the follwoing to create file:
#df.to_csv('prod_similarity.csv', index = False)