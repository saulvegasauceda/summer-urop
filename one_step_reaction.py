"""
Created on Wed Jun 23 2021

@author: neuro
"""
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

rules_df = pd.read_excel("enz_rules.xlsx", sheet_name = "Method A")
start_df = pd.read_excel("Molecule Database.xlsx", sheet_name = "Start")


mol_smiles = start_df["SMILES"]
raw_cof_smiles = rules_df["SMILES"]

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

prod, rxn = one_step_reaction(mol_smiles, enz_rxns, 1000)

### Uncomment the following to create a file containing reaction information.
# #creating a dataframe for the file
# df = pd.DataFrame()

# #initializing list for columns
# start = []
# products = []
# reaction = []

# for i in rxn:
#     #iterating through each edge (rxn) from start to product
#     #adding the appropriate values to each list
#     start.append(i[0])
#     products.append(i[1])
#     reaction.append(i[2])

# # Creates the columns with the corresponding Titles
# df['Reactant SMILES'] = start
# df['Product SMILES'] = products
# df['Reaction SMARTS'] = reaction

# df.to_excel('Start Expansion.xlsx', index = False)
# #Check file for changes
