import rdkit
from rdkit import Chem
import pandas as pd
import rdkit.Chem as Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole




residues = Chem.SDMolSupplier('possible_amines_backbone.sdf')
smiles = [Chem.MolToSmiles(mol) for mol in residues if mol]



rosetta_residues = Chem.SDMolSupplier('rosetta_residues_backbone.sdf')
#for mol in rosetta_residues:
    #if mol:
        #print(mol.GetPropNames())  # List all properties
        #print(mol.GetProp('ResidueName'))


smiles_rosetta = [Chem.MolToSmiles(mol) for mol in rosetta_residues if mol]
mols = [Chem.MolFromSmiles(smiles) for smiles in smiles_rosetta]
Draw.MolsToGridImage(mols)



df = pd.DataFrame({'SMILES': smiles_rosetta})


import random
import string

def generate_random_code():
    return ''.join(random.choices(string.ascii_uppercase + string.digits, k=3))
df['3 letter code']=None
df['3 letter code'] = df['3 letter code'].apply(lambda x: generate_random_code() if pd.isna(x) or x == '' else x)


ligprep_input = df[['SMILES','3 letter code']]
ligprep_input.rename(columns={'3 letter code': 'code'}, inplace=True)
ligprep_input.to_csv('rosetta_ligprep_input.csv',index=False)


# ### reaction smarts for primary amine with one substitution and no other primary amine

df['Full_SMILES']=None



def transform_smiles(row):
    if pd.isna(row['SMILES']): 
        return row['Full_SMILES']
    
    if pd.isna(row['Full_SMILES']):  
        mol = Chem.MolFromSmiles(row['SMILES'])  
        
        if mol:
            # Run the reaction
            products = reaction_1.RunReactants((mol,))
            
            # Check if products list is non-empty and contains valid molecules
            if products and products[0]:  
                product_mol = products[0][0]  
                
                if product_mol:  
                    return Chem.MolToSmiles(product_mol)  
            
      
        return row['SMILES']

    return row['Full_SMILES']


df['Full_SMILES'] = df.apply(transform_smiles, axis=1)



ligprep_input = df[['Full_SMILES','3 letter code']]
ligprep_input.rename(columns={'Full_SMILES': 'SMILES', '3 letter code': 'code'}, inplace=True)


ligprep_input.to_csv('possible_amines_from_review_for_ligprep2.csv',index=False)


# # generate images of molecules to add to datasheet




import os
ligprep_input = ligprep_input.dropna()
ligprep_input['mol'] = ligprep_input['SMILES'].apply(lambda x: Chem.MolFromSmiles(x))

os.makedirs("molecule_images", exist_ok=True)


ligprep_input['image_path'] = ligprep_input['code'].apply(
    lambda name: f"molecule_images/{name}.png"
)
for _, row in ligprep_input.iterrows():
    if row['mol']:
        Draw.MolToFile(row['mol'], row['image_path'])





