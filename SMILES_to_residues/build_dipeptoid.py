# # build capped dipeptoid residues suitable for geometry optimization from SMILES of primary amines

# visualize smarts https://smarts.plus/view/446ad5e1-dd1b-44ec-9053-e26597626596
# https://drzinph.com/learning-reaction-smarts-a-practical-guide-to-reaction-based-patterns/
# https://daylight.com/meetings/summerschool00/course/basics/smirks.html
# https://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html

import pandas as pd
import rdkit.Chem as Chem
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import rdBase
rdBase.DisableLog('rdApp.error')
rdBase.DisableLog('rdApp.warning')
rdBase.DisableLog('rdApp.info')


df = pd.read_csv('../../Downloads/peptoid_sidechains - Sheet1.csv')



#only use the SMILES column
smiles = df[['SMILES']]

clean_smiles = smiles.dropna()

#remove the astrisk
clean_smiles['SMILES'] = clean_smiles['SMILES'].str.replace(r'\*', '', regex=True)


reaction_smarts_1 = '[N;H2:1]>>[C][C](=[O])[N:1][C:2][C:3](=[O:4])[N:1][C:6][C:7](=[O:8])[N:9][C:10]'
reaction_smarts_2 = '[O][C](=[O])[C][N:1]>>[C][C](=[O])[N:1][C:2][C:3](=[O:4])[N:1][C:6][C:7](=[O:8])[N:9][C:10]' 
reaction_smarts_3 = '[C](=[O])[C][N:1]>>[C][C](=[O])[N:1][C:2][C:3](=[O:4])[N:1][C:6][C:7](=[O:8])[N:9][C:10]' #if no oxygen after carbonyl

reaction_1 = AllChem.ReactionFromSmarts(reaction_smarts_1)
reaction_2 = AllChem.ReactionFromSmarts(reaction_smarts_2)
reaction_3 = AllChem.ReactionFromSmarts(reaction_smarts_3)


def apply_reactions(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    for reaction in [reaction_1, reaction_2, reaction_3]:
        products = reaction.RunReactants((mol,))
        if products and products[0]:
            return Chem.MolToSmiles(products[0][0], canonical=True)

    return None



clean_smiles['product'] = clean_smiles['SMILES'].apply(apply_reactions)



from rdkit import Chem
from rdkit.Chem import Draw
product_mols = clean_smiles['product'].dropna().apply(Chem.MolFromSmiles)


product_mols = [mol for mol in product_mols if mol is not None]


Draw.MolsToGridImage(product_mols, molsPerRow=4, subImgSize=(200, 200))


