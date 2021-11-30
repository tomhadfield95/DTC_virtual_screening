import argparse 
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem


def remove_dummy_atom(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol2 = AllChem.ReplaceSubstructs(mol, Chem.MolFromSmiles('*'), Chem.MolFromSmiles('[H]'), True)[0]
    mol3 = Chem.RemoveHs(mol2)

    return Chem.MolToSmiles(mol3)


def join_elab_to_fragment(frag, elab):

    try:
    #remove dummy atom from elab
        elab = remove_dummy_atom(elab)

        combined_mol = AllChem.ReplaceSubstructs(Chem.MolFromSmiles(frag), Chem.MolFromSmiles('*'), Chem.MolFromSmiles(elab), True)[0]

        return Chem.MolToSmiles(combined_mol)

    except:
        return None


def add_full_mols_to_deepfrag_df(df, smiles):
    df['full'] = [join_elab_to_fragment(smiles, row['Fragment SMILES']) for idx, row in df.iterrows()]
    df['frag'] = [smiles for idx, row in df.iterrows()]

    df = df.loc[df['full'] is not None]

    return df


def main(args):

    with open(args.fragment_smiles, 'r') as f:
        fragment_smiles = f.read()

    deepfrag_df = pd.read_csv(args.deepfrag_df)

    deepfrag_df['full'] = [join_elab_to_fragment(fragment_smiles, row['Fragment SMILES']) for idx, row in deepfrag_df.iterrows()]

    #Now generate 3d conformers
    mols3D = []

    for s in deepfrag_df['full'].head(250):
        try:
            m = Chem.MolFromSmiles(s)
            mHs = Chem.AddHs(m)
            AllChem.EmbedMolecule(mHs)
            mols3D.append(mHs)
        except:
            print(s)

    w = Chem.SDWriter(args.output_SDF)
    for m in mols3D:
        w.write(m)
    w.close()

    print(f'Wrote {len(mols3D)} molecules to file at {args.output_SDF}')

    return 0

    

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('deepfrag_df', type=str, help='CSV file containing DeepFrag output')
    parser.add_argument('output_SDF', type = str, help = 'Location to save output SDFs')
    parser.add_argument('--fragment_smiles', '-f', type = str,
                        help = 'Location of text file containing fragment smiles')
   
    

    arguments = parser.parse_args()
    
    main(arguments)

    main(arguments)







