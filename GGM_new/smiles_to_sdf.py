#Script for converting GGM output to an SDF file (for subsequent docking)

import argparse
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd


def main(args):
    
    df = pd.read_csv(args.GGM_df, sep = '\t')
    df.columns = ['gen', 'idx', 'na']
    print(df.head())
    
    mols3D = []
    
    for smiles in df['gen']:
        try:
            mol = Chem.MolFromSmiles(smiles)
            molHs = Chem.AddHs(mol)
            AllChem.EmbedMolecule(molHs)
            mols3D.append(molHs)
        except:
            print(smiles)
    
    
    #Now want to write to file
    w = Chem.SDWriter(args.output_file)
    for m in mols3D:
        w.write(m)
        
    w.close()
    
    print(f'Wrote {len(mols3D)} molecules to file at {args.output_file}')
    
    return 0
    
if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('GGM_df', type=str, help='txt file containing elaborations suggested by GGM')
    parser.add_argument('output_file', type=str, help='Output SDF file name')

    arguments = parser.parse_args()
    
    main(arguments)