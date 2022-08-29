import os
import pandas as pd
from subprocess import PIPE, run

# a function to generate conformations into desired directory
def generate_conf(ident, smi, is_14_isomer, path, n_conf=3):
    
    if is_14_isomer == 1: # marker for 1,4-isomers
        isomer = '14'
    else:
        isomer = '15'
        
    # prepare ligand conformation
    for j in range(1, n_conf+1):
    
        command1 = ['obabel', f'-:{smi}', '-O', f'{path}/mol2/{ident}_{isomer}_{j}.mol2', '--gen3d', 'slow']
                    
        command2 = ['obabel', '-imol2', f'{path}/mol2/{ident}_{isomer}_{j}.mol2', '-opdbqt', '-O', f'{path}/pdbqt/{ident}_{isomer}_{j}.pdbqt']
        
        result = run(command1, stdout=PIPE, stderr=PIPE, universal_newlines=True)
        result = run(command2, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    
    return


if __name__ == "__main__":

    # create dir for molecules mol2 and pdbqt files
    os.makedirs('molecules', exist_ok=True)
    os.makedirs('molecules/mol2', exist_ok=True)
    os.makedirs('molecules/pdbqt', exist_ok=True)
    
    # load compounds to dock
    data_compunds = pd.read_csv('../data/click_compunds_dataset.csv', index_col=0)
    data_compunds = data_compunds[data_compunds['group_rishton'] == 1] # filter rishton compounds
    
    # apply function over molecules in dataframe
    data_compunds.apply(lambda row:
                  generate_conf(row.id, row.smiles, row.is_14isomer, path='./molecules', n_conf=3),
                  axis=1)

