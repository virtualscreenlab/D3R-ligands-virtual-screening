import os
import pandas as pd
from subprocess import PIPE, run

# returns a single docking position for a single conformation
def dock_single_conformation(path, file_name):
    path_to_idock = "path/to/idock-2.2.3/bin/idock" # specify path to idock
    
    name = file_name[:-6]
    
    # dock the molecule
    command = [path_to_idock, "--receptor", "d3.pdbqt",
               "--ligand", f"{path}{file_name}",
               "--center_x", "4.71443478",
               "--center_y", "13.39904348",
               "--center_z", "-9.40869565",
               "--size_x", "20", "--size_y", "20", "--size_z", "20",
               "--out", f"./docking_idock_results/{name}_out"]

    result = run(command, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    output = result.stdout.split("\n")
    
    # parse output
    try:
        score = float(output[-2].split()[3])
        rf_score = float(output[-2].split()[4])
    except:
        score, rf_score = 0, 0
        print(f'Problems with docking molecule {name}')
    
    return score, rf_score


if __name__ == "__main__":
    compound_dir = os.fsencode('./molecules/pdbqt')
    docking_data = [] # list to save docking data
    
    # create dirs for docking outputs
    os.makedirs('./docking_idock_results', exist_ok=True)

    for file in os.listdir(compound_dir):
        # get name of the file
        filename = os.fsdecode(file)
        
        tags = filename.split('_') # get molecule id, isomer, and conf id
        ident, isomer, conf = tags[0], tags[1], tags[2][:1]
        
        idock_score, rf_idock_score = dock_single_conformation('./molecules/pdbqt/', filename)
        docking_data.append([ident, isomer, conf, idock_score, rf_idock_score])
        
    df = pd.DataFrame(docking_data, columns =['id', 'isomer', 'conf_id', 'idock_score', 'rf_score'])
    df.to_csv('idock_data.csv')
