import os
import pandas as pd
from tqdm import tqdm
from subprocess import PIPE, run

# returns a single docking position for a single conformation
def dock_single_conformation(path, file_name):
    path_to_vina = "path/to/vina/bin/vina" # specify path to vina
    
    name = file_name[:-6]
    
    # dock the molecule
    command = [path_to_vina, "--receptor", "d3.pdbqt",
               "--ligand", f"{path}{file_name}",
               "--center_x", "4.71443478",
               "--center_y", "13.39904348",
               "--center_z", "-9.40869565",
               "--size_x", "20", "--size_y", "20", "--size_z", "20",
               "--cpu", "20",
               "--exhaustiveness", "25", "--num_modes", "1",
               "--out", f"./docking_vina_results/pdbqt/{name}_d3.pdbqt"]

    result = run(command, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    output = result.stdout.split("\n")
    
    # parse output
    try:
        affinity = float(output[-3].split()[1])
        
        # convert structure to pdb
        command_convert = ['obabel', '-ipdbqt', f'./docking_vina/pdbqt/{name}_d3.pdbqt', '-opdb', '-O', f'./docking_vina/pdb/{name}_d3.pdb']
        result = run(command_convert, stdout=PIPE, stderr=PIPE, universal_newlines=True)
        
    except:
        affinity = 0
        print(f'Could not retrieve affinity for molecule {file_name}')
        print(result.stdout)

    return affinity


if __name__ == "__main__":
    compound_dir = os.fsencode('./molecules/pdbqt')
    docking_data = [] # list to save docking data
    
    # create dirs for docking outputs
    os.makedirs('./docking_vina_results/pdbqt', exist_ok=True)
    os.makedirs('./docking_vina_results/pdb', exist_ok=True)

    for file in tqdm(os.listdir(compound_dir), desc="Running Vina"):
        # get name of the file
        filename = os.fsdecode(file)
        if 'pdbqt' in filename:
            tags = filename.split('_') # get molecule id, isomer, and conf id
            ident, isomer, conf = tags[0], tags[1], tags[2][:1]
            
            aff = dock_single_conformation('./molecules/pdbqt/', filename)
            docking_data.append([ident, isomer, conf, aff])
        
    df = pd.DataFrame(docking_data, columns =['id', 'isomer', 'conf_id', 'affinity'])
    df.to_csv('vina_data.csv')
