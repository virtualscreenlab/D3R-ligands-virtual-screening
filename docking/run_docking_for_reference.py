import os
import pandas as pd
from subprocess import PIPE, run


# a function to generate conformations into desired directory
def generate_conf(smi, name, path, n_conf=3):

    # prepare ligand conformations
    for j in range(1, n_conf+1):

        command1 = ['obabel', f'-:{smi}',
                    '-O', f'{path}/{name}_{j}.mol2', '--gen3d', 'slow']

        command2 = ['obabel', '-imol2', f'{path}/{name}_{j}.mol2',
                    '-opdbqt', '-O', f'{path}/{name}_{j}.pdbqt']

        result = run(command1, stdout=PIPE, stderr=PIPE, universal_newlines=True)
        result = run(command2, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    return

def docking_vina(path, filename):
    path_to_vina = "path/to/vina/bin/vina" # specify path to vina
    
    name = filename[:-6]
    
    # dock the molecule
    command = [path_to_vina, "--receptor", "d3.pdbqt",
               "--ligand", f"{path}{filename}",
               "--center_x", "4.71443478",
               "--center_y", "13.39904348",
               "--center_z", "-9.40869565",
               "--size_x", "20", "--size_y", "20", "--size_z", "20",
               "--cpu", "8",
               "--exhaustiveness", "25", "--num_modes", "1",
               "--out", f"{path}docking/{name}_d3.pdbqt"]

    result = run(command, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    output = result.stdout.split("\n")
    
    # parse output
    try:
        affinity = float(output[-3].split()[1])
    except:
        affinity = 0
    
    return affinity


def docking_idock(path, filename):
    path_to_idock = "path/to/idock-2.2.3/bin/Linux/idock" # specify path to idock

    name = filename[:-6]

    # dock the molecule
    command = [path_to_idock, "--receptor", "d3.pdbqt",
               "--ligand", f"{path}{filename}",
               "--center_x", "4.71443478",
               "--center_y", "13.39904348",
               "--center_z", "-9.40869565",
               "--size_x", "20", "--size_y", "20", "--size_z", "20",
               "--out", f"{path}docking/{name}_idock_out"]

    result = run(command, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    output = result.stdout.split("\n")

    # parse output
    try:
        score = float(output[-2].split()[3])
        rf_score = float(output[-2].split()[4])
    except:
        score, rf_score = 0, 0
        print(f'Problems with docking molecule {name}')

    return score

if __name__ == "__main__":
    
    # cmoles for the eticlopride and scaffold
    epq = 'CCC1=CC(=C(C(=C1O)C(=O)NCC2CCCN2CC)OC)Cl'
    scaffold = 'N1N=NC=C1CCCCN1CCN(CC1)C1=CC=CC=C1'
    
    # create dir for molecules mol2 and pdbqt files
    os.makedirs('reference/epq', exist_ok=True)
    os.makedirs('reference/scaffold', exist_ok=True)
    
    # generate input conformations for docking
    generate_conf(epq, 'epq', './reference/epq', n_conf=3)
    generate_conf(scaffold, 'scaffold', './reference/scaffold', n_conf=3)
    
    # dock with vina
    print('Vina:')
    
    # create dirs for docking outputs
    os.makedirs('reference/epq/docking', exist_ok=True)
    os.makedirs('reference/scaffold/docking', exist_ok=True)
    
    for filename in os.listdir('./reference/epq'):
        f = os.path.join('./reference/epq', filename)
        if os.path.isfile(f) and filename[-5:] == 'pdbqt': # checking if it is a .pdbqt file:
            aff = docking_vina('./reference/epq/', filename)
            print(f'Affinity for {filename} is {aff}')
            
    for filename in os.listdir('./reference/scaffold'):
        f = os.path.join('./reference/scaffold', filename)
        if os.path.isfile(f) and filename[-5:] == 'pdbqt': # checking if it is a .pdbqt file:
            aff = docking_vina('./reference/scaffold/', filename)
            print(f'Affinity for {filename} is {aff}')
        
    print()
        
    # dock with idock
    print('iDock:')
    
    for filename in os.listdir('./reference/epq'):
        f = os.path.join('./reference/epq', filename)
        if os.path.isfile(f) and filename[-5:] == 'pdbqt': # checking if it is a .pdbqt file:
            aff = docking_idock('./reference/epq/', filename)
            print(f'Affinity for {filename} is {aff}')
            
    for filename in os.listdir('./reference/scaffold'):
        f = os.path.join('./reference/scaffold', filename)
        if os.path.isfile(f) and filename[-5:] == 'pdbqt': # checking if it is a .pdbqt file:
            aff = docking_idock('./reference/scaffold/', filename)
            print(f'Affinity for {filename} is {aff}')
