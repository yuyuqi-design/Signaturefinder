import os
from pathlib import Path
import subprocess as sp
import shutil


# create a job folder under the database folder
def make_new_folder(id):
    print('Creating new folder at:')
    path_to_folder = f"{Path(__file__).parents[0]}/database/{id}"
    print(path_to_folder)
    os.mkdir(path_to_folder)

# delete a job folder
def delete_folder(folder_id):
    path = f"{Path(__file__).parents[0]}/database/{folder_id}"
    print(f'Deleting folder - {path}')
    shutil.rmtree(path)

# copy the reference structure files from database folder to job folder
def copy_references(data_folder):
    structures = ['singleT.pdb', 'dimerT.pdb', 'ref_complex.pdb', 'ref_lig.pdb', 'ref_cof.pdb']
    #copy the structures
    for structure in structures:
        database_folder = f'{Path(__file__).parents[0]}/database'
        sp.run(f'cp {database_folder}/{structure} {data_folder}', shell=True)

# copy the prepared proteins from database folder to the job folder
def prepared_recs(iredfisher_database: str, database_folder:str, data_folder: str):
    # copy the structures
    sp.run(f'cp {database_folder}/{iredfisher_database}/* {data_folder}', shell=True)
    with open(f'{database_folder}/{iredfisher_database}/receptors.txt') as fin:
        proteins = fin.readlines()
        codes = [code.split('.')[0].strip() for code in proteins]

    return codes

