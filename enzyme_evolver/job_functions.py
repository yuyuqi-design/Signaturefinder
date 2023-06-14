import uuid
from enzyme_evolver import database_functions
from enzyme_evolver.mongo.workflow_models import Job
from enzyme_evolver.mongo.user_models import User
from flask_security import current_user
import subprocess as sp
from typing import List


def create_new_job(name, type, initial_status='In Queue', no_user=False):
    if no_user is False:
        user = User.objects(id=current_user.id)[0]
    else:
        user = None
    new_id = str(uuid.uuid4())
    database_functions.make_new_folder(new_id)
    job = Job(folder_id=new_id, owner=user, name=name, type=type, status=initial_status)
    job.save()

    return new_id

def delete_job(folder_id):
    job = Job.objects(folder_id=folder_id)[0]
    job.delete()
    try:
        database_functions.delete_folder(folder_id)
    except:
        print('Could not delete folder')

def update_a_job_status(job: Job, status: str, log=True):
    job.update_status(status, print_log=log)

def update_notes(job: Job):
    job.update_notes('job submitted')

def add_3dmodel_files(job: Job, rec_names: str):
    for rec_name in rec_names:
        rec_name = rec_name.strip()
        model = f'{rec_name}.pdb'
        job.add_file(model)

def add_binding_mode_structures(job: Job, rec_names: List[str], lig_structures: List[str]):
    lig_names = [lig_structure.split('.')[0] for lig_structure in lig_structures]
    for lig_name in lig_names:
        lig_name = lig_name.strip()
        for rec_name in rec_names:
            rec_name = rec_name.strip()
            binding_mode_structure = f'{rec_name}_{lig_name}_binding_poses.pdb'
            job.add_file(binding_mode_structure)

def add_ranking_file(job: Job, lig_structures: List[str]):
    lig_names = [lig_structure.split('.')[0] for lig_structure in lig_structures]
    for lig_name in lig_names:
        lig_name = lig_name.strip()
        ranking_file = f'{lig_name}_IREDFisher_ranking.csv'
        job.add_file(ranking_file)

def add_zipfile(job: Job, folder_id: str):
    job.update_status('Making zip file', print_log=True)
    zipfile = folder_id + '.zip'
    sp.run(f'zip -r {zipfile} *', shell=True)
    job.add_file({zipfile})
