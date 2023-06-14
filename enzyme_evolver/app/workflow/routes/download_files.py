from enzyme_evolver.app.workflow import bp
from flask import redirect, send_file
from enzyme_evolver.mongo.workflow_models import Job
from enzyme_evolver.mongo.user_models import User
from flask_security import current_user
from pathlib import Path
working_dir = Path(__file__).parents[4]

@bp.route('/download_file/<folder_id>/<file_name>', methods=['GET'])
def download_file(folder_id, file_name):
    job = Job.objects(folder_id=folder_id)[0]
    user = User.objects(id=current_user.id)[0]

    if job.owner != user:
        return redirect('/')
    elif file_name not in job.files:
        return redirect('/')

    path_to_file = f"{Path(__file__).parents[4]}/enzyme_evolver/database/{folder_id}/{file_name}"

    return send_file(path_to_file)


if __name__ == "__main__":
    working_dir = Path(__file__).parents[4]
    print(working_dir)
