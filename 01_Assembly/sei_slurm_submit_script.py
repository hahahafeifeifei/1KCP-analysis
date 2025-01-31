from pathlib import Path

import submitit

# Define the work directory and job directory
workdir = Path('/sei_segment_length_impact')
yaml_dir = workdir / 'pending_jobs'
yaml_dir.mkdir(parents=True, exist_ok=True)

finished_jobs_dir = workdir / 'finished_jobs'
finished_jobs_dir.mkdir(parents=True, exist_ok=True)

sub_jobs_need_run = list(yaml_dir.glob('config_*.yaml'))
print('Total sub jobs:', len(sub_jobs_need_run))


# Set up the Submitit executor
executor = submitit.AutoExecutor(folder=workdir / "submitit_logs")

executor.update_parameters(
    slurm_mem='50G',
    gpus_per_node=1,
    cpus_per_task=10,
    timeout_min=60 * 24 * 8,  # Maximum execution time in minutes
    slurm_partition='a100-40g,a40-quad,a40-tmp,v100',
    slurm_qos="gpu-huge",
)

# Submit jobs
jobs = []
with executor.batch():
    for config_file in sub_jobs_need_run:
        sei_command_workflow = submitit.helpers.CommandFunction(
            command=[
            '/path/to/venv/python',
            '/path/to/sei_inference.py',
            config_file.as_posix(),
]
            ,
            verbose=True,
            env={'PATH': '/path/to/venv',}
        )
        job = executor.submit(sei_command_workflow)
        jobs.append(job)

# Optionally, you can wait for the jobs to finish and check the results
for i,job in enumerate(jobs):
    try:
        job.result()  # This will block until the job is finished
        print(f"Job {job.job_id} finished successfully")
        # create a done file
        success_yaml_config = sub_jobs_need_run[i]
        done_file = finished_jobs_dir/f"{success_yaml_config.stem}.done"
        done_file.touch()
    except Exception as e:
        print(f"Job {job.job_id} failed: {e}")