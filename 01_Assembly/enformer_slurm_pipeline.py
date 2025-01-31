# %%
import argparse
import logging
import random
import time

import h5py
import pyBigWig
import submitit
import torch
import torch.nn.functional as F
from enformer_pytorch import Enformer
from pyfaidx import Fasta
from tqdm import tqdm

# %%

'''
This is a pytorch implementation of Enformer from:
'https://github.com/lucidrains/enformer-pytorch'
'''

SEGMENT_LENGTH = 1536 * 128
TARGET_LENGTH = 896 * 128


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta_file', type=Path, required=True)
    parser.add_argument('--sample_save_dir', type=Path, required=True)
    parser.add_argument('--cuda_device', type=str, default='cuda:1')
    parser.add_argument('--embeddings_save_path', type=Path, required=True)
    parser.add_argument('--human_head_save_path', type=Path, required=True)
    parser.add_argument('--mouse_head_save_path', type=Path, required=True)

    args = parser.parse_args()
    return args


def infer_one_contig(seq, enformer, cuda_device, contig_name: str, save_mouse_head=False):
    seq = seq.type(torch.LongTensor).to(cuda_device)
    CHROM_LENGTH = len(seq)
    # segment_start = 0
    # target_start = segment_start + trim_length // 2
    target_start = 0
    target_end = target_start + TARGET_LENGTH
    extra_seq = SEGMENT_LENGTH - TARGET_LENGTH
    extra_left_seq = extra_seq // 2
    extra_right_seq = extra_seq - extra_left_seq
    # segment_end = segment_start + segment_length
    embeddings_chrome = []
    human_head_chrome = []
    mouse_head_chrome = []
    with torch.no_grad():
        while target_start < CHROM_LENGTH:
            left_padding = right_padding = 0
            segment_start = target_start - extra_left_seq
            segment_end = target_end + extra_right_seq
            if segment_start < 0:
                left_padding = -segment_start
                segment_start = 0

            if segment_end > CHROM_LENGTH:
                right_padding = segment_end - CHROM_LENGTH
                segment_end = CHROM_LENGTH
            input_seq = seq[segment_start:segment_end]
            input_seq = F.pad(input_seq, (left_padding, right_padding), value=-1)
            out, x = enformer(input_seq,
                              return_embeddings=True,
                              )
            embeddings_chrome.append(x.cpu().half())
            human_head_chrome.append(out['human'].cpu().half())
            if save_mouse_head:
                mouse_head_chrome.append(out['mouse'].cpu().half())
            target_start = target_end
            target_end = target_start + TARGET_LENGTH
            progress_bar.update(TARGET_LENGTH)
    target_chromosome_length = CHROM_LENGTH // 128
    embeddings_chrome_cat = torch.cat(embeddings_chrome, dim=0).cpu().numpy()[:target_chromosome_length]
    human_head_chrome_cat = torch.cat(human_head_chrome, dim=0).cpu().numpy()[:target_chromosome_length]
    if save_mouse_head:
        mouse_head_chrome_cat = torch.cat(mouse_head_chrome, dim=0).cpu().numpy()[:target_chromosome_length]
        logging.info(f'trimed length: {CHROM_LENGTH - embeddings_chrome_cat.shape[0] * 128}')
        return embeddings_chrome_cat, human_head_chrome_cat, mouse_head_chrome_cat
    else:
        logging.info(f'trimed length: {CHROM_LENGTH - embeddings_chrome_cat.shape[0] * 128}')
        return embeddings_chrome_cat, human_head_chrome_cat


def enformer_inferece(fasta_file, embeddings_save_path, human_head_save_path, enformer,
                      cuda_device, mouse_head_save_path=None):
    # to float16
    enformer = enformer.to(cuda_device)
    enformer.eval()

    vocab = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    fasta_seqs = Fasta(str(fasta_file), )
    total_length = sum([len(fasta_seqs[contig_name]) for contig_name in fasta_seqs.keys()])
    global progress_bar
    progress_bar = tqdm(total=total_length, desc=f'{sample_name}')

    for contig_name in fasta_seqs.keys():
        seq = fasta_seqs[contig_name][:].seq.upper()
        seq = torch.tensor([vocab[i] for i in seq], dtype=torch.long)
        if mouse_head_save_path is None:
            embeddings_chrome_cat, human_head_chrome_cat = infer_one_contig(seq, enformer,
                                                                            cuda_device, contig_name)
        else:
            embeddings_chrome_cat, human_head_chrome_cat, mouse_head_chrome_cat = infer_one_contig(seq, enformer,
                                                                                                   cuda_device,
                                                                                                   contig_name,
                                                                                                   save_mouse_head=True)
        with h5py.File(embeddings_save_path, 'a') as f:
            f.create_dataset(contig_name, data=embeddings_chrome_cat, dtype='float16')
        with h5py.File(human_head_save_path, 'a') as f:
            f.create_dataset(contig_name, data=human_head_chrome_cat, dtype='float16')
        if mouse_head_save_path is not None:
            with h5py.File(mouse_head_save_path, 'a') as f:
                f.create_dataset(contig_name, data=mouse_head_chrome_cat, dtype='float16')


def convert_h5_to_bigwig(hdf5_file, bigwig_save_path, fasta_file, data_dims: int = 5313):
    def get_bigwig_file_name_by_dim(bigwig_save_path, data_dim: int):
        return str(bigwig_save_path / f'{data_dim}.bw')

    bigwig_save_path.mkdir(exist_ok=True, parents=True)
    if isinstance(fasta_file, str):
        fasta_file = Fasta(fasta_file)
    chrom_sizes = [(contig_name, len(fasta_file[contig_name])) for contig_name in fasta_file.keys()]
    # sort by contig size
    chrom_sizes = sorted(chrom_sizes, key=lambda x: x[1], reverse=True)

    bigwig_output_files = [pyBigWig.open(str(get_bigwig_file_name_by_dim(bigwig_save_path, i)), 'w') for i in
                           range(data_dims)]
    for bw in bigwig_output_files:
        bw.addHeader(chrom_sizes)

    with h5py.File(hdf5_file, 'r') as f:
        for contig_name in tqdm([x[0] for x in chrom_sizes]):
            data = f[contig_name][:]
            # add entries to bigwig file
            for i in range(data_dims):
                bigwig_output_files[i].addEntries(
                    contig_name,
                    0,
                    values=data[:, i],
                    span=128,
                    step=128,
                    validate=False,
                )

    for bw in bigwig_output_files:
        bw.close()


def run_enformer_submitit_job(sample_config,
                              model_name_or_path='/Model_weights/enformer-official-rough',
                              # download from https://huggingface.co/EleutherAI/enformer-official-rough
                              ):
    fasta_file = sample_config['fasta_file']
    sample_name = sample_config['sample_name']
    embedding_hdf5_path = sample_config['embedding_hdf5_path']
    human_head_hdf5_path = sample_config['human_head_hdf5_path']
    num_of_gpu = 1

    mouse_head_hdf5_path = None
    cuda_device = f'cuda:{random.randint(0, num_of_gpu - 1)}'
    # print all received params
    print(f'''
           sample_name: {sample_name}
           cuda_device: {cuda_device}
           ''')

    # get the enformer model
    enformer = Enformer.from_pretrained(
        model_name_or_path,
    )

    # make dirs for output
    Path(embedding_hdf5_path).parent.mkdir(parents=True, exist_ok=True)
    Path(human_head_hdf5_path).parent.mkdir(parents=True, exist_ok=True)
    enformer_inferece(
        fasta_file=fasta_file,
        embeddings_save_path=embedding_hdf5_path,
        human_head_save_path=human_head_hdf5_path,
        mouse_head_save_path=mouse_head_hdf5_path,
        enformer=enformer,
        cuda_device=cuda_device,
    )


def generate_yaml_content(sample_name, hap, fasta_path):
    embedding_hdf5_path = workdir / f'hdf5/embedding_hdf5/{sample_name}.hdf5'
    human_head_hdf5_path = workdir / f'hdf5/human_head_hdf5/{sample_name}.hdf5'

    if not embedding_hdf5_path.exists():
        return {
            'embedding_hdf5_path': str(embedding_hdf5_path),
            'human_head_hdf5_path': str(human_head_hdf5_path),
            'fasta_file': str(fasta_path),
            'cuda_device': cuda_device,
            'hap': hap,
            'sample_name': sample_name,

        }


    else:
        raise ValueError(f"Embedding file already exists for sample {sample_name}")


if __name__ == '__main__':
    # get the fasta file
    start_time = time.time()
    from pathlib import Path
    from collections import OrderedDict

    yaml_content = OrderedDict()
    # Paths and constants
    input_dirs = {
        'sample_category_1': Path(
            '/path/to/sample_category_1/'),
        'sample_category_2': Path(
            '/path/to/sample_category_2/'),
    }
    yaml_dir = Path(
        './Enformer/pending_jobs')
    yaml_dir.mkdir(parents=True, exist_ok=True)
    workdir = Path('./Enformer')

    finished_jobs_dir = workdir / 'finished_jobs'
    finished_jobs_dir.mkdir(parents=True, exist_ok=True)

    cuda_device = 0

    # Specify the samples to be processed
    selected_samples = ['sample_1', 'sample_2']
    symbol_link_save_fasta_dir = Path('./input_fasta_dir')

    # Processing files in each directory
    for category, input_dir in input_dirs.items():
        for fasta_file in input_dir.glob('**/*.fasta'):
            sample_name = fasta_file.stem
            if any([sample_name.startswith(selected_sample) for selected_sample in selected_samples]):
                symbol_link_fasta_path = symbol_link_save_fasta_dir / fasta_file.name
                if not symbol_link_fasta_path.exists():
                    symbol_link_fasta_path.symlink_to(fasta_file)

                hap = fasta_file.stem.split('.')[1] if 'hap' in fasta_file.stem else "NA"
                yaml_content[f'{sample_name}'] = generate_yaml_content(sample_name, hap, symbol_link_fasta_path)

    # Set up the Submitit executor
    executor = submitit.DebugExecutor(folder=workdir / "submitit_logs")

    executor.update_parameters(
        slurm_mem='50G',
        gpus_per_node=1,
        cpus_per_task=6,
        timeout_min=60 * 24,  # Maximum execution time in minutes
        slurm_partition='a100-40g,a40-quad,a40-tmp,v100',
        slurm_qos="gpu-huge",  # SLURM QoS
    )

    # Submit jobs
    jobs = []
    with executor.batch():
        for sample_config in yaml_content.values():
            job = executor.submit(run_enformer_submitit_job, sample_config)
            jobs.append(job)

    for i, job in enumerate(jobs):
        try:
            job.result()
            print(f"Job {job.job_id} finished successfully")
            # create a done file
            success_yaml_config = list(yaml_content.values())[i]
            sample_name = success_yaml_config['sample_name']

            done_file = finished_jobs_dir / f"{sample_name}.done"
            done_file.touch()
        except Exception as e:
            print(f"Job {job.job_id} failed: {e}")
