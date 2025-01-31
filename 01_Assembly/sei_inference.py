'''
Script for generating embeddings and sequence class scores using SEI model.
SEI: https://github.com/FunctionLab/sei-framework
selene: https://github.com/FunctionLab/selene
'''

import os
import sys

import numpy as np
import torch.nn as nn
from scipy.interpolate import splev
from multiprocessing import Process, Queue
# %%
import argparse
import logging
import time
# %%
from pathlib import Path

import h5py
import torch
import torch.nn.functional as F
from pyfaidx import Fasta
from tqdm import tqdm

SEGMENT_LENGTH = 4096
BATCH_SIZE = 256 * 2

THRESHOLD = 0.0


def bs(x, df=None, knots=None, degree=3, intercept=False):
    """
    df : int
        The number of degrees of freedom to use for this spline. The
        return value will have this many columns. You must specify at least
        one of `df` and `knots`.
    knots : list(float)
        The interior knots of the spline. If unspecified, then equally
        spaced quantiles of the input data are used. You must specify at least
        one of `df` and `knots`.
    degree : int
        The degree of the piecewise polynomial. Default is 3 for cubic splines.
    intercept : bool
        If `True`, the resulting spline basis will span the intercept term
        (i.e. the constant function). If `False` (the default) then this
        will not be the case, which is useful for avoiding overspecification
        in models that include multiple spline terms and/or an intercept term.

    """

    order = degree + 1
    inner_knots = []
    if df is not None and knots is None:
        n_inner_knots = df - order + (1 - intercept)
        if n_inner_knots < 0:
            n_inner_knots = 0
            print("df was too small; have used %d"
                  % (order - (1 - intercept)))

        if n_inner_knots > 0:
            inner_knots = np.percentile(
                x, 100 * np.linspace(0, 1, n_inner_knots + 2)[1:-1])

    elif knots is not None:
        inner_knots = knots

    all_knots = np.concatenate(
        ([np.min(x), np.max(x)] * order, inner_knots))

    all_knots.sort()

    n_basis = len(all_knots) - (degree + 1)
    basis = np.empty((x.shape[0], n_basis), dtype=float)

    for i in range(n_basis):
        coefs = np.zeros((n_basis,))
        coefs[i] = 1
        basis[:, i] = splev(x, (all_knots, coefs, degree))

    if not intercept:
        basis = basis[:, 1:]
    return basis


def spline_factory(n, df, log=False):
    if log:
        dist = np.array(np.arange(n) - n / 2.0)
        dist = np.log(np.abs(dist) + 1) * (2 * (dist > 0) - 1)
        n_knots = df - 4
        knots = np.linspace(np.min(dist), np.max(dist), n_knots + 2)[1:-1]
        return torch.from_numpy(bs(
            dist, knots=knots, intercept=True)).float()
    else:
        dist = np.arange(n)
        return torch.from_numpy(bs(
            dist, df=df, intercept=True)).float()


class BSplineTransformation(nn.Module):

    def __init__(self, degrees_of_freedom, log=False, scaled=False):
        super(BSplineTransformation, self).__init__()
        self._spline_tr = None
        self._log = log
        self._scaled = scaled
        self._df = degrees_of_freedom

    def forward(self, input):
        if self._spline_tr is None:
            spatial_dim = input.size()[-1]
            self._spline_tr = spline_factory(spatial_dim, self._df, log=self._log)
            if self._scaled:
                self._spline_tr = self._spline_tr / spatial_dim
            if input.is_cuda:
                self._spline_tr = self._spline_tr.to(input.device)
            # Convert self._spline_tr to the same type as input
            self._spline_tr = self._spline_tr.type(input.dtype)

        return torch.matmul(input, self._spline_tr)


class BSplineTransformationNew(nn.Module):

    def __init__(self, degrees_of_freedom, log=False, scaled=False):
        super(BSplineTransformationNew, self).__init__()
        self._spline_tr = None
        self._log = log
        self._scaled = scaled
        self._df = degrees_of_freedom
        _spline_tr = spline_factory(256, self._df, log=self._log)
        # convert to nn.Parameter
        self._spline_tr = nn.Parameter(_spline_tr)

    def forward(self, input):
        return torch.matmul(input, self._spline_tr)


class BSplineConv1D(nn.Module):

    def __init__(self, in_channels, out_channels, kernel_size, degrees_of_freedom, stride=1,
                 padding=0, dilation=1, groups=1, bias=True, log=False, scaled=True):
        super(BSplineConv1D, self).__init__()
        self._df = degrees_of_freedom
        self._log = log
        self._scaled = scaled

        self.spline = nn.Conv1d(1, degrees_of_freedom, kernel_size, stride, padding, dilation,
                                bias=False)
        self.spline.weight = spline_factory(kernel_size, self._df, log=log).view(self._df, 1, kernel_size)
        if scaled:
            self.spline.weight = self.spline.weight / kernel_size
        self.spline.weight = nn.Parameter(self.spline.weight)
        self.spline.weight.requires_grad = False
        self.conv1d = nn.Conv1d(in_channels * degrees_of_freedom, out_channels, 1,
                                groups=groups, bias=bias)

    def forward(self, input):
        batch_size, n_channels, length = input.size()
        spline_out = self.spline(input.view(batch_size * n_channels, 1, length))
        conv1d_out = self.conv1d(spline_out.view(batch_size, n_channels * self._df, length))
        return conv1d_out


class Sei(nn.Module):
    def __init__(self, sequence_length=4096, n_genomic_features=21907):
        """
        Parameters
        ----------
        sequence_length : int
        n_genomic_features : int
        """
        super(Sei, self).__init__()

        self.lconv1 = nn.Sequential(
            nn.Conv1d(4, 480, kernel_size=9, padding=4),
            nn.Conv1d(480, 480, kernel_size=9, padding=4))

        self.conv1 = nn.Sequential(
            nn.Conv1d(480, 480, kernel_size=9, padding=4),
            nn.ReLU(inplace=True),
            nn.Conv1d(480, 480, kernel_size=9, padding=4),
            nn.ReLU(inplace=True))

        self.lconv2 = nn.Sequential(
            nn.MaxPool1d(kernel_size=4, stride=4),
            nn.Dropout(p=0.2),
            nn.Conv1d(480, 640, kernel_size=9, padding=4),
            nn.Conv1d(640, 640, kernel_size=9, padding=4))

        self.conv2 = nn.Sequential(
            nn.Dropout(p=0.2),
            nn.Conv1d(640, 640, kernel_size=9, padding=4),
            nn.ReLU(inplace=True),
            nn.Conv1d(640, 640, kernel_size=9, padding=4),
            nn.ReLU(inplace=True))

        self.lconv3 = nn.Sequential(
            nn.MaxPool1d(kernel_size=4, stride=4),
            nn.Dropout(p=0.2),
            nn.Conv1d(640, 960, kernel_size=9, padding=4),
            nn.Conv1d(960, 960, kernel_size=9, padding=4))

        self.conv3 = nn.Sequential(
            nn.Dropout(p=0.2),
            nn.Conv1d(960, 960, kernel_size=9, padding=4),
            nn.ReLU(inplace=True),
            nn.Conv1d(960, 960, kernel_size=9, padding=4),
            nn.ReLU(inplace=True))

        self.dconv1 = nn.Sequential(
            nn.Dropout(p=0.10),
            nn.Conv1d(960, 960, kernel_size=5, dilation=2, padding=4),
            nn.ReLU(inplace=True))
        self.dconv2 = nn.Sequential(
            nn.Dropout(p=0.10),
            nn.Conv1d(960, 960, kernel_size=5, dilation=4, padding=8),
            nn.ReLU(inplace=True))
        self.dconv3 = nn.Sequential(
            nn.Dropout(p=0.10),
            nn.Conv1d(960, 960, kernel_size=5, dilation=8, padding=16),
            nn.ReLU(inplace=True))
        self.dconv4 = nn.Sequential(
            nn.Dropout(p=0.10),
            nn.Conv1d(960, 960, kernel_size=5, dilation=16, padding=32),
            nn.ReLU(inplace=True))
        self.dconv5 = nn.Sequential(
            nn.Dropout(p=0.10),
            nn.Conv1d(960, 960, kernel_size=5, dilation=25, padding=50),
            nn.ReLU(inplace=True))

        self._spline_df = int(128 / 8)
        self.spline_tr = nn.Sequential(
            nn.Dropout(p=0.5),
            # BSplineTransformation(self._spline_df, scaled=False)
            BSplineTransformationNew(self._spline_df, scaled=False)
        )

        self.classifier = nn.Sequential(
            nn.Linear(960 * self._spline_df, n_genomic_features),
            nn.ReLU(inplace=True),
            nn.Linear(n_genomic_features, n_genomic_features),
            nn.Sigmoid())

    def forward(self, x):
        """Forward propagation of a batch.
        """
        lout1 = self.lconv1(x)
        out1 = self.conv1(lout1)

        lout2 = self.lconv2(out1 + lout1)
        out2 = self.conv2(lout2)

        lout3 = self.lconv3(out2 + lout2)
        out3 = self.conv3(lout3)

        dconv_out1 = self.dconv1(out3 + lout3)
        cat_out1 = out3 + dconv_out1
        dconv_out2 = self.dconv2(cat_out1)
        cat_out2 = cat_out1 + dconv_out2
        dconv_out3 = self.dconv3(cat_out2)
        cat_out3 = cat_out2 + dconv_out3
        dconv_out4 = self.dconv4(cat_out3)
        cat_out4 = cat_out3 + dconv_out4
        dconv_out5 = self.dconv5(cat_out4)
        out = cat_out4 + dconv_out5

        spline_out = self.spline_tr(out)
        reshape_out = spline_out.view(spline_out.size(0), 960 * self._spline_df)
        predict = self.classifier(reshape_out)
        return predict


def criterion():
    """
    The criterion the model aims to minimize.
    """
    return nn.BCELoss()


def get_optimizer(lr):
    """
    The optimizer and the parameters with which to initialize the optimizer.
    At a later time, we initialize the optimizer by also passing in the model
    parameters (`model.parameters()`). We cannot initialize the optimizer
    until the model has been initialized.
    """
    return (torch.optim.SGD,
            {"lr": lr, "weight_decay": 1e-7, "momentum": 0.9})


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


def seq_indices_to_one_hot(t, padding=-1):
    is_padding = t == padding
    t = t.clamp(min=0)
    one_hot = F.one_hot(t, num_classes=5)
    out = one_hot[..., :4].float()
    out = out.masked_fill(is_padding[..., None], 0.25)
    return out


def convert_numpy_seq_to_model_input(sequences, cuda_device):
    sequences = sequences.transpose(1, 2)
    sequences = sequences.to(cuda_device).half()
    return sequences


def infer_one_contig(seq, enformer, cuda_device, contig_name: str, embeddings_save_path):

    max_len = len(seq) // TARGET_LENGTH
    with h5py.File(embeddings_save_path, 'a') as f:
        f.create_dataset(contig_name, shape=(max_len, 21907), dtype='float16', compression='lzf', chunks=(1, 21907))
        f.create_dataset(f'{contig_name}_seq_class_scores', shape=(max_len, len(clustervfeat))
                         , dtype='float16',
                         compression='lzf', chunks=(1, len(clustervfeat)))

    save_queue = Queue()
    save_process = Process(target=save_worker, args=(save_queue, embeddings_save_path))
    save_process.start()

    seq = seq.type(torch.LongTensor).to(cuda_device)
    CHROM_LENGTH = len(seq)
    target_start = 0
    target_end = target_start + TARGET_LENGTH
    extra_seq = SEGMENT_LENGTH - TARGET_LENGTH
    extra_left_seq = extra_seq // 2
    extra_right_seq = extra_seq - extra_left_seq

    # embeddings_chrome = []
    # seq_class_chrome = []
    sequences = torch.zeros((BATCH_SIZE, SEGMENT_LENGTH, 4), device=cuda_device)  # Allocate on GPU directly
    batch_counter = 0

    with torch.no_grad():
        save_start = 0
        while target_end <= CHROM_LENGTH:
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
            input_seq = seq_indices_to_one_hot(input_seq)  # Optimize this function for GPU if possible

            sequences[batch_counter, :, :] = input_seq
            batch_counter += 1
            if batch_counter == BATCH_SIZE:
                preds = enformer(sequences.transpose(1, 2).half())
                embeddings_batch = preds.cpu().half()
                seq_class_batch = sc_projection(preds)
                save_end = save_start + BATCH_SIZE
                logger.debug(f"Saving {contig_name} from {save_start} to {save_end}")
                save_queue.put((contig_name, save_start, save_end, embeddings_batch, seq_class_batch))
                # embeddings_chrome.append(preds.cpu().half())
                # seq_class_chrome.append(sc_projection(preds))

                batch_counter = 0  # Reset counter after processing a batch
                progress_bar.update(BATCH_SIZE * TARGET_LENGTH)
                save_start += BATCH_SIZE

            target_start = target_end
            target_end = target_start + TARGET_LENGTH


        if batch_counter > 0:  # Process the remaining sequences
            preds = enformer(sequences[:batch_counter].transpose(1, 2).half())
            embeddings_batch = preds.cpu().half()
            seq_class_batch = sc_projection(preds)
            save_end = save_start + batch_counter
            logger.debug(f"Saving {contig_name} from {save_start} to {save_end}")
            save_queue.put((contig_name, save_start, save_end, embeddings_batch, seq_class_batch))
            progress_bar.update(batch_counter * TARGET_LENGTH)

    save_queue.put(None)
    save_process.join()
    # embeddings_chrome_cat = torch.cat(embeddings_chrome, dim=0).cpu().numpy()
    # seq_class_chrome_cat = torch.cat(seq_class_chrome, dim=0).cpu().numpy()
    # return embeddings_chrome_cat, seq_class_chrome_cat


def sc_projection(chromatin_profile_preds):
    dot_product = torch.mm(chromatin_profile_preds, clustervfeat.T)
    seq_class_scores = dot_product / clustervfeat_norm
    return seq_class_scores.cpu().half()


def save_to_hdf5(embeddings_chrome_cat, seq_class_chrome_cat, embeddings_save_path, contig_name):
    with h5py.File(embeddings_save_path, 'a') as f:
        embeddings_chrome_cat[embeddings_chrome_cat < THRESHOLD] = 0
        f.create_dataset(contig_name, data=embeddings_chrome_cat, dtype='float16', compression='lzf')
        f.create_dataset(f'{contig_name}_seq_class_scores', data=seq_class_chrome_cat, dtype='float16',
                         compression='lzf')


def save_worker(save_queue, embeddings_save_path):
    f= h5py.File(embeddings_save_path, 'a')
    while True:
        task = save_queue.get()
        if task is None:  # Termination signal
            break
        contig_name, start, stop, embeddings_chrome_batch, seq_class_chrome_batch = task
        # with h5py.File(embeddings_save_path, 'a') as f:
            # Apply threshold
        embeddings_chrome_batch[embeddings_chrome_batch < THRESHOLD] = 0

        f[contig_name][start:stop] = embeddings_chrome_batch
        f[f'{contig_name}_seq_class_scores'][start:stop] = seq_class_chrome_batch
    f.close()

def enformer_inferece(fasta_file, embeddings_save_path, enformer,
                      cuda_device, ):
    # to float16
    enformer = enformer.to(cuda_device).half()
    enformer.eval()

    vocab = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}
    fasta_seqs = Fasta(str(fasta_file), )
    total_length = sum([len(fasta_seqs[contig_name]) for contig_name in fasta_seqs.keys()])

    global progress_bar
    start_time = time.time()
    logger.info(f"Processing {sample_name} {haplotype}")
    progress_bar = tqdm(total=total_length, desc=f'{sample_name} {haplotype}')

    # save_queue = Queue()
    # save_process = Process(target=save_worker, args=(save_queue, embeddings_save_path))
    # save_process.start()

    # processes = []
    for contig_name in fasta_seqs.keys():
        seq = fasta_seqs[contig_name][:].seq.upper()
        logger.info(f"Processing {contig_name} with length {len(seq)}")
        seq = torch.tensor([vocab.get(i,4) for i in seq], dtype=torch.long)
        infer_one_contig(seq, enformer, cuda_device, contig_name, embeddings_save_path)
    end_time = time.time()
    logger.info(f"--- Using H:M:S --- {time.strftime('%H:%M:%S', time.gmtime(end_time - start_time))}")

def read_config(config_path):
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
    return config


def process_samples(config):
    global TARGET_LENGTH, sample_name, haplotype
    for sample_name_hap, details in config.items():
        fasta_file = details['fasta_path']
        haplotype = details['hap']
        cuda_device = f'cuda:{details["cuda_device"]}'
        embedding_hdf5_path = Path(details['embedding_hdf5_path'])
        TARGET_LENGTH = int(details['target_length'])
        sample_name = details['sample_name']

        # Process each sample here
        print(f"Processing {sample_name} with fasta file {fasta_file} and haplotype {haplotype}")

        print(f'''
        sample_name: {sample_name}
        haplotype: {haplotype}
        cuda_device: {cuda_device}
        TARGET_LENGTH: {TARGET_LENGTH}
        ''')
        # make dirs for output
        Path(embedding_hdf5_path).parent.mkdir(parents=True, exist_ok=True)
        if embedding_hdf5_path.exists():
            print(f"The file {embedding_hdf5_path} already exists. Skipping...")
            print(f"For re-generating the embeddings, please delete the file and re-run the script.")
        else:
            enformer_inferece(
                fasta_file=fasta_file,
                embeddings_save_path=embedding_hdf5_path,
                enformer=model,
                cuda_device=cuda_device,
            )

# %%
if __name__ == '__main__':
    # get the fasta file
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    logger.addHandler(logging.StreamHandler())
    logger.handlers[0].setFormatter(formatter)

    start_time = time.time()
    sei_dir = "sei-framework/model" # https://github.com/FunctionLab/sei-framework
    clustervfeat = np.load(os.path.join(sei_dir, 'projvec_targets.npy'))
    clustervfeat = torch.from_numpy(clustervfeat).cuda().half()
    clustervfeat_norm = torch.linalg.norm(clustervfeat, dim=1)

    model_name_or_path = 'sei-framework/model/sei.pth'

    model_configs = {
        'path': 'sei-framework/model/sei.py',
        'class': 'Sei',
        'class_args': {},
        # 'non_strand_specific': 'mean'
    }
    model = Sei(**model_configs['class_args'])
    if "non_strand_specific" in model_configs:
        print("Using non-strand specific model.")
        from selene_sdk.utils import NonStrandSpecific # https://github.com/FunctionLab/selene
        model = NonStrandSpecific(model, mode=model_configs["non_strand_specific"])
        model_paramter_dict = (torch.load(model_name_or_path))
        new_model_paramter_dict = {}
        for k, v in model_paramter_dict.items():
            if k.startswith('module'):
                new_model_paramter_dict[k[7:]] = v
            else:
                new_model_paramter_dict[k] = v
        model.load_state_dict(new_model_paramter_dict, strict=False)

    else:
        print("Using strand-specific model.")
        model_paramter_dict = (torch.load(model_name_or_path))
        new_model_paramter_dict = {}
        for k, v in model_paramter_dict.items():
            if k.startswith('module'):
                new_model_paramter_dict[k[len('module.model.'):]] = v
            else:
                new_model_paramter_dict[k] = v
        model.load_state_dict(new_model_paramter_dict, strict=False)

    import yaml
    config_file_path = sys.argv[1]
    config = read_config(config_file_path)
    process_samples(config)
