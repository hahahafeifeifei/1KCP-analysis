#!/usr/bin/env python3
import sys
sys.path.append("scripts")
import gzip
import re
import numpy as np
import py_kmc_api as pka

def openfile(filename):
    if filename == "-":
        return sys.stdin
    elif filename.endswith('.gz'):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")

def rev_edge_info(edge):
    direction_dict = {"+" : "-", "-" : "+"}
    regex = r'(\d+)([-+])(\d+)([-+])'
    match = re.match(regex, edge)
    node1, direction1, node2, direction2 = match.groups()    
    rev_edge = str(node2) + direction_dict[str(direction2)] + str(node1) + direction_dict[str(direction1)]
    return rev_edge

def rev_comp(seq):
    rev_dict = {'A':'T', 'G':'C', 'C':'G', 'T':'A', 'N':'N'}
    rev_seq = ''.join([rev_dict[base] for base in seq[::-1]])
    return rev_seq

def mismatch_validate_seq(seq, kmer_data_base, min_supp, klen, conf_range):
    lower_conf = max(1, conf_range[0])
    upper_conf = conf_range[1]
    bases = ['A', 'C', 'G', 'T']
    res = pka.CountVec()
    kmer_data_base.GetCountersForRead(seq, res)
    supp_len = len([x for x in res.value if x >= lower_conf and x <= upper_conf])
    conf_len = len([x for x in res.value if x <= upper_conf])
    if conf_len == 0 or supp_len/conf_len >= min_supp:
        supp_flag = True
    else:
        supp_flag = False
    return supp_flag

gfa_file = sys.argv[1]
kmer_file = sys.argv[2]
ref_sample_list = sys.argv[3].split(",")
min_supp = float(sys.argv[6])

####Load the kmer file list
sample_kmer_dict = {}
for line in openfile(kmer_file):
    sample_kmer_dict[line.strip().split()[0]] = line.strip().split()[1]

####Load the node and edge dit
node_seq_dict = {}
edge_dict = {}
for line in openfile(gfa_file):
    if line[0] == "S":
        node = line.strip().split()[1]
        seq = line.strip().split()[2].upper()
        node_seq_dict[node] = seq
    if line[0] == "L":
        edge_info = line.strip().split()
        edge = edge_info[1] + edge_info[2] + edge_info[3] + edge_info[4]
        edge_dict[edge] = True

####Load the reference node and edge
ref_node_dict = {}
ref_edge_dict = {}
for line in openfile(gfa_file):
    if line[0] == "W":
        sample = line.split()[1]
        if sample in ref_sample_list:
            path_info = line.strip().split()[6]
            path_node = re.split('<|>', path_info)[1:]
            path_direct = re.split('\d+', path_info)[:-1]
            pre_node_info = "null"
            for i in range(len(path_node)):
                node = path_node[i]
                direct = path_direct[i]
                node_info = node + direct.replace(">","+").replace("<","-")
                ref_node_dict[node] = True
                if pre_node_info != "null":
                    edge = pre_node_info + node_info
                    rev_edge = rev_edge_info(edge)
                    if edge in edge_dict:
                        ref_edge_dict[edge] = True
                    elif rev_edge in edge_dict:
                        ref_edge_dict[rev_edge] = True
                pre_node_info = node_info

####Non ref node and edge and their position
pre_sample = "null"
nonref_node_stat = {}
nonref_edge_stat = {}
sample_nonref_node_stat = {}
sample_nonref_edge_stat = {}
for line in openfile(gfa_file):
    if line[0] == "W":
        sample = line.split()[1]
        if sample not in sample_kmer_dict:
            continue
        if sample != pre_sample:
            ####Load kmer database
            kmer_data_base = pka.KMCFile()
            kmer_data_base.OpenForRA(sample_kmer_dict[sample])
            kmer_size = kmer_data_base.KmerLength()
            ####kmer confident range
            range_file = openfile(sample_kmer_dict[sample] + ".conf")
            conf_range = [int(float(cov)) for cov in range_file.read().strip().split()]
            ####Update the nonref statistics
            if pre_sample != "null":
                for node, error in sample_nonref_node_stat.items():
                    nonref_node_stat[node] = nonref_node_stat.get(node, [0, 0])
                    nonref_node_stat[node][0] += 1
                    if error:
                        nonref_node_stat[node][1] += 1
                for edge, error in sample_nonref_edge_stat.items():
                    nonref_edge_stat[edge] = nonref_edge_stat.get(edge, [0, 0])
                    nonref_edge_stat[edge][0] += 1
                    if error:
                        nonref_edge_stat[edge][1] += 1
                sample_nonref_node_stat = {}
                sample_nonref_edge_stat = {}
            pre_sample = sample
        ####Kmer set extesion: kmer of node (k-1 bp extention + node) and edge (2bp node connection + k-2 bp extention)
        path_seq = ""
        path_size = int(line.split()[5]) - int(line.split()[4])
        path_info = line.strip().split()[6]
        path_node = re.split('<|>', path_info)[1:]
        path_direct = re.split('\d+', path_info)[:-1]
        pre_node_pos = 0
        node_pos = 0
        pre_node_info = "null"
        nonref_node_region = {}
        nonref_edge_region = {}
        for i in range(len(path_node)):
            node = path_node[i]
            direct = path_direct[i]
            node_seq = node_seq_dict[node] if direct == ">" else rev_comp(node_seq_dict[node])
            path_seq += node_seq
            node_pos += len(node_seq)
            node_info = node + direct.replace(">","+").replace("<","-")
            if node not in ref_node_dict:
                nonref_node_region[node] = [max(pre_node_pos - kmer_size + 1, 0), min(node_pos + kmer_size - 1, path_size)]
            if pre_node_info != "null":
                pre_node = pre_node_info[:-1]
                edge = pre_node_info + node_info
                rev_edge = rev_edge_info(edge)
                if edge in edge_dict:
                    edge = edge
                elif rev_edge in edge_dict:
                    edge = rev_edge
                if edge not in ref_edge_dict:
                    nonref_edge_region[edge] = [max(pre_node_pos - kmer_size + 1, 0), min(pre_node_pos + kmer_size - 1, path_size)]
            pre_node_pos = node_pos
            pre_node_info = node_info
        ####Search the kmer count, add and retain the node and edge with all error k-mer to the error dict, delete them when they are correct in another path position
        for node, region in nonref_node_region.items():
            seq = path_seq[region[0]:region[1]]
            if len(seq) < kmer_size:
                sample_nonref_node_stat[node] = True
            else:
                if mismatch_validate_seq(seq, kmer_data_base, min_supp, kmer_size, conf_range):
                    sample_nonref_node_stat[node] = True
                elif node not in sample_nonref_node_stat:
                    sample_nonref_node_stat[node] = False
        for edge, region in nonref_edge_region.items():
            seq = path_seq[region[0]:region[1]]
            if len(seq) < kmer_size:
                sample_nonref_edge_stat[edge] = True
            else:
                if mismatch_validate_seq(seq, kmer_data_base, min_supp, kmer_size, conf_range):
                    sample_nonref_edge_stat[edge] = True
                elif edge not in sample_nonref_edge_stat:
                    sample_nonref_edge_stat[edge] = False   

for node, valid in sample_nonref_node_stat.items():
    nonref_node_stat[node] = nonref_node_stat.get(node, [0, 0])
    nonref_node_stat[node][0] += 1
    if valid:
        nonref_node_stat[node][1] += 1
for edge, valid in sample_nonref_edge_stat.items():
    nonref_edge_stat[edge] = nonref_edge_stat.get(edge, [0, 0])
    nonref_edge_stat[edge][0] += 1
    if valid:
        nonref_edge_stat[edge][1] += 1

####output error node and edge
node_out_file = open(sys.argv[4],"w")
for node, stat in nonref_node_stat.items():
    node_out_file.write(node + "\t" + str(stat[1]/stat[0]) + "\n")
    #if stat[1]/stat[0] < min_ratio:
        #fo.write(node + "\n")
edge_out_file = open(sys.argv[5],"w")
for edge, stat in nonref_edge_stat.items():
    edge_out_file.write(edge + "\t" + str(stat[1]/stat[0]) + "\n")
    #if stat[1]/stat[0] < min_ratio:
        #fo.write(edge + "\n")
