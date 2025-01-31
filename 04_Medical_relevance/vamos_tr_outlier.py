#!/usr/bin/env python3
import sys
import gzip
import math
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn.neighbors import NearestNeighbors
from scipy import stats

#Determine the normal and outlier length threshold of TR sites
#python3 vamos_tr_outlier.py input.vcf output.vcf

def openfile(filename):
    if filename == "-":
        return sys.stdin
    elif filename.endswith('.gz'):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")

input_file = sys.argv[1]
output_len_file = open(sys.argv[2], "w")

for line in openfile(input_file):
    if line[0:2] == "##":
        output_len_file.write(line)
    elif line[0] == "#":
        output_len_file.write('##FORMAT=<ID=CN,Number=1,Type=String,Description="The TR copy number">' + '\n')
        output_len_file.write('##INFO=<ID=LEN,Number=A,Type=Integer,Description="TR length">' + '\n')
        output_len_file.write('##INFO=<ID=MIN_PCT,Number=1,Type=Integer,Description="Repeat size of 99.5th percentile">' + '\n')
        output_len_file.write('##INFO=<ID=MAX_PCT,Number=1,Type=Integer,Description="Repeat size of 0.5th percentile">' + '\n')
        output_len_file.write('##INFO=<ID=MAX_NORMAL,Number=1,Type=Integer,Description="Maximum repeat size not belong to DBSCAN outlier">' + '\n')
        output_len_file.write(line)
    else:
        #initialize
        line_info = line.strip().split()
        info = line_info[7]
        info_dict = {info.split("=")[0] : info.split("=")[1] for info in line_info[7].split(";")}
        ru_list = info_dict["RU"].split(",")
        allele_len_list = []
        #parse allele
        allele_list = info_dict["ALTANNO"].split(",")
        for allele in allele_list:
            allele_len = len(allele.split("-"))
            allele_len_list.append(allele_len)
        #parse gt
        gt_list = line_info[9:]
        gt_len_list = np.array([])
        for gt in gt_list:
            for gt_hap in gt.replace("|","/").split("/"):
                if gt_hap != ".":
                    gt_len_list = np.append(gt_len_list, allele_len_list[int(gt_hap) - 1])
                else:
                    gt_len_list = np.append(gt_len_list, np.nan)
        #len statistics
        gt_clean_len_list = gt_len_list[~np.isnan(gt_len_list)]
        max_len = max(gt_clean_len_list)
        eps_value = 2 * stats.mode(gt_clean_len_list, keepdims=True).mode[0]
        sample_value = math.floor(np.log2(len(gt_clean_len_list)))
        dbscan = DBSCAN(eps=eps_value, min_samples=sample_value)
        dbscan.fit(gt_clean_len_list.reshape(-1,1))
        core_sample_len = gt_clean_len_list[dbscan.core_sample_indices_]
        max_core_sample_len = max(core_sample_len)
        nbrs = NearestNeighbors(radius=eps_value)
        nbrs.fit(core_sample_len.reshape(-1,1))
        distances, indices = nbrs.radius_neighbors(np.array(range(0,10000)).reshape(-1, 1))
        #the boundary with maximum cluster
        non_outlier_list = []
        for i,j in enumerate(indices):
            if len(j) != 0 and i>=max_core_sample_len:
                non_outlier_list.append(i)
        if non_outlier_list == []:
            max_normal = "."
        else:
            max_normal = str(max(non_outlier_list))
        max_pct = str(int(stats.scoreatpercentile(gt_clean_len_list, interpolation_method="lower",per=99.5)))
        min_pct = str(int(stats.scoreatpercentile(gt_clean_len_list, interpolation_method="higher",per=0.5)))

        #output len
        len_out_list = []
        out_gt_list = []
        for j in range(len(gt_len_list)):
            if j % 2 == 0:
                if (not np.isnan(gt_len_list[j])) and (gt_len_list[j] not in len_out_list):
                    len_out_list.append(gt_len_list[j])
                if (not np.isnan(gt_len_list[j+1])) and (gt_len_list[j] not in len_out_list):
                    len_out_list.append(gt_len_list[j+1])
                gt_list = []
                cn_list = []
                if np.isnan(gt_len_list[j]):
                    gt_list.append(".")
                    cn_list.append(".")
                else:
                    gt_list.append(str(len_out_list.index(gt_len_list[j])+1))
                    cn_list.append(str(int(gt_len_list[j])))
                if np.isnan(gt_len_list[j+1]):
                    gt_list.append(".")
                    cn_list.append(".")
                else:
                    gt_list.append(str(len_out_list.index(gt_len_list[j+1])+1))
                    cn_list.append(str(int(gt_len_list[j+1])))
                gt = gt_list[0] + "|" + gt_list[1]
                cn = cn_list[0] + "|" + cn_list[1] 
                out_gt_list.append(gt + ":" + cn)
        out_line_info = line_info[:9].copy()
        out_line_info[7] += ";LEN=" + ",".join([str(int(len_out)) for len_out in len_out_list])
        out_line_info[7] +=  ";MIN_PCT=" + min_pct
        out_line_info[7] +=  ";MAX_PCT=" + max_pct
        out_line_info[7] +=  ";MAX_NORMAL=" + max_normal
        out_line_info[8] = "GT:CN"
        output_len_file.write("\t".join(out_line_info + out_gt_list) + "\n")
output_len_file.close()
