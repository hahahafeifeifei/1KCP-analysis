#!/usr/bin/env python3
import sys
import gzip
import numpy as np

def openfile(filename):
    if filename == "-":
        return sys.stdin
    elif filename.endswith('.gz'):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")

#Filter the outlier allele and genrate the TR length and TR variation dosage for TR site
#python3 vamos_tr_dosage.py input.vcf output.length output.motif

input_file = sys.argv[1]
output_len_file = open(sys.argv[2], "w")
output_motif_file = open(sys.argv[3], "w")
outlier_sd_range = 3
min_allele = 3

for line in openfile(input_file):
    if line[0:2] == "##":
        output_len_file.write(line)
        output_motif_file.write(line)
    elif line[0] == "#":
        output_len_file.write('##FORMAT=<ID=CN,Number=1,Type=String,Description="The TR copy number">' + '\n')
        output_len_file.write('##FORMAT=<ID=DS,Number=1,Type=Float,Description="TR normalized dosage format">' + '\n')
        output_len_file.write('##INFO=<ID=LEN,Number=A,Type=Integer,Description="TR length">' + '\n')
        output_len_file.write(line)
        output_motif_file.write('##FORMAT=<ID=CN,Number=1,Type=String,Description="The TR copy number">' + '\n')
        output_motif_file.write('##FORMAT=<ID=DS,Number=1,Type=Float,Description="TR normalized dosage format">' + '\n')
        output_motif_file.write('##INFO=<ID=LEN,Number=A,Type=Integer,Description="TR length">' + '\n')
        output_motif_file.write('##INFO=<ID=MOTIF,Number=1,Type=String,Description="TR motif">' + '\n')
        output_motif_file.write(line)
    else:
        #initialize
        line_info = line.strip().split()
        info = line_info[7]
        info_dict = {info.split("=")[0] : info.split("=")[1] for info in line_info[7].split(";")}
        ru_list = info_dict["RU"].split(",")
        allele_len_list = []
        allele_ru_list = []
        for ru in ru_list:
            allele_ru_list.append([])
        #parse allele
        allele_list = info_dict["ALTANNO"].split(",")
        for allele in allele_list:
            allele_len = len(allele.split("-"))
            allele_ru = [0] * len(ru_list)
            for unit in allele.split("-"):
                allele_ru[int(unit)] += 1
            for i, ru_len in enumerate(allele_ru):
                allele_ru_list[i].append(ru_len)
            allele_len_list.append(allele_len)
        #parse gt
        gt_list = line_info[9:]
        gt_len_list = np.array([])
        gt_ru_list = []
        for i in range(len(ru_list)):
            gt_ru_list.append(np.array([]))
        for gt in gt_list:
            for gt_hap in gt.replace("|","/").split("/"):
                if gt_hap != ".":
                    gt_len_list = np.append(gt_len_list, allele_len_list[int(gt_hap) - 1])
                    for i in range(len(ru_list)):
                        gt_ru_list[i] = np.append(gt_ru_list[i],allele_ru_list[i][int(gt_hap) - 1])
                else:
                    gt_len_list = np.append(gt_len_list, np.nan)
                    for i in range(len(ru_list)):
                        gt_ru_list[i] = np.append(gt_ru_list[i], np.nan)
        #len statistics
        len_mean = np.nanmean(gt_len_list)
        len_sd = np.nanstd(gt_len_list)
        len_range = (gt_len_list > len_mean + outlier_sd_range * len_sd) | (gt_len_list < len_mean - outlier_sd_range * len_sd)
        gt_len_list[len_range] = np.nan
        #len normalize
        len_norm_list = gt_len_list.copy()
        len_norm_max = np.nanmax(len_norm_list)
        len_norm_min = np.nanmin(len_norm_list)
        if len_norm_max == len_norm_min:
            len_norm_list = np.array([0.0] * len(len_norm_list))
            len_norm_list[np.isnan(gt_len_list)] = np.nan
        else:
            len_norm_list = (len_norm_list - len_norm_min) / (len_norm_max - len_norm_min)
        #output len
        len_out_list = []
        out_gt_list = []
        for j in range(len(gt_len_list)):
            if j % 2 == 0:
                if np.isnan(len_norm_list[j]) or np.isnan(len_norm_list[j+1]):
                    gt = ".|."
                    cn = ".|."
                    ds = "."
                else:
                    if gt_len_list[j] not in len_out_list:
                        len_out_list.append(gt_len_list[j])
                    if gt_len_list[j+1] not in len_out_list:
                        len_out_list.append(gt_len_list[j+1])
                    gt = str(len_out_list.index(gt_len_list[j])+1) + "|" + str(len_out_list.index(gt_len_list[j+1])+1)
                    cn = str(int(gt_len_list[j])) + "|" + str(int(gt_len_list[j+1]))
                    ds = str("{:.2f}".format(len_norm_list[j] + len_norm_list[j+1]))
                out_gt_list.append(gt + ":" + cn + ":" + ds)
        if len(len_out_list) >= min_allele:
            out_line_info = line_info[:9].copy()
            out_line_info[7] += ";LEN=" + ",".join([str(int(len_out)) for len_out in len_out_list])
            out_line_info[8] = "GT:CN:DS"
            output_len_file.write("\t".join(out_line_info + out_gt_list) + "\n")


        for i in range(len(ru_list)):
            ru_motif = ru_list[i]
            #len statistics
            ru_mean = np.nanmean(gt_ru_list[i])
            ru_sd = np.nanstd(gt_ru_list[i])
            ru_range = (gt_ru_list[i] > ru_mean + outlier_sd_range * ru_sd) | (gt_ru_list[i] < ru_mean - outlier_sd_range * ru_sd)
            gt_ru_list[i][ru_range] = np.nan
            #len normalize
            ru_norm_list = gt_ru_list[i].copy()
            ru_norm_max = np.nanmax(ru_norm_list)
            ru_norm_min = np.nanmin(ru_norm_list)
            if ru_norm_max == ru_norm_min:
                ru_norm_list = np.array([0.0] * len(ru_norm_list))
                ru_norm_list[np.isnan(gt_ru_list[i])] = np.nan
            else:
                ru_norm_list = (ru_norm_list - ru_norm_min) / (ru_norm_max - ru_norm_min)
            #output len
            ru_out_list = []
            out_gt_list = []
            for j in range(len(gt_ru_list[i])):
                if j % 2 == 0:
                    if np.isnan(ru_norm_list[j]) or np.isnan(ru_norm_list[j+1]):
                        gt = ".|."
                        cn = ".|."
                        ds = "."
                    else:
                        if gt_ru_list[i][j] not in ru_out_list:
                            ru_out_list.append(gt_ru_list[i][j])
                        if gt_ru_list[i][j+1] not in ru_out_list:
                            ru_out_list.append(gt_ru_list[i][j+1])
                        gt = str(ru_out_list.index(gt_ru_list[i][j])+1) + "|" + str(ru_out_list.index(gt_ru_list[i][j+1])+1)
                        cn = str(int(gt_ru_list[i][j])) + "|" + str(int(gt_ru_list[i][j+1]))
                        ds = str("{:.2f}".format(ru_norm_list[j] + ru_norm_list[j+1]))
                    out_gt_list.append(gt + ":" + cn + ":" + ds)
            if len(ru_out_list) >= min_allele:
                out_line_info = line_info[:9].copy()
                out_line_info[7] += ";MOTIF=" + ru_list[i] + ";LEN=" + ",".join([str(int(ru_out)) for ru_out in ru_out_list])
                out_line_info[8] = "GT:CN:DS"
                output_motif_file.write("\t".join(out_line_info + out_gt_list) + "\n")
output_len_file.close()
output_motif_file.close()

