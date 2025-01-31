#!/usr/bin/env python3
import sys
import gzip

#Normalize the graph decomposed vcf based on sample list
#python3 graph_vcf_norm.py input.vcf sample.list > output.vcf

def openfile(filename):
    if filename == "-":
        return sys.stdin
    elif filename.endswith('.gz'):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")

input_file = sys.argv[1]
sample_file = sys.argv[2]

final_sample_list = []
final_sample_sex_list = []
for line in openfile(sample_file):
    sample = line.strip().split()[0]
    sex = line.strip().split()[1]
    final_sample_list.append(sample)
    final_sample_sex_list.append(sex)

line_number = 0
sample_index = []
for line in openfile(input_file):
    if line.startswith("##"):
        print(line.strip())
    elif line.startswith("#"):
        line_number = len(line.strip().split())
        sample_list = line.strip().split()[9:]
        for sample in sample_list:
            if sample in final_sample_list:
                sample_index.append(final_sample_list.index(sample))
            else:
                sample_index.append("NA")
        print('\t'.join(line.strip().split()[0:9] + final_sample_list))
    else:
        if len(line.strip().split()) != line_number:
            continue
        if len(line)>=10000000000:
            continue
        chr = line.strip().split()[0]
        sample_gt_list = line.strip().split()[9:]
        sample_gt_add_list = []
        final_sample_gt_list = [".|."] * len(final_sample_list)
        for i, sample_gt in enumerate(sample_gt_list):
            if sample_index[i] == "NA":
                continue
            else:
                final_sample_gt_list[sample_index[i]] = sample_gt
        for i, sample_gt in enumerate(final_sample_gt_list):
            hap_gt_list = sample_gt.replace("/","|").split("|")
            if len(hap_gt_list) == 1:
                hap_gt_list = hap_gt_list + hap_gt_list
            if final_sample_sex_list[i] == "male" and ("chrX" in chr or "chrY" in chr):
                if hap_gt_list[0] == "." and hap_gt_list[1] != ".":
                    hap_gt_list[0] = hap_gt_list[1]
                elif hap_gt_list[0] != "." and hap_gt_list[1] == ".":
                    hap_gt_list[1] = hap_gt_list[0]
            elif final_sample_sex_list[i] == "female" and "chrY" in chr:
                hap_gt_list = [".", "."]
            sample_gt_add_list.append(hap_gt_list[0] + "|" + hap_gt_list[1])
        print('\t'.join(line.strip().split()[0:9] + sample_gt_add_list))