#!/usr/bin/env python3
import sys
import gzip

def openfile(filename):
    if filename == "-":
        return sys.stdin
    elif filename.endswith('.gz'):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")
    
mode = sys.argv[1]

if mode == "split":
    for line in sys.stdin:
        if line.startswith("##"):
            print(line.strip())
        elif line.startswith("#"):
            sample_list = line.strip().split()[9:]
            final_sample_list = []
            for sample in sample_list:
                final_sample_list.append(sample + ".hap1")
                final_sample_list.append(sample + ".hap2")
            print('\t'.join(line.strip().split()[0:9] + final_sample_list))
        else:
            sample_gt_list = line.strip().split()[9:]
            final_sample_gt_list = []
            for sample_gt in sample_gt_list:
                hap_gt = sample_gt.replace("/","|").split("|")
                final_sample_gt_list.append(hap_gt[0] + "|" + hap_gt[0])
                final_sample_gt_list.append(hap_gt[1] + "|" + hap_gt[1])
            print('\t'.join(line.strip().split()[0:9] + final_sample_gt_list))

if mode == "merge":
    for line in sys.stdin:
        if line.startswith("##"):
            print(line.strip())
        elif line.startswith("#"):
            final_sample_gt = ""
            hap_list = line.strip().split()[9:]
            final_sample_list = []
            for hap in hap_list:
                sample = ".".join(hap.split(".")[:-1])
                if sample not in final_sample_list:
                    final_sample_list.append(sample)
            print('\t'.join(line.strip().split()[0:9] + final_sample_list))
        else:
            sample_gt_list = line.strip().split()[9:]
            final_sample_gt_list = []
            for i in range(len(sample_gt_list)):
                hap_gt = sample_gt_list[i].replace("/","|").split("|")[0]
                if i % 2 == 0:
                    final_sample_gt += hap_gt
                    final_sample_gt += "|"
                elif i % 2 == 1:
                    final_sample_gt += hap_gt
                    final_sample_gt_list.append(final_sample_gt)
                    final_sample_gt = ""
            print('\t'.join(line.strip().split()[0:9] + final_sample_gt_list))

