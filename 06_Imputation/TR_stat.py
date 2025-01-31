#!/usr/bin/env python3
import sys
import gzip
import numpy as np
from scipy import stats

#Calculate R2 for TR variant
#python3 TR_stat.py input output

def openfile(filename):
    if filename == "-":
        return sys.stdin
    elif filename.endswith('.gz'):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")

input_file = sys.argv[1]
variant_dict = {}
for line in openfile(input_file):
    line_info = line.strip().split()
    variant = line_info[1]
    if variant not in variant_dict:
        variant_dict[variant] = [[], [], []]
    truth = float(line_info[2])
    gt = float(line_info[3])
    dosage = float(line_info[4])
    variant_dict[variant][0].append(truth)
    variant_dict[variant][1].append(dosage)
    variant_dict[variant][2].append(gt)

for variant, vector in variant_dict.items():
    if len(vector[0]) != 1:
        r_dosage = stats.pearsonr(vector[0], vector[1])[0]
        r2_dosage = r_dosage * r_dosage
        r_gt = stats.pearsonr(vector[0], vector[2])[0]
        r2_gt = r_gt * r_gt
    else:
        r2_dosage = "nan"
        r2_gt = "nan"
    print("\t".join([variant, str(r2_dosage), str(r2_gt)]))