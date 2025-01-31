#!/usr/bin/env python3
import sys
import gzip
import numpy as np
from scipy import stats

#Calculate R2 for non-TR variant
#python3 nonTR_stat.py input output

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
        variant_dict[variant] = [[], [], [0, 0]]
    truth = int(line_info[2])
    gt = int(line_info[3])
    dosage = float(line_info[4])
    variant_dict[variant][0].append(truth)
    variant_dict[variant][1].append(dosage)
    if truth >= 1:
        variant_dict[variant][2][0] += 1
        if truth == gt:
            variant_dict[variant][2][1] += 1

for variant, vector in variant_dict.items():
    if len(vector[0]) != 1:
        r = stats.pearsonr(vector[0], vector[1])[0]
        r2 = r * r
    else:
        r2 = "nan"
    if vector[2][0] != 0:
        concordance = vector[2][1]/vector[2][0]
    else:
        concordance = "nan"
    print("\t".join([variant, str(r2), str(concordance)]))