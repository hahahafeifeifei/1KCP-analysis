import sys
import subprocess
import h5py
from sklearn.metrics import roc_auc_score
from sklearn.metrics import auc
from sklearn.metrics import precision_recall_curve
import numpy as np

# Compare the epigenome prediction performance based on AUROC with truth
# python3 epi_benchmark.py input_full.hdf compare_track.list truth.peak.bed truth.nonref_peck.bed prefix 128

input_file = sys.argv[1]
atac_file = sys.argv[2]
truth_bed = sys.argv[3]
nonref_bed = sys.argv[4]
prefix = sys.argv[5]

step = int(sys.argv[6])
atac_index = []
for line in open(atac_file):
    atac_index.append(int(line.strip())-1)

profile_array_list = []
fo=open(prefix + ".pos.bed","w")
f = h5py.File(input_file)
for sei_key, sei_value in f.items():
    if "seq_class_scores" not in sei_key:
        contig = sei_key
        region_start = 0
        region_end = step
        region_value = ""
        for i in range(sei_value.shape[0]):
            profile_array_list.append(sei_value[i][atac_index])
            fo.write("\t".join([contig,str(region_start),str(region_end)]) + "\n")
            region_start += step
            region_end += step
fo.close()
profile_array = np.stack(profile_array_list, axis=0)
profile_array_list = []

label_result = subprocess.run("bedtools intersect -a {} -b {} -loj | awk -v OFS=\'\\t\'  \'{{if($4==\".\") print $1,$2,$3,\"0\";else print $1,$2,$3,\"1\"}}\' | uniq | awk \'{{print $4}}\' ".format(prefix + ".pos.bed", truth_bed), shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
nonref_result = subprocess.run("bedtools intersect -a {} -b {} -loj | awk -v OFS=\'\\t\' \'{{if($4==\".\") print $1,$2,$3,\"0\";else print $1,$2,$3,\"1\"}}\' | uniq | awk \'{{print $4}}\' ".format(prefix + ".pos.bed", nonref_bed), shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
label_array = np.array(list(map(int,label_result.stdout.decode('UTF-8').strip().split("\n"))))
nonref_array = np.array(list(map(int,nonref_result.stdout.decode('UTF-8').strip().split("\n"))))

for i in range(len(atac_index)):
    all_auroc = roc_auc_score(label_array,profile_array[:,i])
    nonref_auroc = roc_auc_score(label_array[nonref_array==1],profile_array[nonref_array==1,i])
    print(str(atac_index[i]+1) + "\t" + str(all_auroc) + "\t" + str(nonref_auroc) )