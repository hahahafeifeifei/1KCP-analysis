import re
import sys
import gzip

#Add unique id for each variant
#python3 graph_vcf_id_add_biallelic.py input.vcf output.vcf

def openfile(filename):
    if filename == "-":
        return sys.stdin
    elif filename.endswith('.gz'):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")

vcf_file = openfile(sys.argv[1])
id_dict = {}
########Process vcf file
for line in vcf_file:
    if line[0] == "#":
        print(line, end="")
    else:
        line_list = line.strip().split()
        id = line_list[2]
        id_dict[id] = id_dict.get(id, 0) + 1
        new_id = id + "-" + str(id_dict[id])
        line_list[2] = new_id
        print('\t'.join(line_list))