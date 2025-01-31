import re
import sys
import gzip

#Annotate the variant type based on allele traversal
#python3 graph_vcf_vt_annotate.py input.vcf input.gfa > output.vcf

def openfile(filename):
    if filename == "-":
        return sys.stdin
    elif filename.endswith('.gz'):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")

vcf_file = openfile(sys.argv[1])
gfa_file = openfile(sys.argv[2])

########Read graph data
node_len_dict = {}
for line in gfa_file:
    if line[0] == "S":
        node_id = line.strip().split()[1]
        node_len_dict[node_id] = len(line.strip().split()[2])

########Process vcf file
for line in vcf_file:
    if line[0:2] == "##":
        print(line.strip())
    elif line[0] == "#":
        print('##INFO=<ID=TL,Number=R,Type=Integer,Description="Length of allele Traversal">')
        print('##INFO=<ID=VT,Number=1,Type=String,Description="Variant type">')
        print(line.strip())
    else:
        line_list = line.strip().split()
        info_dict = {info.split("=")[0] : info.split("=")[1] for info in line_list[7].split(";")}

        at_list = info_dict["AT"].split(",")
        tl_list = []
        for at in at_list:
            if at == 'None' or at == ".":
                tl_list.append('None')
            else:
                at_node = re.split('<|>', at)[2:-1]
                tl = 0
                for node in at_node:
                    tl += node_len_dict[node]
                tl_list.append(tl)
        tl = ','.join([str(x) for x in tl_list])
        tl_clean_list = []
        for length in tl_list:
            if length != "None":
                tl_clean_list.append(length)
        vt = None
        if len(tl_clean_list) == 2:
            if max(tl_clean_list) - min(tl_clean_list) == 0:
                vt = 'Biallelic_SNV'
            elif max(tl_clean_list) - min(tl_clean_list) >= 1 and max(tl_clean_list) - min(tl_clean_list) < 50:
                vt = 'Biallelic_INDEL'
            elif max(tl_clean_list) - min(tl_clean_list) >= 50:
                vt = 'Biallelic_SV'
        elif len(tl_clean_list) >= 3:
            if max(tl_clean_list) - min(tl_clean_list) == 0:
                vt = 'Multiallelic_SNV'
            elif max(tl_clean_list) - min(tl_clean_list) >= 1 and max(tl_clean_list) - min(tl_clean_list) < 50:
                vt = 'Multiallelic_INDEL'
            elif max(tl_clean_list) - min(tl_clean_list) >= 50:
                vt = 'Multiallelic_SV'
        if vt == None:
            continue
        info_dict["TL"] = tl
        info_dict["VT"] = vt
        line_list[7] = ";".join([key + "=" + value for key, value in info_dict.items()])
        print('\t'.join(line_list))
