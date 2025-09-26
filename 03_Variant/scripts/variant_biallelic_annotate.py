import sys
import gzip

def openfile(filename):
    if filename == "-":
        return sys.stdin
    elif filename.endswith('.gz'):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")


vcf_file = openfile(sys.argv[1])
annotate_file = openfile(sys.argv[2])
out_nc_file = open(sys.argv[3], "w")
out_nf_file = open(sys.argv[4], "w")

anno_class = ["CDS","five_prime_UTR","three_prime_UTR","nc_exon","intron",
              "CTCF","E","HET", "P","PC","TF",
              "LINE","SINE","LTR","DNA","Retroposon","Low_complexity","TR"]

########Read node data
node_len_dict = {}
node_anno_dict = {}
for line in annotate_file:
    node_id = line.strip().split()[0]
    node_len = int(line.strip().split()[1])
    node_anno = line.strip().split()[2]
    node_len_dict[node_id] = node_len
    if node_id not in node_anno_dict:
        node_anno_dict[node_id] = [False] * len(anno_class)
    node_anno_dict[node_id][anno_class.index(node_anno)] = True

########Process vcf file
for line in vcf_file:
    if line[0:1] == "#":
        continue
    else:
        line_field = line.strip().split()
        info_dict = {info.split("=")[0] : info.split("=")[1] for info in line_field[7].split(";")}
        if "AT" not in info_dict:
            continue
        if "ID" in info_dict:
            bubble_id = info_dict["ID"]
        elif "BUBBLE_ID" in info_dict:
            bubble_id = info_dict["BUBBLE_ID"]
        else:
            continue
        variant_id = line_field[2]
        allele_count = len(line_field[4].split(",")) + 1
        if allele_count != 2:
            print(variant_id + "is not biallelic variant.")

        nf_list = [0] * len(anno_class)
        bubble_node = bubble_id.replace('<','>').split(">")[1:]
        for node in bubble_node:
            if node in node_anno_dict:
                for j, class_type in enumerate(node_anno_dict[node]):
                    if class_type:
                        nf_list[j] += 1
        for i, nf in enumerate(nf_list):
            nc = anno_class[i]
            if nf > 0:
                out_nf_file.write("\t".join([variant_id, nc, str(nf)]) + '\n')

        nl_list = [[0] * allele_count for i in range(len(anno_class))]
        at_list = info_dict["AT"].split(",")
        for i, at in enumerate(at_list):
            if at == "None" or at == ".":
                for j in range(len(anno_class)):
                    nl_list[j][i] = "None"
            else:
                at_node = at.replace('<','>').split(">")[2:-1]
                for node in at_node:
                    if node in node_anno_dict:
                        for j, class_type in enumerate(node_anno_dict[node]):
                            if class_type:
                                nl_list[j][i] += node_len_dict[node]
                            
        for i, nl_len in enumerate(nl_list):
            nl_clean_len = [len for len in nl_len if len != "None" ]
            nc = anno_class[i]
            if max(nl_clean_len) - min(nl_clean_len) > 0:
                out_nc_file.write("\t".join([variant_id, nc, str(nl_len[0]), str(nl_len[1])]) + '\n')
out_nc_file.close()
out_nf_file.close()
