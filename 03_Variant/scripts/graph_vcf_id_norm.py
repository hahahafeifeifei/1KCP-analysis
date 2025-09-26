import re
import sys
import gzip

def openfile(filename):
    if filename == "-":
        return sys.stdin
    elif filename.endswith('.gz'):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")

def id_norm(id):
    direct_dict = {">" : "<", "<" : ">"}
    id_node = re.split('<|>', id)[1:]
    id_direct = re.split('\d+', id)[:-1]
    if id_direct[0] == "<":
        id = direct_dict[id_direct[1]] + id_node[1] + direct_dict[id_direct[0]] + id_node[0]
    return id

vcf_file = openfile(sys.argv[1])

########Process vcf file
for line in vcf_file:
    if line[0:2] == "##":
        print(line, end="")
    elif line[0] == "#":
        print('##INFO=<ID=BUBBLE_ID,Number=1,Type=String,Description="ID of pangenome bubble">')
        print(line, end="")
    else:
        line_list = line.strip().split()
        info_dict = {info.split("=")[0] : info.split("=")[1] for info in line_list[7].split(";")}
        chr = line_list[0]
        pos = line_list[1]
        id = id_norm(line_list[2])
        info_dict["BUBBLE_ID"] = id
        ref = line_list[3]
        alt = line_list[4]
        vt = info_dict["VT"]
        allele_type = "Multiallelic" if vt.split("_")[0] == "Multiallelic" else "Biallelic"
        variant = vt.split("_")[1]
        new_id = "-".join([chr, pos, allele_type, variant, id])
        line_list[2] = new_id
        line_list[7] = ";".join([key + "=" + value for key, value in info_dict.items()])
        print('\t'.join(line_list))
