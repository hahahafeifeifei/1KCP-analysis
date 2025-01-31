import sys
import gzip

#Annotate the repeat content for variant site based on allele tranversal
#python3 graph_vcf_repeat_annotate.py input.id_at.list input.anno input.hap.count > output.repeat.anno

def openfile(filename):
    if filename == "-":
        return sys.stdin
    elif filename.endswith('.gz'):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")


list_file = openfile(sys.argv[1])
annotate_file = openfile(sys.argv[2])
hap_count_file = openfile(sys.argv[3])
anno_class = ["LINE","SINE","LTR","DNA","Retroposon","Low_complexity","TR"]
########Read len data
node_len_dict = {}
for line in hap_count_file:
    node_id = line.strip().split()[0]
    node_len = int(line.strip().split()[1])
    node_len_dict[node_id] = node_len

########Read node data
node_anno_dict = {}
for line in annotate_file:
    node_id = line.strip().split()[0]
    node_anno = line.strip().split()[2]
    if node_anno in anno_class:
        if node_id not in node_anno_dict:
            node_anno_dict[node_id] = [False] * len(anno_class)
        node_anno_dict[node_id][anno_class.index(node_anno)] = True

########Process vcf file
for line in list_file:
    line_field = line.strip().split()
    id = line_field[0]
    at = line_field[1]
    at_len = 0
    non_repeat_len = 0
    nl_list = [0] * len(anno_class)

    at_node = at.replace('<','>').split(">")[2:-1]
    for node in at_node:
        at_len += node_len_dict[node]
        if node in node_anno_dict:
            for i, class_type in enumerate(node_anno_dict[node]):
                if class_type:
                    nl_list[i] += node_len_dict[node]
        else:
            non_repeat_len += node_len_dict[node]

    if non_repeat_len/at_len > 0.7:
        type = "Low repeat"
    elif max(nl_list)/at_len >= 0.3:
        type = anno_class[nl_list.index(max(nl_list))]
    else:
        type = "Mixed repeat"

    print("\t".join([id, str(at_len), type]))