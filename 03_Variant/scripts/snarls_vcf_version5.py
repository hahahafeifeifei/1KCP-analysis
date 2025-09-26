from numpy import delete
from treelib import Node, Tree
import re
import sys
import json
import gzip

def openfile(filename):
    if filename == "-":
        return sys.stdin
    elif filename.endswith('.gz'):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")

def get_seq(seq, direct):
    if direct == ">":
        return seq
    elif direct == "<":
        rev_dict = {'A':'T', 'G':'C', 'C':'G', 'T':'A', 'N':'N'}
        rev_seq = ''.join([rev_dict[base] for base in seq[::-1]])
        return rev_seq

vcf_file = openfile(sys.argv[1])
gfa_file = openfile(sys.argv[2])
snarls_file = openfile(sys.argv[3])
out_file = open(sys.argv[4], "w")

########Read graph data
node_len_dict = {}
node_seq_dict = {}
for line in gfa_file:
    if line[0] == "S":
        node_id = line.strip().split()[1]
        node_len_dict[node_id] = len(line.strip().split()[2])
        node_seq_dict[node_id] = line.strip().split()[2]

########Read bubble data
snarls_tree_dict = {}
snarls_node_dict = {}
last_tree = None
for line in snarls_file:
    snarls_json = json.loads(line)
    id = '_'.join([str(x) for x in sorted([int(snarls_json['start']['node_id']), int(snarls_json['end']['node_id'])])])
    if "parent" not in snarls_json.keys():
        tree = Tree()
        tree.create_node(id, id, data = 0)
        if last_tree != None:
            root_id = last_tree.root
            snarls_tree_dict[root_id] = last_tree
    else:
        parent_id = '_'.join([str(x) for x in sorted([int(snarls_json['parent']['start']['node_id']), int(snarls_json['parent']['end']['node_id'])])])
        tree.create_node(id, id, parent = parent_id, data = 0)
    snarls_node_dict[id] = tree.root
    last_tree = tree
root_id = last_tree.root
snarls_tree_dict[root_id] = last_tree

########Process vcf file
for line in vcf_file:
    if line[0:2] == "##":
        out_file.write(line)
    elif line[0] == "#":
        out_file.write('##INFO=<ID=BD,Number=1,Type=Integer,Description="Depth of Bubble located in snarl tree">' + '\n')
        out_file.write('##INFO=<ID=PB,Number=1,Type=String,Description="Parent Bubble">' + '\n')
        out_file.write('##INFO=<ID=RT,Number=1,Type=String,Description="Whether the bubble contains Reference Traversal">' + '\n')
        out_file.write(line)
    else:
        line_list = line.strip().split()
        chr = line_list[0]
        pos = line_list[1]
        qual = line_list[5]
        filter = line_list[6]
        format = line_list[8]
        gt_list = line_list[9:]

        ####Read info column
        info_dict = {info.split("=")[0] : info.split("=")[1] for info in line_list[7].split(";")}
        an = info_dict["AN"]
        ns = info_dict["NS"]
        lv = info_dict["LV"]
        at_list = info_dict["AT"].split(",")
        ac_list = info_dict["AC"].split(",")
        af_list = info_dict["AF"].split(",")
        filtered_at_index = [0] + [i + 1 for i, ac in enumerate(info_dict["AC"].split(",")) if int(ac) != 0]

        ####Construct traversal list
        at_node_list = []
        at_path_list = []
        for at in at_list:
            at_node = re.split('<|>', at)[1:]
            at_direct = re.split('\d+', at)[:-1]
            at_path = []
            for i in range(len(at_node)):
                at_path.append(at_direct[i] + at_node[i])
            at_node_list.append(at_node)
            at_path_list.append(at_path)

        ####Search the snarls tree of variant
        root_id = '_'.join([str(x) for x in sorted([int(at_node_list[0][0]), int(at_node_list[0][-1])])])
        try:
            snarls_tree = snarls_tree_dict[root_id]
        except:
            try:
                snarls_tree = snarls_tree_dict[snarls_node_dict[root_id]].subtree(root_id)
            except:
                continue

        ####Refine snarls tree and filter bubbles with < 2 transversals that have diferent nodes go through
        for bubble in snarls_tree.expand_tree(mode = Tree.DEPTH):
            at_list = []
            for i in filtered_at_index:
                if bubble.split("_")[0] in at_node_list[i] and bubble.split("_")[1] in at_node_list[i]:
                    bubble_index = sorted([at_node_list[i].index(bubble.split("_")[0]), at_node_list[i].index(bubble.split("_")[1])])
                    at = "".join(at_path_list[i][bubble_index[0] + 1 : bubble_index[1]])
                    if at not in at_list:
                        at_list.append(at)
                        snarls_tree[bubble].data += 1
        bubble_list = [bubble.identifier for bubble in snarls_tree.filter_nodes(lambda x: x.data<2)]

        root_remove = False
        for bubble in bubble_list:
            if bubble == snarls_tree.root:
                root_remove = True
            else:
                snarls_tree.link_past_node(bubble)

        ####Write record for each bubble
        bn_list = []
        for bubble in snarls_tree.expand_tree(mode = Tree.DEPTH):
            if bubble == snarls_tree.root and root_remove:
                continue
            tl_list = []
            alt_list = []
            at_list = []
            rt = None
            at_map_list = [0] * len(at_node_list)
            bubble_map_list = [0] * len(at_node_list)
            del_bubble = False
            bubble_at_index_list = []
            for i in filtered_at_index:
                if bubble.split("_")[0] in at_node_list[i] and bubble.split("_")[1] in at_node_list[i]:
                    bubble_index = sorted([at_node_list[i].index(bubble.split("_")[0]), at_node_list[i].index(bubble.split("_")[1])])
                    bubble_id = at_path_list[i][bubble_index[0]] + at_path_list[i][bubble_index[1]]
                    at = "".join(at_path_list[i][bubble_index[0] : bubble_index[1] + 1])
                    bubble_at_index_list.append(i)
                    if i == 0:
                        del_bubble = True
                        bubble_map_list = [1] * len(at_node_list)
                    if del_bubble:
                        bubble_map_list[i] = 0
                    else:
                        bubble_map_list[i] = 1
                    ####Merge alleles
                    if at not in at_list:
                        at_map_list[i] = len(at_list)
                        at_list.append(at)
                        tl = 0
                        for node in at_node_list[i][bubble_index[0] + 1 : bubble_index[1]]:
                            tl += node_len_dict[node]
                        tl_list.append(tl)
                        prefix_path = at_path_list[i][bubble_index[0]]
                        prefix_seq = node_seq_dict[re.split('<|>', prefix_path)[-1]]
                        prefix_direct = re.split('\d+', prefix_path)[0]
                        seq = get_seq(prefix_seq, prefix_direct)[-1]
                        for path in at_path_list[i][bubble_index[0] + 1 : bubble_index[1]]:
                            path_seq = node_seq_dict[re.split('<|>', path)[-1]]
                            path_direct = re.split('\d+', path)[0]
                            seq += get_seq(path_seq, path_direct)
                        if i == 0:
                            ref = seq
                            rt = 'True'
                        else:
                            alt_list.append(seq)
                    else:
                        at_map_list[i] = at_list.index(at)
                elif i == 0:
                    at_map_list[i] = 0
                    at = '.'
                    at_list.append(at)
                    tl = '.'
                    tl_list.append(tl)
                    ref = '*'
                    rt = 'False'
            
            bubble_write = len([True for i in filtered_at_index if bubble_map_list[i] == 1]) >= 2
            if del_bubble:
                bubble_ref = '<BUBBLE>'
                bubble_alt = line_list[3][0]
                bt = "DEL"
            else:
                bubble_ref = line_list[3][0]
                bubble_alt = '<BUBBLE>'
                bt = "INS"

            bubble_node = set()
            if root_remove:
                bd = snarls_tree.level(bubble)
            else:
                bd = snarls_tree.level(bubble) + 1
            if bd != 1:
                pb = None
                parent_bubble = snarls_tree.parent(bubble).tag
                for i in bubble_at_index_list:
                    if parent_bubble.split("_")[0] in at_node_list[i] and parent_bubble.split("_")[1] in at_node_list[i]:
                        parent_bubble_index = sorted([int(at_node_list[i].index(parent_bubble.split("_")[0])), int(at_node_list[i].index(parent_bubble.split("_")[1]))])
                        parent_bubble_node_at = set(at_path_list[i][parent_bubble_index[0] + 1 : parent_bubble_index[1]])
                        if bubble_node == set():
                            bubble_node = parent_bubble_node_at
                        else:
                            bubble_node = bubble_node & parent_bubble_node_at
                        pb = at_path_list[i][parent_bubble_index[0]] + at_path_list[i][parent_bubble_index[1]]
                    else:
                        pb = "Unknown"
            bn = ",".join([node for node in bubble_node])
            if bn in bn_list:
                bubble_write = False
            else:
                bn_list.append(bn)

            ####Bubble type classification
            tl_clean_list = []
            if '.' in tl_list:
                tl_clean_list = tl_list[1:]
            else:
                tl_clean_list = tl_list
            
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
            
            if vt == 'Biallelic_SNV' or vt == 'Multiallelic_SNV':
                alt_list = [seq[1:] for seq in alt_list]
                if ref != '*':
                    ref = ref[1:] 
            alt = ','.join(alt_list)

            if vt == None:
                continue

            ####INFO column
            at_ac_merge_list = [0] * max(at_map_list)
            at_af_merge_list = [0] * max(at_map_list)
            for origin in range(len(at_map_list)):
                at_merge = at_map_list[origin]
                if at_merge != 0:
                    at_ac_merge_list[at_merge - 1] = int(at_ac_merge_list[at_merge - 1]) + int(ac_list[origin - 1])
                    at_af_merge_list[at_merge - 1] = float(at_af_merge_list[at_merge - 1]) + float(af_list[origin - 1])
            at_ac_merge_list = [str(x) for x in at_ac_merge_list]
            at_af_merge_list = ['%.6g'%x for x in at_af_merge_list]

            bubble_ac_merge = 0
            bubble_af_merge = 0
            for origin in range(len(bubble_map_list)):
                if bubble_map_list[origin] != 0:
                    bubble_ac_merge = int(bubble_ac_merge) + int(ac_list[origin - 1])
                    bubble_af_merge = float(bubble_af_merge) + float(af_list[origin - 1])
            bubble_ac_merge = str(bubble_ac_merge)
            bubble_af_merge = '%.6g'%bubble_af_merge

            if bd == 1:
                at_info = ";".join(["AC=" + ','.join(at_ac_merge_list), "AF=" + ','.join(at_af_merge_list), "AN=" + an, "AT=" + ','.join(at_list),
                "BD=" + str(bd), "LV=" + lv, "NS=" + ns, "RT=" + rt])
                bubble_info = ";".join(["AC=" + bubble_ac_merge, "AF=" + bubble_af_merge, "AN=" + an,
                "BD=" + str(bd), "BT=" + bt, "LV=" + lv, "NS=" + ns, "RT=" + rt, "ST=" + ','.join(at_list)])
            else:
                at_info = ";".join(["AC=" + ','.join(at_ac_merge_list), "AF=" + ','.join(at_af_merge_list), "AN=" + an, "AT=" + ','.join(at_list),
                "BD=" + str(bd), "LV=" + lv, "NS=" + ns, "PB=" + pb , "RT=" + rt])
                bubble_info = ";".join(["AC=" + bubble_ac_merge, "AF=" + bubble_af_merge, "AN=" + an,
                "BD=" + str(bd), "BN=" + bn, "BT=" + bt, "LV=" + lv, "NS=" + ns, "PB=" + pb, "RT=" + rt, "ST=" + ','.join(at_list)])

            ####GT transform
            at_gt_trans_list = []
            for gt in gt_list:
                gt_trans = []
                for gt_hap in gt.split("|"):
                    if gt_hap == ".":
                        gt_trans.append(gt_hap)
                    else:
                        gt_trans.append(str(at_map_list[int(gt_hap)]))
                at_gt_trans_list.append("|".join(gt_trans))

            bubble_gt_trans_list = []
            for gt in gt_list:
                gt_trans = []
                for gt_hap in gt.split("|"):
                    if gt_hap == ".":
                        gt_trans.append(gt_hap)
                    else:
                        gt_trans.append(str(bubble_map_list[int(gt_hap)]))
                bubble_gt_trans_list.append("|".join(gt_trans))
            if rt == "False":
                out_file.write('\t'.join([chr, pos, bubble_id, ref, alt, qual, filter, at_info, format] + at_gt_trans_list) + '\n')
            #if bubble_write:
                #out_file.write('\t'.join([chr, pos, bubble_id, bubble_ref, bubble_alt, qual, filter, bubble_info, format] + bubble_gt_trans_list) + '\n')
out_file.close()

