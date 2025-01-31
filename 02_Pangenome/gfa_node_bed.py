import sys
import gzip
import re

#Generate the bed file of nodes traversed by the given assembly
#python3 gfa_node_bed.py input.gfa output.bed sample haplotype

def openfile(filename):
    if filename == "-":
        return sys.stdin
    elif filename.endswith('.gz'):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")

def main():
    gfa_file = sys.argv[1]
    out_bed_file = sys.argv[2]
    target_sample = sys.argv[3]
    target_hap = sys.argv[4]

    node_len_dict = {}
    for line in openfile(gfa_file):
        if line[0] == "S":
            node = line.strip().split()[1]
            node_len_dict[node] = len(line.strip().split()[2])

    fo = open(out_bed_file, "w")
    for line in openfile(gfa_file):
        if line[0] == "W":
            sample = line.strip().split()[1]
            hap = line.strip().split()[2]
            if sample == target_sample and hap == target_hap:
                contig = line.strip().split()[3]
                contig_start = int(line.strip().split()[4])
                path_info = line.strip().split()[6]
                path_node = re.split('<|>', path_info)[1:]
                node_start = contig_start
                for node in path_node:
                    node_end = node_start + node_len_dict[node]
                    fo.write('\t'.join([contig, str(node_start), str(node_end), "s" + node]) + '\n')
                    node_start = node_end
    fo.close()
                    
if __name__ == "__main__":
    main()