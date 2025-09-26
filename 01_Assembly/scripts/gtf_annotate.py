import sys
import gzip

def openfile(filename):
    if filename == "-":
        return sys.stdin
    elif filename.endswith('.gz'):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")
    
input_file = sys.argv[1]
type_file = sys.argv[2]
out_file = sys.argv[3]

type_dict = {}
for line in openfile(type_file):
    type_dict[line.strip().split()[0]] = line.strip().split()[1]

fo = open(out_file, "w")
for line in openfile(input_file):
    if line.startswith("#"):
        fo.write(line)
    else:
        fields = line.strip().split("\t")
        feature_type = fields[2]
        attributes_dict = {attribute.strip().split(" ")[0]:attribute.strip().split(" ")[1].strip('"') for attribute in fields[8].split(";")[:-1]}
        if attributes_dict["transcript_id"] in type_dict:
            fields[8] += ' gene_type "' + type_dict[attributes_dict["transcript_id"]] + '";'
            fo.write("\t".join(fields) + "\n")
fo.close()
