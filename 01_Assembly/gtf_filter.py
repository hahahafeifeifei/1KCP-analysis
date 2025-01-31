import sys
import gzip

# Filter transcripts based on TPM and exon number
# python3 gtf_filter.py input.gtf output.filter.gtf

def openfile(filename):
    if filename == "-":
        return sys.stdin
    elif filename.endswith('.gz'):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")

input_file = sys.argv[1]
out_file = sys.argv[2]
tpm_threshold = 0.5
exon_count_threshold = 1

transcripts_tpm = {}  # Dictionary to store TPM values for transcripts
transcripts_exon = {}

f = openfile(input_file)
for line in f:
    if not line.startswith("#"):
        fields = line.strip().split("\t")
        feature_type = fields[2]
        attributes_dict = {attribute.strip().split(" ")[0]:attribute.strip().split(" ")[1].strip('"') for attribute in fields[8].split(";")[:-1]}
        if feature_type == "transcript":
            transcripts_tpm[attributes_dict["transcript_id"]] = float(attributes_dict["TPM"])
        if feature_type == "exon":
            transcripts_exon[attributes_dict["transcript_id"]] = transcripts_exon.get(attributes_dict["transcript_id"], 0) + 1

f = openfile(input_file)
fo = open(out_file, "w")
for line in f:
    if line.startswith("#"):
        fo.write(line)
    else:
        fields = line.strip().split("\t")
        attributes_dict = {attribute.strip().split(" ")[0]:attribute.strip().split(" ")[1].strip('"') for attribute in fields[8].split(";")[:-1]}
        if transcripts_exon[attributes_dict["transcript_id"]] > exon_count_threshold and transcripts_tpm[attributes_dict["transcript_id"]] >= tpm_threshold:
            fo.write(line)
fo.close()