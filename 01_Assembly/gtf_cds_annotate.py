import sys
import gzip

# Add the CDS annotation from one gff fill to another gtf file
# python3 gtf_filter.py input.no_cds.gtf input.cds.gff input.add_cds.gtf

def openfile(filename):
    if filename == "-":
        return sys.stdin
    elif filename.endswith('.gz'):
        return gzip.open(filename, "rt")
    else:
        return open(filename, "r")

def interval_overlap(start1, end1, start2, end2):
    max_start = start1 if start1 >= start2 else start2
    min_end = end1 if end1 <= end2 else end2
    overlap_len = max(min_end - max_start + 1, 0)
    return overlap_len

input_gtf1_file = sys.argv[1]
input_gtf2_file = sys.argv[2]
output_file = sys.argv[3]

cds_coords = {}
for line in openfile(input_gtf2_file):
    fields = line.strip().split('\t')
    if len(fields) == 9 and fields[2] == "CDS":
        transcript_id = fields[0]
        start = int(fields[3])
        end = int(fields[4])
        cds_coords[transcript_id] = [start, end]

cds_len_dict = {}
for line in openfile(input_gtf1_file):
    if not line.startswith("#"):
        fields = line.strip().split('\t')
        feature_type = fields[2]
        attributes_dict = {attribute.strip().split(" ")[0]:attribute.strip().split(" ")[1].strip('"') for attribute in fields[8].split(";")[:-1]}
        if feature_type == 'exon':
            transcript_id = attributes_dict["transcript_id"]
            exon_len = int(fields[4]) - int(fields[3]) + 1
            cds_len_dict[transcript_id] = cds_len_dict.get(transcript_id, 0) + exon_len

fo = open(output_file, "w")
for line in openfile(input_gtf1_file):
    if line.startswith("#"):
        fo.write(line)
    else:
        fo.write(line)
        fields = line.strip().split('\t')
        feature_type = fields[2]
        attributes_dict = {attribute.strip().split(" ")[0]:attribute.strip().split(" ")[1].strip('"') for attribute in fields[8].split(";")[:-1]}
        if feature_type == 'transcript':
            transcripts_attributes = attributes_dict
            cds_len = 0
        if feature_type == 'exon':
            transcript_id = attributes_dict["transcript_id"]
            if transcript_id in cds_coords:
                exon_start = int(fields[3])
                exon_end = int(fields[4])
                strand = fields[6]
                if strand == "+" or strand == ".":
                    cds_start = cds_coords[transcript_id][0]
                    cds_end = cds_coords[transcript_id][1]
                elif strand == "-":
                    cds_start = cds_len_dict[transcript_id] - cds_coords[transcript_id][1] + 1
                    cds_end = cds_len_dict[transcript_id] - cds_coords[transcript_id][0] + 1
                cds_genome_start = exon_start - 1 + cds_start - cds_len
                cds_genome_end = exon_start - 1 + cds_end - cds_len
                fo.write("\t".join([str(exon_start), str(exon_end), str(cds_genome_start), str(cds_genome_end), str(cds_len)])+"\n")
                if strand == "+":
                    ###Exon start downstream CDS 5-UTR
                    if cds_genome_start > exon_start:
                        fo.write("\t".join(fields[0:2] + ['five_prime_UTR', str(exon_start), str(min(cds_genome_start - 1, exon_end))] + fields[5:]) + "\n")
                    ###Exon end upstream CDS 3-UTR
                    if cds_genome_end < exon_end:
                        fo.write("\t".join(fields[0:2] + ['three_prime_UTR', str(max(cds_genome_end + 1, exon_start)), str(exon_end)] + fields[5:]) + "\n")
                    ###Exon CDS overlap stop codon CDS
                    overlen = interval_overlap(exon_start, exon_end, cds_genome_start, cds_genome_end)
                    if overlen > 0:
                        fo.write("\t".join(fields[0:2] + ['CDS', str(max(cds_genome_start, exon_start)), str(min(cds_genome_end, exon_end))] + fields[5:]) + "\n")
                elif strand == "-":
                    ###Exon start downstream CDS 5-UTR
                    if cds_genome_start > exon_start:
                        fo.write("\t".join(fields[0:2] + ['three_prime_UTR', str(exon_start), str(min(cds_genome_start - 1, exon_end))] + fields[5:]) + "\n")
                    ###Exon end upstream CDS 3-UTR
                    if cds_genome_end < exon_end:
                        fo.write("\t".join(fields[0:2] + ['five_prime_UTR', str(max(cds_genome_end + 1, exon_start)), str(exon_end)] + fields[5:]) + "\n")
                    ###Exon CDS overlap stop codon CDS
                    overlen = interval_overlap(exon_start, exon_end, cds_genome_start, cds_genome_end)
                    if overlen > 0:
                        fo.write("\t".join(fields[0:2] + ['CDS', str(max(cds_genome_start, exon_start)), str(min(cds_genome_end, exon_end))] + fields[5:]) + "\n")
                cds_len = cds_len + exon_end - exon_start + 1
fo.close()