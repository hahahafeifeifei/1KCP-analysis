import edlib
import sys

input_vcf_file = sys.argv[1]
output_vcf_file = sys.argv[2]
info_flag = False
format_flag = False
fo = open(output_vcf_file, "w")
for line in open(input_vcf_file):
    if line.startswith("#"):
        if line[:6] == "##INFO":
            if not info_flag:
                fo.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">' + "\n")
                fo.write('##INFO=<ID=RU,Number=1,Type=String,Description="Comma separated motif sequences list in the reference orientation">' + "\n")
                fo.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">' + "\n")
                fo.write('##INFO=<ID=ALTANNO_H1,Number=1,Type=String,Description="Motif representation for the h1 alternate allele>"' + "\n")
                fo.write('##INFO=<ID=ALTANNO_H2,Number=1,Type=String,Description="Motif representation for the h2 alternate allele>"' + "\n")
                fo.write('##INFO=<ID=LEN_H1,Number=1,Type=Integer,Description="Length of the motif annotation for the h1 alternate allele>"' + "\n")
                fo.write('##INFO=<ID=LEN_H2,Number=1,Type=Integer,Description="Length of the motif annotation for the h2 alternate allele>"' + "\n")
                info_flag = True
        elif line[:8] == "##FORMAT":
            if not format_flag:
                fo.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' + "\n")
                fo.write('##FORMAT=<ID=SS,Number=1,Type=Float,Description="Sequence similarity (edit distance)">' + "\n")
                format_flag = True
        else:
            fo.write(line)
    else:
        line_list = line.strip().split()
        info_dict = {info.split("=")[0]:info.split("=")[1] for info in line_list[7].split(";")[:-1]}
        motif_list = info_dict['RU'].split(",")
        allele_path_list = info_dict['ALTANNO_H1'].split("-")
        target_seq = ""
        for motif_index in allele_path_list: 
            target_seq += motif_list[int(motif_index)]
        query_seq = line_list[9].split(":")[1]
        ss = 1 - edlib.align(target_seq, query_seq, mode = "NW", task = "distance")["editDistance"]/max(len(target_seq), len(query_seq))
        line_list[8] = "GT:SS"
        line_list[9] = "1/1:" + "{:.4f}".format(ss)
        fo.write("\t".join(line_list) + "\n")
fo.close()        


