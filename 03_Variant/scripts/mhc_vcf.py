import sys

# Paths to your files (update these paths accordingly)
input_data_file = sys.argv[1]      # File containing AL lines
sample_list_file = sys.argv[2]    # File containing the list of sample names

# Read the sample list from the sample list file
with open(sample_list_file, 'r') as f:
    samples = [line.strip().split()[0] for line in f if line.strip()]
all_haps = [s + ".hap1" for s in samples] + [s + ".hap2" for s in samples]
# Create a mapping of sample names to indices for easier lookup
sample_indices = {sample: idx for idx, sample in enumerate(samples)}
num_samples = len(samples)

gene_hap_dict = {}
with open(input_data_file, 'r') as f:
    for line in f:
        line_info = line.strip().split()
        name = line_info[2] + "_" + line_info[4]
        if name not in gene_hap_dict:
            gene_hap_dict[name] = line_info[6].split(',')
        else:
            gene_hap_dict[name] += line_info[6].split(',')
gene_missing_dict = {}
for name, haps in gene_hap_dict.items():
    gene_missing_dict[name] = list(set(all_haps) - set(haps))


print("\n".join(['##fileformat=VCFv4.2',
        '##INFO=<ID=Gene,Number=1,Type=String,Description="HLA gene">',        
        '##INFO=<ID=Allele,Number=1,Type=String,Description="HLA allele">',
        '##INFO=<ID=Field,Number=1,Type=Integer,Description="Resolution of HLA allele">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##contig=<ID=chr6>']))
print('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'] + samples))
# Read the input data from the AL lines file
with open(input_data_file, 'r') as f:
    for line in f:
        line_info = line.strip().split()
        chrom = line_info[0]
        pos = line_info[1]
        gene = line_info[2]
        allele = line_info[3]
        field = line_info[4]
        supporting_samples = line_info[6].split(',')
        id = chrom + "-" + str(pos) + "-" + allele + "-" + field + "Field"
        # Store the genotypes for this allele
        genotypes = ['0|0'] * num_samples
        name = gene + "_" + field
        for sample in gene_missing_dict[name]:
            sample_name = sample.split('.hap')[0]
            haplotype = sample.split('.hap')[1].split(".")[0]
            if sample_name in sample_indices:
                sample_idx = sample_indices[sample_name]
                if haplotype == '1':  # hap1
                    gt = genotypes[sample_idx].split('|')
                    gt[0] = '.'
                    genotypes[sample_idx] = '|'.join(gt)
                elif haplotype == '2':  # hap2
                    gt = genotypes[sample_idx].split('|')
                    gt[1] = '.'
                    genotypes[sample_idx] = '|'.join(gt)
        for sample in supporting_samples:
            # Extract sample name and haplotype info
            sample_name = sample.split('.hap')[0]
            haplotype = sample.split('.hap')[1].split(".")[0]
            # Get the index of the sample
            if sample_name in sample_indices:
                sample_idx = sample_indices[sample_name]
                # Update genotype based on haplotype
                if haplotype == '1':  # hap1
                    gt = genotypes[sample_idx].split('|')
                    gt[0] = '1'
                    genotypes[sample_idx] = '|'.join(gt)
                elif haplotype == '2':  # hap2
                    gt = genotypes[sample_idx].split('|')
                    gt[1] = '1'
                    genotypes[sample_idx] = '|'.join(gt)
        info = 'Gene=' + gene + ";" + 'Allele=' + allele + ";" + 'Field=' + field
        print("\t".join([chrom, str(pos), id, 'N', '<Allele>', '.', '.', info , 'GT'] + genotypes))
