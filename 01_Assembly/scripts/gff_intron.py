import sys

def parse_gtf(gtf_file):
    """Parse GTF file and return a dictionary of gene IDs with their corresponding exons."""
    gene_exons = {}
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if fields[2] == 'exon':
                attributes = dict(item.strip().split("=") for item in fields[8].split(';') if item.strip())
                gene_id = attributes.get('gene_id').strip('"')
                exon_chr = fields[0]
                exon_start = int(fields[3])
                exon_end = int(fields[4])
                if gene_id not in gene_exons:
                    gene_exons[gene_id] = []
                gene_exons[gene_id].append([exon_chr, exon_start, exon_end])
    return gene_exons

def get_introns(gene_exons):
    """Compute intron regions based on exon regions."""
    gene_introns = {}
    for gene_id, exons in gene_exons.items():
        intron_chr = exons[0][0]
        exons = [(exon[1], exon[2]) for exon in exons]
        exons.sort()
        introns = []
        prev_end = 0
        for exon_start, exon_end in exons:
            if prev_end != 0 and prev_end < exon_start - 1:
                introns.append([intron_chr, prev_end, exon_start - 1])
            prev_end = max(exon_end, prev_end)
        gene_introns[gene_id] = introns
    return gene_introns

def write_gff(gene_introns, output_file):
    """Write intron regions to a GFF file."""
    with open(output_file, 'w') as f:
        for gene_id, introns in gene_introns.items():
            for index, intron in enumerate(introns, start=1):
                f.write(f"{intron[0]}\t{intron[1]}\t{intron[2]}\tintron\n")


# Parse GTF file to get exon regions
gene_exons = parse_gtf(sys.argv[1])

# Compute intron regions based on exon regions
gene_introns = get_introns(gene_exons)

# Write intron regions to a GFF file
write_gff(gene_introns, sys.argv[2])
