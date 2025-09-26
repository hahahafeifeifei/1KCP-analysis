import sys
import re

def parse_info(info_str):
    info_dict = {}
    for item in info_str.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info_dict[key] = value
        else:
            info_dict[item] = True
    return info_dict

def reconstruct_sequence(altanno_str, ru_list):
    indices = [int(x) for x in altanno_str.split('-')]
    sequence = ''.join([ru_list[i] for i in indices])
    return sequence

def main():
    if len(sys.argv) != 4:
        print("Usage: python script.py <vcf_file> <ref_altanno_file> <output_file>")
        sys.exit(1)
    
    vcf_file = sys.argv[1]
    ref_altanno_file = sys.argv[2]
    output_file = sys.argv[3]
    
    ref_altanno_dict = {}
    with open(ref_altanno_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            chrom, pos, ref_altanno = parts[0], parts[1], parts[2]
            ref_altanno_dict[(chrom, pos)] = ref_altanno
    
    with open(vcf_file, 'r') as vcf, open(output_file, 'w') as out:
        for line in vcf:
            if line.startswith('#'):
                out.write(line)
                continue
            
            fields = line.strip().split('\t')
            chrom, pos = fields[0], fields[1]
            info = parse_info(fields[7])
            
            if 'RU' not in info or 'ALTANNO' not in info:
                continue
            ru_list = info['RU'].split(',')
            altanno_list = info['ALTANNO'].split(',')
            
            key = (chrom, pos)
            if key not in ref_altanno_dict:
                continue
            ref_altanno = ref_altanno_dict[key]
            
            if ref_altanno not in altanno_list:
                continue
            ref_index = altanno_list.index(ref_altanno)
            
            allele_sequences = []
            for altanno in altanno_list:
                allele_sequences.append(reconstruct_sequence(altanno, ru_list))
            
            ref_sequence = allele_sequences[ref_index]
            alt_sequences = [seq for i, seq in enumerate(allele_sequences) if i != ref_index]
            
            fields[3] = ref_sequence
            fields[4] = ','.join(alt_sequences) if alt_sequences else '.'
            
            new_altanno_list = [altanno_list[ref_index]] + [altanno for i, altanno in enumerate(altanno_list) if i != ref_index]
            info['ALTANNO'] = ','.join(new_altanno_list)

            info_str = []
            for key, value in info.items():
                if value is True:
                    info_str.append(key)
                else:
                    info_str.append(f"{key}={value}")
            fields[7] = ';'.join(info_str)
            
            if len(fields) > 8:
                format_field = fields[8].split(':')
                gt_index = format_field.index('GT') if 'GT' in format_field else -1
                
                if gt_index >= 0:
                    for i in range(9, len(fields)):
                        sample_fields = fields[i].split(':')
                        if len(sample_fields) > gt_index:
                            gt_value = sample_fields[gt_index]
                            
                            new_gt_parts = []
                            for part in re.split(r'[/|]', gt_value):
                                if part == '.':
                                    new_gt_parts.append(part)
                                else:
                                    allele_num = int(part)
                                    if allele_num == ref_index + 1:
                                        new_gt_parts.append('0')
                                    elif allele_num > ref_index + 1:
                                        new_gt_parts.append(str(allele_num - 1))
                                    else:
                                        new_gt_parts.append(part)
                            
                            if '/' in gt_value:
                                new_gt = '/'.join(new_gt_parts)
                            else:
                                new_gt = '|'.join(new_gt_parts)
                            
                            sample_fields[gt_index] = new_gt
                            fields[i] = ':'.join(sample_fields)
            
            out.write('\t'.join(fields) + '\n')

if __name__ == '__main__':
    main()
