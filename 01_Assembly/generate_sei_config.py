from pathlib import Path
import yaml
from pathlib import Path
from collections import OrderedDict

# %%
from pathlib import Path
import yaml
from collections import OrderedDict

yaml_content = {}
input_dirs = {
    'sample_category_1': Path(
        '/path/to/sample_category_1/'),
    'sample_category_2': Path(
        '/path/to/sample_category_2/'),
}

workdir = Path('/SEI_SEGMENT_LENGTH_IMPACT')

yaml_dir = workdir / 'pending_jobs'
yaml_dir.mkdir(parents=True, exist_ok=True)

TARGET_LENGTH_LIST = [50, 100, 200, 300, 400, 500, 1000]
cuda_device = 0

selected_samples = ['sample_1', 'sample_2']
symbol_link_save_fasta_dir = Path('./input_fasta_dir')


# Function to generate YAML content
def generate_yaml_content(sample_name, hap, fasta_path, target_length):
    embedding_path = workdir / f'hdf5/{sample_name}.{target_length}.hdf5'
    if not embedding_path.exists():
        return {
            'embedding_hdf5_path': str(embedding_path),
            'fasta_path': str(fasta_path),
            'target_length': target_length,
            'cuda_device': cuda_device,
            'hap': hap,
            'sample_name': sample_name
        }
    else:
        raise ValueError(f"Embedding file already exists for sample {sample_name}")
    
# Processing files in each directory
for category, input_dir in input_dirs.items():
    for fasta_file in input_dir.glob('**/*.fasta'):
        sample_name = fasta_file.stem
        if any([sample_name.startswith(selected_sample) for selected_sample in selected_samples]):
            hap = fasta_file.stem.split('.')[1] if 'hap' in fasta_file.stem else "NA"
            symbol_link_fasta_path = symbol_link_save_fasta_dir / fasta_file.name
            if not symbol_link_fasta_path.exists():
                symbol_link_fasta_path.symlink_to(fasta_file)
            for target_length in TARGET_LENGTH_LIST:
                yaml_content[f'{sample_name}_{target_length}'] = generate_yaml_content(sample_name, hap, symbol_link_fasta_path, target_length)


# Write to YAML file
with open(yaml_dir / 'config.yaml', 'w') as file:
    yaml.dump(yaml_content, file, sort_keys=False)
print("High coverage YAML file created successfully.")

# Optional: Create ordered dict, split and save as different parts
yaml_content_ordered = OrderedDict(sorted(yaml_content.items(), key=lambda t: t[0]))
n = 28
for i in range(n):
    start = i * len(yaml_content_ordered) // n
    end = (i + 1) * len(yaml_content_ordered) // n
    print(f"Writing part {i} from {start} to {end} total {end - start}")
    with open(yaml_dir / f'config_{i}.yaml', 'w') as file:
        yaml.dump(dict(list(yaml_content_ordered.items())[start:end]), file, sort_keys=False)