'''
Read the JSON files and find the overlapping regions for each protein ID.
Compare the overlapping keys (indices) and values (amino acid codes) between the two files for each protein.
Convert the amino acid codes to their single-letter format using the provided dictionary.
Group consecutive indices and format them as per the requirements.
Output the final structure in JSON format.
'''

import json

# Load JSON data from files
with open("filtered_protein_data.json", "r") as file1, open("beta_sheet_output.json", "r") as file2:
    file1_data = json.load(file1)
    file2_data = json.load(file2)


# Amino acid mapping
amino_acid_mapping = {
    'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU',
    'F': 'PHE', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
    'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N': 'ASN',
    'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER',
    'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
}

# Reverse the mapping for easier lookup
single_letter_mapping = {v: k for k, v in amino_acid_mapping.items()}

# Find overlap and process
output_data = {}

for protein_id in file1_data:
    if protein_id in file2_data:
        overlap = {}
        
        # Find overlap in indices and identical amino acids
        for index in file1_data[protein_id]:
            if (index in file2_data[protein_id] and
                    file1_data[protein_id][index] == file2_data[protein_id][index]):
                overlap[int(index)] = single_letter_mapping[file1_data[protein_id][index]]
        
        # Group consecutive indices
        if overlap:
            sorted_indices = sorted(overlap.keys())
            sequences = []
            current_seq = overlap[sorted_indices[0]]
            
            for i in range(1, len(sorted_indices)):
                if sorted_indices[i] == sorted_indices[i - 1] + 1:
                    current_seq += overlap[sorted_indices[i]]
                else:
                    sequences.append(current_seq)
                    current_seq = overlap[sorted_indices[i]]
            
            sequences.append(current_seq)  # Add the last sequence
            output_data[protein_id] = sequences

# Save the output to a JSON file
with open("combined_protein_data.json", "w") as outfile:
    json.dump(output_data, outfile, indent=4)

print("Output saved to combined_protein_data.json")
