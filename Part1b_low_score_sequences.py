import json

# Load JSON data from files
with open("filtered_protein_data.json", "r") as file1:
    file1_data = json.load(file1)

# Amino acid mapping
amino_acid_mapping = {
    "A": "ALA",
    "C": "CYS",
    "D": "ASP",
    "E": "GLU",
    "F": "PHE",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "K": "LYS",
    "L": "LEU",
    "M": "MET",
    "N": "ASN",
    "P": "PRO",
    "Q": "GLN",
    "R": "ARG",
    "S": "SER",
    "T": "THR",
    "V": "VAL",
    "W": "TRP",
    "Y": "TYR",
}

# Reverse the mapping for easier lookup
single_letter_mapping = {v: k for k, v in amino_acid_mapping.items()}

# Find overlap and process
output_data = {}

for protein_id in file1_data:
    output = {}
    # Find overlap in indices and identical amino acids
    for index in file1_data[protein_id]:
        output[int(index)] = single_letter_mapping[file1_data[protein_id][index]]
    
    if output:
        # Group consecutive indices
        sorted_indices = sorted(output.keys())
        sequences = []
        current_seq = output[sorted_indices[0]]

        for i in range(1, len(sorted_indices)):
            if sorted_indices[i] == sorted_indices[i - 1] + 1:
                current_seq += output[sorted_indices[i]]
            else:
                sequences.append(current_seq)
                current_seq = output[sorted_indices[i]]

        sequences.append(current_seq)  # Add the last sequence
        output_data[protein_id] = sequences

# Save the output to a JSON file
with open("lowscore_protein_data.json", "w") as outfile:
    json.dump(output_data, outfile, indent=4)

print("Output saved to combined_protein_data.json")
