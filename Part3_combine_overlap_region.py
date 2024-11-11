import json

# Load JSON data from files
with open("filtered_protein_data.json", "r") as file1, open("beta_sheet_output.json", "r") as file2:
    json1 = json.load(file1)
    json2 = json.load(file2)

# Prepare output list
output_lines = []

# Find common protein keys and overlapping regions
for key in json1.keys() & json2.keys():
    # Find common residues
    common_residues = {k: json1[key][k] for k in json1[key] if k in json2[key] and json1[key][k] == json2[key][k]}
    
    # Concatenate values and add to output if there are common residues
    if common_residues:
        output_lines.append(f">{key}")
        output_lines.append("".join(common_residues.values()))

# Write output to a text file
with open("combined_protein_data", "w") as outfile:
    outfile.write("\n".join(output_lines))

print("Output saved to combined_protein_data")
