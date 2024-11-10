
"""
Part 1
Yuyan Zhang

Reading pdb files to identify the low-score regions.
"""
import os
import json

input_folder = "/Users/zhangyuyan/Desktop/CBB520/HW/Assignment3/pdb_data" # the path of your pdb_data saved
output_file = "/Users/zhangyuyan/Desktop/CBB520/HW/Assignment3/filtered_protein_data.json" # the path of your pdb_data after filtering score less than 35
result = {}

for filename in os.listdir(input_folder):

    if filename.endswith(".pdb"):
        protein_id = filename.split(".")[0]  
        result[protein_id] = {}
        
        with open(os.path.join(input_folder, filename), "r") as file:
            for line in file:
            
                if line.startswith("ATOM"):
                    # extract useful information
                    index = int(line[22:26].strip())  # Amino Acid's position
                    amino_acid_letter = line[17:20].strip()  #Amino Acid name's position
                    score = float(line[60:66].strip())  # socre's position
                    
                    # filtering conditions
                    if score < 35 and index not in result[protein_id]:
                        result[protein_id][index] = amino_acid_letter


with open(output_file, "w") as outfile:
    json.dump(result, outfile, indent=4)

print("Finished! Result saved at:", output_file)