"""
find amino acid sequence by going through the stride output

"""
import os
import json


def extract_beta_sheet_from_file(file_path):
    beta_sheet_dict = {}

    with open(file_path, 'r') as file:
        for line in file:
            # Skip lines that do not start with "ASG"
            if not line.startswith("ASG"):
                continue

            # Split the line into columns
            columns = line.split()

            # Check if the region is a beta-sheet
            if "Strand" in columns[6]:  # Assuming "BetaSheet" is in the 7th column
                position = columns[4]  # Position is in the 5th column (1-based index)
                amino_acid = columns[1]  # Amino acid is in the 2nd column

                # Store the data in the dictionary
                beta_sheet_dict[position] = amino_acid

    return beta_sheet_dict


def process_all_files(input_dir, output_json_file):
    all_beta_sheets = {}

    for file_name in os.listdir(input_dir):
        if file_name.endswith("_beta_sheet.txt"):
            file_path = os.path.join(input_dir, file_name)
            beta_sheet_data = extract_beta_sheet_from_file(file_path)

            # Add the beta sheet data to the dictionary
            if beta_sheet_data:
                # Remove "_beta_sheet" from the file name for the key
                base_name = file_name.replace("_beta_sheet.txt", "")
                all_beta_sheets[base_name] = beta_sheet_data
            else:
                # Add an empty dictionary if no beta-sheet region is found
                base_name = file_name.replace("_beta_sheet.txt", "")
                all_beta_sheets[base_name] = {}

    # Write the data to a JSON file
    with open(output_json_file, 'w') as json_file:
        json.dump(all_beta_sheets, json_file, indent=4)


# Main function
def main():
    input_dir = "stride_output"  # Corrected path for stride output folder
    output_json_file = "beta_sheet_output.json"  # Output JSON file

    # Process all files and generate the output JSON
    process_all_files(input_dir, output_json_file)

    print(f"All files processed! Output saved in {output_json_file}")


if __name__ == "__main__":
    main()

