"""
download stride and create stride out

Sarah Hui
"""
import os
import subprocess
import requests
import tarfile


# Step 1: Download and extract STRIDE (if not already done)
def download_and_extract_stride():
    url = "https://webclu.bio.wzw.tum.de/stride/stride.tar.gz"
    tar_path = "stride.tar.gz"

    # Download STRIDE tarball
    response = requests.get(url)
    with open(tar_path, "wb") as file:
        file.write(response.content)

    print(f"File downloaded as {tar_path}")

    # Extract STRIDE tarball
    extracted_path = "stride_extracted"
    with tarfile.open(tar_path, "r:gz") as tar:
        tar.extractall(path=extracted_path)
        print(f"File extracted to {extracted_path}/")

    # Return the path to the compiled STRIDE binary
    stride_path = os.path.join(extracted_path, "stride")
    return stride_path


# Step 2: Run STRIDE on each PDB file and extract the beta-sheet regions
def process_pdb_files(input_dir, stride_path, output_dir):
    os.makedirs(output_dir, exist_ok=True)  # Ensure output directory exists

    for pdb_file in os.listdir(input_dir):
        if pdb_file.endswith(".pdb"):
            pdb_path = os.path.join(input_dir, pdb_file)
            output_file = os.path.join(output_dir, pdb_file[:-4] + '_beta_sheet.txt')

            # Run STRIDE and capture the output
            with open(output_file, 'w') as outfile:
                subprocess.run([stride_path, pdb_path], stdout=outfile)
            print(f"Processed: {pdb_file}")

            # Parse STRIDE output to extract beta-sheet regions
            extract_beta_sheet_sequences(output_file)


# Step 3: Extract beta-sheet amino acid sequences
def extract_beta_sheet_sequences(output_file):
    with open(output_file, 'r') as file:
        lines = file.readlines()

    beta_sheet_sequence = []

    for line in lines:
        # STRIDE output lines that define secondary structure typically contain 'B' for beta-sheet
        if "B" in line:
            # Extract the amino acid sequence for beta-sheet regions
            columns = line.split()
            amino_acid = columns[1]  # Assuming the amino acid is in the second column
            beta_sheet_sequence.append(amino_acid)

    if beta_sheet_sequence:
        print(f"Beta-sheet sequence: {''.join(beta_sheet_sequence)}")
    else:
        print("No beta-sheet region found.")

    return beta_sheet_sequence


# Main function
def main():
    input_dir = "/Users/tsz/Documents/School/CBB520/assignment3"
    output_dir = "stride_output"

    # Check if the STRIDE binary is available
    stride_path = download_and_extract_stride()

    # Process each PDB file in the input directory
    process_pdb_files(input_dir, stride_path, output_dir)

    print("All files processed!")


if __name__ == "__main__":
    main()


