Assignment:

That final column of numbers in the alphafold pdb files, ranging from 0 to 100, is an estimate of the confidence of the structure around that amino acid. For low confidence regions of S. cerevisiae alphafold protein predicted structures (score <=35), what amino acid sequences (1-6 amino acids) can be found over-represented or under-represented in beta sheet regions of proteins? How do you identify over-represented and under-represented patterns?

Step 1: Download Yeast pdb files from Alphafold

Run this to download the files.
```
curl -o Yeast_pdb.tar https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000002311_559292_YEAST_v4.tar
```
Unzip them and note down the location.
```
tar -xvf Yeast_pdb.tar -C <folder_name>
```
Unfortunately, the output is once more a bunch of compressed files. To unzip, run
```
gzip -d <folder_name>/*.pdb.gz
```

Step 2: Reading pdb files to identify the low-score regions.

Run this to identify regions whose scores are lower than 35.
```
python Part1_identify_low_score.py
```
And the result will save as filtered_protein_data.json, you can also find this in the repository.

Step 3: Use Stride to find beta-sheet region. 
```
python Part2_beta_sheet_region.py and Part2_amino_acid_in_beta_region.py
```
Result is save as beta_sheet_output.json.

Step 4: Combine the results from Part1 and Part2. 
```
python Part3_combine_overlap_region.py
```
Result is save as combined_protein_data.json.