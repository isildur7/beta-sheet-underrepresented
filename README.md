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

