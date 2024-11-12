# we should compare the sequences in the low score data to sequences in the beta sheet
# data to find outliers. For each length, we will first count how many unique motifs
# are present in the beta sheet data. Motif extraction should happen for each individual
# sequence. Then we will will randomize every sequence in the low score data 100 times,
# and count how many times each motif appears in the randomized sequences. We will then
# compare the number of times each motif appears in the original sequence to the mean
# and standard deviation of the counts in the randomized sequences. If the count of a
# motif in the original sequence is more than 3 standard deviations away from the mean,
# we will consider it an outlier.

import json

from seqgen import find_amino_acid_distribution, generate_random_sequence

# Load JSON data from files

with open("lowscore_protein_data.json", "r") as file1:
    file1_data = json.load(file1)

with open("combined_protein_data.json", "r") as file2:
    file2_data = json.load(file2)

# For length 1, we can combine all sequences in the beta sheet data into one sequence
# and simply find the amino acid distribution. We will then combine all sequences in the
# low score data into one sequence and randomize it 100 times to find the distribution
# of amino acids in the randomized sequences. We will then compare the amino acid
# distribution in the original sequence to the mean and standard deviation of the
# distributions in the randomized sequences. If the count of an amino acid in the
# original sequence is more than 3 standard deviations away from the mean, we will
# consider it an outlier.

# Combine all sequences in the beta sheet data and find amino acid distribution
beta_sheet_sequence = "".join(
    [sequence for protein_id in file2_data for sequence in file2_data[protein_id]]
)
beta_sheet_amino_acids, beta_sheet_probs = find_amino_acid_distribution(
    beta_sheet_sequence
)

# # Combine all sequences in the beta sheet data and find amino acid distribution in low score data
low_score_sequence = "".join(
    [sequence for protein_id in file1_data for sequence in file1_data[protein_id]]
)
low_score_amino_acids, low_score_probs = find_amino_acid_distribution(
    low_score_sequence
)

# Randomize the low score sequence 100 times
random_sequences = [
    generate_random_sequence(
        len(low_score_sequence), low_score_amino_acids, low_score_probs
    )
    for _ in range(100)
]

# Find amino acid distribution in the randomized sequences
random_amino_acids = [
    find_amino_acid_distribution(sequence) for sequence in random_sequences
]

# Find mean and standard deviation for each amino acid in the randomized sequences
random_means = [
    sum(probs) / len(random_sequences) for amino_acids, probs in random_amino_acids
]
random_stds = [
    (sum([(prob - mean) ** 2 for prob in probs]) / len(random_sequences)) ** 0.5
    for mean, (amino_acids, probs) in zip(random_means, random_amino_acids)
]

# Find outliers in the beta sheet data compared to the low score data
outliers = {
    amino_acid: prob
    for amino_acid, prob, mean, std in zip(
        beta_sheet_amino_acids, beta_sheet_probs, random_means, random_stds
    )
    if prob > mean + 3 * std
}
# also store the mean and std for each amino acid in a different dictionary per amino acid
outliers_mean_std = {
    amino_acid: (mean, std)
    for amino_acid, mean, std in zip(beta_sheet_amino_acids, random_means, random_stds)
}

# store the outliers and their mean and std in a json file
with open("outliers_len1.json", "w") as outfile:
    json.dump(outliers, outfile)

with open("outliers_mean_std_len1.json", "w") as outfile:
    json.dump(outliers_mean_std, outfile)

# For length 2, we will find all unique two letter motifs in the beta sheet data by
# processing each sequence individually. We will ignore all sequences with length less
# than 2. We will similarly find all unique two letter motifs in the low score data by
# ignoring all sequences with length less than 2. We will then randomize each sequence
# 100 times and count how many times each motif appears in the randomized sequences. We
# will then compare the number of times each motif appears in the original sequence to
# the mean and standard deviation of the counts in the randomized sequences. If the count
# of a motif in the original sequence is more than 3 standard deviations away from the
# mean, we will consider it an outlier.
