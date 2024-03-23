#!/usr/bin/env python3
import numpy as np
import math

input_file = input("Input File: ")
output_file = input("Output File: ")

# Open MSA
headers = []
with open(input_file, 'r') as f:
    # Obtain all line information, 
    # read the second line for length of each alignment sequence (should be the same), 
    # count the number of lines and divide by 2 to get the number of sequences in the fasta
    lines = f.read().splitlines()
    nrows, ncols = (int(len(lines)/2), len(lines[1]))
    # Initiate an empty np array to store information
    all_seqs = np.empty(shape=(nrows,ncols), dtype=np.str_)

    for index, line in enumerate(lines):
        if '>' in line:
            # Pass lines that are headers
            headers.append(line)
        else:
            # Push list form of sequence into np array at correct index
            all_seqs[math.floor(index/2)] = list(line)

# Track columns we want to delete
cols_to_remove = []
for col in range(ncols):
    # Check to see if there's only gaps
    remove_col = True
    for site in all_seqs[:, col]:
        if site != "-":
            remove_col = False
            break
    if remove_col:
        cols_to_remove.append(col)


# Delete all columns marked as gap-only
cut_cols = np.delete(all_seqs, cols_to_remove, axis=1)

#print(cut_cols.shape)
# Write out output
with open(output_file, "w") as of:
    for index, row in enumerate(cut_cols):
        of.writelines(headers[index])
        of.write("\n")
        of.writelines(''.join(row))
        of.write("\n")
            
