#!/usr/bin/env python3
import numpy as np
import math

input_file = input("Input File: ")
output_file = input("Output File: ")
start_col = int(input("Start column to truncate from (inclusive): "))
end_col = int(input("Final column to truncate to (exclusive): "))

# Open MSA
headers = []
with open(input_file, 'r') as f:
    # Obtain all line information, 
    # read the second line for length of each alignment sequence (should be the same), 
    # count the number of lines and divide by 2 to get the number of sequences in the fasta
    lines = f.readlines()
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


cols_to_remove = np.arange(start_col, end_col)
cut_cols = np.delete(all_seqs, cols_to_remove, axis=1)


with open(output_file, "w") as of:
    for index, row in enumerate(cut_cols):
        of.writelines(headers[index])
        of.writelines(''.join(row))