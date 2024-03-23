#!/usr/bin/env python3
#Combine a sequence in MSA in one single line

# Define input and output file paths
input_file = input("Enter the input MSA file: ")
output_file = input("Enter the output MSA file: ")

# Open the input and output files
with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    sequence = []  # To store the sequence lines
    for line in infile:
        line = line.strip()
        if line.startswith(">"):
            # If it's a header line, write the previous sequence lines and then the header
            if sequence:
                outfile.write("".join(sequence) + "\n")
                sequence = []
            outfile.write(line + "\n")
        else:
            # If it's a sequence line, append it to the sequence list
            sequence.append(line)
    # Write the last sequence after the loop
    if sequence:
        outfile.write("".join(sequence) + "\n")

print(f"Output written to {output_file}")