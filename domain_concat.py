#!/usr/bin/env python3

import argparse
from header_templates import Entry, Pfam_Header, UniProt_Header, Header


def concatenate_all_domains(concat_dicts, curr_concat, key_attr, writer):
    # Recursive scheme: repeatedly grab the first entry in concat dicts,
    # receive a domain dictionary, and concatenate every present sequence to the current concatenation
    # base case where all other files have been addressed for its sequence and concatenated to curr_concat
    # write out the concatenated sequence and the header.
    if len(concat_dicts) == 1:
        for seq in concat_dicts[0][key_attr]:
            writer.write(key_attr + "\n")
            writer.write(curr_concat + seq + "\n")
    else:
        for seq in concat_dicts[0][key_attr]:
            concatenate_all_domains(concat_dicts[1:], curr_concat + seq, key_attr, writer)


# Create a parser to get input file names as space separated values (multiple domains), output parameter to specify filename
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inputs", dest="input_files", help="Space-separated input files", nargs="+")
parser.add_argument("-o", "--output", dest="output_filename", help="Output filename")
args = parser.parse_args()

# Specify the type of header 
header_type = Pfam_Header

concat_dicts = []

for filename in args.input_files:
    # Open all input files and read in the fasta data
    with open(filename) as f:
        fasta_data = f.read()

    # Create a dictionary to track sequences with the same header accession ID (the grouping term)
    domain_dict = dict()

    # Every entry is split across > characters. Entries consist of headers and sequences.
    # e.g.
    # >sp|P0C7P0|P0C7P0_PANTR/120-130 Description of protein name OS=Pan Troglodytes OX=12345 GN=NDUFV2 PE=1 SV=1
    # GALFV ....

    for entry_text in fasta_data.split(">")[1:]:
        # Make an entry object to track header information
        entry = Entry(">"+entry_text, header_type)
        # Get the current sequences associated with that header
        seq_list_of_ID = domain_dict.get(entry.header.ID, [])
        # Add our new found sequence to it.
        seq_list_of_ID.append(entry.sequence.strip())
        # Update the dictionary accordingly
        domain_dict[entry.header.ID] = seq_list_of_ID
    # After finishing a file, add it to the list that tracks all domain_dicts
    concat_dicts.append(domain_dict)

with open(args.output_filename, 'w') as of:
    for ID in concat_dicts[0].keys():
        # If the ID (grouping term) is found in all domain files, 
        # we know there's sequences we can concatenate across all domains.
        if False not in [(ID in domain_dict) for domain_dict in concat_dicts]:
            concatenate_all_domains(concat_dicts, "", ID, of)

