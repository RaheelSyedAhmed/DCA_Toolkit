#!/usr/bin/env python3

import re
import argparse

class Header:
    def __init__(self, header_text: str):
        self.text = header_text

class UniProt_Header(Header):
    def __init__(self, header_text):
        super().__init__(header_text)
        # Finds UniProt Database (either Swiss-Prot/sp or TrEMBL/tr)
        self.db = re.search(r">(sp|tr)", self.text).group(1)
        # UniProt ID regex from their website (slightly modified)
        self.ID = re.search(r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z,0-9][A-Z,0-9][0-9]([A-Z][A-Z,0-9][A-Z,0-9][0-9]){0,1})", self.text).group(0)
        # Finds entry name and protein name right after the ID. 
        # Entry name may also encode a five character species name after an underscore
        post_ID = header_text.split("|")[2].split("OS=")
        self.entry_name, self.protein_name = re.search(r"(\S+)(?:\s(.+))?", post_ID[0]).groups()
        # This finds the subsection of an entry name (e.g. P0C7P0_PANTR -> PANTR, where PANTR is Pan Troglodytes)
        self.alt_species = re.search(r"_(\w{5})", self.entry_name).group(1)
        # Ensures secondary options are available. Inactive entries only have db, ID, and entry name.
        if len(post_ID) > 1:
            # Finds OS, a taxonomic identifier
            self.OS = re.search(r"OS=\S+\s", header_text).group(0).strip()
            # Finds OX, a taxonomic identifier
            self.OX = re.search(r"OX=\d+\s", header_text).group(0).strip()
            # Search for gene name. If it exists, we can add it as an attribute
            GN_search = re.search(r"GN=\S+\s", header_text)
            if GN_search:
                self.GN = GN_search.group(0).strip()
            else:
                self.GN = None
            # Search for Protein Existence tag
            self.PE = re.search(r"PE=[1-5]", header_text).group(0)
            # Search for Sequence Version tag
            self.SV = re.search(r"SV=\d+", header_text).group(0)
            

class Pfam_Header(Header):
    def __init__(self, header_text):
        super().__init__(header_text)
        # For headers directly off of Pfam's alignment page. 
        # It has the format ACCESSION_CODE
        self.ID, self.res_range = self.text.split("/")
        self.ID = re.sub(r"\.\d", "", self.ID)

class Entry:
    def __init__(self, entry_text, header_type):
        # For parsing header text and sequence text. 
        header_text = entry_text.split("\n")[0]
        sequence = "\n".join(entry_text.split("\n")[1:])
        self.header = header_type(header_text)
        self.sequence = sequence


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

