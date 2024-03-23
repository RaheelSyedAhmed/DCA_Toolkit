#!/usr/bin/env python3

import re
import argparse

class Header:
    def __init__(self, header_text: str):
        self.text = header_text
class UniProt_Header(Header):
    def __init__(self, header_text):
        super().__init__(header_text)
        self.db = re.search(r">(sp|tr)", self.text).group(1)
        print("DB:", self.db)
        self.ID = re.search(r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z,0-9][A-Z,0-9][0-9]([A-Z][A-Z,0-9][A-Z,0-9][0-9]){0,1})", self.text).group(0)
        print("ID:", self.ID)
        self.entry_name, self.protein_name = re.search(r"(\S+)\s(.+)", header_text.split("|")[2].split("OS=")[0]).groups()
        print("Entry Name:", self.entry_name)
        print("Protein Name:", self.protein_name)
        self.OS = re.search(r"OS=\S+\s", header_text).group(0).strip()
        print("OS:", self.OS)
        self.OX = re.search(r"OX=\d+\s", header_text).group(0).strip()
        print("OX:", self.OX)
        GN_search = re.search(r"GN=\S+\s", header_text)
        if GN_search:
            self.GN = GN_search.group(0).strip()
        else:
            self.GN = None
        print("Gene Name:", self.GN)
        self.PE = re.search(r"PE=[1-5]", header_text).group(0)
        print("Protein Existence:", self.PE)
        self.SV = re.search(r"SV=\d+", header_text).group(0)
        print("Sequence Version:", self.SV)
class Pfam_Header(Header):

    def __init__(self, header_text):
        super().__init__(header_text)
        self.ID, self.res_range = self.text.split("/")

class Entry:
    def __init__(self, entry_text, header_type):
        header_text = entry_text.split("\n")[0]
        sequence = "\n".join(entry_text.split("\n")[1:])
        self.header = header_type(header_text)
        self.sequence = sequence


def concatenate_all_domains(concat_dicts, curr_concat, key_attr, writer):
    if len(concat_dicts) == 1:
        for seq in concat_dicts[0][key_attr]:
            writer.write(ID + "\n")
            writer.write(curr_concat + seq + "\n")
    else:
        for seq in concat_dicts[0][key_attr]:
            concatenate_all_domains(concat_dicts[1:], curr_concat + seq, key_attr, writer)


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inputs", dest="input_files", help="Space-separated input files", nargs="+")
parser.add_argument("-o", "--output", dest="output_filename", help="Output filename")
args = parser.parse_args()


database_type = Pfam_Header
concat_dicts = []
for filename in args.input_files:
    with open(filename) as f:
        fasta_data = f.read()
    domain_dict = dict()
    for entry_text in fasta_data.split(">")[1:]:
        entry = Entry(">"+entry_text, database_type)
        seq_list_of_ID = domain_dict.get(entry.header.ID, [])
        seq_list_of_ID.append(entry.sequence.strip())
        domain_dict[entry.header.ID] = seq_list_of_ID
    concat_dicts.append(domain_dict)

with open(args.output_filename, 'w') as of:
    for ID in concat_dicts[0].keys():
        if False not in [(ID in domain_dict) for domain_dict in concat_dicts]:
            concatenate_all_domains(concat_dicts, "", ID, of)

