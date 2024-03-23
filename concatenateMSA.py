#!/usr/bin/env python3
from Bio import SeqIO
from scipy.spatial import distance
import pandas as pd
import numpy as np

def strtobool(val):
    val = val.lower()
    if val in ('y', 'yes', 't', 'true', 'on', '1'):
        return 1
    elif val in ('n', 'no', 'f', 'false', 'off', '0'):
        return 0
    else:
        return ValueError("invalid truth value %r" % (val,))

def parse_desc(record):
    # Obtain description, then split it into general information and parameter information by indexing our first parameter, OS
    desc = record.description
    parameter_start_index = desc.index("OS=")
    general_desc = desc[:parameter_start_index]
    parameter_part = desc[parameter_start_index:]
    
    # Add general description attribute to the record
    setattr(record, "general_desc", general_desc)


    # Keep track of parameter info
    params = []
    # Set up variable to keep words until next parameter is hit (noted by = sign)
    curr_param = ""
    for param_comp in parameter_part.split(" "):
        # If word does not have =, we're not dealing with a new parameter yet
        if "=" not in param_comp:
            curr_param += " " + param_comp
        else:
            # We've encountered a new parameter and need to reset the tracking variable
            if curr_param != "":
                params.append(curr_param)
            curr_param = param_comp
    # Final append for lagging tracking variable
    params.append(curr_param)

    # For every pair of parameter information (e.g. OX=135677), break it down into attribute and value
    for param_pair in params:
        param, val = param_pair.split("=")
        setattr(record, param, val)
    return record


# Specify threshold for hamming distance
h_dist_thresh = float(input("Specify threshold for hamming distance (0-1): "))

# Flags for requirements that may be needed
same_species_req = True
filter_entries = strtobool(input("Would you like to filter out entries?: "))

if filter_entries:
    keep_other_homologs = strtobool(input("Would you like to retain homologs instead of filtering them out?: "))

# Provide file names here
a_records = list(SeqIO.parse(input("First MSA file name: "), "fasta"))
b_records = list(SeqIO.parse(input("Second MSA file name: "), "fasta"))

# Choose reference proteins
reference_A = list(str(a_records[0].seq))
reference_B = list(str(b_records[0].seq))

# Parse the description and assign parameters as attributes
a_records = list(map(parse_desc, a_records))
b_records = list(map(parse_desc, b_records))

# Store hamming distance data (each record to its reference) in dictionaries
dist_a = dict()
for a in a_records:
    dist_a[a.seq] = distance.hamming(a, reference_A)

dist_b = dict()
for b in b_records:
    dist_b[b.seq] = distance.hamming(b, reference_B)

# Variable to store outputs for eventual writing
outputs = []

pairs = []


"""
#a_OXs = []
#b_OXs = []
#for a in a_records:
#    a_OXs.append(a.OX)
#for b in b_records:
#    b_OXs.append(b.OX)
    
# Unique organisms in a, b, then their intersection over taxonomic identifier
#print(len(set(a_OXs)))
#print(len(set(b_OXs)))
#print(len(set(a_OXs).intersection(set(b_OXs))))
"""

### Generate all possible pairs
# Compare every record in the first file to the second
for a in a_records:
    for b in b_records:
        # If same species is required, ensure OX parameter is the same for both records
        if (not same_species_req) or (a.OX == b.OX):
            # If the dissimilarity between both records in respect to their reference proteins is close, allow concatenation
            h_dist_a = dist_a[a.seq]
            h_dist_b = dist_b[b.seq]
            if(abs(h_dist_a - h_dist_b) <= h_dist_thresh):
                # Make a header with information from both records
                header = ">{} -- {}".format(a.id, b.id)
                header += ", OS1={} -- OX2={}".format(a.OS, b.OS) if (not same_species_req) else ", OS={}".format(a.OS)
                header += ", OX1={} -- OX2={}".format(a.OX, b.OX) if (not same_species_req) else ", OX={}".format(a.OX)
                header += ", GN1={} -- GN2={}".format(a.GN, b.GN) if (hasattr(a, "GN") and hasattr(b, "GN")) else ""
                header += ", dA={} -- dB={}".format(h_dist_a, h_dist_b)
                header += "\n"
                pairs.append({"OX": a.OX, 
                              "a_ID": a.id, 
                              "b_ID": b.id, 
                              "cust_header": header, 
                              "h_dist_a": h_dist_a, 
                              "h_dist_b": h_dist_b, 
                              "h_dist_diff": abs(h_dist_a - h_dist_b),
                              "concat_seq": str(a.seq) + str(b.seq) + "\n"})

if not filter_entries:
    concat_fname = "Concat_{}.fasta".format(str(h_dist_thresh).replace(".", "_"))
    with open(concat_fname, "w") as of:
        for pair in pairs:
            of.writelines(pair['cust_header'])
            of.writelines(pair['concat_seq'])
else:
    final_pairs = []
    pairs_df = pd.DataFrame(pairs)

    while(not pairs_df.empty):
        kept_pair = pairs_df.sort_values(by=['h_dist_diff']).head(1)
        final_pairs.append(kept_pair.values)
        if keep_other_homologs:
            pairs_df = pairs_df[pairs_df['a_ID'] != kept_pair.iloc[0]['a_ID']]
            pairs_df = pairs_df[pairs_df['b_ID'] != kept_pair.iloc[0]['b_ID']]
        else:
            pairs_df = pairs_df[pairs_df['OX'] != kept_pair.iloc[0]['OX']]

    concat_fname = "organism_select_concat.fasta" if keep_other_homologs == False else "homolog_select_concat.fasta"
    with open(concat_fname, 'w') as of:
        for pair in final_pairs:
            pair_arr = pair[0]
            # Field for cust_header and concat_seq
            of.writelines(pair_arr[3])
            of.writelines(pair_arr[7])

print("concatenateMSA.py has finished running.")
