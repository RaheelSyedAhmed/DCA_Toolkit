#!/usr/bin/env python3
import plotly.express as px
import pandas as pd
import numpy as np
import sys
import re

MSA_fname = sys.argv[1]
count_dict = dict()

with open(MSA_fname, 'r') as MSA_fs:
    data = MSA_fs.read().splitlines()
# From the second element onwards, skip every other element.
# This skips all of the headers.
seqs = data[1::2]

for seq in seqs:
    dash_match = re.findall(r"-+", seq)
    gap_counts = [len(x) for x in dash_match]
    if len(gap_counts) > 0:
        max_continuous_gaps = max(gap_counts)
        count_dict[max_continuous_gaps] = count_dict.get(max_continuous_gaps, 0) + 1
    else:
        count_dict[0] = count_dict.get(0, 0) + 1

sorted_count_dict = sorted(count_dict.items())
cumul_count = []
curr_count = 0
for key, value in sorted_count_dict:
    curr_count += value
    cumul_count.append(curr_count)
perc_retained = [round(x/cumul_count[-1], 4) for x in cumul_count]


count_df = pd.DataFrame(sorted_count_dict, columns =['num_cont_gaps', 'count'])
count_df['percent_retained'] = perc_retained
overall_largest_num_cont_gaps = max(count_dict.keys())
fig = px.histogram(count_df, x="num_cont_gaps", y="count", nbins=overall_largest_num_cont_gaps + 1)
fig.show()

fig2 = px.line(count_df, x="num_cont_gaps", y="percent_retained", markers=True)
fig2.show()