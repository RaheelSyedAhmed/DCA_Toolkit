#!/usr/bin/env python3

import sys

if len(sys.argv) != 2:
    sys.exit("Supply only one filename.")
else:
    with open(sys.argv[1]) as f:
        text = f.read()
    print("Number of sequences:", (len(text.split(">")) - 1))
    