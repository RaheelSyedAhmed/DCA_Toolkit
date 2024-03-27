import re

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
        sequence = entry_text.split("\n")[1]
        self.header = header_type(header_text)
        self.sequence = sequence