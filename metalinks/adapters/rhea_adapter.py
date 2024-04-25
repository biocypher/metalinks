#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
BioCypher - rhea adapter
"""

from enum import Enum
from typing import Optional
import pandas as pd
import numpy as np
import hashlib
from tqdm import tqdm
import numpy as np
import scipy.io as sio

from biocypher._logger import logger

logger.debug(f"Loading module {__name__}.")


class RheaEdgeType(Enum):
    """
    RECON edge types.
    """
    PD_rhea = "PD_rhea"

class RheaMetaboliteToProteinEdgeField(Enum):
    """
    RECON metabolite to protein edge fields.
    """
    
    METABOLITE_ID = "METABOLITE_ID"
    PROTEIN_ID = "PROTEIN_ID"
    REACTION_ID = "REACTION_ID"

    _PRIMARY_SOURCE_ID = 'METABOLITE_ID'
    _PRIMARY_TARGET_ID = 'PROTEIN_ID'
    _PRIMARY_REACTION_ID = 'REACTION_ID'

    DIRECTION = "direction"
    STATUS = "status"
    SUBSYSTEM = "subsystem"
    TRANSPORT = "transport"
    TRANSPORT_DIRECTION = 'transport_direction'
    REV = 'rev'


class RheaAdapter:

    def __init__(
        self, 
        id_batch_size: int = int(1e6),
        edge_types: Optional[list] = None,
        edge_fields: Optional[list] = None,
        test_mode: bool = False,
    ):
        self.id_batch_size = id_batch_size
        self.data_source = "RECON"
        self.data_source_version = "3D"
        self.data_license = 'None'
        self.test_mode = test_mode

    def get_edges(self):
        """
        Get edges from RECON.
        """

        rhea_rxn_path = 'data/rhea/rhea-reactions.txt'
        rhea_protein_path = 'data/rhea/rhea2uniprot_human.tsv'
        hmdb_mapping_path = 'data/mapping_tables/hmdb_mapping.csv'

        entries = parse_text_file(rhea_rxn_path)
        df = create_dataframe(entries)

        df['RHEA_ID'] = [x.split(':')[1] for x in df['ENTRY']]
        df['CHEBI_ID'] = df['CHEBI_ID'].str.split(',').str[0]

        rhea_uniprot = pd.read_csv(rhea_protein_path, sep=',')
        rhea_uniprot['RHEA_ID'] = rhea_uniprot['RHEA_ID'].astype(str)

        rhea = df.merge(rhea_uniprot, on='RHEA_ID', how='inner')
        rhea.dropna(subset=['ID'], inplace=True)
        rhea.drop_duplicates(subset=['ID', 'CHEBI_ID'], inplace=True)

        hmdb = pd.read_csv(hmdb_mapping_path, sep=',')
        hmdb.dropna(subset=['chebi_id'], inplace=True)
        hmdb_dict = dict(zip( hmdb['chebi_id'], hmdb['accession']))

        rhea['HMDB'] = rhea['CHEBI_ID'].map(hmdb_dict)
        rhea.dropna(subset=['HMDB'], inplace=True)

        df['POSITION'][df['EQUATION_SYMBOL'] == '<='] = ['left' if x == 'right' else 'right' for x in df['POSITION']]

        rhea['direction'] = rhea['POSITION'].replace({'right': 'producing', 'left': 'degrading'})

        rhea['uniprot'] = 'uniprot:' + rhea['ID']

        rhea.drop_duplicates(subset=['HMDB', 'uniprot'], inplace=True)


        for row in tqdm(rhea.iterrows()):
            attributes  = {
                'status': 'rhea',
                'direction': row[1]['direction'],
            }
            r = row[1].astype(str)
            h = hashlib.md5(''.join(r).encode('utf-8')).hexdigest()
            yield h, row[1]['HMDB'], row[1]['uniprot'], 'PD_rhea', attributes

# Function to parse the text file and extract the relevant information
def parse_text_file(file_path):
    entries = []
    with open(file_path, 'r') as file:
        entry = {}
        for line in file:
            line = line.strip()
            if line.startswith("ENTRY"):
                if entry:
                    entries.append(entry)
                entry = {}
                entry["ENTRY"] = line.split()[1]
            elif line.startswith("DEFINITION"):
                entry["DEFINITION"] = line.split(" ", 1)[1]
            elif line.startswith("EQUATION"):
                equation_parts = line.split()[1:]
                equation_string = ' '.join(equation_parts)
                equation_symbols = ['<=', '<=>', '=', '=>']
                equation_symbol = None
                for symbol in equation_symbols:
                    if symbol in equation_string:
                        equation_symbol = symbol
                        break
                if equation_symbol:
                    left, right = equation_string.split(equation_symbol)
                    left_chebis = [chebi for chebi in left.split(" ") if "CHEBI" in chebi]
                    right_chebis = [chebi for chebi in right.split(" ") if "CHEBI" in chebi]
                    entry["LEFT_CHEBI"] = left_chebis
                    entry["EQUATION_SYMBOL"] = equation_symbol
                    entry["RIGHT_CHEBI"] = right_chebis
        if entry:
            entries.append(entry)
    return entries

# Function to process entries and create DataFrame
def create_dataframe(entries):
    data = []
    for entry in entries:
        left_chebis = entry.get("LEFT_CHEBI", [])
        right_chebis = entry.get("RIGHT_CHEBI", [])
        equation_symbol = entry.get("EQUATION_SYMBOL", None)
        entry_value = entry.get("ENTRY", None)
        
        # Create separate entries for each left CHEBI ID
        for left_chebi in left_chebis:
            data.append((entry_value, left_chebi, equation_symbol, 'left'))
        
        # Create separate entries for each right CHEBI ID
        for right_chebi in right_chebis:
            data.append((entry_value, right_chebi, equation_symbol, 'right'))
            
    df = pd.DataFrame(data, columns=['ENTRY', 'CHEBI_ID', 'EQUATION_SYMBOL', 'POSITION'])
    return df