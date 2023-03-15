#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
BioCypher - stitch adapter prototype
"""

from enum import Enum
from typing import Optional
import pandas as pd
import numpy as np
import hashlib
from pypath.utils import mapping
from tqdm import tqdm

from biocypher._logger import logger

logger.debug(f"Loading module {__name__}.")

class STITCHEdgeType(Enum):
    """
    STITCH edge types.
    """
    MR = "MR"

class STITCHMetaboliteToProteinEdgeField(Enum):
    """
    STITCH metabolite to protein edge fields.
    """
    
    METABOLITE_ID = "METABOLITE_ID"
    PROTEIN_ID = "PROTEIN_ID"
    REACTION_ID = "REACTION_ID"

    _PRIMARY_SOURCE_ID = 'METABOLITE_ID'
    _PRIMARY_TARGET_ID = 'PROTEIN_ID'
    _PRIMARY_REACTION_ID = 'REACTION_ID'

    MODE = "mode"
    DATABASE = "database"
    EXPERIMENT = "experiment"
    PREDICTION = "prediction"
    TEXTMINING = "textmining"
    COMBINED_SCORE = "combined_score"


class STITCHAdapter:

    def __init__(
        self, 
        id_batch_size: int = int(1e6),
        edge_types: Optional[list] = None,
        edge_fields: Optional[list] = None,
        test_mode: bool = False,
    ):
        self.id_batch_size = id_batch_size
        self.data_source = "STITCH"
        self.data_source_version = "v5.0"
        self.data_license = 'None'
        self.test_mode = test_mode

    def get_edges(self):
        """
        Get edges from STITCH.
        """

        details_path = '/home/efarr/Documents/metalinks/Data/Source/Stitch/9606.protein_chemical.links.detailed.v5.0.tsv'
        actions_path = '/home/efarr/Documents/metalinks/Data/Source/Stitch/9606.actions.v5.0.tsv'
        metabolite_map_path = '/home/efarr/Documents/metalinks/Data/Intermediate/Mapping/hmdb_metabolites.csv'


        for hash, source, target, label, attributes in tqdm(details_generator(details_path, actions_path, metabolite_map_path)):
            yield hash, source, target, label, attributes

def details_generator(details_file, actions_file, metabolite_map_path ):
    actions = pd.read_csv(actions_file, sep='\t', dtype=str, usecols=['item_id_a', 'item_id_b', 'mode'])
    actions['chemical'], actions['protein'] = np.where(actions['item_id_a'].str.startswith('9606'),
                                                    (actions['item_id_b'], actions['item_id_a']),
                                                    (actions['item_id_a'], actions['item_id_b']))
    actions.drop(columns=['item_id_a', 'item_id_b'], inplace=True)

    details_chunks = pd.read_csv(details_file, sep='\t', dtype=str, chunksize=1000000)
    map_table = pd.read_csv(metabolite_map_path, sep=',', dtype=str)
    map_table.dropna(subset= ['pubchem_id'], inplace=True)
    map_dict = dict(zip(map_table['pubchem_id'].astype(int) , map_table['accession']))

    counter = 0
    for details_chunk in details_chunks:
        merged = pd.merge(details_chunk, actions, on=['chemical', 'protein'], how='left')
        merged.drop_duplicates(subset=['chemical', 'protein', 'mode'], inplace=True)
        merged['protein'] = merged['protein'].str.replace('9606.', '')
        merged['protein'] = ensp_to_genesymbol(merged['protein'])
        merged['chemical'] = merged['chemical'].str[4:].astype(int)
        merged['chemical'] = merged['chemical'].map(map_dict)
        merged.dropna(subset=['chemical'], inplace=True)
       
        # change dtype of merge to str
        merged = merged.astype(str)


        counter += 1
        if counter > 1:
            break

        for _, row in merged.iterrows():
            r = row.astype(str)
            h = hashlib.md5(''.join(r).encode('utf-8')).hexdigest()
            if pd.notna(row['mode']):
                attributes = { 'mode': row['mode'], **row[2:].to_dict()}
                yield h, row['chemical'], row['protein'], 'MR', attributes
            else:
                attributes = { **row[2:].to_dict()}
                yield h, row['chemical'], row['protein'], 'MR', attributes

            
def ensp_to_genesymbol(ensp_list):
    gene_symbol_list = []
    for element in ensp_list:
        symbol = mapping.map_name(element, 'ensp_biomart', 'uniprot')
        if symbol != set():
            gene_symbol_list.append(symbol.pop())
        else:
            symbol = mapping.map_name(element, 'ensp', 'uniprot')
            if symbol != set():
                gene_symbol_list.append(symbol.pop())
            else:
                gene_symbol_list.append('NA')
    return gene_symbol_list










