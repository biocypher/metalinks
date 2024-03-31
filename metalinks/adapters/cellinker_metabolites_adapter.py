#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
BioCypher - recon adapter prototype
"""

from enum import Enum
from typing import Optional
import pandas as pd
import numpy as np
import hashlib
from tqdm import tqdm
import numpy as np
from pypath.utils import mapping

from biocypher._logger import logger

logger.debug(f"Loading module {__name__}.")

class CellinkerEdgeType(Enum):
    """
    Cellinker edges.
    """
    CL = "CL"

class CellinkerMetaboliteToProteinEdgeField(Enum):
    """
    Cellinker metabolite to protein edge fields.
    """
    
    METABOLITE_ID = "METABOLITE_ID"
    PROTEIN_ID = "PROTEIN_ID"

    _PRIMARY_SOURCE_ID = 'METABOLITE_ID'
    _PRIMARY_TARGET_ID = 'PROTEIN_ID'
    _PRIMARY_REACTION_ID = 'REACTION_ID'

    MODE = 'mode'


class CellinkerAdapter:

    def __init__(
        self, 
        id_batch_size: int = int(1e6),
        edge_types: Optional[list] = None,
        edge_fields: Optional[list] = None,
        test_mode: bool = False,
    ):
        self.id_batch_size = id_batch_size
        self.data_source = "Cellinker"
        self.data_source_version = ""
        self.data_license = 'None'
        self.test_mode = test_mode

    def get_edges(self):
        """
        Get edges from Cellinker (curated file)
        """
        cellinker_path = 'data/Cellinker/human-sMOL.txt'
        hmdb_mapping = 'data/mapping_tables/hmdb_mapping.csv'
        
        cellinker = pd.read_csv(cellinker_path, sep='\t')
        hmdb = pd.read_csv(hmdb_mapping, sep=',')
        
        cellinker = cellinker.dropna(subset=['ligand_pubchem_cid'])
        hmdb = hmdb.dropna(subset=['pubchem_id'])
        hmdb_dict = dict(zip(hmdb['pubchem_id'].astype(int), hmdb['accession']))

        cellinker['ligand_pubchem_cid'] = cellinker['ligand_pubchem_cid'].astype(int)
        cellinker['HMDB'] = cellinker['ligand_pubchem_cid'].map(hmdb_dict)
        cellinker.dropna(subset=['HMDB'], inplace=True)
        cellinker.dropna(subset=['Receptor_symbol'], inplace=True)
        cellinker.drop_duplicates(inplace=True)
        
        cellinker.rename(columns={'Receptor_uniprot': 'uniprot'}, inplace=True)
        
        for row in cellinker.iterrows():
            if row[1]['uniprot'] == set():
                print(row[1]['Receptor_symbol'])
                continue
            attributes  = {
                'mode': 'activation'
            }
            id = 'uniprot:' + row[1]['uniprot']
            r = row[1].astype(str)
            h = hashlib.md5(''.join(r).encode('utf-8')).hexdigest()
            yield h, row[1]['HMDB'], id, 'CL', attributes
