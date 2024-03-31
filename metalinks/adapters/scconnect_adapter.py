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

class ScconnectEdgeType(Enum):
    """
    scConnect edge types.
    """
    SCC = "SCC"

class ScconnectMetaboliteToProteinEdgeField(Enum):
    """
    RECON metabolite to protein edge fields.
    """
    
    METABOLITE_ID = "METABOLITE_ID"
    PROTEIN_ID = "PROTEIN_ID"

    _PRIMARY_SOURCE_ID = 'METABOLITE_ID'
    _PRIMARY_TARGET_ID = 'PROTEIN_ID'
    _PRIMARY_REACTION_ID = 'REACTION_ID'

    MODE = 'mode'


class ScconnectAdapter:

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
        Get edges from Scconnect (curated file)
        """
        scc_int_path = 'data/scConnect/interactions.csv'
        scc_lig_path = 'data/scConnect/ligands.csv'
        hmdb_path = 'data/mapping_tables/hmdb_mapping.csv'

        scconnect = pd.read_csv(scc_lig_path)
        scconnect.dropna(subset=['PubChem CID'], inplace=True)
        scconnect = scconnect[scconnect['Type'].isin(['Metabolite', 'Inorganic'])]
        interactions = pd.read_csv(scc_int_path)
        interactions = interactions[interactions['ligand'].isin(scconnect['Name'])]
        interactions = interactions[['ligand', 'target', 'target_uniprot', 'type']]
        interactions = interactions.merge(scconnect[['Name', 'PubChem CID']], left_on='ligand', right_on='Name')
        interactions.drop(columns=['Name'], inplace=True)
        interactions.drop_duplicates(inplace=True)
        interactions.dropna(subset= ['PubChem CID', 'target_uniprot'], inplace=True)
        interactions.rename(columns={'target_uniprot': 'uniprot'}, inplace=True)
        interactions = interactions.assign(uniprot=interactions['uniprot'].str.split('|')).explode('uniprot')

        interactions['type'] = interactions['type'].replace('Agonist', 'activation')
        interactions['type'] = interactions['type'].replace('Antagonist', 'inhibition')

        hmdb = pd.read_csv(hmdb_path)
        hmdb.dropna(subset=['pubchem_id'], inplace=True)
        hmdb_dict = dict(zip( hmdb['pubchem_id'].astype(int), hmdb['accession']))
        interactions['hmdb'] = interactions['PubChem CID'].astype(int).map(hmdb_dict)
        interactions.dropna(subset=['hmdb'], inplace=True)

        # missing_symbols = pd.read_csv('data/mapping_tables/missing_symbols.csv')
        # missing_symbols_dict = dict(zip(missing_symbols['symbol'], missing_symbols['uniprot']))

        for row in interactions.iterrows():
            attributes  = {
                'mode': row[1]['type']
            }
            id = 'uniprot:' + row[1]['uniprot']
            r = row[1].astype(str)
            h = hashlib.md5(''.join(r).encode('utf-8')).hexdigest()
            yield h, row[1]['hmdb'], id, 'SCC', attributes
