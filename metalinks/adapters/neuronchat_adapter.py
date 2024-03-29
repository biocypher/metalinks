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

class NeuronchatEdgeType(Enum):
    """
    RECON edge types.
    """
    NC = "NC"

class NeuronchatMetaboliteToProteinEdgeField(Enum):
    """
    RECON metabolite to protein edge fields.
    """
    
    METABOLITE_ID = "METABOLITE_ID"
    PROTEIN_ID = "PROTEIN_ID"

    _PRIMARY_SOURCE_ID = 'METABOLITE_ID'
    _PRIMARY_TARGET_ID = 'PROTEIN_ID'
    _PRIMARY_REACTION_ID = 'REACTION_ID'

    MODE = 'mode'


class NeuronchatAdapter:

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
        Get edges from Neuronchat (curated file)
        """
        neuronchat_path = 'data/NeuronChatDB_human.csv'
        neuronchat_table_path = 'data/mapping_tables/Neuronchat_table.csv'

        ncdb                = pd.read_csv(neuronchat_table_path, sep=',')

        ncdb_cut            = pd.read_csv(neuronchat_path, sep=',')
        ncdb_cut['Sensor']  = ncdb_cut['interaction_name'].str.split('_').str[1]

        ncdb_dict           = dict(zip(ncdb['Query'], ncdb['HMDB']))
        ncdb_cut['Query']   = ncdb_cut['interaction_name'].str.split('_').str[0]
        ncdb_cut['HMDB']    = ncdb_cut['Query'].map(ncdb_dict)

        ncdb_cut['gene']    = ncdb_cut['interaction_name'].str.split('_').str[1]
        ncdb_cut['uniprot'] = ncdb_cut['gene'].apply(lambda x: mapping.map_name(x, 'genesymbol', 'uniprot'))


        missing_symbols = pd.read_csv('data/mapping_tables/missing_symbols.csv')
        missing_symbols_dict = dict(zip(missing_symbols['symbol'], missing_symbols['uniprot']))


        for row in ncdb_cut.iterrows():
            if row[1]['uniprot'] == set():
                row[1]['uniprot'] = [missing_symbols_dict[row[1]['gene']]]
                if row[1]['uniprot'] == 'nan':
                    print('missing uniprot for', row[1]['symbol'])
            attributes  = {
                'mode': 'activation'
            }
            id = 'uniprot:' + row[1]['uniprot'].pop()
            r = row[1].astype(str)
            h = hashlib.md5(''.join(r).encode('utf-8')).hexdigest()
            yield h, row[1]['HMDB'], id, 'NC', attributes
