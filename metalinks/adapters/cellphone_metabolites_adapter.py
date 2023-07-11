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

class CellphoneEdgeType(Enum):
    """
    RECON edge types.
    """
    CP = "CP"

class CellphoneMetaboliteToProteinEdgeField(Enum):
    """
    RECON metabolite to protein edge fields.
    """
    
    METABOLITE_ID = "METABOLITE_ID"
    PROTEIN_ID = "PROTEIN_ID"

    _PRIMARY_SOURCE_ID = 'METABOLITE_ID'
    _PRIMARY_TARGET_ID = 'PROTEIN_ID'
    _PRIMARY_REACTION_ID = 'REACTION_ID'

    MODE = 'mode'


class CellphoneAdapter:

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
        Get edges from Cellphone (curated file)
        """

        cpdb = pd.read_excel('data/Cellphone_suptab4_curated.xlsx')
        cpdb['symbol'] = cpdb['protein_name_b'].str.split('_').str[0]
        cpdb['uniprot'] = cpdb['symbol'].apply(lambda x: mapping.map_name(x, 'genesymbol', 'uniprot'))

        for row in cpdb.iterrows():
            if row[1]['uniprot'] == set():
                print(row)
                continue
            attributes  = {
                'mode': 'activation'
            }
            id = 'uniprot:' + row[1]['uniprot'].pop()
            r = row[1].astype(str)
            h = hashlib.md5(''.join(r).encode('utf-8')).hexdigest()
            yield h, row[1]['Chebi/HMDB_name_a'], id, 'CP', attributes
