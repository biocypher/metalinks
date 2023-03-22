#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
BioCypher - stitch adapter prototype
"""

from enum import Enum
from typing import Optional
from pypath.utils import mapping
from pypath.inputs import hmdb
import polars as pl

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

        print( 'Getting MR connections from STITCH... ')
              
        # prepare stitch data
        details_path = '/home/efarr/Documents/metalinks/Data/Source/Stitch/9606.protein_chemical.links.detailed.v5.0.tsv'
        actions_path = '/home/efarr/Documents/metalinks/Data/Source/Stitch/9606.actions.v5.0.tsv'

        actions = pl.scan_csv(actions_path, sep='\t')
        flipped = actions.filter(pl.col('item_id_a').str.contains('9606')).select(['item_id_b', 'item_id_a', 'mode'])
        df1 = flipped.collect()
        df1 = df1.rename({'item_id_a': 'protein', 'item_id_b': 'chemical', 'mode': 'mode'})

        normal = actions.filter(pl.col('item_id_b').str.contains('9606')).select(['item_id_a', 'item_id_b', 'mode'])
        df2 = normal.collect()
        df2 = df2.rename({'item_id_a': 'chemical', 'item_id_b': 'protein', 'mode': 'mode'})
        df = pl.concat([df1, df2], how='vertical')

        details = pl.scan_csv(details_path, sep='\t')
        data = details.collect()

        interactions = data.join(df, on=['chemical', 'protein'], how='left')
        del data
        del df
        interactions = interactions.filter((pl.col('combined_score') > 700) & (pl.col('mode').is_not_null()) ) # change to lower cutoff later

        #aesthetics
        interactions = interactions.with_columns(pl.col('protein').str.replace('9606.', ''))
        interactions = interactions.with_columns(pl.col('chemical').str.slice(4, None).cast(pl.Int64))

        # mapping
        map_dict = hmdb.hmdb_mapping('pubchem_compound_id', 'accession', head = 10)
        for k, v in map_dict.items():
            map_dict[k] = v.pop()
        interactions = interactions.with_columns(pl.col('chemical').cast(pl.Utf8).map_dict(map_dict).alias('metabolite'))
        interactions = interactions.filter(pl.col('metabolite').is_not_null())

        # add as hash per reaction
        reaction_id  = interactions.hash_rows(seed=42) # change to md5 later
        interactions = interactions.with_columns(reaction_id)
        interactions = interactions.rename({'': 'reaction_id'})

        counter = 0

        for row in interactions.iterrows():
            attributes = {'mode': row[7], 
                          'database': row[2], 
                          'experiment': row[3], 
                          'prediction': row[4], 
                          'textmining': row[5], 
                          'combined_score': row[6],
                        }
            
            uniprot = mapping.map_name(row[1], 'ensp_biomart', 'uniprot')
            if uniprot is None:
                uniprot = mapping.map_name(row[1], 'ensp', 'uniprot')
            if uniprot != set():
                uniprot =  'uniprot:' + uniprot.pop()

            else:
                continue

            counter += 1
            if counter > 100:
                break

            yield row[9], row[8], uniprot, 'MR', attributes