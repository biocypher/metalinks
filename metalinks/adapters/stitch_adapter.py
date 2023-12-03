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
        details_path = '/Users/ef6/Documents/Saez/Data/Source/Stitch/9606.protein_chemical.links.detailed.v5.0.tsv' # can be dowloaded here: https://zenodo.org/records/10200150 or the stitch website: http://stitch.embl.de/cgi/download.pl?UserId=n5QdzJfmvzSj&sessionId=5nq4BYfHQNmU
        actions_path = '/Users/ef6/Documents/Saez/Data/Source/Stitch/9606.actions.v5.0.tsv'

        actions = pl.scan_csv(actions_path, separator='\t')
        flipped = actions.filter(pl.col('item_id_a').str.contains('9606')).select(['item_id_b', 'item_id_a', 'mode'])
        df1 = flipped.collect()
        df1 = df1.rename({'item_id_a': 'protein', 'item_id_b': 'chemical', 'mode': 'mode'})

        normal = actions.filter(pl.col('item_id_b').str.contains('9606')).select(['item_id_a', 'item_id_b', 'mode'])
        df2 = normal.collect()
        df2 = df2.rename({'item_id_a': 'chemical', 'item_id_b': 'protein', 'mode': 'mode'})
        df = pl.concat([df1, df2], how='vertical')

        details = pl.scan_csv(details_path, separator='\t')
        data = details.collect()

        interactions = data.join(df, on=['chemical', 'protein'], how='left')
        del data
        del df
        interactions = interactions.filter((pl.col('combined_score') > 150) & (pl.col('mode').is_not_null()) ) # change to lower cutoff later

        interactions = interactions.with_columns(pl.col('protein').str.replace('9606.', ''))
        interactions = interactions.with_columns(pl.col('chemical').str.slice(4, None).cast(pl.Int64))

        met_table = hmdb.metabolites_table('accession',  'pubchem_compound_id')

        map_dict = dict(zip(met_table.iloc[:,1], met_table.iloc[:,0]))
        interactions = interactions.with_columns(pl.col('chemical').cast(pl.Utf8).map_dict(map_dict).alias('metabolite'))
        interactions = interactions.filter(pl.col('metabolite').is_not_null())

        order = {
            'activation': 1,
            'inhibition': 2,
            'binding': 3,
            'pred_bind': 4,
            'reaction': 5,
            'expression': 6,
            'catalysis': 7
        }

        interactions = interactions.with_columns(pl.col("mode").map_dict(order).alias("mode_order"))
        interactions = interactions.sort("mode_order").drop("mode_order")       

        reaction_id  = interactions.hash_rows(seed=42) # change to md5 later
        interactions = interactions.with_columns(reaction_id)
        interactions = interactions.rename({'': 'reaction_id'})

        uniprot_df = mapping.translation_df('ensp_biomart', 'uniprot')
        uniprot_df.rename(columns={'ensp_biomart': 'protein'}, inplace=True)

        uniprot_dict = dict(zip(uniprot_df['protein'], uniprot_df['uniprot']))

        interactions = interactions.with_columns(pl.col('protein').cast(pl.Utf8).map_dict(uniprot_dict).alias('uniprot'))
        interactions = interactions.with_columns(pl.col('uniprot').cast(pl.Utf8).apply(lambda x: 'uniprot:' + x).alias('uniprot'))

        for row in interactions.iter_rows():

            attributes = {
                'mode': row[7],
                'database': row[2],
                'experiment': row[3],
                'prediction': row[4],
                'textmining': row[5],
                'combined_score': row[6],
            }

            yield str(row[9]), row[8], row[10], 'MR', attributes # only temporary fix 
