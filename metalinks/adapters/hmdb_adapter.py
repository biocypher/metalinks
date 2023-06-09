#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
BioCypher - HMDB adapter prototype
"""

from enum import Enum
from typing import Optional
import hashlib
from pypath.inputs import hmdb
from biocypher._logger import logger
from pandas import read_csv

logger.debug(f"Loading module {__name__}.")


class HMDBNodeType(Enum):
    """
    HMDB nodes.
    """

    METABOLITE = "hmdb_metabolite"
    PROTEIN = "hmdb_protein"


class HMDBMetaboliteNodeField(Enum):
    """
    Fields available for HMDB metabolites.
    """

    METABOLITE_ID = "metaboliteHmdbId:ID(Metabolite-ID)"
    _PRIMARY_ID = METABOLITE_ID

    METABOLITE_NAME         = "name"
    METABOLITE_KEGG_ID      = "kegg_id"
    METABOLITE_CHEBI_ID     = "chebiId"
    METABOLITE_PUBCHEM_ID   = "pubchemId"
    METABOLITE_PROTEINS     = "proteins"
    METABOLITE_PATHWAYS     = "pathways"
    METABOLITE_CELLULAR_LOCATIONS = "cellular_locations"
    METABOLITE_BIOSPECIMEN_LOCATIONS = "biospecimen_locations"
    METABOLITE_TISSUE_LOCATIONS = "tissue_locations"
    METABOLITE_DISEASES = "diseases"
    METABOLITE_KINGDOM = "kingdom"
    METABOLITE_CLASS = "class"
    METABOLITE_SUB_CLASS = "sub_class"
    METABOLITE_MOLECULAR_FRAMEWORK = "molecular_framework"

    

    


# class HMDBProteinNodeField(Enum):
#     """
#     Fields available for DepMap compounds.
#     """

#     PROTEIN_ID = "uniprot"
#     _PRIMARY_ID = PROTEIN_ID

#     PROTEIN_SYMBOL = "symbol"
#     PROTEIN_HMDBP_ID = "hmdbp_id"
#     PROTEIN_PATHWAYS = "pathways"
#     PROTEIN_METABOLITES = "metabolites"




class HMDBEdgeType(Enum):
    """
    HMDB edges.
    """

    PD = "PD"


class HMDBMetaboliteToProteinEdgeField(Enum):
    """
    Fields available for DepMap gene to gene edges.
    """
    REACTION_ID = "reactionId"
    _PRIMARY_REACTION_ID = REACTION_ID

    METABOLITE_ID = "metaboliteHmdbId:START_ID(Metabolite-ID)"
    _PRIMARY_SOURCE_ID = METABOLITE_ID

    PROTEIN_ID = "proteinHmdbId:END_ID(Protein-ID)"
    _PRIMARY_TARGET_ID = PROTEIN_ID

    SOURCE_DATABASES = "source"
    DIRECTION = "direction"
    MET_NAME = "met_name"
    STATUS = "status"



class HMDBAdapter:
    def __init__(
        self,
        id_batch_size: int = int(1e6),
        node_types: Optional[list] = None,
        node_fields: Optional[list] = None,
        edge_types: Optional[list] = None,
        edge_fields: Optional[list] = None,
        test_mode: bool = False,
    ):

        self.id_batch_size = id_batch_size
        self.test_mode = test_mode

        self.data_source = "HMDB"
        self.data_version = "v5.0"
        self.data_licence = "None"

    def get_nodes(self):
        """
        Get nodes from .XML and yield them to the batch writer.

        Args:
            label: input label of nodes to be read

        Returns:
            generator of tuples representing nodes
        """

        print(  "Getting metabolites"  )

        data = hmdb.hmdb_table(
            'accession',
            'kegg_id',
            'chebi_id',
            'pubchem_compound_id', 
            'name',
            'biological_properties',
            'diseases',
            'taxonomy',
            #  head = 100
        )

        data = data[data['pubchem_compound_id'] != '']
        data['cellular_locations'] = data['biological_properties'].apply(lambda x: x['cellular_locations'])
        data['biospecimen_locations'] = data['biological_properties'].apply(lambda x: x['biospecimen_locations'])
        data['tissue_locations'] = data['biological_properties'].apply(lambda x: x['tissue_locations'])
        # extract pathways name information from pathways key in biological_properties column
        data['pathways'] = data['biological_properties'].apply(lambda x: [pathway['name'] for pathway in x['pathways']])
        data['diseases'] = data['diseases'].apply(lambda x: [disease['name'] for disease in x])
        data['kingdom'] = data['taxonomy'].apply(lambda x: x['kingdom'])
        data['class'] = data['taxonomy'].apply(lambda x: x['class'])
        data['sub_class'] = data['taxonomy'].apply(lambda x: x['sub_class'])
        data['molecular_framework'] = data['taxonomy'].apply(lambda x: x['molecular_framework'])
        


        for index in range(len(data)):
            attributes = data.iloc[index, 1:].to_dict()
            yield data.iloc[index, 0], 'hmdb_metabolite', attributes





    def get_edges(self):
        """
        Get edges from web and yield them to the batch writer.

        Args:
            label: input label of edges to be read

        Returns:
            generator of tuples representing edges
        """
        
        print(  "Getting mappings"  )

        protein_mapping_path = 'data/mapping_tables/hmdb_protein_mapping.csv'
        protein_mapping = read_csv(protein_mapping_path, sep=',')
        id_conversion = dict(zip(protein_mapping['hmdbp_id'], protein_mapping['uniprot']))

        print(  "Getting edges"  )

        reactions_path = 'data/hmdb_reactions_full_status.csv'
        reactions = read_csv(reactions_path, sep=',')
        reactions['HMDBP'] = reactions['HMDBP'].apply(lambda x: id_conversion[x] if x in id_conversion else None)
        reactions.rename(columns={'HMDBP': 'uniprot'}, inplace=True)
        reactions['uniprot'] = reactions['uniprot'].apply(lambda x: 'uniprot:' + x if x is not None else None)
        reactions['reaction_id'] = reactions.apply(lambda x: hashlib.md5(str(x).encode('utf-8')).hexdigest(), axis=1)
        reactions = reactions[['reaction_id'] + [col for col in reactions.columns if col != 'reaction_id']]

        
        for index in range(len(reactions)):
            attributes = reactions.iloc[index, 3:].to_dict()
            yield reactions["reaction_id"][index], reactions['Metabolite'][index], reactions['uniprot'][index], 'PD', attributes
