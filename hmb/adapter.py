#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
BioCypher - HMDB adapter prototype
"""

from enum import Enum
from typing import Optional

import sys
sys.path.append("/home/efarr/Documents/BC/BioCypher/")
from biocypher import *

from biocypher._logger import logger

logger.debug(f"Loading module {__name__}.")


class HMDBNodeType(Enum):
    """
    HMDB nodes.
    """

    METABOLITE = "metabolite"
    PROTEIN = "protein"


class HMDBMetaboliteNodeField(Enum):
    """
    Fields available for HMDB metabolites.
    """

    METABOLITE_ID = "metaboliteHmdbId:ID(Metabolite-ID)"
    _PRIMARY_ID = METABOLITE_ID

    METABOLITE_NAME         = "name"
    METABOLITE_KEGG_ID      = "keggId"
    METABOLITE_CHEBI_ID     = "chebiId"
    METABOLITE_PUBCHEM_ID   = "pubchemId"
    METABOLITE_PROTEINS     = "proteins"
    METABOLITE_PATHWAYS     = "keggId"
    


class HMDBProteinNodeField(Enum):
    """
    Fields available for DepMap compounds.
    """

    PROTEIN_ID = "proteinHmdbId:ID(Hmdb-ID)"
    _PRIMARY_ID = PROTEIN_ID

    PROTEIN_SYMBOL = "symbol"
    PROTEIN_UNIPROT_ID = "uniprotId"
    PROTEIN_PATHWAYS = "pathways"
    PROTEIN_METABOLITES = "metabolites"




class HMDBEdgeType(Enum):
    """
    HMDB edges.
    """

    METABOLITE_TO_PROTEIN = "metabolite_to_protein"


class HMDBMetaboliteToProteinEdgeField(Enum):
    """
    Fields available for DepMap gene to gene edges.
    """

    METABOLITE_ID = "metaboliteHmdbId:START_ID(Metabolite-ID)"
    _PRIMARY_SOURCE_ID = METABOLITE_ID

    PROTEIN_ID = "proteinHmdbId:END_ID(Protein-ID)"
    _PRIMARY_TARGET_ID = PROTEIN_ID

    SOURCE_DATABASES = "source"
    TYPE = "type"
    DIRECTION = "direction"



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

        self._set_up_types_and_fields(
            node_types, node_fields, edge_types, edge_fields
        )

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

        loc_dict = {
            HMDBNodeType.METABOLITE.value: "/home/efarr/Documents/metalinks/Data/Source/HMDB/hmdb_metabolites.xml",
            HMDBNodeType.PROTEIN.value: "/home/efarr/Documents/metalinks/Data/Source/HMDB/hmdb_proteins.xml",
        }

        ### get in the fun !!!! ###

        ### outptu format: 3 tuple

        # for label in self.node_types:
        #     # read csv for each label

        #     with (open(loc_dict[label], "r")) as f:

        #         reader = csv.reader(f)
        #         prop_items = next(reader)

        #         if self.test_mode:
        #             reader = islice(reader, 0, 100)

        #         for row in reader:
        #             _id = _process_node_id(row[0], label)
        #             _label = label
        #             _props = self._process_properties(
        #                 dict(zip(prop_items[1:], row[1:]))
        #             )
        #             yield _id, _label, _props

    def get_edges(self):
        """
        Get edges from web and yield them to the batch writer.

        Args:
            label: input label of edges to be read

        Returns:
            generator of tuples representing edges
        """

       # get in the fun !!!!

        # output has to be a 5 tuple ?


        # for label in self.edge_types:

        #     # read csv for each label
        #     with (open(loc_dict[label], "r")) as f:

        #         reader = csv.reader(f)
        #         prop_items = next(reader)

        #         if self.test_mode:
        #             reader = islice(reader, 0, 100)

        #         for row in reader:
        #             _src = self._process_source_id(row[0], label)
        #             _tar = _process_target_id(row[1], label)
        #             _label = label
        #             _props = self._process_properties(
        #                 dict(zip(prop_items[2:], row[2:]))
        #             )

        #             if not _src and _tar:
        #                 continue

        #             yield _src, _tar, _label, _props

    


# multi-line fields: only due to line 832 in cellModels_all.csv?


class HMDBHandler(xml.sax.handler.ContentHandler):
    def __init__(self):
        # ...
    
    def startElement(self, name, attrs):
        # ...

    def endElement(self, name):
        # ...

        if name == "metabolite" and self.current_accession:
            yield (self.current_accession, self.current_kegg, self.current_pubchem, self.current_chebi, 
                   self.current_met_name, self.current_inchi, self.protein_dict.get(self.current_accession, []),
                   self.pathway_dict.get(self.current_pathway, []))

            self.current_accession = ""
            self.current_kegg = ""
            self.current_pubchem = ""
            self.current_chebi = ""
            self.current_name = ""
            self.current_met_name = ""
            self.current_inchi = ""
            self.protein_dict = {}
            self.pathway_dict = {}

    def characters(self, content):
        self.current_element_text += content
