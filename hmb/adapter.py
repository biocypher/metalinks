#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
BioCypher - HMDB adapter prototype
"""

from enum import Enum
from typing import Optional
import requests
import re
import xml.sax.handler
import xml.etree.ElementTree as ET


from bs4 import BeautifulSoup


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
    METABOLITE_PATHWAYS     = "pathways"
    


class HMDBProteinNodeField(Enum):
    """
    Fields available for DepMap compounds.
    """

    PROTEIN_ID = "uniprot_id"
    _PRIMARY_ID = PROTEIN_ID

    PROTEIN_SYMBOL = "symbol"
    HMDBP_ID = "hmdbp_id"
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
    MET_NAME = "met_name"



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

        # self._set_up_types_and_fields(
        #     node_types, node_fields, edge_types, edge_fields
        # )

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
        # loc_dict = {
        #     HMDBNodeType.METABOLITE.value: "/home/efarr/Documents/metalinks/Data/Source/HMDB/hmdb_metabolites_testing.xml",
        #     HMDBNodeType.PROTEIN.value: "/home/efarr/Documents/metalinks/Data/Source/HMDB/hmdb_proteins.xml",
        # }

        print(  "Getting metabolites"  )
        # handler = HMDBMetaboliteHandler()
        # xml.sax.parse("file:///home/efarr/Documents/metalinks/Data/Source/HMDB/hmdb_metabolites_testing.xml", handler)

        print(  "Getting proteins"  )

        ProteinDataExtractor('/home/efarr/Documents/metalinks/Data/Source/HMDB/hmdb_proteins.xml')

    def get_edges(self):
        """
        Get edges from web and yield them to the batch writer.

        Args:
            label: input label of edges to be read

        Returns:
            generator of tuples representing edges
        """
        print(  "Getting edges"  )
        for reaction_id in range(1, 14):
            # specify the HMDB webfile URL to be scraped
            hmdb_url = 'https://hmdb.ca/reactions/' + str(reaction_id)

            # send a GET request to the HMDB webfile URL and store the response
            response = requests.get(hmdb_url)

            
            try: 
                # parse the html using beautiful soup
                soup = BeautifulSoup(response.text, 'html.parser')

                # find all the reaction panels
                reaction = soup.find(class_='reaction-panel')

                # extract all the metabolite ids
                metabolite_ids = re.findall(r'/metabolites/(HMDB\d+)', str(reaction))

                status = reaction.text.split('Status')[1].strip()
                status = status.split(' ')[0]


                enzyme_id = re.search(r'/proteins/(HMDBP\d+)', str(reaction)).group(1)

                reaction_str = soup.find(class_='panel-heading').text
                # split reaction string by either + or = and count how many object were before the =
                reactands, products = reaction_str.split('=')
                reactands = reactands.split('+')
                products = products.split('+')

                reactand_ids = metabolite_ids[:len(reactands)]

                #make dictionary of metabolite ids and reactands and products by concatenating the reactand and product lists and then use them as value
                names = reactands + products
                name_dict = dict(zip(metabolite_ids, names))    

                # for every id in the metabolite_ids list, check if it is a reactand or product and yield the appropriate edge
                for id in metabolite_ids:
                    if id in reactand_ids:
                        attributes = {'met_name' :  name_dict[id],
                                      'type': 'reactand', 
                                      'status': status, 
                                      'reaction_id': reaction_id}
                        print(id, enzyme_id, 'reactand', attributes)
                        yield (id, enzyme_id, 'reactand', attributes)
                    else:    
                        attributes = {'met_name' :  name_dict[id],
                                      'type': 'product', 
                                      'status': status, 
                                      'reaction_id': reaction_id}
                        print(id, enzyme_id, 'product', attributes)
                        yield (id, enzyme_id, 'product', attributes)




            except:
                print(f'Could not parse reaction {reaction_id}')
                continue
    


# multi-line fields: only due to line 832 in cellModels_all.csv?


class HMDBMetaboliteHandler(xml.sax.handler.ContentHandler):
    def __init__(self):
        self.current_element = ""
        self.current_element_text = ""
        self.current_accession = ""
        self.current_kegg = ""
        self.current_pubchem = ""
        self.current_chebi = ""
        self.current_met_name = ""
        self.current_name = ""
        self.current_pathway = ""
        self.current_pathways = ""
        self.current_disease = ""
        self.current_inchi = ""
        self.protein_dict = {}
        self.pathway_dict = {}
        self.disease_dict = {}


    def startElement(self, name, attrs):
        self.current_element = name

    def endElement(self, name):
        print(name)
        if name == "accession" and not self.current_accession:
            self.current_accession = self.current_element_text.strip()
            self.current_pathways = ""
        elif name == "kegg_id":
            self.current_kegg = self.current_element_text.strip()
        elif name == "inchi":
            self.current_inchi = self.current_element_text.strip()
        elif name == "pubchem_compound_id":
            self.current_pubchem = self.current_element_text.strip()
        elif name == "chebi_id":
            self.current_chebi = self.current_element_text.strip()
        elif name == "name" and not self.current_met_name:
            self.current_met_name = self.current_element_text.strip()
        elif name == "protein_accession":
            protein_accession = self.current_element_text.strip()
            if protein_accession:
                self.protein_dict.setdefault(self.current_accession, []).append(protein_accession)
        elif name == "name" and not self.current_pathways == "pathways":
            pathway_name = self.current_element_text.strip()
            if pathway_name:
                self.pathway_dict.setdefault(self.current_pathway, []).append(pathway_name)
        elif name == "pathways":
            self.current_pathways = "pathways"

        self.current_element_text = ""
        self.current_element = ""
        

        if name == "metabolite" and self.current_accession:
            print(self.current_accession)
            attributes = {'kegg_id': self.current_kegg, 'pubchem_compound_id': self.current_pubchem,
                            'chebi_id': self.current_chebi, 'name': self.current_met_name, 'inchi': self.current_inchi,
                            'protein_accession': self.protein_dict.get(self.current_accession, []), 
                            'pathways': self.pathway_dict.get(self.current_accession, [])
                            }
            yield (self.current_accession, self.current_name, attributes)


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
        print(content)
        self.current_element_text += content





# def extract_protein_data(path):
#     print('heloo')
#     tree = ET.parse(path)
#     root = tree.getroot()
#     print('dadidu')
#     for protein in root.findall('.//{http://www.hmdb.ca}protein'):
#         print('dada')
#         accession = protein.find('{http://www.hmdb.ca}accession').text
#         gene_name = protein.find('{http://www.hmdb.ca}gene_name').text
#         uniprot_id = protein.find('{http://www.hmdb.ca}uniprot_id').text
#         metabolites = protein.findall('.//{http://www.hmdb.ca}metabolite_associations/')
#         pathways = protein.findall('.//{http://www.hmdb.ca}pathways/')
#         metabolite_array = []
#         pathway_array = []
#         for m in metabolites:
#             accession = m.find('{http://www.hmdb.ca}accession').text
#             metabolite_array.append(accession)
#         for p in pathways:
#             name = p.find('{http://www.hmdb.ca}name').text
#             pathway_array.append(name)
#         attributes = {'metabolites': metabolite_array, 
#                       'pathways': pathway_array, 
#                       'accession': accession,}
#         yield (uniprot_id, gene_name, attributes)

class ProteinDataExtractor:
    def __init__(self, path):
        self.path = path
        self.tree = None
        self.root = None
        
    def __iter__(self):
        self.tree = ET.parse(self.path)
        self.root = self.tree.getroot()
        return self

    def __next__(self):
        for protein in self.root.findall('.//{http://www.hmdb.ca}protein'):
            print('dada')
            accession = protein.find('{http://www.hmdb.ca}accession').text
            gene_name = protein.find('{http://www.hmdb.ca}gene_name').text
            uniprot_id = protein.find('{http://www.hmdb.ca}uniprot_id').text
            metabolites = protein.findall('.//{http://www.hmdb.ca}metabolite_associations/')
            pathways = protein.findall('.//{http://www.hmdb.ca}pathways/')
            metabolite_array = []
            pathway_array = []
            for m in metabolites:
                accession = m.find('{http://www.hmdb.ca}accession').text
                metabolite_array.append(accession)
            for p in pathways:
                name = p.find('{http://www.hmdb.ca}name').text
                pathway_array.append(name)
            attributes = {'metabolites': metabolite_array, 
                          'pathways': pathway_array, 
                          'accession': accession,}
            yield (uniprot_id, gene_name, attributes)
        else:
            raise StopIteration
