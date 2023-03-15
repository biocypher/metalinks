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
import hashlib
from tqdm import tqdm

from bs4 import BeautifulSoup
from biocypher._logger import logger

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
    


class HMDBProteinNodeField(Enum):
    """
    Fields available for DepMap compounds.
    """

    PROTEIN_ID = "uniprot"
    _PRIMARY_ID = PROTEIN_ID

    PROTEIN_SYMBOL = "symbol"
    PROTEIN_HMDBP_ID = "hmdbp_id"
    PROTEIN_PATHWAYS = "pathways"
    PROTEIN_METABOLITES = "metabolites"




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

        filename = "/home/efarr/Documents/metalinks/Data/Source/HMDB/hmdb_metabolites_testing.xml"
        handler = HMDBMetaboliteHandler()
        parser = xml.sax.make_parser()
        parser.setContentHandler(handler)

        with open(filename, 'r') as f:
            parser.parse(f)

        for accession, attributes in tqdm(handler.metabolites()):
            yield accession, 'hmdb_metabolite', attributes


    def get_edges(self):
        """
        Get edges from web and yield them to the batch writer.

        Args:
            label: input label of edges to be read

        Returns:
            generator of tuples representing edges
        """
        
        print(  "Getting mappings"  )

        path = '/home/efarr/Documents/metalinks/Data/Source/HMDB/hmdb_proteins.xml'

        tree = ET.parse(path)
        root = tree.getroot()
        
        accession_array = []
        uniprot_array = []
        
        for protein in tqdm(root.findall('.//{http://www.hmdb.ca}protein')):
            accession = protein.find('{http://www.hmdb.ca}accession').text
            uniprot = protein.find('{http://www.hmdb.ca}uniprot_id').text
            accession_array.append(accession)
            uniprot_array.append(uniprot)
 
        id_conversion = dict(zip(accession_array, uniprot_array))

        print(  "Getting edges"  )

        a = []

        for reaction_id in tqdm(range(1, 18)):
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
                
                if 'External Links' in reaction.text:
                    external_links = reaction.text.split('External Links')[1].split('Status')[0].strip()
                    if external_links.startswith('Kegg Reaction ID'):
                        external_links = external_links.split(':')[1].strip()
                else:
                    external_links = 'Not available'

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
                                      'reaction_id': reaction_id,
                                      'external_links': external_links,
                                      'hmdbp_id': enzyme_id,}
                        if (id, enzyme_id) not in a:
                            a.append((id, enzyme_id))
                            uniprot = id_conversion[enzyme_id]
                            edge_id = hashlib.md5((id + uniprot + 'PD' + str(attributes)).encode('utf-8')).hexdigest()                       
                            yield edge_id,  id, uniprot, 'PD', attributes
                    else:    
                        attributes = {'met_name' :  name_dict[id],
                                      'type': 'product', 
                                      'status': status, 
                                      'reaction_id': reaction_id, 
                                      'external_links': external_links,
                                      'hmdbp_id': enzyme_id,}
                        if (id, enzyme_id) not in a:
                            a.append((id, enzyme_id))
                            uniprot = id_conversion[enzyme_id]
                            edge_id = hashlib.md5((id + uniprot + 'PD' + str(attributes)).encode('utf-8')).hexdigest()
                            yield edge_id, id, uniprot, 'PD', attributes
            except:
                print(f'Could not parse reaction {reaction_id}')
                continue




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
        self._metabolites = []

    def startElement(self, name, attrs):
        self.current_element = name

    def endElement(self, name):
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
            attributes = {'kegg_id': self.current_kegg, 'pubchem_compound_id': self.current_pubchem,
                            'chebi_id': self.current_chebi, 'name': self.current_met_name, 'inchi': self.current_inchi,
                            'protein_accession': self.protein_dict.get(self.current_accession, []), 
                            'pathways': self.pathway_dict.get(self.current_pathway, [])
                            }

            self._metabolites.append((self.current_accession, attributes))

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

    def metabolites(self):
        return self._metabolites
    






