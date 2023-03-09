import xml.sax.handler
import pandas as pd


in_path = 'file:///home/efarr/Documents/metalinks/Data/Source/HMDB/hmdb_metabolites.xml'
test_path = 'file:///home/efarr/Documents/metalinks/Data/Source/HMDB/hmdb_metabolites_testing.xml'
out_path = '/home/efarr/Documents/metalinks/Data/Intermediate/HMDB/hmdb_metabolites_test.csv'


class HMDBHandler(xml.sax.handler.ContentHandler):
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

        self.accession_array = []
        self.kegg_array = []
        self.pubchem_array = []
        self.chebi_array = []
        self.protein_array = []
        self.name_array = []
        self.pathway_array = []
        self.disease_array = []
        self.inchi_array = []


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
        

        if name == "metabolite" and self.current_accession:#and self.current_kegg and self.current_pubchem and self.current_chebi:
            self.accession_array.append(self.current_accession)
            self.kegg_array.append(self.current_kegg)
            self.pubchem_array.append(self.current_pubchem)
            self.chebi_array.append(self.current_chebi)
            self.name_array.append(self.current_met_name)
            self.inchi_array.append(self.current_inchi)
            self.protein_array.append(self.protein_dict.get(self.current_accession, []))
            self.pathway_array.append(self.pathway_dict.get(self.current_pathway, []))

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

handler = HMDBHandler()
xml.sax.parse(test_path, handler)

df = pd.DataFrame({
    "accession": handler.accession_array,
    "kegg_id": handler.kegg_array,
    "pubchem_id": handler.pubchem_array,
    "chebi_id": handler.chebi_array,
    "name": handler.name_array,
    "inchi": handler.inchi_array,
    "protein_accession": handler.protein_array,
    "pathways": handler.pathway_array,

})

df.to_csv(out_path, index=False)