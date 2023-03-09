import xml.etree.ElementTree as ET
import pandas as pd


in_path = '/home/efarr/Documents/metalinks/Data/Source/HMDB/hmdb_proteins.xml'
out_path = '/home/efarr/Documents/metalinks/Data/Intermediate/HMDB/hmdb_proteins_test.csv'

tree = ET.parse(in_path)
root = tree.getroot()
accession_array = []
uniprot_array = []
gene_name_array = []
metabolite_array = []
pathway_array = []
disease_array = []

for protein in root.findall('.//{http://www.hmdb.ca}protein'):
    accession = protein.find('{http://www.hmdb.ca}accession').text
    gene_name = protein.find('{http://www.hmdb.ca}gene_name').text
    uniprot = protein.find('{http://www.hmdb.ca}uniprot_id').text
    metabolites = protein.findall('.//{http://www.hmdb.ca}metabolite_associations/')
    pathways = protein.findall('.//{http://www.hmdb.ca}pathways/')
    accession_array.append(accession)
    l = []
    s = []
    for m in metabolites:
        accession = m.find('{http://www.hmdb.ca}accession').text
        l.append(accession)
    for p in pathways:
        name = p.find('{http://www.hmdb.ca}name').text
        s.append(name)
    uniprot_array.append(uniprot)
    gene_name_array.append(gene_name)
    metabolite_array.append(l)
    pathway_array.append(s)


df = pd.DataFrame({
    "accession": accession_array,
    "gene_name": gene_name_array,
    "metabolites": metabolite_array,
    "pathways": pathway_array,
    "uniprot": uniprot_array
})

df.to_csv(out_path, index=False)
