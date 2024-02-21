from os import path
import sqlite3
import pandas as pd
import numpy as np
from ast import literal_eval

def expand_list_column(df, column_name, pk='hmdb'):
    """
    Expand a list-containing column into a separate DataFrame. Drops the original column from the original DataFrame.

    Parameters:
    - df: The original DataFrame.
    - column_name: The name of the column to expand.
    - pk: The primary key column name to include in the expanded DataFrame for reference.

    Returns:
    - A DataFrame with each item of the list in a new row, along with the primary key.
    """
    # Ensure the column is in list form
    df[column_name] = df[column_name].apply(lambda x: literal_eval(x) if isinstance(x, str) else x)
    
    # Explode the column into a new DataFrame
    expanded_df = df[[pk, column_name]].explode(column_name)
    expanded_df = expanded_df.dropna().reset_index(drop=True).drop_duplicates()
    
    df.drop(column_name, axis=1, inplace=True)
    
    return expanded_df


def create_table_query(annotation):
    """
    Generates a SQL query to create a table with a foreign key reference to the 'metabolites' table,
    where the table name and the annotation column name are the same.

    Parameters:
    - annotation (str): The name of the new table and the annotation column in the new table.

    Returns:
    - A string containing the SQL command to create the table.
    """
    query = f"""
    CREATE TABLE IF NOT EXISTS {annotation} (
        hmdb TEXT,
        {annotation} TEXT,
        PRIMARY KEY (hmdb, {annotation}),
        FOREIGN KEY (hmdb) REFERENCES metabolites(hmdb)
    );
    """
    return query

# Connect to the SQLite database (or create it if it doesn't exist)
conn = sqlite3.connect(path.join('data', 'metalinks.db'))
conn.execute("PRAGMA foreign_keys = ON;")

# Load the data
## Metabolites
mets = pd.read_csv(path.join('data', 'MetaboliteTable.csv'))
# TODO: Fix this issue in the Cypher query
mets['hmdb'] = mets['hmdb'].replace(to_replace='"', value='', regex=True)
mets['metabolite'] = mets['metabolite'].replace(to_replace='"', value='', regex=True)
mets['pubchem'] = mets['pubchem'].apply(lambda x: '' if np.isnan(x) else str(int(x)))
for column in mets.columns:
    if column not in ['hmdb', 'metabolite', 'pubchem']:
        mets[column] = mets[column].apply(lambda x: literal_eval(x) if pd.notnull(x) else x)

# Create Metabolite Annotation DataFrames
columns_of_interest = ['cell_location', 'tissue_location', 'biospecimen_location', 'disease', 'pathway']
expanded_dataframes = {}

for column_name in columns_of_interest:
    df = expand_list_column(mets, column_name)
    expanded_dataframes[column_name] = df

# Proteins
prots = pd.read_csv(path.join('data', 'ProteinTable.csv'))

# TODO: Fix this issue in the Cypher query
prots['uniprot'] = prots['uniprot'].replace(to_replace='"', value='', regex=True)
prots['gene_symbol'] = prots['gene_symbol'].replace(to_replace='"', value='', regex=True)

# Metabolite-Protein Edges
edges = pd.read_csv(path.join('data', 'EdgeTable.csv'))
edges['mor'] = edges['mor'].apply(lambda x: literal_eval(x) if isinstance(x, str) else x)
# mor of length 1 then just keep that element, mor of length > 1 then 0
edges['mor'] = edges['mor'].apply(lambda x: x[0] if len(x) == 1 else 0)
edges['source'] = edges['source'].apply(lambda x: literal_eval(x) if isinstance(x, str) else x)
edges = edges.replace(to_replace='"', value='', regex=True)

# Create a DataFrame for Sources
source_expanded = edges[['hmdb', 'uniprot', 'source']].explode('source')
source_expanded = source_expanded.drop_duplicates().reset_index(drop=True)

edges.drop('source', axis=1, inplace=True) # Drop the original source column


### Create the SQL Tables
create_metabolites_table = """
CREATE TABLE IF NOT EXISTS metabolites (
    hmdb TEXT PRIMARY KEY,
    metabolite TEXT,
    pubchem TEXT,
    metabolite_subclass TEXT
);
"""

create_proteins_table = """
CREATE TABLE IF NOT EXISTS proteins (
    uniprot TEXT PRIMARY KEY,
    gene_symbol TEXT,
    protein_type TEXT
);
"""

create_edges_table = """
CREATE TABLE IF NOT EXISTS edges (
    hmdb TEXT,
    uniprot TEXT,
    db_score REAL,
    experiment_score REAL,
    combined_score REAL,
    interaction_mode TEXT,
    mor INTEGER,
    -- source TEXT,
    PRIMARY KEY (hmdb, uniprot),
    FOREIGN KEY (hmdb) REFERENCES metabolites(hmdb),
    FOREIGN KEY (uniprot) REFERENCES proteins(uniprot)
);
"""

source_query = """
CREATE TABLE source (
    hmdb VARCHAR(255),
    uniprot VARCHAR(255),
    source VARCHAR(255),
    PRIMARY KEY (hmdb, uniprot, source),
    FOREIGN KEY (hmdb, uniprot) REFERENCES edges(hmdb, uniprot)
    ); """

# Execute the SQL commands
with conn:
    conn.execute(create_metabolites_table)
    conn.execute(create_proteins_table)
    conn.execute(create_edges_table)
    conn.execute(source_query)

# Populate the tables
mets.to_sql('metabolites', conn, if_exists='append', index=False)
prots.to_sql('proteins', conn, if_exists='append', index=False)
edges.to_sql('edges', conn, if_exists='append', index=False)
source_expanded.to_sql('source', conn, if_exists='append', index=False)

# Create Annotation Queries and Tables
for key in expanded_dataframes:
    # Create the table
    query = create_table_query(key)
    with conn:
        conn.execute(query)
    # Populate
    expanded_dataframes[key].to_sql(key, conn, if_exists='append', index=False)

""" Cypher Queries used to generate the Input Data for the SQLite Database,
from a Metalinks Neo4j dump, available at https://zenodo.org/records/10200150"""
"""
//EdgeTable
MATCH (m)-[a]->(p:Protein)
WHERE 
  (type(a) IN ['StitchMetaboliteReceptor', 'NeuronchatMetaboliteReceptor', 'CellphoneMetaboliteReceptor'] AND 
  ((type(a) = 'StitchMetaboliteReceptor' AND (a.database > 500 OR a.experiment > 500 OR a.combined_score > 700)) OR
  (type(a) <> 'StitchMetaboliteReceptor')))
AND ANY(value in m.cellular_locations WHERE value = 'Extracellular')
AND NOT a.mode in ['reaction', 'catalysis', 'expression', 'pred_binding']
WITH DISTINCT m.id AS hmdb, REPLACE(p.id, "uniprot:", "") AS uniprot, a
RETURN 
  hmdb,
  uniprot,
  COLLECT(DISTINCT CASE WHEN type(a) = 'StitchMetaboliteReceptor' THEN 'Stitch'
                        WHEN type(a) = 'NeuronchatMetaboliteReceptor' THEN 'NeuronChat'
                        WHEN type(a) = 'CellphoneMetaboliteReceptor' THEN 'CellPhoneDB'
                        ELSE 'Other'
            END) AS sources,
  MAX(a.database) AS db_score,
  MAX(a.experiment) AS experiment_score,
  MAX(a.combined_score) AS combined_score,
  COLLECT(DISTINCT CASE a.mode
                    WHEN 'activation' THEN 1
                    WHEN 'inhibition' THEN -1
                    ELSE 0
            END) AS mor_modes
"""

"""
//MetaboliteTable
MATCH (m)-[a]->(p:Protein)
WHERE 
  (type(a) IN ['StitchMetaboliteReceptor', 'NeuronchatMetaboliteReceptor', 'CellphoneMetaboliteReceptor'] AND 
  ((type(a) = 'StitchMetaboliteReceptor' AND (a.database > 500 OR a.experiment > 500 OR a.combined_score > 700)) OR
  (type(a) <> 'StitchMetaboliteReceptor')))
AND ANY(value in m.cellular_locations WHERE value = 'Extracellular')
AND NOT a.mode in ['reaction', 'catalysis', 'expression', 'pred_binding']
RETURN 
  DISTINCT m.id as hmdb,
  m.name as metabolite,
  m.pubchem_compound_id as pubchem,
  m.cellular_locations as cell_location,
  m.tissue_locations as tissue_location,
  m.biospecimen_locations as biospecimen_location,
  m.sub_class as metabolite_subclass,
  m.diseases as disease,
  m.pathways as pathway
"""

"""
//ProteinTable
MATCH (m)-[a]->(p:Protein)
WHERE 
 (type(a) IN ['StitchMetaboliteReceptor', 'NeuronchatMetaboliteReceptor', 'CellphoneMetaboliteReceptor'] AND 
 ((type(a) = 'StitchMetaboliteReceptor' AND (a.database > 500 OR a.experiment > 500 OR a.combined_score > 700)) OR
 (type(a) <> 'StitchMetaboliteReceptor')))
AND NOT a.mode in ['reaction', 'catalysis', 'expression', 'pred_binding']
AND ANY(value in m.cellular_locations WHERE value = 'Extracellular')
RETURN DISTINCT
 REPLACE(p.id, "uniprot:", "") as uniprot,
 p.symbol as gene_symbol,
CASE WHEN p.receptor_type = "NA" THEN null ELSE p.receptor_type END as protein_type
"""