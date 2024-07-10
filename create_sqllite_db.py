from os import path, remove
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
    expanded_df = df[[pk, column_name]].explode(column_name)
    expanded_df[column_name] = expanded_df[column_name].replace(to_replace='"', value='', regex=True)
    expanded_df[column_name] = expanded_df[column_name].replace(to_replace="'", value='', regex=True)
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

# Metabolite-Protein Edges (Ligand-Receptor)
edges = pd.read_csv(path.join('data', 'EdgeTable.csv'))
edges['type'] = 'lr'
# Metabolite-Protein Edges (Production-Degradation)
prod = pd.read_csv(path.join('data', 'ProductionTable.csv'))
prod['type'] = 'pd'
prod['transport_direction'] = prod['transport_direction'].replace(to_replace='"unknown"', value=np.nan)

edges = pd.concat([edges, prod], ignore_index=True)
edges['source'] = edges['source'].apply(lambda x: literal_eval(x) if isinstance(x, str) else x)
edges['mor'] = edges['mor'].apply(lambda x: literal_eval(x) if isinstance(x, str) else x)

# Create a DataFrame for Sources
edges = edges.replace(to_replace='"', value='', regex=True)
edges = edges.explode('source').explode('mor').drop_duplicates()

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
    source TEXT,
    db_score REAL,
    experiment_score REAL,
    combined_score REAL,
    mor INTEGER,
    type TEXT,
    transport_direction TEXT,
    PRIMARY KEY (hmdb, uniprot, source, mor),
    FOREIGN KEY (hmdb) REFERENCES metabolites(hmdb),
    FOREIGN KEY (uniprot) REFERENCES proteins(uniprot)
);
"""

# Execute the SQL commands
with conn:
    conn.execute(create_metabolites_table)
    conn.execute(create_proteins_table)
    conn.execute(create_edges_table)
    
# Populate the tables
try:
    mets.to_sql('metabolites', conn, if_exists='replace', index=False)
    prots.to_sql('proteins', conn, if_exists='replace', index=False)
    edges.to_sql('edges', conn, if_exists='replace', index=False)
    
    # Create Annotation Queries and Tables
    for key in expanded_dataframes:
        # Create the table
        query = create_table_query(key)
        with conn:
            conn.execute(query)
        # Populate
        expanded_dataframes[key].to_sql(key, conn, if_exists='replace', index=False)
        
except Exception as e:
    print(e)
    remove('metalinks.db')
finally:
    conn.close()
    
    




""" Cypher Queries used to generate the Input Data for the SQLite Database,
from a Metalinks Neo4j dump, available at https://zenodo.org/records/10200150"""
"""
//EdgeTable
MATCH (m)-[a]->(p:Protein)
WHERE 
type(a) IN ['CellinkerMetaboliteReceptor', 'ScconnectMetaboliteReceptor', 'StitchMetaboliteReceptor', 'NeuronchatMetaboliteReceptor', 'CellphoneMetaboliteReceptor'] 
AND ((a.database >= 200 OR a.experiment >= 300 OR a.predicted >= 700 OR a.combined_score >= 900) OR
  (type(a) <> 'StitchMetaboliteReceptor'))
AND ANY(value in m.cellular_locations WHERE value = 'Extracellular')
AND ((p.receptor_type in ['catalytic_receptor', 'gpcr', 'nhr']) OR ((p.receptor_type in ['lgic',  'other_ic', 'transporter', 'vgic'] AND a.mode in ['activation', 'inhibition'])))
AND NOT a.mode in ['reaction', 'catalysis', 'expression']
WITH DISTINCT m.id AS hmdb, REPLACE(p.id, "uniprot:", "") AS uniprot, a
RETURN 
  hmdb,
  uniprot,
  COLLECT(CASE WHEN type(a) = 'StitchMetaboliteReceptor' THEN 'Stitch'
                        WHEN type(a) = 'NeuronchatMetaboliteReceptor' THEN 'NeuronChat'
                        WHEN type(a) = 'CellphoneMetaboliteReceptor' THEN 'CellPhoneDB'
                        WHEN type(a) = 'ScconnectMetaboliteReceptor' THEN 'scConnect'
                        WHEN type(a) = 'CellinkerMetaboliteReceptor' THEN 'Cellinker'
                        ELSE 'Other'
            END) AS source,
  MAX(a.database) AS db_score,
  MAX(a.experiment) AS experiment_score,
  MAX(a.combined_score) AS combined_score,
  COLLECT(CASE a.mode
                    WHEN 'activation' THEN 1
                    WHEN 'inhibition' THEN -1
                    ELSE 0
            END) AS mor
"""

"""
//ProductionTable
MATCH (m)-[a]->(p:Protein)
WHERE type(a) in ['ReconProductionDegradation','HMDBProductionDegradation', 'HmrProductionDegradation', 'RheaProductionDegradation']  
AND ANY(value in m.cellular_locations WHERE value = 'Extracellular')  
AND NOT (a.transport_direction = 'out' AND a.direction = 'degradation')
RETURN DISTINCT 
m.id as hmdb,
REPLACE(p.id, "uniprot:", "") as uniprot,
COLLECT(CASE a.direction
              WHEN 'producing' THEN 1
              WHEN 'degrading' THEN -1
              ELSE 0
       END) as mor,
       a.transport_direction as transport_direction,
  COLLECT(a.status) AS source
"""

"""
//MetaboliteTable
MATCH (m)-[a]->(p:Protein)
WHERE 
type(a) IN ['CellinkerMetaboliteReceptor', 'ScconnectMetaboliteReceptor', 'StitchMetaboliteReceptor', 'NeuronchatMetaboliteReceptor', 'CellphoneMetaboliteReceptor'] 
AND ((a.database >= 200 OR a.experiment >= 300 OR a.predicted >= 700 OR a.combined_score >= 900) OR
  (type(a) <> 'StitchMetaboliteReceptor'))
AND ANY(value in m.cellular_locations WHERE value = 'Extracellular')
AND ((p.receptor_type in ['catalytic_receptor', 'gpcr', 'nhr']) OR ((p.receptor_type in ['lgic',  'other_ic', 'transporter', 'vgic'] AND a.mode in ['activation', 'inhibition'])))
AND NOT a.mode in ['reaction', 'catalysis', 'expression']
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
  type(a) IN ['StitchMetaboliteReceptor', 'NeuronchatMetaboliteReceptor', 'CellphoneMetaboliteReceptor']
  AND (
    (type(a) = 'StitchMetaboliteReceptor' AND 
     (a.database >= 200 OR a.experiment >= 300 OR a.predicted >= 700 OR a.combined_score >= 900)
    )
    OR type(a) <> 'StitchMetaboliteReceptor'
  )
  AND NOT a.mode IN ['reaction', 'catalysis', 'expression', 'pred_binding']
  AND ANY(value IN m.cellular_locations WHERE value = 'Extracellular')
RETURN DISTINCT
  REPLACE(p.id, "uniprot:", "") as uniprot,
  p.symbol as gene_symbol,
  CASE WHEN p.receptor_type = "NA" THEN null ELSE p.receptor_type END as protein_type
"""

