# metalinks-biocypher 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10200150.svg)](https://doi.org/10.5281/zenodo.10200150)

This is the Repo for the Metalinks database, which collects and curates
knowledge around metabolite-protein interactions. You can query our dataset in
our web app here: [MetalinksDB webpage](https://metalinks.omnipathdb.org/)

## Reproducing Metalinks

1. Clone this repository.

```bash
git clone https://github.com/biocypher/metalinks.git
cd metalinks
```

2. Install the required packages.

```bash
poetry install
```

or

```bash 
conda env create -f env.yml
conda activate metalinks
```

To create the knowledge graph, run the following within the created environment:
```bash 
python create_knowledge_graph.py
```

<!-- TODO rest -->

All additional files and a webpage dump can be found here:
[zenodo](https://zenodo.org/records/10200150)

## Simple MetalinksDB Access / Generating an SQL Database

To generate a sqlite3 database from the knowledge graph, run the following within the created environment:

1. Load a neo4j dump using a neo4j instance following running `create_knowledge_graph.py`.
Alternatively, a dump of the database used in the paper is available here https://zenodo.org/records/11066196, preferably use neo4j v4.4.32 if loading this dump.
2. Export the knowledge graph to a csv files, see the queries in create_sqlite_db.py
3. Run the create_sqlite_db.py script to generate the sqlite3 database.

We also provide MetalinkDB as csv files: https://github.com/saezlab/MetalinksDB/tree/main/metalinksDB

MetalinksDB is also queryable via LIANA+: https://liana-py.readthedocs.io/en/latest/notebooks/prior_knowledge.html#Metabolite-Receptor-Interactions