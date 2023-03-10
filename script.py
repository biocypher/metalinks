import cProfile
import io
import pstats

from hmdb.adapter import (
    HMDBAdapter,
    HMDBEdgeType,
    HMDBNodeType,
    HMDBMetaboliteNodeField,
    HMDBProteinNodeField,
    HMDBMetaboliteToProteinEdgeField,
)

import biocypher

PROFILE = False

# Configure node types and fields
node_types = [
    HMDBNodeType.METABOLITE,
    HMDBNodeType.PROTEIN,
]

node_fields = [
    HMDBMetaboliteNodeField._PRIMARY_ID,
    HMDBMetaboliteNodeField.METABOLITE_NAME,
    HMDBMetaboliteNodeField.METABOLITE_KEGG_ID,
    HMDBMetaboliteNodeField.METABOLITE_CHEBI_ID,
    HMDBMetaboliteNodeField.METABOLITE_PUBCHEM_ID,
    HMDBMetaboliteNodeField.METABOLITE_PROTEINS,
    HMDBMetaboliteNodeField.METABOLITE_PATHWAYS,
    HMDBProteinNodeField._PRIMARY_ID,
    HMDBProteinNodeField.PROTEIN_SYMBOL,
    HMDBProteinNodeField.PROTEIN_HMDBP_ID,
    HMDBProteinNodeField.PROTEIN_PATHWAYS,
    HMDBProteinNodeField.PROTEIN_METABOLITES,

]


edge_types = [
    HMDBEdgeType.METABOLITE_TO_PROTEIN,

]

edge_fields = [
    HMDBMetaboliteToProteinEdgeField._PRIMARY_SOURCE_ID,
    HMDBMetaboliteToProteinEdgeField._PRIMARY_TARGET_ID,
    HMDBMetaboliteToProteinEdgeField.TYPE,
    HMDBMetaboliteToProteinEdgeField.SOURCE_DATABASES,
    HMDBMetaboliteToProteinEdgeField.DIRECTION,
    HMDBMetaboliteToProteinEdgeField.MET_NAME,
]


def main():
    """
    Connect BioCypher to HMDB adapter to import data into Neo4j.

    Optionally, run with profiling.
    """
    if PROFILE:
        profile = cProfile.Profile()
        profile.enable()

    ###############
    # ACTUAL CODE #
    ###############

    # start biocypher
    driver = biocypher.Driver(
        offline=True,
        db_name="neo4j",
        user_schema_config_path="/home/efarr/Documents/GitHub/metalinks-biocypher/config/schema_config.yaml",
        quote_char='"',
        skip_duplicate_nodes=True,
        skip_bad_relationships=True,
        strict_mode=False,
        delimiter=","
    )

    # check schema
    driver.show_ontology_structure()

    # create adapter
    HMDB = HMDBAdapter(
        node_types=node_types,
        node_fields=node_fields,
        edge_types=edge_types,
        edge_fields=edge_fields,
        test_mode=True,
    )

    # write nodes and edges to csv
    driver.write_nodes(HMDB.get_nodes())
    driver.write_edges(HMDB.get_edges())

    # convenience and stats
    driver.write_import_call()
    driver.log_missing_bl_types()
    driver.log_duplicates()

    ######################
    # END OF ACTUAL CODE #
    ######################

    if PROFILE:
        profile.disable()

        s = io.StringIO()
        sortby = pstats.SortKey.CUMULATIVE
        ps = pstats.Stats(profile, stream=s).sort_stats(sortby)
        ps.print_stats()

        ps.dump_stats("adapter.prof")
        # look at stats using snakeviz


if __name__ == "__main__":
    main()
