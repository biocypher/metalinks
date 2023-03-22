import cProfile
import io
import pstats
import sys
from biocypher import BioCypher
sys.path.append("/home/efarr/Documents/BC/CROssBAR-BioCypher-Migration/bccb/")

# set working directory
import os
os.chdir("/home/efarr/Documents/GitHub/metalinks-biocypher")

from metalinks.adapters.hmdb_adapter import (
    HMDBAdapter,
    HMDBEdgeType,
    HMDBNodeType,
    HMDBMetaboliteNodeField,
    HMDBMetaboliteToProteinEdgeField,
)

from metalinks.adapters.stitch_adapter import (
    STITCHAdapter,
    STITCHEdgeType,
    STITCHMetaboliteToProteinEdgeField,
)

from uniprot_adapter import (
    Uniprot,
    UniprotNodeType,
    UniprotNodeField,
)

from metalinks.adapters.recon_adapter import (
    ReconAdapter,
    ReconEdgeType,
    ReconMetaboliteToProteinEdgeField,
)


import biocypher

PROFILE = False

# Configure node types and fields
hmdb_node_types = [
    HMDBNodeType.METABOLITE,
]

uniprot_node_types = [
    UniprotNodeType.PROTEIN,
]

hmdb_node_fields = [
    HMDBMetaboliteNodeField._PRIMARY_ID,
    HMDBMetaboliteNodeField.METABOLITE_NAME,
    HMDBMetaboliteNodeField.METABOLITE_KEGG_ID,
    HMDBMetaboliteNodeField.METABOLITE_CHEBI_ID,
    HMDBMetaboliteNodeField.METABOLITE_PUBCHEM_ID,
    HMDBMetaboliteNodeField.METABOLITE_PROTEINS,
    HMDBMetaboliteNodeField.METABOLITE_PATHWAYS,
]

uniprot_node_fields = [
    UniprotNodeField.PROTEIN_SECONDARY_IDS,
    UniprotNodeField.PROTEIN_LENGTH,
    UniprotNodeField.PROTEIN_MASS,
    UniprotNodeField.PROTEIN_ORGANISM,
    UniprotNodeField.PROTEIN_ORGANISM_ID,
    UniprotNodeField.PROTEIN_NAMES,
    UniprotNodeField.PROTEIN_PROTEOME,
    UniprotNodeField.PROTEIN_EC,
    UniprotNodeField.PROTEIN_GENE_NAMES,
    UniprotNodeField.PROTEIN_ENSEMBL_TRANSCRIPT_IDS,
    UniprotNodeField.PROTEIN_ENSEMBL_GENE_IDS,
    UniprotNodeField.PROTEIN_ENTREZ_GENE_IDS,
    UniprotNodeField.PROTEIN_VIRUS_HOSTS,
    UniprotNodeField.PROTEIN_KEGG_IDS,

]

hmdb_edge_types = [
    HMDBEdgeType.PD,

]

stitch_edge_types = [
    STITCHEdgeType.MR,
]

recon_edge_types = [
    ReconEdgeType.PD,
]

hmdb_edge_fields = [
    HMDBMetaboliteToProteinEdgeField._PRIMARY_SOURCE_ID,
    HMDBMetaboliteToProteinEdgeField._PRIMARY_TARGET_ID,
    HMDBMetaboliteToProteinEdgeField._PRIMARY_REACTION_ID,
    HMDBMetaboliteToProteinEdgeField.SOURCE_DATABASES,
    HMDBMetaboliteToProteinEdgeField.DIRECTION,
    HMDBMetaboliteToProteinEdgeField.MET_NAME,
    HMDBMetaboliteToProteinEdgeField.STATUS,

]

stitch_edge_fields = [
    STITCHMetaboliteToProteinEdgeField._PRIMARY_SOURCE_ID,
    STITCHMetaboliteToProteinEdgeField._PRIMARY_TARGET_ID,
    STITCHMetaboliteToProteinEdgeField._PRIMARY_REACTION_ID,
    STITCHMetaboliteToProteinEdgeField.MODE,
    STITCHMetaboliteToProteinEdgeField.DATABASE,
    STITCHMetaboliteToProteinEdgeField.EXPERIMENT,
    STITCHMetaboliteToProteinEdgeField.PREDICTION,
    STITCHMetaboliteToProteinEdgeField.TEXTMINING,
    STITCHMetaboliteToProteinEdgeField.COMBINED_SCORE,

]


recon_edge_fields = [
    ReconMetaboliteToProteinEdgeField._PRIMARY_SOURCE_ID,
    ReconMetaboliteToProteinEdgeField._PRIMARY_TARGET_ID,
    ReconMetaboliteToProteinEdgeField._PRIMARY_REACTION_ID,
    ReconMetaboliteToProteinEdgeField.STATUS,
    ReconMetaboliteToProteinEdgeField.DIRECTION,
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

    bc = BioCypher(
    biocypher_config_path="config/biocypher_config.yaml",
    )
    
    # check schema
    bc.show_ontology_structure()
    

    # create adapter
    HMDB = HMDBAdapter(
        node_types=hmdb_node_types,
        node_fields=hmdb_node_fields,
        edge_types=hmdb_edge_types,
        edge_fields=hmdb_edge_fields,
        test_mode=True,
    )

    uniprot_adapter = Uniprot(
            organism="9606",
            node_types=uniprot_node_types,
            node_fields=uniprot_node_fields,
            test_mode=True,
        )
    
    uniprot_adapter.download_uniprot_data(
        cache=True,
        retries=5,
    )


    STITCH = STITCHAdapter(
        edge_types=stitch_edge_types,
        edge_fields=stitch_edge_fields,
        test_mode=True,
    )

    RECON = ReconAdapter(
        edge_types=recon_edge_types,
        edge_fields=recon_edge_fields,
        test_mode=True,
    )

    # write nodes and edges to csv
    bc.write_nodes(HMDB.get_nodes())
    bc.write_nodes(uniprot_adapter.get_nodes())
    bc.write_edges(HMDB.get_edges())
    bc.write_edges(STITCH.get_edges())
    bc.write_edges(RECON.get_edges())

    # convenience and stats
    bc.write_import_call()
    bc.log_missing_bl_types()
    bc.log_duplicates()

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
