import cProfile
import io
import pstats
from biocypher import BioCypher

import os
import requests

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
    ACTIONS_PATH,
    DETAILS_PATH,
)

from metalinks.adapters.uniprot_metalinks import (
    Uniprot,
    UniprotNodeType,
    UniprotNodeField,
)

from metalinks.adapters.recon_adapter import (
    ReconAdapter,
    ReconEdgeType,
    ReconMetaboliteToProteinEdgeField,
    METMAP_PATH,
)

from metalinks.adapters.hmr_adapter import (
    HmrAdapter,
    HmrEdgeType,
    HmrMetaboliteToProteinEdgeField,
)

from metalinks.adapters.rhea_adapter import (
    RheaAdapter,
    RheaEdgeType,
    RheaMetaboliteToProteinEdgeField,
)

from metalinks.adapters.cellphone_metabolites_adapter import (
    CellphoneAdapter,
    CellphoneEdgeType,
    CellphoneMetaboliteToProteinEdgeField,
)


from metalinks.adapters.neuronchat_adapter import (
    NeuronchatAdapter,
    NeuronchatEdgeType,
    NeuronchatMetaboliteToProteinEdgeField,
)

from metalinks.adapters.cellinker_metabolites_adapter import (
    CellinkerAdapter,
    CellinkerEdgeType,
    CellinkerMetaboliteToProteinEdgeField,
)

PROFILE = False

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
    HMDBMetaboliteNodeField.METABOLITE_CELLULAR_LOCATIONS,
    HMDBMetaboliteNodeField.METABOLITE_BIOSPECIMEN_LOCATIONS,
    HMDBMetaboliteNodeField.METABOLITE_TISSUE_LOCATIONS,
    HMDBMetaboliteNodeField.METABOLITE_DISEASES,
    HMDBMetaboliteNodeField.METABOLITE_KINGDOM,
    HMDBMetaboliteNodeField.METABOLITE_CLASS,
    HMDBMetaboliteNodeField.METABOLITE_SUB_CLASS,
    HMDBMetaboliteNodeField.METABOLITE_MOLECULAR_FRAMEWORK,
]

uniprot_node_fields = [
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
    UniprotNodeField.PROTEIN_SYMBOL,
    UniprotNodeField.PROTEIN_RECEPTOR_TYPE,
    # UniprotNodeField.PROTEIN_SUBCELLULAR_LOCATION,
]

hmdb_edge_types = [
    HMDBEdgeType.PD_hmdb,
]

stitch_edge_types = [
    STITCHEdgeType.MR,
]

recon_edge_types = [
    ReconEdgeType.PD_recon,
]

hmr_edge_types = [
    HmrEdgeType.PD_hmr,
]

rhea_edge_types = [
    RheaEdgeType.PD_rhea,
]

cellphone_edge_types = [
    CellphoneEdgeType.CP,
]

neuronchat_edge_types = [
    NeuronchatEdgeType.NC,
]

cellinker_edge_types = [
    CellinkerEdgeType.CL,
]


hmdb_edge_fields = [
    HMDBMetaboliteToProteinEdgeField._PRIMARY_SOURCE_ID,
    HMDBMetaboliteToProteinEdgeField._PRIMARY_TARGET_ID,
    HMDBMetaboliteToProteinEdgeField._PRIMARY_REACTION_ID,
    HMDBMetaboliteToProteinEdgeField.SOURCE_DATABASES,
    HMDBMetaboliteToProteinEdgeField.DIRECTION,
    HMDBMetaboliteToProteinEdgeField.MET_NAME,
    HMDBMetaboliteToProteinEdgeField.STATUS,
    HMDBMetaboliteToProteinEdgeField.SUBSYSTEM,
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
    ReconMetaboliteToProteinEdgeField.SUBSYSTEM,
    ReconMetaboliteToProteinEdgeField.TRANSPORT,
    ReconMetaboliteToProteinEdgeField.TRANSPORT_DIRECTION,
    ReconMetaboliteToProteinEdgeField.REV,
]

hmr_edge_fields = [
    HmrMetaboliteToProteinEdgeField._PRIMARY_SOURCE_ID,
    HmrMetaboliteToProteinEdgeField._PRIMARY_TARGET_ID,
    HmrMetaboliteToProteinEdgeField._PRIMARY_REACTION_ID,
    HmrMetaboliteToProteinEdgeField.STATUS,
    HmrMetaboliteToProteinEdgeField.DIRECTION,
    HmrMetaboliteToProteinEdgeField.SUBSYSTEM,
    HmrMetaboliteToProteinEdgeField.TRANSPORT,
    HmrMetaboliteToProteinEdgeField.TRANSPORT_DIRECTION,
    HmrMetaboliteToProteinEdgeField.REV,
]

rhea_edge_fields = [
    RheaMetaboliteToProteinEdgeField._PRIMARY_SOURCE_ID,
    RheaMetaboliteToProteinEdgeField._PRIMARY_TARGET_ID,
    RheaMetaboliteToProteinEdgeField._PRIMARY_REACTION_ID,
    RheaMetaboliteToProteinEdgeField.DIRECTION,
    # RheaMetaboliteToProteinEdgeField.STATUS,
    # RheaMetaboliteToProteinEdgeField.SUBSYSTEM,
    # RheaMetaboliteToProteinEdgeField.TRANSPORT,
    # RheaMetaboliteToProteinEdgeField.TRANSPORT_DIRECTION,
    # RheaMetaboliteToProteinEdgeField.REV,
]

cellphone_edge_fields = [
    CellphoneMetaboliteToProteinEdgeField._PRIMARY_SOURCE_ID,
    CellphoneMetaboliteToProteinEdgeField._PRIMARY_TARGET_ID,
    CellphoneMetaboliteToProteinEdgeField._PRIMARY_REACTION_ID,
    CellphoneMetaboliteToProteinEdgeField.MODE,
]

neuronchat_edge_fields = [
    NeuronchatMetaboliteToProteinEdgeField._PRIMARY_SOURCE_ID,
    NeuronchatMetaboliteToProteinEdgeField._PRIMARY_TARGET_ID,
    NeuronchatMetaboliteToProteinEdgeField._PRIMARY_REACTION_ID,
    NeuronchatMetaboliteToProteinEdgeField.MODE,
]

cellinker_edge_fields = [
    CellinkerMetaboliteToProteinEdgeField._PRIMARY_SOURCE_ID,
    CellinkerMetaboliteToProteinEdgeField._PRIMARY_TARGET_ID,
    CellinkerMetaboliteToProteinEdgeField._PRIMARY_REACTION_ID,
    CellinkerMetaboliteToProteinEdgeField.MODE,
]


def download_files(file_mappings):
    """
    Download cached files from Zenodo and store them in the given paths.

    Parameters
        file_mappings: A dictionary where keys are URLs and values are the local file paths to store the downloaded files.
    """
    for url, path in file_mappings.items():
        # Ensure the directory exists
        os.makedirs(os.path.dirname(path), exist_ok=True)

        # Skip downloading if the file already exists
        if os.path.exists(path):
            print(f"File already exists: {path}")
            continue

        # Download the file
        response = requests.get(url, allow_redirects=True)
        if response.status_code == 200:
            with open(path, "wb") as file:
                file.write(response.content)
            print(f"Downloaded and saved: {path}")
        else:
            print(f"Failed to download {url}")


# Define the mappings of URLs to local storage paths
# TODO replace with BioCypher Resource classes (need to implement requests with parameters)
file_mappings = {
    "https://zenodo.org/records/10200150/files/9606.actions.v5.0.tsv?download=1": ACTIONS_PATH,
    "https://zenodo.org/records/10200150/files/9606.protein_chemical.links.detailed.v5.0.tsv?download=1": DETAILS_PATH,
    "https://zenodo.org/records/10200150/files/metmap_curated.csv?download=1": METMAP_PATH,
}


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

    # download cached files
    download_files(file_mappings)

    bc = BioCypher(
        biocypher_config_path="config/biocypher_config.yaml",
    )

    # check schema
    bc.show_ontology_structure()

    # create adapter
    HMDB = HMDBAdapter(
        node_types=hmdb_node_types,
        node_fields=hmdb_node_fields,
        test_mode=True,
    )

    UNIPROT = Uniprot(
        organism="9606",
        node_types=uniprot_node_types,
        node_fields=uniprot_node_fields,
        test_mode=False,
    )

    UNIPROT.download_uniprot_data(
        cache=True,
        retries=5,
    )

    STITCH = STITCHAdapter(
        edge_types=stitch_edge_types,
        edge_fields=stitch_edge_fields,
        test_mode=False,
    )

    RECON = ReconAdapter(
        edge_types=recon_edge_types,
        edge_fields=recon_edge_fields,
        test_mode=True,
    )

    HMR = HmrAdapter(
        edge_types=hmr_edge_types,
        edge_fields=hmr_edge_fields,
        test_mode=True,
    )

    RHEA = RheaAdapter(
        edge_types=rhea_edge_types,
        edge_fields=rhea_edge_fields,
        test_mode=True,
    )

    CELLPHONE = CellphoneAdapter(
        edge_types=cellphone_edge_types,
        edge_fields=cellphone_edge_fields,
        test_mode=True,
    )

    NEURONCHAT = NeuronchatAdapter(
        edge_types=neuronchat_edge_types,
        edge_fields=neuronchat_edge_fields,
        test_mode=True,
    )

    CELLINKER = CellinkerAdapter(
        edge_types=cellinker_edge_types,
        edge_fields=cellinker_edge_fields,
        test_mode=True,
    )

    # write nodes and edges to csv
    bc.write_nodes(HMDB.get_nodes())
    bc.write_edges(CELLPHONE.get_edges())
    bc.write_edges(NEURONCHAT.get_edges())
    bc.write_edges(CELLINKER.get_edges())
    bc.write_edges(STITCH.get_edges())  # high RAM, thus attention at the beginning
    bc.write_edges(RECON.get_edges())
    bc.write_edges(HMR.get_edges())
    bc.write_edges(RHEA.get_edges())
    bc.write_edges(HMDB.get_edges())
    bc.write_nodes(UNIPROT.get_nodes())


    # convenience and stats
    bc.write_import_call()
    bc.summary()

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
