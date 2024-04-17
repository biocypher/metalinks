# adapted from crossbar project https://github.com/HUBioDataLab/CROssBAR-BioCypher-Migration

from time import time
from typing import Optional
from enum import Enum, auto
from functools import lru_cache
import pandas as pd
import re

from tqdm import tqdm  # progress bar
from pypath.share import curl, settings
from pypath.utils import mapping
from pypath.inputs import uniprot
from biocypher._logger import logger
from contextlib import ExitStack
from bioregistry import normalize_curie

logger.debug(f"Loading module {__name__}.")


class UniprotNodeType(Enum):
    """
    Node types of the UniProt API represented in this adapter.
    """

    PROTEIN = auto()
    GENE = auto()
    ORGANISM = auto()
    CELLULAR_COMPARTMENT = auto()


class UniprotNodeField(Enum):
    """
    Fields of nodes the UniProt API represented in this adapter. Overview of
    uniprot fields: https://www.uniprot.org/help/return_fields
    """

    # core attributes
    PROTEIN_LENGTH = "length"
    PROTEIN_SUBCELLULAR_LOCATION = "subcellular_location"
    PROTEIN_MASS = "mass"
    PROTEIN_ORGANISM = "organism_name"
    PROTEIN_ORGANISM_ID = "organism_id"
    PROTEIN_NAMES = "protein_name"
    PROTEIN_EC = "ec"
    PROTEIN_GENE_NAMES = "gene_names"
    PRIMARY_GENE_NAME = "gene_primary"
    PROTEIN_CC_DISEASE = "cc_disease"
    
    # xref attributes
    PROTEIN_ENSEMBL_TRANSCRIPT_IDS = "xref_ensembl"
    PROTEIN_PROTEOME = "xref_proteomes"
    PROTEIN_ENTREZ_GENE_IDS = "xref_geneid"
    PROTEIN_VIRUS_HOSTS = "virus_hosts"
    PROTEIN_KEGG_IDS = "xref_kegg"

    # not from uniprot REST
    # we provide these by mapping ENSTs via pypath
    PROTEIN_ENSEMBL_GENE_IDS = "ensembl_gene_ids"

    PROTEIN_SYMBOL = "symbol"
    PROTEIN_RECEPTOR_TYPE = "receptor_type"


class Uniprot:
    """
    Class that downloads uniprot data using pypath and reformats it to be ready
    for import into a BioCypher database.

    Args:
        organism: organism code in NCBI taxid format, e.g. "9606" for human.

        rev: if True, it downloads reviewed entries only.
    """

    def __init__(
        self,
        organism="*",
        rev=True,
        node_types: Optional[list] = None,
        node_fields: Optional[list] = None,
        normalise_curies: bool = True,
        test_mode: bool = False,
    ):

        # params
        self.organism = organism
        self.rev = rev
        self.normalise_curies = normalise_curies
        self.test_mode = test_mode

        # provenance
        self.data_source = "uniprot"
        self.data_version = "2022_04"  # TODO get version from pypath
        self.data_licence = "CC BY 4.0"

        self._configure_fields()

        self._set_node_and_edge_fields(
            node_types=node_types,
            node_fields=node_fields,
        )

        # loading of subcellular locations set
        self.locations = set()

    def download_uniprot_data(
        self,
        cache=False,
        debug=False,
        retries=3,
    ):
        """
        Wrapper function to download uniprot data using pypath; used to access
        settings.

        Args:
            cache: if True, it uses the cached version of the data, otherwise
            forces download.

            debug: if True, turns on debug mode in pypath.

            retries: number of retries in case of download error.
        """

        # stack pypath context managers
        with ExitStack() as stack:

            stack.enter_context(settings.context(retries=retries))

            if debug:
                stack.enter_context(curl.debug_on())

            if not cache:
                stack.enter_context(curl.cache_off())

            self._download_uniprot_data()

            # preprocess data
            self._preprocess_uniprot_data()

    def _download_uniprot_data(self):
        """
        Download uniprot data from uniprot.org through pypath.

        Here is an overview of uniprot return fields:
        https://www.uniprot.org/help/return_fields

        TODO make use of multi-field query
        """

        logger.info("Downloading uniprot data...")

        t0 = time()

        # download all swissprot ids
        self.uniprot_ids = list(uniprot._all_uniprots(self.organism, self.rev))

        # limit to 100 for testing
        if self.test_mode:
            self.uniprot_ids = self.uniprot_ids[:100]

        # download attribute dicts
        self.data = {}
        for query_key in tqdm(self.node_fields):
            if query_key in [
                UniprotNodeField.PROTEIN_ENSEMBL_GENE_IDS.value,
                UniprotNodeField.PROTEIN_RECEPTOR_TYPE.value,
                "symbol",
            ]:
                continue

            elif query_key == UniprotNodeField.PROTEIN_SUBCELLULAR_LOCATION.value:
                self.data[query_key] = uniprot.uniprot_locations(
                    self.organism, self.rev
                )
            else:
                self.data[query_key] = uniprot.uniprot_data(
                    query_key, self.organism, self.rev
                )

            logger.debug(f"{query_key} field is downloaded")

        # add ensembl gene ids
        self.data[UniprotNodeField.PROTEIN_ENSEMBL_GENE_IDS.value] = {}

        t1 = time()
        msg = f"Acquired UniProt data in {round((t1-t0) / 60, 2)} mins."
        logger.info(msg)

    def _preprocess_uniprot_data(self):
        """
        Preprocess uniprot data to make it ready for import. First, three types
        of processing are applied:
        - nothing is done (for ensembl gene ids, which come from pypath)
        - simple string replacement
        - replace separators in integers and convert to int
        - field splitting

        Then, special treatment is applied to some fields:
        - ensg ids are extracted from the ensembl transcript ids
        - protein names and virus hosts have dedicated normalisation functions
        """

        logger.info("Preprocessing UniProt data.")

        for arg in tqdm(self.node_fields):

            # do not process ensembl gene ids (we will get them from pypath)
            if arg in [
                UniprotNodeField.PROTEIN_ENSEMBL_GENE_IDS.value,
                UniprotNodeField.PROTEIN_RECEPTOR_TYPE.value,
                "symbol",
            ]:
                pass

            # Integers
            elif arg in [
                UniprotNodeField.PROTEIN_LENGTH.value,
                UniprotNodeField.PROTEIN_MASS.value,
                UniprotNodeField.PROTEIN_ORGANISM_ID.value,
            ]:
                for protein, attribute_value in self.data.get(arg).items():
                    self.data[arg][protein] = int(str(attribute_value).replace(",", ""))
            # Simple replace
            elif arg not in self.split_fields:
                if not arg in [
                    UniprotNodeField.PROTEIN_SUBCELLULAR_LOCATION.value,
                    "subcellular_location",
                ]:
                    for protein, attribute_value in self.data.get(arg).items():

                        self.data[arg][protein] = (
                            attribute_value.replace("|", ",").replace("'", "^").strip()
                        )

            # Split fields
            else:

                for protein, attribute_value in self.data.get(arg).items():
                    # Field splitting
                    self.data[arg][protein] = self._split_fields(arg, attribute_value)

            # Special treatment
            # ENST and ENSG ids
            if arg == UniprotNodeField.PROTEIN_ENSEMBL_TRANSCRIPT_IDS.value:

                for protein, attribute_value in self.data.get(arg).items():

                    attribute_value, ensg_ids = self._find_ensg_from_enst(
                        attribute_value
                    )

                    # update enst in data dict
                    self.data[UniprotNodeField.PROTEIN_ENSEMBL_TRANSCRIPT_IDS.value][
                        protein
                    ] = attribute_value

                    if ensg_ids:
                        # add ensgs to data dict
                        self.data[UniprotNodeField.PROTEIN_ENSEMBL_GENE_IDS.value][
                            protein
                        ] = ensg_ids

            # Protein names
            elif arg == UniprotNodeField.PROTEIN_NAMES.value:

                for protein, attribute_value in self.data.get(arg).items():

                    self.data[arg][protein] = self._split_protein_names_field(
                        attribute_value
                    )

            elif arg == UniprotNodeField.PROTEIN_VIRUS_HOSTS.value:

                for protein, attribute_value in self.data.get(arg).items():
                    self.data[arg][protein] = self._split_virus_hosts_field(
                        attribute_value
                    )

            elif arg == UniprotNodeField.PROTEIN_SUBCELLULAR_LOCATION.value:
                for protein, attribute_value in self.data.get(arg).items():
                    individual_protein_locations = []
                    for element in attribute_value:
                        loc = (
                            str(element.location)
                            .replace("'", "")
                            .replace("[", "")
                            .replace("]", "")
                            .strip()
                        )
                        individual_protein_locations.append(loc)
                        self.locations.add(loc)

                    self.data[arg][protein] = individual_protein_locations
            elif arg == UniprotNodeField.PROTEIN_CC_DISEASE.value:
                for protein, attribute_value in self.data.get(arg).items():
                    self.data[arg][protein] = self._get_mim_ids(attribute_value)
    def get_nodes(self):
        """
        Yield nodes (protein, gene, organism) from UniProt data.
        """

        logger.info(
            "Preparing UniProt nodes of the types "
            f"{[type.name for type in self.node_types]}."
        )

        GtP = pd.read_csv("data/targets_and_families.csv", sep=",", skiprows=1)
        target_dict = dict(zip(GtP["Human SwissProt"], GtP["Type"]))

        for uniprot_entity in self._reformat_and_filter_proteins():

            protein_id, all_props = uniprot_entity

            protein_props = self._get_protein_properties(all_props)

            symbol = mapping.map_name(protein_id.split(":")[1], "uniprot", "genesymbol")

            receptor_type = target_dict.get(protein_id.split(":")[1])

            if receptor_type:
                protein_props[UniprotNodeField.PROTEIN_RECEPTOR_TYPE.value] = (
                    receptor_type
                )
            else:
                protein_props[UniprotNodeField.PROTEIN_RECEPTOR_TYPE.value] = "NA"
            if symbol:
                protein_props[UniprotNodeField.PROTEIN_SYMBOL.value] = symbol.pop()
            else:
                protein_props[UniprotNodeField.PROTEIN_SYMBOL.value] = "NA"

            yield (protein_id, "protein", protein_props)

            # append gene node to output if desired
            if UniprotNodeType.GENE in self.node_types:

                gene_list = self._get_gene(all_props)

                for gene_id, gene_props in gene_list:

                    if gene_id:
                        yield (gene_id, "gene", gene_props)

            # append organism node to output if desired
            if UniprotNodeType.ORGANISM in self.node_types:

                organism_id, organism_props = self._get_organism(all_props)

                if organism_id:
                    yield (
                        organism_id,
                        "organism",
                        organism_props,
                    )

    def _reformat_and_filter_proteins(self):
        """
        For each uniprot id, select desired fields and reformat to give a tuple
        containing id and properties. Yield a tuple for each protein.
        """

        for protein in tqdm(self.uniprot_ids):

            protein_id = self._normalise_curie_cached("uniprot", protein)

            _props = {}

            for arg in self.node_fields:
                if arg in [
                    "symbol",
                    "receptor_type",
                ]:
                    continue

                else:

                    _props[arg] = self.data.get(arg).get(protein)

            yield protein_id, _props

    def _get_gene(self, all_props: dict) -> list:
        """
        Get gene node representation from UniProt data per protein. Since one
        protein can have multiple genes, return a list of tuples.
        """

        # if genes and database(GeneID) fields exist, define gene_properties
        if not (
            UniprotNodeField.PROTEIN_GENE_NAMES.value in all_props.keys()
            and UniprotNodeField.PROTEIN_ENTREZ_GENE_IDS.value in all_props.keys()
        ):
            return []

        # Find preferred identifier for gene and check if it exists
        if UniprotEdgeField.GENE_ENTREZ_ID in self.edge_fields:

            id_type = UniprotNodeField.PROTEIN_ENTREZ_GENE_IDS.value

        elif UniprotEdgeField.GENE_ENSEMBL_GENE_ID in self.edge_fields:

            id_type = UniprotNodeField.PROTEIN_ENSEMBL_GENE_IDS.value

        gene_raw = all_props.pop(id_type)

        if not gene_raw:
            return []

        type_dict = {
            UniprotNodeField.PROTEIN_ENTREZ_GENE_IDS.value: "ncbigene",
            UniprotNodeField.PROTEIN_ENSEMBL_GENE_IDS.value: "ensembl",
        }

        gene_props = dict()

        for k in all_props.keys():

            if k not in self.gene_properties:
                continue

            if "(" in k:
                k_new = k.split("(")[1].split(")")[0].lower()
            else:
                k_new = k.lower()

            # select parenthesis content in field names and make lowercase
            gene_props[k_new] = all_props[k]

        # source, licence, and version fields
        gene_props["source"] = self.data_source
        gene_props["licence"] = self.data_licence
        gene_props["version"] = self.data_version

        gene_list = []

        genes = self._ensure_iterable(gene_raw)

        for gene in genes:

            gene_id = self._normalise_curie_cached(
                type_dict[id_type],
                gene,
            )

            gene_list.append((gene_id, gene_props))

        return gene_list

    def _get_organism(self, all_props: dict):

        organism_props = dict()

        organism_id = self._normalise_curie_cached(
            "ncbitaxon",
            all_props.pop(UniprotNodeField.PROTEIN_ORGANISM_ID.value),
        )

        for k in all_props.keys():

            if k in self.organism_properties:
                organism_props[k] = all_props[k]

        # source, licence, and version fields
        organism_props["source"] = self.data_source
        organism_props["licence"] = self.data_licence
        organism_props["version"] = self.data_version

        return organism_id, organism_props

    def _get_protein_properties(self, all_props: dict) -> dict:

        protein_props = dict()

        for k in all_props.keys():

            # define protein_properties
            if k not in self.protein_properties:
                continue

            if k == UniprotNodeField.PROTEIN_NAMES.value:
                protein_props["primary_protein_name"] = (
                    self._ensure_iterable(all_props[k])[0] if all_props[k] else None
                )

            # replace hyphens and spaces with underscore
            protein_props[k.replace(" ", "_").replace("-", "_")] = all_props[k]

        # source, licence, and version fields
        protein_props["source"] = self.data_source
        protein_props["licence"] = self.data_licence
        protein_props["version"] = self.data_version

        return protein_props
    
    def _get_mim_ids(self, text):
        """
        Split and extract MIM IDs from the text field in UniProt.
        
        Args:
            text: text field in UniProt
        Example:
            "Defects in ABCA4 are the cause of Stargardt disease (STGD1) [MIM:248200]." -> ['248200']
        """
        mim_ids = re.findall(r'\[MIM:(\d+)\]', text)
        return mim_ids

    def _split_fields(self, field_key, field_value):
        """
        Split fields with multiple entries in uniprot
        Args:
            field_key: field name
            field_value: entry of the field
        """
        if field_value:
            # replace sensitive elements for admin-import
            field_value = field_value.replace("|", ",").replace("'", "^").strip()

            # define fields that will not be splitted by semicolon
            split_dict = {
                UniprotNodeField.PROTEIN_PROTEOME.value: ",",
                UniprotNodeField.PROTEIN_GENE_NAMES.value: " ",
            }

            # if field in split_dict split accordingly
            if field_key in split_dict.keys():
                field_value = field_value.split(split_dict[field_key])
                # if field has just one element in the list make it string
                if len(field_value) == 1:
                    field_value = field_value[0]

            # split semicolons (;)
            else:
                field_value = field_value.strip().strip(";").split(";")

                # split colons (":") in kegg field
                if field_key == UniprotNodeField.PROTEIN_KEGG_IDS.value:
                    _list = []
                    for e in field_value:
                        _list.append(e.split(":")[1].strip())
                    field_value = _list

                # take first element in database(GeneID) field
                if field_key == UniprotNodeField.PROTEIN_ENTREZ_GENE_IDS.value:
                    field_value = field_value[0]

                # if field has just one element in the list make it string
                if isinstance(field_value, list) and len(field_value) == 1:
                    field_value = field_value[0]

            return field_value

        else:
            return None

    def _split_protein_names_field(self, field_value):
        """
        Split protein names field in uniprot
        Args:
            field_value: entry of the protein names field
        Example:
            "Acetate kinase (EC 2.7.2.1) (Acetokinase)" -> ["Acetate kinase", "Acetokinase"]
        """
        field_value = field_value.replace("|", ",").replace(
            "'", "^"
        )  # replace sensitive elements

        if "[Cleaved" in field_value:
            # discarding part after the "[Cleaved"
            clip_index = field_value.index("[Cleaved")
            protein_names = field_value[:clip_index].replace("(Fragment)", "").strip()

            # handling multiple protein names
            if "(EC" in protein_names[0]:
                splitted = protein_names[0].split(" (")
                protein_names = []

                for name in splitted:
                    if not name.strip().startswith("EC"):
                        if not name.strip().startswith("Fragm"):
                            protein_names.append(name.rstrip(")").strip())

            elif " (" in protein_names[0]:
                splitted = protein_names[0].split(" (")
                protein_names = []
                for name in splitted:
                    if not name.strip().startswith("Fragm"):
                        protein_names.append(name.rstrip(")").strip())

        elif "[Includes" in field_value:
            # discarding part after the "[Includes"
            clip_index = field_value.index("[Includes")
            protein_names = field_value[:clip_index].replace("(Fragment)", "").strip()
            # handling multiple protein names
            if "(EC" in protein_names[0]:

                splitted = protein_names[0].split(" (")
                protein_names = []

                for name in splitted:
                    if not name.strip().startswith("EC"):
                        if not name.strip().startswith("Fragm"):
                            protein_names.append(name.rstrip(")").strip())

            elif " (" in protein_names[0]:
                splitted = protein_names[0].split(" (")
                protein_names = []
                for name in splitted:
                    if not name.strip().startswith("Fragm"):
                        protein_names.append(name.rstrip(")").strip())

        # handling multiple protein names
        elif "(EC" in field_value.replace("(Fragment)", ""):
            splitted = field_value.split(" (")
            protein_names = []

            for name in splitted:
                if not name.strip().startswith("EC"):
                    if not name.strip().startswith("Fragm"):
                        protein_names.append(name.rstrip(")").strip())

        elif " (" in field_value.replace("(Fragment)", ""):
            splitted = field_value.split(" (")
            protein_names = []
            for name in splitted:
                if not name.strip().startswith("Fragm"):
                    protein_names.append(name.rstrip(")").strip())

        else:
            protein_names = field_value.replace("(Fragment)", "").strip()

        return protein_names

    def _split_virus_hosts_field(self, field_value):
        """
        Split virus hosts fields in uniprot

        Args:
            field_value: entry of the virus hosts field

        Example:
            "Pyrobaculum arsenaticum [TaxID: 121277]; Pyrobaculum oguniense [TaxID: 99007]" -> ['121277', '99007']
        """
        if field_value:
            if ";" in field_value:
                splitted = field_value.split(";")
                virus_hosts_tax_ids = []
                for v in splitted:
                    virus_hosts_tax_ids.append(
                        v[v.index("[") + 1 : v.index("]")].split(":")[1].strip()
                    )
            else:
                virus_hosts_tax_ids = (
                    field_value[field_value.index("[") + 1 : field_value.index("]")]
                    .split(":")[1]
                    .strip()
                )

            return virus_hosts_tax_ids
        else:
            return None

    def _find_ensg_from_enst(self, enst_list):
        """
        take ensembl transcript ids, return ensembl gene ids by using pypath mapping tool

        Args:
            field_value: ensembl transcript list

        """

        enst_list = self._ensure_iterable(enst_list)

        enst_list = [enst.split(" [")[0] for enst in enst_list]

        ensg_ids = set()
        for enst_id in enst_list:
            ensg_id = list(
                mapping.map_name(enst_id.split(".")[0], "enst_biomart", "ensg_biomart")
            )
            ensg_id = ensg_id[0] if ensg_id else None
            if ensg_id:
                ensg_ids.add(ensg_id)

        ensg_ids = list(ensg_ids)

        if len(ensg_ids) == 1:
            ensg_ids = ensg_ids[0]

        if len(enst_list) == 1:
            enst_list = enst_list[0]

        return enst_list, ensg_ids

    @lru_cache
    def _normalise_curie_cached(
        self, prefix: str, identifier: str, sep: str = ":"
    ) -> Optional[str]:
        """
        Wrapper to call and cache `normalize_curie()` from Bioregistry.
        """

        if not self.normalise_curies:
            return identifier

        return normalize_curie(f"{prefix}{sep}{identifier}", sep=sep)

    def _configure_fields(self):
        # fields that need splitting
        self.split_fields = [
            UniprotNodeField.PROTEIN_PROTEOME.value,
            UniprotNodeField.PROTEIN_GENE_NAMES.value,
            UniprotNodeField.PROTEIN_EC.value,
            UniprotNodeField.PROTEIN_ENTREZ_GENE_IDS.value,
            UniprotNodeField.PROTEIN_ENSEMBL_TRANSCRIPT_IDS.value,
            UniprotNodeField.PROTEIN_KEGG_IDS.value,
        ]

        # properties of nodes
        self.protein_properties = [
            UniprotNodeField.PROTEIN_LENGTH.value,
            UniprotNodeField.PROTEIN_MASS.value,
            UniprotNodeField.PROTEIN_NAMES.value,
            UniprotNodeField.PROTEIN_PROTEOME.value,
            UniprotNodeField.PROTEIN_EC.value,
            UniprotNodeField.PROTEIN_VIRUS_HOSTS.value,
            UniprotNodeField.PROTEIN_ORGANISM_ID.value,
        ]

        self.gene_properties = [
            UniprotNodeField.PROTEIN_GENE_NAMES.value,
            UniprotNodeField.PROTEIN_ENTREZ_GENE_IDS.value,
            UniprotNodeField.PROTEIN_KEGG_IDS.value,
            UniprotNodeField.PROTEIN_ENSEMBL_TRANSCRIPT_IDS.value,
            UniprotNodeField.PROTEIN_ENSEMBL_GENE_IDS.value,
            UniprotNodeField.PRIMARY_GENE_NAME.value,
        ]

        self.organism_properties = [UniprotNodeField.PROTEIN_ORGANISM.value]

    def _set_node_and_edge_fields(
        self, node_types, node_fields  # , edge_types, edge_fields
    ):

        # ensure computation of ENSGs
        if UniprotNodeField.PROTEIN_ENSEMBL_GENE_IDS in node_fields and not (
            UniprotNodeField.PROTEIN_ENSEMBL_TRANSCRIPT_IDS in node_fields
        ):
            node_fields.append(UniprotNodeField.PROTEIN_ENSEMBL_TRANSCRIPT_IDS)

        # check which node types and fields to include
        if node_types:

            self.node_types = node_types

        else:

            self.node_types = [field for field in UniprotNodeType]

        if node_fields:

            self.node_fields = [field.value for field in node_fields]

        else:

            self.node_fields = [field.value for field in UniprotNodeField]

    def _ensure_iterable(self, value):
        if isinstance(value, str):
            return [value]
        else:
            return value
