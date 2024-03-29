#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
BioCypher - hmr adapter prototype
"""

from enum import Enum
from typing import Optional
import pandas as pd
import numpy as np
import hashlib
from pypath.utils import mapping
from tqdm import tqdm
import numpy as np
import scipy.io as sio
from pypath.utils import mapping

from biocypher._logger import logger

logger.debug(f"Loading module {__name__}.")

class HmrEdgeType(Enum):
    """
    HMR edge types.
    """
    PD_hmr = "PD_hmr"

class HmrMetaboliteToProteinEdgeField(Enum):
    """
    HMR metabolite to protein edge fields.
    """
    
    METABOLITE_ID = "METABOLITE_ID"
    PROTEIN_ID = "PROTEIN_ID"
    REACTION_ID = "REACTION_ID"

    _PRIMARY_SOURCE_ID = 'METABOLITE_ID'
    _PRIMARY_TARGET_ID = 'PROTEIN_ID'
    _PRIMARY_REACTION_ID = 'REACTION_ID'

    DIRECTION = "direction"
    STATUS = "status"
    SUBSYSTEM = "subsystem"
    TRANSPORT = "transport"
    TRANSPORT_DIRECTION = 'transport_direction'
    REV = 'rev'


class HmrAdapter:

    def __init__(
        self, 
        id_batch_size: int = int(1e6),
        edge_types: Optional[list] = None,
        edge_fields: Optional[list] = None,
        test_mode: bool = False,
    ):
        self.id_batch_size = id_batch_size
        self.data_source = "HMR"
        self.data_source_version = "1.5.0"
        self.data_license = 'None'
        self.test_mode = test_mode

    def get_edges(self):
        """
        Get edges from HMR.
        """

        hmr = sio.loadmat('data/HMR/Human-GEM.mat')
        hmr = hmr['ihuman']
        genes = pd.read_csv('data/HMR/genes.tsv', sep='\t', index_col=0)
        reactions = pd.read_csv('data/HMR/reactions.tsv', sep='\t', index_col=0)
        metabolites = pd.read_csv('data/HMR/metabolites.tsv', sep='\t', index_col=0)

        data = hmr

        reaction_ids = np.array(reactions.index)
        mets = np.array(metabolites.index)
        rxn_gene_df = pd.DataFrame(data['rxnGeneMat'][0][0].toarray(), index=reaction_ids, columns=genes['geneSymbols'])
        S = pd.DataFrame(data['S'][0][0].toarray(), index=mets, columns=reaction_ids)
        lb_ub = pd.DataFrame(data['lb'][0][0].flatten(), index=reaction_ids, columns=['lb'])
        lb_ub['ub'] = data['ub'][0][0].flatten()
        lb_ub['rev'] = lb_ub.apply(lambda x: 'reversible' if x['lb'] < 0 and x['ub'] > 0 else 'irreversible', axis=1)
        lb_ub['direction'] = lb_ub.apply(lambda x: 'forward' if x['ub'] > 0 else 'backward', axis=1)
        subsystem = pd.DataFrame(data['subSystems'][0][0].flatten(), index=reaction_ids, columns=['subsystem'])
        subsystem = [x[0][0][0] for x in subsystem['subsystem']]

        reaction_to_genes = get_gene_symbols(rxn_gene_df)

        reaction_to_metabolites_prod = get_metabolites(S, d = 1)
        reaction_to_metabolites_deg = get_metabolites(S, d = -1)

        reaction_to_metabolites_prod['transport'] = get_comp_dir(reaction_to_metabolites_prod, S.columns)
        reaction_to_metabolites_deg['transport'] = get_comp_dir(reaction_to_metabolites_deg, S.columns)

        metabolite_to_gene = get_metabolite_to_gene(reaction_to_metabolites_prod, reaction_to_metabolites_deg, reaction_to_genes, lb_ub)

        ss_dict = dict(zip(reaction_ids, subsystem))
        metabolite_to_gene['subsystem'] = metabolite_to_gene['reaction_id'].map(ss_dict)    
        metabolite_to_gene['compartment'] = metabolite_to_gene['metabolite_id'].str[-1]

        metabolite_to_gene['transport_direction'] = 'unknown'
        metabolite_to_gene.loc[metabolite_to_gene['subsystem'].str.contains('Transport'), 'transport_direction'] = 'out'
        metabolite_to_gene.loc[metabolite_to_gene['transport'].str.startswith('c'), 'transport_direction'] = 'in'
        metabolite_to_gene.loc[metabolite_to_gene['transport'] == 'c->e', 'transport_direction'] = 'out'
        metabolite_to_gene.loc[metabolite_to_gene['transport'] == 'e->c', 'transport_direction'] = 'in'

        print(f'collapsed metabolites to genes, now have {len(metabolite_to_gene)} metabolite to gene links')

        df = metabolites.iloc[:,2:6]

        df.rename(columns={'metPubChemID': 'pubchem_id', 
                           'metHMDBID': 'hmdb_id',
                           'metKEGGID' : 'kegg_id',
                           'metChEBIID':'chebi_id'}, inplace=True)

        print(f'loaded metabolite mapping files')

        print(f'filled missing values in metabolite mapping files')

        met_dict = dict(zip(mets, df['hmdb_id']))
        metabolite_to_gene['hmdb_id'] = metabolite_to_gene['metabolite_id'].apply(lambda x: met_dict[x])

        metabolite_to_gene.drop(['metabolite_id', 'reaction_id'], axis=1, inplace=True)
        metabolite_to_gene.drop_duplicates(inplace=True)
        metabolite_to_gene.dropna(subset=['hmdb_id'], inplace=True)
        metabolite_to_gene['status'] = 'hmr'
        uniprot_df = mapping.translation_df('uniprot', 'genesymbol')
        if 'RORA' not in uniprot_df['genesymbol'].values: #solves weird error that sometimes pypath gives the wrong column names
            uniprot_df = uniprot_df.rename(columns={'genesymbol': 'uniprot', 'uniprot': 'genesymbol'})

        uniprot_dict = dict(zip(uniprot_df['genesymbol'], uniprot_df['uniprot']))
        metabolite_to_gene['uniprot'] = metabolite_to_gene['gene_id'].map(uniprot_dict)
        print(f'{metabolite_to_gene["uniprot"].isna().sum()} uniprot ids are missing')
        metabolite_to_gene.dropna(subset=['uniprot'], inplace=True)
        metabolite_to_gene['uniprot'] = metabolite_to_gene['uniprot'].apply(lambda x: 'uniprot:' + x if x is not np.nan else x)


        for row in tqdm(metabolite_to_gene.iterrows()):
            attributes  = {
                'status': row[1]['status'],
                'direction': row[1]['direction'],
                'symbol': row[1]['gene_id'],
                'subsystem': row[1]['subsystem'],
                'transport': row[1]['transport'],
                'transport_direction': row[1]['transport_direction'], 
                'rev': row[1]['rev'],
            }
            r = row[1].astype(str)
            h = hashlib.md5(''.join(r).encode('utf-8')).hexdigest()
            yield h, row[1]['hmdb_id'], row[1]['uniprot'], 'PD_hmr', attributes


def get_gene_symbols(rxn_gene_df):
    row_sums = rxn_gene_df.sum(axis=1)
    rxn_genes = rxn_gene_df.index[row_sums > 0]
    rxn_gene_df = rxn_gene_df.loc[rxn_genes]
    rxn_gene_df = rxn_gene_df.stack().reset_index()
    rxn_gene_df.columns = ['reaction_id', 'gene_id', 'value']
    rxn_gene_df = rxn_gene_df[rxn_gene_df['value'] == 1]
    rxn_gene_df = rxn_gene_df.drop('value', axis=1)
    rxn_gene_df.drop_duplicates(inplace=True)
    return rxn_gene_df

def get_comp_dir(reaction_to_metabolites, reactions):
    reaction_to_metabolites['comp_dir'] = 'unknown'
    for i in range(len(reactions)):
        a = reaction_to_metabolites[ reaction_to_metabolites['reaction_id'] == reactions[i]]['comp_out'].unique()   
        if len(a) != 2:
            continue       
        else:
            c1 = a[0] + '->' + a[1]
            c2 = a[1] + '->' + a[0]
            reaction_to_metabolites.loc[(reaction_to_metabolites['reaction_id'] == reactions[i]) & (reaction_to_metabolites['comp_out'] == a[0]), 'comp_dir'] = c1
            reaction_to_metabolites.loc[(reaction_to_metabolites['reaction_id'] == reactions[i]) & (reaction_to_metabolites['comp_out'] == a[1]), 'comp_dir'] = c2
    return reaction_to_metabolites['comp_dir']

def get_metabolites(S, d = 1):
    S = S.copy()
    S[S != d] = 0
    S[S == d] = 1
    S = S.stack().reset_index()
    S.columns = ['metabolite_id', 'reaction_id', 'value']
    S['comp_out'] = S['metabolite_id'].apply(lambda x: x[-1])
    S = S[S['value'] == 1]
    S = S.drop('value', axis=1)
    S.drop_duplicates(inplace=True)
    return S

def get_metabolite_to_gene(reaction_to_metabolites_prod, reaction_to_metabolites_deg, reaction_to_genes, lb_ub):
    metabolite_to_gene = pd.merge(reaction_to_metabolites_prod, reaction_to_genes, on='reaction_id')
    metabolite_to_gene_deg = pd.merge(reaction_to_metabolites_deg, reaction_to_genes, on='reaction_id')
    metabolite_to_gene_deg['direction'] = 'degrading'
    metabolite_to_gene = pd.concat([metabolite_to_gene, metabolite_to_gene_deg])
    metabolite_to_gene['direction'] = metabolite_to_gene['direction'].apply(lambda x: 'producing' if x != 'degrading' else x)
    reversible_reactions = lb_ub[lb_ub['rev'] == 'reversible'].index
    metabolite_to_gene['rev'] = metabolite_to_gene['reaction_id'].apply(lambda x: 'reversible' if x in reversible_reactions else 'irreversible')
    return metabolite_to_gene

def get_hmdb_ids_s(df, metmap3):
    df = df.merge(metmap3, on='pubchem_id', how='left')
    df.drop_duplicates(inplace=True)
    return df

# def symbol_to_uniprot(ensp_list):
#     gene_symbol_list = []
#     for element in ensp_list:
#         symbol = mapping.map_name(str(element), 'genesymbol', 'uniprot')
#         if symbol != set():
#             gene_symbol_list.append(symbol.pop())
#         else:
#             gene_symbol_list.append('NA')
#     return gene_symbol_list
