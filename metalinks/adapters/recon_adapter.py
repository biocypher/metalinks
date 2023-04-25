#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
BioCypher - recon adapter prototype
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

class ReconEdgeType(Enum):
    """
    RECON edge types.
    """
    PD = "PD"

class ReconMetaboliteToProteinEdgeField(Enum):
    """
    RECON metabolite to protein edge fields.
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


class ReconAdapter:

    def __init__(
        self, 
        id_batch_size: int = int(1e6),
        edge_types: Optional[list] = None,
        edge_fields: Optional[list] = None,
        test_mode: bool = False,
    ):
        self.id_batch_size = id_batch_size
        self.data_source = "RECON"
        self.data_source_version = "3D"
        self.data_license = 'None'
        self.test_mode = test_mode

    def get_edges(self):
        """
        Get edges from RECON.
        """

        recon_path = 'data/Recon3D_301.mat'
        recon_symbols_path = 'data/recon_gene_symbols.csv'

        map1_path = 'data/mapping_tables/global_ID_mapping_curated.tsv'
        map2_path = '/home/efarr/Documents/metalinks/Data/Intermediate/Mapping/metmap_curated.csv'
        map3_path = 'data/mapping_tables/hmdb_mapping.csv'
        # map4_path = '/home/efarr/Documents/metalinks/Data/Intermediate/Mapping/recon3D_anno_curated.csv'

        recon = sio.loadmat(recon_path)
        symbols = pd.read_csv(recon_symbols_path, sep=';')
        recon = recon['Recon3D']

        data = recon

        rxn_gene_df = pd.DataFrame(data['rxnGeneMat'][0][0])
        reaction_ids = data['rxns'][0][0].flatten()
        reaction_ids = [x[0] for x in reaction_ids]
        mets = data['mets'][0][0].flatten()
        mets = [x[0] for x in mets]
        rxn_gene_df.columns = symbols['symbols']
        rxn_gene_df.index = reaction_ids
        S = pd.DataFrame(data['S'][0][0].toarray(), index=mets, columns=reaction_ids)
        lb_ub = pd.DataFrame(data['lb'][0][0], index=reaction_ids, columns=['lb'])
        lb_ub['ub'] = data['ub'][0][0]
        lb_ub['rev'] = lb_ub.apply(lambda x: 'reversible' if x['lb'] < 0 and x['ub'] > 0 else 'irreversible', axis=1)
        lb_ub['direction'] = lb_ub.apply(lambda x: 'forward' if x['ub'] > 0 else 'backward', axis=1)
        subsystem = pd.DataFrame(data['subSystems'][0][0].flatten(), index=reaction_ids, columns=['subsystem'])
        subsystem = [x[0][0][0] for x in subsystem['subsystem']]


        reaction_to_genes = get_gene_symbols(rxn_gene_df)

        reaction_to_metabolites_prod = get_metabolites(S, d = 1)
        reaction_to_metabolites_deg = get_metabolites(S, d = -1)

        metabolite_to_gene = get_metabolite_to_gene(reaction_to_metabolites_prod, reaction_to_metabolites_deg, reaction_to_genes, lb_ub)

        ss_dict = dict(zip(reaction_ids, subsystem))
        metabolite_to_gene['subsystem'] = metabolite_to_gene['reaction_id'].map(ss_dict)    
        metabolite_to_gene['compartment'] = metabolite_to_gene['metabolite_id'].apply(lambda x: x.split('[')[1])
        metabolite_to_gene['compartment'] = metabolite_to_gene['compartment'].apply(lambda x: x.split(']')[0])
        metabolite_to_gene['subsystem'].value_counts()
        metabolite_to_gene['transport'] = 'unknown'

        # if compartment is c, the subsystem is 'Transport, Extracellular'then transport is 'e->c'
        metabolite_to_gene.loc[(metabolite_to_gene['compartment'] == 'c') & (metabolite_to_gene['subsystem'] == 'Transport, extracellular'), 'transport'] = 'c->e'
        metabolite_to_gene.loc[(metabolite_to_gene['compartment'] == 'e') & (metabolite_to_gene['subsystem'] == 'Transport, extracellular'), 'transport'] = 'e->c'
        metabolite_to_gene.loc[(metabolite_to_gene['compartment'] == 'm') & (metabolite_to_gene['subsystem'] == 'Transport, mitochondrial'), 'transport'] = 'm->c'
        metabolite_to_gene.loc[(metabolite_to_gene['compartment'] == 'c') & (metabolite_to_gene['subsystem'] == 'Transport, mitochondrial'), 'transport'] = 'c->m'
        metabolite_to_gene.loc[(metabolite_to_gene['compartment'] == 'r') & (metabolite_to_gene['subsystem'] == 'Transport, endoplasmic reticular'), 'transport'] = 'r->c'
        metabolite_to_gene.loc[(metabolite_to_gene['compartment'] == 'c') & (metabolite_to_gene['subsystem'] == 'Transport, endoplasmic reticular'), 'transport'] = 'c->r'
        metabolite_to_gene.loc[(metabolite_to_gene['compartment'] == 'l') & (metabolite_to_gene['subsystem'] == 'Transport, lysosomal'), 'transport'] = 'l->c'
        metabolite_to_gene.loc[(metabolite_to_gene['compartment'] == 'c') & (metabolite_to_gene['subsystem'] == 'Transport, lysosomal'), 'transport'] = 'c->l'
        metabolite_to_gene.loc[(metabolite_to_gene['compartment'] == 'x') & (metabolite_to_gene['subsystem'] == 'Transport, peroxisomal'), 'transport'] = 'x->c'
        metabolite_to_gene.loc[(metabolite_to_gene['compartment'] == 'c') & (metabolite_to_gene['subsystem'] == 'Transport, peroxisomal'), 'transport'] = 'c->x'
        metabolite_to_gene.loc[(metabolite_to_gene['compartment'] == 'g') & (metabolite_to_gene['subsystem'] == 'Transport, golgi apparatus'), 'transport'] = 'g->c'
        metabolite_to_gene.loc[(metabolite_to_gene['compartment'] == 'c') & (metabolite_to_gene['subsystem'] == 'Transport, golgi apparatus'), 'transport'] = 'c->g'
        metabolite_to_gene.loc[(metabolite_to_gene['compartment'] == 'n') & (metabolite_to_gene['subsystem'] == 'Transport, nuclear'), 'transport'] = 'n->c'
        metabolite_to_gene.loc[(metabolite_to_gene['compartment'] == 'c') & (metabolite_to_gene['subsystem'] == 'Transport, nuclear'), 'transport'] = 'c->n'

        metabolite_to_gene['transport_direction'] = 'unknown'
        metabolite_to_gene[metabolite_to_gene['subsystem'].str.contains('Transport')]['transport_direction'] = 'out'
        metabolite_to_gene.loc[metabolite_to_gene['transport'].str.startswith('c'), 'transport_direction'] = 'in'
        metabolite_to_gene.loc[metabolite_to_gene['transport'] == 'c->e', 'transport_direction'] = 'out'
        metabolite_to_gene.loc[metabolite_to_gene['transport'] == 'e->c', 'transport_direction'] = 'in'

        print(f'collapsed metabolites to genes, now have {len(metabolite_to_gene)} metabolite to gene links')

        metmap1 = pd.read_csv(map1_path, sep='\t', dtype=object)
        metmap2 = pd.read_csv(map2_path, sep='\t', dtype=str)
        metmap3 = pd.read_csv(map3_path, sep=',', dtype=str)

        field_names = ['metPubChemID',  'metHMDBID', 'metKEGGID', 'metCHEBIID',  'mets']

        df = pd.DataFrame(columns=field_names)

        for field_name in field_names:
            df[field_name] = data[field_name][0][0].flatten()
            df[field_name] = df[field_name].apply(lambda x: np.nan if x.size == 0 else x[0])   

        df.rename(columns={'metPubChemID': 'pubchem_id', 
                           'metHMDBID': 'hmdb_id',
                           'metKEGGID' : 'kegg_id',
                           'metCHEBIID':'chebi_id'}, inplace=True)

        # df =      pd.read_csv(map4_path, sep=',', dtype=str)

        print(f'loaded metabolite mapping files')

        dfs = preprocess_metmaps(df, metmap1, metmap2, metmap3)

        test = fill_missing_values(df1=dfs[0], df2=dfs[3], df3=dfs[1])
        test = fill_missing_values(test[0], test[1], test[2], 'kegg_id', 'chebi_id', 'hmdb_id', 'pubchem_id')
        test = fill_missing_values(test[0], test[1], test[2], 'hmdb_id', 'chebi_id', 'kegg_id', 'pubchem_id')
        test = fill_missing_values(test[0], test[1], test[2], 'pubchem_id', 'chebi_id', 'kegg_id', 'hmdb_id')

        print(f'filled missing values in metabolite mapping files')

        met_dict = dict(zip(mets, test[0]['hmdb_id']))
        metabolite_to_gene['hmdb_id'] = metabolite_to_gene['metabolite_id'].apply(lambda x: met_dict[x])

        # drop metabolite and reaction ids
        metabolite_to_gene.drop(['metabolite_id', 'reaction_id'], axis=1, inplace=True)
        metabolite_to_gene.drop_duplicates(inplace=True)
        metabolite_to_gene.dropna(subset=['hmdb_id'], inplace=True)
        metabolite_to_gene['status'] = 'recon'
        metabolite_to_gene['uniprot'] = symbol_to_uniprot(metabolite_to_gene['gene_id'])
        # fuse a 'uniprot:' uniprot columns
        metabolite_to_gene['uniprot'] = metabolite_to_gene['uniprot'].apply(lambda x: 'uniprot:' + x if x is not np.nan else x)
        

        for row in tqdm(metabolite_to_gene.iterrows()):
            attributes  = {
                'status': row[1]['status'],
                'direction': row[1]['direction'],
                'symbol': row[1]['gene_id'],
                'subsystem': row[1]['subsystem'],
                'transport': row[1]['transport'],
                'transport_direction': row[1]['transport_direction']
            }
            r = row[1].astype(str)
            h = hashlib.md5(''.join(r).encode('utf-8')).hexdigest()
            yield h, row[1]['hmdb_id'], row[1]['uniprot'], 'PD', attributes



def get_metabolites(S):
    S = S.copy()
    S[S != 0] = 1
    S = S.stack().reset_index()
    S.columns = ['metabolite_id', 'reaction_id', 'value']
    S = S[S['value'] == 1]
    S = S.drop('value', axis=1)
    S.drop_duplicates(inplace=True)
    return S

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


def get_metabolites(S, d = 1):
    S = S.copy()
    S[S != d] = 0
    S[S == d] = 1
    S = S.stack().reset_index()
    S.columns = ['metabolite_id', 'reaction_id', 'value']
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
    rev_df = metabolite_to_gene[metabolite_to_gene['reaction_id'].isin(reversible_reactions)]
    rev_df['direction'] = rev_df['direction'].apply(lambda x: 'degrading' if x == 'producing' else 'producing')
    metabolite_to_gene = pd.concat([metabolite_to_gene, rev_df])
    return metabolite_to_gene


def fill_missing_values(df1, df2, df3, str1 = 'chebi_id', str2 = 'kegg_id', str3 = 'hmdb_id', str4 = 'pubchem_id'):
    df1 = df1.copy()
    # store unique values in each column of df1 in one list
    before = []
    for i in range(0, len(df1.columns)):
        before.append(df1.iloc[:,i].unique())

    dfx = df1.dropna(subset=[str1])
    # order the columns in dfx
    dfx = dfx[[str1, str2, str3, str4]]
    test = pd.merge(dfx, df2, how='left', on= str1).drop_duplicates()
    test = pd.merge(test, df3, how='left', on= str1).drop_duplicates()

    test.iloc[:,[1,4,7]] = check_fill(test.iloc[:,[1,4,7]])
    test.iloc[:,[2,5,8]] = check_fill(test.iloc[:,[2,5,8]])
    test.iloc[:,[3,6,9]] = check_fill(test.iloc[:,[3,6,9]])
    
    dict1 = create_dict(test, [str1, str2 + '_x'])
    dict2 = create_dict(test, [str1, str3 + '_x'])
    dict3 = create_dict(test, [str1, str4 + '_x'])

    df1 = fillna_with_map(df1, str1, str2, str3, str4, dict1, dict2, dict3) # still used the dictionaries but in the loop
    df2 = fillna_with_map(df2, str1, str2, str3, str4, dict1, dict2, dict3)
    df3 = fillna_with_map(df3, str1, str2, str3, str4, dict1, dict2, dict3)

    after = []
    for i in range(0, len(df1.columns)):
        after.append(df1.iloc[:,i].unique())

    for i in range(0, len(before)):
        print(f'before: {len(before[i])}, after: {len(after[i])}, additional: {len(set(after[i]) - set(before[i]))}')

    return df1, df2, df3


def get_hmdb_ids_s(df, metmap3):
    df = df.merge(metmap3, on='pubchem_id', how='left')
    df.drop_duplicates(inplace=True)
    return df


def drop_nan(df, col1, col2, col3):
    df = df.copy()
    df = df.dropna(subset=[col1, col2, col3], how='all')
    df[col1] = df[col1].fillna(df[col2])
    df[col1] = df[col1].fillna(df[col3])
    df = df.drop(col2, axis=1)
    df = df.drop(col3, axis=1)
    # if length is 9 insert two zeros after the first 4 characters
    df[col1] = np.where(df[col1].str.len() == 9, df[col1].str[:4] + '00' + df[col1].str[4:], df[col1])
    return df


def preprocess_metmaps(df, metmap1, metmap2, metmap3):
    metmap2.rename(columns={'CID': 'pubchem_id', 'KEGG' : 'kegg_id', 'HMDB' : 'hmdb_id', 'ChEBI' : 'chebi_id'}, inplace=True)
    metmap2['chebi_id'] = 'CHEBI:' + metmap2['chebi_id']
    df.rename(columns={'metCHEBIID': 'chebi_id', 'metKEGGID': 'kegg_id', 'metHMDBID': 'hmdb_id', 'metPubChemID': 'pubchem_id'}, inplace=True)
    metmap1['pubchem_id'] = metmap1['pubchem_id'].apply(lambda x: str(int(x)) if not pd.isnull(x) else x)
    df['chebi_id'] = df['chebi_id'].apply(lambda x: 'CHEBI:' + x if not pd.isnull(x) and not x.startswith('CHEBI:') else x)

    # remove spaces from the beginning and end of the strings in df['kegg_id'] by extracting a substring starting with a 'C' and having 6 characters
    df['kegg_id'] = df['kegg_id'].apply(lambda x: x[x.find('C'):x.find('C')+6] if not pd.isnull(x) else x)
    # do the same for metmap2['kegg_id']
    metmap2['kegg_id'] = metmap2['kegg_id'].apply(lambda x: x[x.find('C'):x.find('C')+6] if not pd.isnull(x) else x)
    # change dtype of hmdb_id to str
    metmap2['hmdb_id'] = metmap2['hmdb_id'].astype(str)
    # if length of metmap2[str3 + '_id'] is lower that 11 add zeros after the first 4 characters to top up to 11 characters only if the string is not nan
    metmap2['hmdb_id'] = metmap2['hmdb_id'].apply(lambda x: x if len(x) == 11 else x[:4] + '0'*(11-len(x)) + x[4:])
    # set value in hmd_id to nan if it starts with 'nan
    metmap2['hmdb_id'] = metmap2['hmdb_id'].apply(lambda x: np.nan if x.startswith('nan') else x)
    # count length of strings in hmdb_id
    metmap2['hmdb_id'].str.len().value_counts()

    # do the same for df['hmdb_id']
    df['hmdb_id'] = df['hmdb_id'].astype(str)
    df['hmdb_id'] = df['hmdb_id'].apply(lambda x: x if len(x) == 11 else x[:4] + '0'*(11-len(x)) + x[4:])
    df['hmdb_id'] = df['hmdb_id'].apply(lambda x: np.nan if x.startswith('nan') else x)
    df['hmdb_id'].str.len().value_counts()

    # rename metmap3['accession'] to 'hmdb_id'
    metmap3.rename(columns={'accession': 'hmdb_id'}, inplace=True)
    # add a 'CHEBI:' to the beginning of the strings in metmap3['chebi_id'] if they are not nan
    metmap3['chebi_id'] = metmap3['chebi_id'].apply(lambda x: 'CHEBI:' + x if not pd.isnull(x) else x)

    df1 = df[['chebi_id', 'kegg_id', 'hmdb_id', 'pubchem_id']]
    df2 = metmap1[['chebi_id', 'kegg_id', 'hmdb_id', 'pubchem_id']]
    df3 = metmap2[['chebi_id', 'kegg_id', 'hmdb_id', 'pubchem_id']]
    df4 = metmap3[['chebi_id', 'kegg_id', 'hmdb_id', 'pubchem_id']]

    return df1, df2, df3, df4


def fillna_with_map(df, str1, str2, str3, str4, dict1, dict2, dict3):
    df = df.copy()
    dictlist = [dict1, dict2, dict3]
    counter = 0
    for col in [str2, str3, str4]:
        df[col] = df[col].fillna(df[str1].map(dictlist[counter]))
        counter += 1
    return df


def check_fill(df):
    counter = []
    df = df.copy()
    # extract name of columns bi extracting the characters of the first column names until the first '_'
    str1 = df.columns[0][:df.columns[0].find('_')]
    for i in range(0, len(df)):
        row = df.iloc[i]
        # if row has one na value, check if other two rows have the same value, if yes, fill in the na value
        if row.isna().sum() <= 1:
            # if all three values are the same continue to next row
            if row[0] == row[1] and row[0] == row[2]:
                continue
            elif row[0] == row[1]:
                df.iloc[i,2] = row[0]
                continue
            elif row[0] == row[2]:
                df.iloc[i,1] = row[0]
                continue
            elif row[1] == row[2]:
                df.iloc[i,0] = row[1]
                continue
            else:
            # entry that is nan will be set to on of the other two values
                df.iloc[i, df.iloc[i].isna()] = row[~row.isna()][0]
                counter.append(row)
        # if row has two na values set both to the third value
        elif row.isna().sum() == 2:
            df.iloc[i, df.iloc[i].isna()] = row[~row.isna()]
        # if row has no na values, do nothing
    print(f' found {len(counter)} entries that had a conflict in type: {str1}')
    return df



# write function that takes in a dataframe and a list of columns and fills creates a dictionary with the values of the columns
# keys are the first columns and values are the values of the other columns, if the value is nan, it is not added to the dictionary
def create_dict(df, list_of_columns):
    df1 = df.dropna(subset=[list_of_columns[1]])
    keys = df1[list_of_columns[0]]
    values = df1[list_of_columns[1]]
    out_dict = dict(zip(keys, values))
    return out_dict



def symbol_to_uniprot(ensp_list):
    gene_symbol_list = []
    for element in ensp_list:
        symbol = mapping.map_name(str(element), 'genesymbol', 'uniprot')
        if symbol != set():
            gene_symbol_list.append(symbol.pop())
        else:
            gene_symbol_list.append('NA')
    return gene_symbol_list



# solution to exclude local files ( very slow)

# from pypath.inputs import hmdb
# map_table = hmdb.hmdb_table('pubchem_compound_id', 'accession', 'kegg_id', 'chebi_id')
# names = ['pubchem', 'hmdb', 'kegg', 'chebi']
# map_table['chebi_id'] = map_table['chebi_id'].apply(lambda x: 'CHEBI:' + x if not pd.isnull(x) else x)

# map_dict = {}
# for i in range(4):
#     for j in range(4):
#         if i != j:
#             map_dict[(i, j)] = dict(zip(map_table.iloc[:, i], map_table.iloc[:, j]))
#             map_dict[names[i] + '_' + names[j]] = map_dict.pop((i, j))

# for key in map_dict.keys():
#     map_dict[key] = {k: v for k, v in map_dict[key].items() if k is not None and v is not None}

# for i in names:
#     names_remain = names.copy()
#     names_remain.remove(i)
#     for j in df[i]:
#         if j is not np.nan:
#             for k in names_remain:
#                 if j in map_dict[i + '_' + k].keys():
#                     df[k][df[i] == j] = map_dict[i + '_' + k][j]

# # update map_dict to include the mapping from df
# for i in names:
#     for j in names:
#         if i != j:
#             update = dict(zip(df[i], df[j]))
#             # remove all the NA values in key and value
#             update = {k: v for k, v in update.items() if k is not np.nan and v is not np.nan}
#             map_dict[i + '_' + j].update(update)


