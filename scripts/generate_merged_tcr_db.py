# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 10:13:08 2021
@author: sschattg

This script is used to parse the mcpas and vdj-db databases and append them to the current paired TCR database package with conga.
These database files can be used for single chain matching.
Two tsv files are generated as output, one for mouse and for with human.
"""

import numpy as np
import pandas as pd
import math
import re
import argparse
import sys
from os.path import exists

parser = argparse.ArgumentParser(description="Regenerate human and mouse databases for single chain matching starting with the paired db file")

#type is str by default
parser.add_argument('--paired_db', help='Paired TCR db file new_paired_tcr_db_for_matching_nr.tsv packaged with conga')
parser.add_argument('--mcpas', help='McPAS db file. Usually named McPAS-TCR.csv by default')
parser.add_argument('--vdjdb', help='vdjdb db file. Either vdjdb.txt or vdjdb_full.txt')
parser.add_argument('--output_dir', help='Where the output files should go')
parser.add_argument('--remove_10x_200k_from_vdjdb', action='store_true', help='Filter out 10X_200K dataset from vdjdb')
args = parser.parse_args()

assert exists(args.paired_db)
assert exists(args.mcpas)
assert exists(args.vdjdb)

if args.output_dir is None:
    output_dir = './'
    print('No output directory indicated. Using current directory.')
else:
    output_dir = args.output_dir

# start with the original phil one and build from there
conga_paired = pd.read_csv(args.paired_db, delimiter='\t')
conga_paired['match_cdr'] = conga_paired['cdr3a'] + ';' + conga_paired['cdr3b']
conga_paired['species'] = 'human'

# dictionaries for remapping
mcpas_column_map = {'CDR3.alpha.aa' : 'cdr3a',
 'CDR3.beta.aa' : 'cdr3b',
 'TRAV' : 'va',
 'TRAJ' : 'ja',
 'TRBV' : 'vb',
 'TRBJ' : 'jb',
 'Pathology': 'epitope_species',
 'Antigen.protein': 'epitope_gene',
 'Epitope.peptide': 'epitope',
 'Species' :'species',
 'MHC' : 'mhc',
 'PubMed.ID': 'dataset'}

vdjdb_column_map = {'cdr3.alpha' : 'cdr3a',
 'cdr3.beta' : 'cdr3b',
 'v.alpha' : 'va',
 'j.alpha' : 'ja',
 'v.beta' : 'vb',
 'j.beta' : 'jb',
 'antigen.species': 'epitope_species',
 'antigen.gene': 'epitope_gene',
 'antigen.epitope': 'epitope',
 'reference.id' : 'dataset'}

add_cols = ('donor', 'pmhc', 'is_second_alpha_chain')

# reformat the mcpas to match paired db
mcpas = pd.read_csv(args.mcpas, encoding = 'mbcs', low_memory=False)
mcpas = mcpas.rename(columns= mcpas_column_map)
mcpas = mcpas[ mcpas['species'].isin( ['Mouse','Human'])]
mcpas['db'] = 'mcpas'
mcpas['match_cdr'] = mcpas['cdr3a'] + ';' + mcpas['cdr3b']

for col in add_cols:
    mcpas[col] = np.nan

hla_trim = []
for l in mcpas.itertuples():
    if isinstance(l.mhc, str):
        if re.match('-', l.mhc):
            trim = l.mhc.split('-')[1].split(':')[0]
        else:
            trim = l.mhc.split(':')[0]    
    else:
        trim = l.mhc
    hla_trim.append(trim)
    
mcpas['mhc_trim'] = hla_trim

mcpas_filter = mcpas[ ~mcpas['match_cdr'].isin( conga_paired['match_cdr']) ].copy() #remove matches in conga_db

mcpas_clean = mcpas_filter[['cdr3a',
 'cdr3b',
 'dataset',
 'epitope',
 'epitope_gene',
 'epitope_species',
 'ja',
 'jb',
 'mhc',
 'va',
 'vb',
 'pmhc',
 'donor',
 'is_second_alpha_chain',
 'db',
 'mhc_trim',
 'match_cdr',
 'species']] #can't pass a list of columns in pandas...


# reformat the vdjdb to match ====
vdjdb = pd.read_csv(args.vdjdb, delimiter='\t', low_memory=False)
vdjdb = vdjdb.rename(columns= vdjdb_column_map)
vdjdb = vdjdb[ vdjdb['species'].isin( ['MusMusculus','HomoSapiens'])]
vdjdb['db'] = 'vdjdb'
vdjdb['mhc'] = vdjdb['mhc.a'] + ',' + vdjdb['mhc.b']
vdjdb['mhca'] = vdjdb['mhc.a']
vdjdb['mhcb'] = vdjdb['mhc.b']
vdjdb['match_cdr'] = vdjdb['cdr3a'] + ';' + vdjdb['cdr3b']

for col in add_cols:
    vdjdb[col] = np.nan

if args.remove_10x_200k_from_vdjdb:
    #remove 10X data
    vdjdb = vdjdb[ ~vdjdb['dataset'].isin( ['https://www.10xgenomics.com/resources/application-notes/a-new-way-of-exploring-immunity-linking-highly-multiplexed-antigen-recognition-to-immune-repertoire-and-phenotype/#'])] 

hla_trim_vdj = []
for l in vdjdb.itertuples():
    if l.mhcb == 'B2M':
        trim1 = l.mhca.split('-')[1].split(':')[0]
        hla_trim_vdj.append(trim1)
    else:
        trim2 = l.mhcb.split('-')[1].split(':')[0]
        hla_trim_vdj.append(trim1 + ',' + trim2)
  
vdjdb['mhc_trim'] = hla_trim_vdj

vdjdb_filter = vdjdb[ ~vdjdb['match_cdr'].isin( conga_paired['match_cdr']) ].copy() #remove matches in conga db

vdjdb_clean = vdjdb_filter[['cdr3a',
 'cdr3b',
 'dataset',
 'epitope',
 'epitope_gene',
 'epitope_species',
 'ja',
 'jb',
 'mhc',
 'va',
 'vb',
 'pmhc',
 'donor',
 'is_second_alpha_chain',
 'db',
 'mhc_trim',
 'match_cdr',
 'species']]

#bind the 3 df together
big_db_df = pd.concat([conga_paired, mcpas_clean, vdjdb_clean], ignore_index=True, sort=False)
big_db_df = big_db_df.drop(columns= 'match_cdr')
big_db_df = big_db_df.drop(columns= 'mhc_trim') #drop trim for now

#split into human and mouse
human_df = big_db_df[big_db_df['species'].isin(['Human','HomoSapiens'])].copy()
mouse_df = big_db_df[big_db_df['species'].isin(['Mouse','MusMusculus'])].copy()

#annotate the mhc class type
human_mhc_class = []
for l in human_df.itertuples():
    if isinstance(l.mhc, str):
        if re.match('D', l.mhc):
            mhc_class = 'class_2'
        else:
            mhc_class = 'class_1' 
    else:
        mhc_class = l.mhc
    human_mhc_class.append(mhc_class)

human_df['mhc_class'] = human_mhc_class

mouse_mhc_class = []
for l in mouse_df.itertuples():
    if isinstance(l.mhc, str):
        if re.match('DR', l.mhc) or re.match('DQ', l.mhc) or re.match('I-Ab', l.mhc) or re.match('Eb', l.mhc):
            mhc_class = 'class_2'
        else:
            mhc_class = 'class_1' 
    else:
        mhc_class = l.mhc
    mouse_mhc_class.append(mhc_class)
mouse_df['mhc_class'] = mouse_mhc_class

#write out tsv files
human_df.to_csv(output_dir + 'human_tcr_db_for_matching.tsv', sep = '\t', index= False)
mouse_df.to_csv(output_dir + 'mouse_tcr_db_for_matching.tsv', sep = '\t', index= False)

print('Successfully regenerated database files:')
print(output_dir + 'human_tcr_db_for_matching.tsv')
print(output_dir + 'mouse_tcr_db_for_matching.tsv')
