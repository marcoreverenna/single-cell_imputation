#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Module containing the functions for testing Beagle imputation algorithm on single cell SNP array data

__author__ = Marco Reverenna
__copyright__ = Copyright 2022-2023
__version__ = 1.0
__maintainer__ = Marco Reverenna
__email__ = marcoreverenna@gmail.com
__status__ = Dev
"""

# load the packages
import os
import logging
import pandas as pd
import random as rd
from sklearn.metrics import jaccard_score
from sklearn.metrics import recall_score

def similarity_recall_after(dictionary_imp, dictionary_gen, seeds, path_output, path_processed, path_result, sc_samples):
    j_list_post = []
    recall_list_post = []


    for sample in sc_samples:
        dict_post = []
        print(f'Analysing the pair: {sample}')
        for chrom in list(range(1,23)):
            for seed in seeds:

                SC_merg_chr_imputed = pd.read_csv(os.path.join(path_output, f'SC_cons_chr{chrom}_seed{seed}_imputed.vcf.gz'), header=8, sep='\t', dtype='object')
                gDNA_merg_chr_ctrl = pd.read_csv(os.path.join(path_processed, f'gDNA_control/gDNA_cons_chr{chrom}_seed{seed}_ctrl.vcf.gz'), header=0, sep='\t', dtype='object')   

                SC_merg_chr_imputed['FORMAT'] = 'GT'
                SC_merg_chr_imputed.iloc[:,9:] = SC_merg_chr_imputed.iloc[:,9:].apply(lambda x : x.str.split(':').str.get(0))               

                # first cleaning step: remove extra alleles (perhaps due to indels or other non common SNP)
                # SC_merg_chr_imputed = SC_merg_chr_imputed[SC_merg_chr_imputed.iloc[:, 9:].isin(['0|0','0|1','1|0','1|1','1|2','2|1','2|2','2|0', '0|2']).any(axis=1)]
                # second cleaning step: remove duplicates, rarely there are some positions that are repeated due to rare INDELS
                # third cleaning step: remove all the ungenotyped positions in bulk because we can not compare an imputed genotype with a gDNA missing value
                SC_merg_chr_imputed = SC_merg_chr_imputed[~SC_merg_chr_imputed.iloc[:, 9:].isin(['3|0', '0|3', '3|3']).any(axis=1)]
                SC_merg_chr_imputed.drop_duplicates(subset=['POS'], inplace=True)

                gDNA_merg_chr_ctrl_nomiss = gDNA_merg_chr_ctrl[~(gDNA_merg_chr_ctrl.iloc[:, 9:] == './.').any(axis=1)] 

                # create the final datasets
                SC_merg_chr_imputed_filtered = SC_merg_chr_imputed[SC_merg_chr_imputed.POS.isin(gDNA_merg_chr_ctrl_nomiss['POS'].to_list())==True]
                gDNA_merg_chr_ctrl_nomiss_filtered = gDNA_merg_chr_ctrl_nomiss[gDNA_merg_chr_ctrl_nomiss.POS.isin(SC_merg_chr_imputed_filtered['POS'].to_list())==True]    
                print(f'analising singlecell chromosome:{chrom} seed:{seed} imputed')

                # define two vectors to compare and then calculate the jaccard score considering 'micro' option
                vec_gDNA_final = gDNA_merg_chr_ctrl_nomiss_filtered['gDNA_consensus'].map(dictionary_gen).tolist()
                vec_SC_final = SC_merg_chr_imputed_filtered[sample].map(dictionary_imp).tolist()

                j_value = jaccard_score(vec_gDNA_final, vec_SC_final, average='micro')
                recall_value = recall_score(vec_gDNA_final, vec_SC_final, average='micro')

                recall_list_post.append(recall_value)
                j_list_post.append(j_value)

                dict_js_re = {'sample_name':sample,
                           'gDNA': 'consensus',
                           'chromosome':chrom,
                           'seed':seed,
                           'j_score':j_value,
                           'recall':recall_value,
                           'tot_snps':len(vec_SC_final)}

                dict_post.append(dict_js_re)
        # collect all the results inside a dataframe
        pd.DataFrame(dict_post).to_excel(os.path.join(path_result, f'tables/similarity_recall_{sample}_after.xlsx'), index=False)




def concat_tables(path_result, sc_samples):
    # tables before imputation
    file_names_pre = [os.path.join(path_result,f'tables/similarity_recall_{sample}_before.xlsx') for sample in sc_samples]
    dfs_pre = [pd.read_excel(file) for file in file_names_pre]
    result_pre = pd.concat(dfs_pre)

    # tables after imputation
    file_names_post = [os.path.join(path_result,f'tables/similarity_recall_{sample}_after.xlsx') for sample in sc_samples]
    dfs_post = [pd.read_excel(file) for file in file_names_post] 
    result_post = pd.concat(dfs_post)

    # prepare a concatenate table
    result_merged = result_post.copy()
    result_merged['j_score_pre'] = result_pre['j_score']
    result_merged['tot_snps_pre'] = result_pre['tot_snps']
    result_merged['recall_pre'] = result_pre['recall']

    result_merged.rename(columns={'j_score': 'j_score_post',
                                  'tot_snps':'tot_snps_post',
                                  'recall': 'recall_post'
                                  },
                         inplace=True)

    result_merged = result_merged[['sample_name',
                                   'chromosome',
                                   'seed',
                                   'j_score_pre',
                                   'j_score_post',
                                   'recall_pre',
                                   'recall_post',
                                   'tot_snps_pre',
                                   'tot_snps_post']
                                  ]
    result_merged.to_excel(os.path.join(path_result, f'tables/similarity_recall_concatenated.xlsx'), index=False)


if __name__ == "__main__":

    # define the paths
    data_dir = 'data/'
    processed_dir = 'data/processed/'
    output_dir = 'data/output/'
    plots_dir = 'plots/'
    results_dir = 'results/'
    map_dir = 'map/'
    reference_dir = 'reference/'
    ref_pos_dir = 'reference_positions/'

    # define global variables
    chromosomes_int = list(range(1,23))
    chromosomes_str = [str(chrom) for chrom in chromosomes_int]
    seeds = range(1, 3)
    percentages = ['10']

    sc_samples = ['201666260058_R01C01','201285550024_R06C01','201285550024_R04C01', '201285550024_R02C01',
                  '201666260058_R02C01', '201285550024_R06C02','201285550024_R04C02', '201285550024_R02C02',
                  '201666260058_R03C01', '201662330032_R06C01','201662330032_R04C01', '201662330032_R02C01',
                  '201666260058_R04C01', '201662330032_R06C02','201662330032_R04C02', '201662330032_R02C02',
                  '201666260058_R05C01', '201662330098_R06C01','201662330098_R04C01', '201662330098_R02C01']
    # create other step-folders inside processed-folder
    directories_processed = ['gDNA_chroms', 'SC_chroms', 'gDNA_control', 'SC_toimp', 'SC_control']
    for directory in directories_processed:
        dir_path = os.path.join('data/processed/', directory)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)

    try:
        logging.info("Running module: similarity_recall_after")
        similarity_recall_after(dictionary_imp= {'0|0': 0, '0|1': 1, '1|0': 1, '1|1': 2, '1|2': 3, '2|1': 3, '2|2': 4,'2|0': 5, '0|2': 5},
                        dictionary_gen= {'0/0': 0, '0/1': 1, '1/0': 1, '1/1': 2, '1/2': 3, '2/1': 3, '2/2': 4,'2/0': 5, '0/2': 5},
                        seeds= range(1, 3),
                        path_output= output_dir,
                        path_processed= processed_dir,
                        path_result =results_dir,
                        sc_samples=sc_samples

        )
    except Exception as e:
        logging.error(f"Error in similarity_recall_after: {e}")

    try:
        logging.info("Running module: concat_tables")
        concat_tables(path_result=results_dir, sc_samples = sc_samples)
    except Exception as e:
        logging.error(f"Error in concat_tables: {e}")
