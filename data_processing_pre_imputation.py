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



pd.options.mode.chained_assignment = None

def concatenate_SC(files, path_data, output_filename):
    vcf_list = []

    for i, file in enumerate(files):
        file_path = os.path.join(path_data, f'raw/GM12878_{file}.vcf.gz')
        vcf_data = pd.read_csv(file_path, header=30, sep='\t', dtype='object')

        if i == 0:
            # if it is the first file then include all the columns
            vcf_list.append(vcf_data)
        else:
            # or include only the last column
            vcf_list.append(vcf_data.iloc[:, -1])

    file_merged = pd.concat(vcf_list, axis=1)

    output_path = os.path.join(path_data, 'processed', output_filename)
    file_merged.to_csv(output_path, sep='\t', index=False)



def filtering_all_SC(path_processed):
    chromosomes_str = [str(chrom) for chrom in list(range(1,23))]
    df_consensus_filtered = pd.read_csv(os.path.join(path_processed,f'gDNA_consensus_dataframe.vcf.gz'), header=0, sep='\t', dtype='object')
    df_sc = pd.read_csv(os.path.join(path_processed, f'GM12878_SC_merged.vcf.gz'), header=0, sep='\t', dtype='object')
    df_sc_autosomes = df_sc[df_sc['#CHROM'].isin(chromosomes_str)]
    print(df_consensus_filtered.shape[0], df_sc_autosomes.shape[0])
    id_list = df_consensus_filtered['ID'].tolist()
    df_sc_filtered = df_sc_autosomes[df_sc_autosomes['ID'].isin(id_list)]
    print(df_consensus_filtered.shape[0], df_sc_filtered.shape[0])
    df_sc_filtered.to_csv(os.path.join(path_processed, f'SC_consensus_dataframe.vcf.gz'), sep='\t', index=False)



def splitting_chromosomes(kind, chrom, path_processed):
    vcf = pd.read_csv(os.path.join(path_processed, f'{kind}_consensus_dataframe.vcf.gz'), header=0, sep='\t', dtype='object')
    vcf_chr = vcf[vcf['#CHROM'] == chrom]
    vcf_chr['FORMAT'] = 'GT'
    vcf_chr.loc[:, vcf_chr.columns[9:]] = vcf_chr.loc[:, vcf_chr.columns[9:]].apply(lambda x: x.str.split(':').str.get(0))
    vcf_chr.to_csv(os.path.join(path_processed, f'{kind}_chroms/{kind}_consensus_chr{chrom}.vcf.gz'), sep='\t', index=False)



def dataset_arrangement(chrom, seed, path, path_result, perc='10'):
    vcf_gDNA = pd.read_csv(os.path.join(path, f'gDNA_chroms/gDNA_consensus_chr{chrom}.vcf.gz'), header=0, sep='\t', dtype='object')
    vcf_SC = pd.read_csv(os.path.join(path, f'SC_chroms/SC_consensus_chr{chrom}.vcf.gz'), header=0, sep='\t', dtype='object')

    rd.seed(seed)
    id = vcf_SC['ID'].to_list()
    id_seed = rd.sample(id, round((len(id)*(int(perc)/100))))

    # save the snps considered in a text file
    id_masked = os.path.join(path_result, f'lists/ids_chr{chrom}_seed{seed}.txt')
    with open(id_masked,'w') as output:
        output.write(str(id_seed))

    vcf_SC_toimp = vcf_SC[~vcf_SC['ID'].isin(id_seed)]
    vcf_gDNA_ctrl = vcf_gDNA[vcf_gDNA['ID'].isin(id_seed)]
    vcf_SC_ctrl = vcf_SC[vcf_SC['ID'].isin(id_seed)]

    print(vcf_SC_ctrl.shape[0], vcf_gDNA_ctrl.shape[0])

    # save the new dataset to impute and test before imputation
    vcf_gDNA_ctrl = vcf_gDNA_ctrl.to_csv(os.path.join(path, f'gDNA_control/gDNA_cons_chr{chrom}_seed{seed}_ctrl.vcf.gz'), sep='\t', index=False)
    vcf_SC_ctrl = vcf_SC_ctrl.to_csv(os.path.join(path, f'SC_control/SC_cons_chr{chrom}_seed{seed}_ctrl.vcf.gz'), sep='\t', index=False)
    vcf_SC_toimp = vcf_SC_toimp.to_csv(os.path.join(path, f'SC_toimp/SC_cons_chr{chrom}_seed{seed}_toimp.vcf.gz'), sep='\t', index=False)




def similarity_recall_before(dictionary, seeds, path_processed, path_result, sc_samples):
    j_list_pre = []
    recall_list_pre = []
    for sample in sc_samples:
        dict_pre = []
        for chrom in [str(chrom) for chrom in list(range(1,23))]:
            for seed in seeds:
                gDNA_cons_chr_ctrl = pd.read_csv(os.path.join(path_processed, f'gDNA_control/gDNA_cons_chr{chrom}_seed{seed}_ctrl.vcf.gz'), header=0, sep='\t', dtype='object')
                SC_cons_chr_ctrl = pd.read_csv(os.path.join(path_processed, f'SC_control/SC_cons_chr{chrom}_seed{seed}_ctrl.vcf.gz'), header=0, sep='\t', dtype='object')

                vec_gDNA = gDNA_cons_chr_ctrl['gDNA_consensus'].map(dictionary).tolist()
                vec_SC = SC_cons_chr_ctrl[sample].map(dictionary).tolist()

                j_value = jaccard_score(vec_gDNA, vec_SC, average='micro')
                recall_value = recall_score(vec_gDNA, vec_SC, average='micro')

                recall_list_pre.append(recall_value)

                j_list_pre.append(j_value)
                dict_js_re = {'sample_name':sample,
                              'gDNA': 'consensus',
                              'chromosome':chrom,
                              'seed':seed,
                              'j_score':j_value,
                              'recall':recall_value,
                              'tot_snps': len(vec_SC)
                               }
                dict_pre.append(dict_js_re)

        pd.DataFrame(dict_pre).to_excel(os.path.join(path_result, f'tables/similarity_recall_{sample}_before.xlsx'), index=False)



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
        logging.info("Running module: concatenate_SC")
        concatenate_SC(files = ["SC_129", "SC_131", "SC_133", "SC_135", "SC_137",
                               "SC_139", "SC_141", "SC_143", "SC_145", "SC_147",
                               "SC_149", "SC_151", "SC_153", "SC_155", "SC_157",
                               "SC_159", "SC_161", "SC_163", "SC_165", "SC_167"],
                               path_data = data_dir,
                               output_filename = "GM12878_SC_merged.vcf.gz"
                               )
    except Exception as e:
        logging.error(f"Error in concatenate_SC: {e}")

    print('Running module: filtering_all_SC')
    filtering_all_SC(path_processed = processed_dir)

    print('Running module: splitting_chromosomes SC')
    for chromosome in chromosomes_str:
        splitting_chromosomes(kind='SC',
                                 chrom = chromosome,
                                 path_processed = processed_dir
                                 )

    print('Running module: splitting_chromosomes gDNA')
    for chromosome in chromosomes_str:
        splitting_chromosomes(kind='gDNA',
                                 chrom = chromosome,
                                 path_processed = processed_dir
                                 )

    print('Running module: dataset_arrangement')
    for chromosome in chromosomes_str:
        print(f'chromosome: {chromosome}')
        for percentage in percentages:
            for sd in seeds:
                dataset_arrangement(chrom=chromosome,
                                       seed=sd,
                                       path = processed_dir,
                                       path_result = results_dir,
                                       perc='10'
                                       )

    print('Running module: similarity_recall_before')
    similarity_recall_before(dictionary={'0/0': 0, '0/1': 1, '1/0': 1, '1/1': 2, '1/2': 3, '2/2': 4, './.': 5},
                         seeds= range(1, 3),
                         path_processed= processed_dir,
                         path_result= results_dir,
                         sc_samples= sc_samples
                         )