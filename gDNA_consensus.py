#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Questo script crea una matrice di similarita' in riferimento al gDNA, inoltre crea un vettore gDNA privo di missing values 

__author__ = Marco Reverenna
__copyright__ = Copyright 2022-2023
__version__ = 1.0
__maintainer__ = Marco Reverenna
__email__ = marcoreverenna@gmail.com
__status__ = Dev
"""

import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import combinations
from sklearn.metrics import jaccard_score



def concatenate_gDNA(file1, file2, file3, file4, file5, path_data):
    # load vcfs
    vcf_1 = pd.read_csv(os.path.join(path_data, f'raw/GM12878_{file1}.vcf.gz'), header=30, sep='\t', dtype='object')
    vcf_2 = pd.read_csv(os.path.join(path_data, f'raw/GM12878_{file2}.vcf.gz'), header=30, sep='\t', dtype='object')
    vcf_3 = pd.read_csv(os.path.join(path_data, f'raw/GM12878_{file3}.vcf.gz'), header=30, sep='\t', dtype='object')
    vcf_4 = pd.read_csv(os.path.join(path_data, f'raw/GM12878_{file4}.vcf.gz'), header=30, sep='\t', dtype='object')
    vcf_5 = pd.read_csv(os.path.join(path_data, f'raw/GM12878_{file5}.vcf.gz'), header=30, sep='\t', dtype='object')
    # merge vcf files
    file_merged = pd.concat([vcf_1, vcf_2.iloc[:,-1], vcf_3.iloc[:,-1], vcf_4.iloc[:,-1], vcf_5.iloc[:,-1]], axis=1)
    # save to processed directory
    file_merged.to_csv(os.path.join(path_data, f'processed/GM12878_merged_gDNA.vcf.gz'), sep='\t', index=False)
    return file_merged
    


def filter_autosomes_and_cleaning(dataframe):
    dataframe.loc[:, dataframe.columns[9:]] = dataframe.loc[:, dataframe.columns[9:]].apply(lambda x: x.str.split(':').str.get(0))
    file_filtered = dataframe[dataframe['#CHROM'].isin(chromosomes_str)]
    dictionary_array ={'0/0': 0, '0/1': 1, '1/1': 2, '1/2': 3, '2/2': 4, './.': 5}
    columns_to_map = file_filtered.columns[9:]  
    file_filtered[columns_to_map] = file_filtered[columns_to_map].applymap(dictionary_array.get)
    return file_filtered


def plot_similarity_matrix(data, labels, min_value, max_value):
    result_matrix = [[0.0 for _ in range(len(labels))] for _ in range(len(labels))]

    for i, j in combinations(range(len(labels)), 2):
        vector1 = data[labels[i]].astype(int)
        vector2 = data[labels[j]].astype(int)
        jaccard = jaccard_score(vector1, vector2, average='micro')
        result_matrix[i][j] = round(jaccard, 4)
        result_matrix[j][i] = round(jaccard, 4)

    for i in range(len(labels)):
        result_matrix[i][i] = 1.0

    result_df = pd.DataFrame(result_matrix, columns=labels, index=labels)

    cmap_custom = sns.diverging_palette(10, 240, as_cmap=True)

    fig_width = 6 
    fig_height = 5  
    plt.figure(figsize=(fig_width, fig_height))

    sns.heatmap(result_df, annot=True, fmt='.4f', cmap=cmap_custom, linewidths=0.75, vmin=min_value, vmax=max_value)
    plt.title('Similarity matrix in gDNA data', fontsize=10)  
    plt.xticks(fontsize=8)  
    plt.yticks(fontsize=8) 

    plt.subplots_adjust(top=0.9, left=0.2, right=0.9, bottom=0.2)
    plt.savefig('plots/similarity_matrix_expID_MaRe_min{min_value}_max_{max_value}.pdf', bbox_inches='tight')


def dataframe_reducing(dataframe):
    new_df_id_pos_samples = dataframe[['ID','POS', '201666260058_R06C01', '201666260058_R06C02', '201662330194_R02C01', '201662330194_R01C01', '201662330194_R01C02']]
    
    return new_df_id_pos_samples



def value_more_frequent(dataframe):
    numeric_columns = dataframe.columns[2:]

    most_frequent_numbers = dataframe[numeric_columns].apply(lambda row: np.argmax(np.bincount(row)), axis=1)

    dataframe['value_more_frequent'] = most_frequent_numbers

    return dataframe


def frequence_majority_value(row, dataframe):
    counts = 0
    number_more_frequent = row['value_more_frequent']
    for col in dataframe.columns[2:-1]:
        if row[col] == number_more_frequent:
            counts += 1
    return counts


def missing_counts(row, dataframe):
    labels = ['201666260058_R06C01', '201666260058_R06C02', '201662330194_R02C01', '201662330194_R01C01', '201662330194_R01C02']
    num_samples = len(labels)+ 2
    counts = 0
    for col in dataframe.columns[2:num_samples]:
        if row[col] == 5:
            counts += 1
    return counts


def include_value(row, dataframe):
    if row['number_no_calls'] == 0 and row['value_more_frequent'] != 5 and row['frequence_majority_value'] == 5:
        return 'yes'
    return 'no'


def convert_back(dataframe):
    dictionary_opp ={0 :'0/0', 1: '0/1', 2:'1/1', 3:'1/2', 4:'2/2', 5:'./.'}

    labels = ['201666260058_R06C01', '201666260058_R06C02', '201662330194_R02C01', '201662330194_R01C01', '201662330194_R01C02']

    num_samples = len(labels)+ 3

    new_df_id_pos_included = dataframe[dataframe['include']== 'yes']

    columns_to_map = new_df_id_pos_included.columns[2 : num_samples] 

    new_df_id_pos_included[columns_to_map] = new_df_id_pos_included[columns_to_map].applymap(dictionary_opp.get)
    
    return new_df_id_pos_included
    
def create_gDNA_vector(dataframe):    
    dataframe = dataframe.iloc[:,:3]
    dataframe.rename(columns={'201666260058_R06C01': 'gDNA_consensus'}, inplace=True)
    return dataframe

def create_consensus_gDNA_dataframe(dataframe):
    id_list_proxy_gDNA = dataframe['ID'].tolist()
    
    df_gDNA = pd.read_csv(f'data/processed/GM12878_merged_gDNA.vcf.gz', header=0, sep='\t', dtype='object')

    chromosomes_int = list(range(1, 23))
    chromosomes_str = [str(chrom) for chrom in chromosomes_int]
    df_gDNA_autosomes = df_gDNA[df_gDNA['#CHROM'].isin(chromosomes_str)]
    
    df_consensus_filtered = df_gDNA_autosomes[df_gDNA_autosomes['ID'].isin(id_list_proxy_gDNA)]
    
    df_consensus_filtered.loc[:, df_consensus_filtered.columns[9:]] = df_consensus_filtered.loc[:, df_consensus_filtered.columns[9:]].apply(lambda x: x.str.split(':').str.get(0))

    df_consensus_filtered = df_consensus_filtered.iloc[:,:10]

    df_consensus_filtered.rename(columns={'201666260058_R06C01': 'gDNA_consensus'}, inplace=True)

    df_consensus_filtered['FORMAT'] = 'GT'  

    #for i in range(19):
    #    df_consensus_filtered[f'gDNA_consensus_{i+1}'] = df_consensus_filtered['gDNA_consensus']
    
    df_consensus_filtered.to_csv(os.path.join(processed_dir,'gDNA_consensus_dataframe.vcf.gz'), sep='\t', index=False)
    
    return df_consensus_filtered



if __name__ == "__main__":
    
    # sets paths
    data_dir = 'data/'
    map_dir = 'map/'
    plots_dir = 'plots/'
    results_dir = 'results/'
    output_dir = 'data/output/'
    processed_dir = 'data/processed/'
    
    # set variables
    chromosomes_int = list(range(1, 23))
    chromosomes_str = [str(chrom) for chrom in chromosomes_int]
    labels = ['201666260058_R06C01', '201666260058_R06C02', '201662330194_R02C01', '201662330194_R01C01', '201662330194_R01C02']

    # create sub-folder inside 
    directories_results = ['lists', 'tables', 'statistics']
    for directory in directories_results:
        dir_path = os.path.join(results_dir, directory)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
    
    # create the directories to store new processed data and final outcomes
    directories_data = ['processed', 'output']
    for directory in directories_data:
        dir_path = os.path.join(data_dir, directory)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
            
    new_folders = ['reference_positions', 'plots', 'reference', 'map']
    for directory in new_folders:
        if not os.path.exists(directory):
            os.makedirs(directory)
    
    # calling the functions
    file_merged = concatenate_gDNA(file1= 'gDNA_1', file2='gDNA_2', file3='gDNA_3', file4='gDNA_4', file5='gDNA_5', path_data= data_dir)    
    
    file_filtered = filter_autosomes_and_cleaning(dataframe=file_merged)
    
    plot_similarity_matrix(data=file_filtered, labels = labels, min_value = 0.986, max_value = 1.00)
    
    reduced_dataframe = dataframe_reducing(dataframe=file_filtered)
    
    dataframe_freq = value_more_frequent(dataframe=reduced_dataframe)
    
    dataframe_freq['frequence_majority_value'] = dataframe_freq.apply(lambda row: frequence_majority_value(row, dataframe_freq), axis=1)
        
    dataframe_freq['number_no_calls'] = dataframe_freq.apply(lambda row: missing_counts(row, dataframe_freq), axis=1)
    
    dataframe_freq['include'] = dataframe_freq.apply(lambda row: include_value(row, dataframe_freq), axis=1)
    
    print("Count of 'no' in 'include' column:", dataframe_freq[dataframe_freq['include'] == 'no'].shape[0])
    print("Count of 'yes' in 'include' column:", dataframe_freq[dataframe_freq['include'] == 'yes'].shape[0])

    dataframe_back = convert_back(dataframe=dataframe_freq)
    
    dataframe_back.to_excel(os.path.join(results_dir, 'tables/gDNA_dataframe_consensus_raw.xlsx'), index=False)

    gDNA_vec = create_gDNA_vector(dataframe_back)
    
    gDNA_vec.to_excel(os.path.join(results_dir, 'tables/gDNA_vector_consensus.xlsx'), index=False)
    
    final_dataframe = create_consensus_gDNA_dataframe(dataframe=dataframe_back)