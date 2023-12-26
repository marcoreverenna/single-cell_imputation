#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" This functions creates a gDNA variant call file used for next steps in the analysis which no contains missing values, furthermore a similarity matrix. 

__author__ = Marco Reverenna
__copyright__ = Copyright 2022-2023
__version__ = 1.0
__maintainer__ = Marco Reverenna
__email__ = marcoreverenna@gmail.com
__status__ = Dev
"""

import os
import logging
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import combinations
from sklearn.metrics import jaccard_score


def concatenate_gDNA(file1, file2, file3, file4, file5, path_data):
    """
    Concatenates five VCF files related to GM12878 genomic DNA. Each file is read, and the last column of the 
    files (except the first one) is concatenated into a single DataFrame. The concatenated DataFrame is then saved 
    to a specified directory as a compressed VCF file. 

    Args:
    file1, file2, file3, file4, file5 (str): Filenames (without path or extension) of the VCF files to be concatenated.
    path_data (str): The directory path where the files are located and where the merged file will be saved.

    Returns:
    pandas.DataFrame: A DataFrame containing the merged data from the five VCF files.
    """
    # load bulk VCFs
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
    """
    Filters and cleans a given DataFrame, assumed to be a VCF (Variant Call Format) file. The function first 
    simplifies the data in columns after the 9th (exclusively genetic data) by keeping only the first part of 
    each string before a colon. Then, it filters the DataFrame to include only those rows where the '#CHROM' 
    column values are within a predefined list of autosomes (chromosomes_str). Finally, it maps specific 
    genetic variant representations to numerical values based on a predefined dictionary.

    Args:
    dataframe (pandas.DataFrame): A DataFrame representing VCF data.

    Returns:
    pandas.DataFrame: The filtered and cleaned DataFrame with simplified genetic data and rows corresponding 
    only to specified autosomes.
    """
    dataframe.loc[:, dataframe.columns[9:]] = dataframe.loc[:, dataframe.columns[9:]].apply(lambda x: x.str.split(':').str.get(0))
    file_filtered = dataframe[dataframe['#CHROM'].isin(chromosomes_str)]
    dictionary_array ={'0/0': 0, '0/1': 1, '1/1': 2, '1/2': 3, '2/2': 4, './.': 5}
    columns_to_map = file_filtered.columns[9:]
    file_filtered[columns_to_map] = file_filtered[columns_to_map].applymap(dictionary_array.get)
    return file_filtered


def plot_similarity_matrix(data, labels, min_value, max_value):
    """
    Generates a similarity matrix plot for genomic DNA data. The function calculates the Jaccard similarity 
    score between pairs of vectors, corresponding to different labels, in the provided data. The similarity 
    scores are calculated and stored in a matrix, which is then visualized using a heatmap. The diagonal 
    of the matrix is set to 1.0, indicating maximum similarity. The heatmap allows for adjustments in color 
    mapping based on minimum and maximum value parameters.

    Args:
    data (pandas.DataFrame): A DataFrame containing genomic data.
    labels (list): A list of labels corresponding to the columns in the DataFrame.
    min_value (float): The minimum value for the color scale in the heatmap.
    max_value (float): The maximum value for the color scale in the heatmap.

    Outputs:
    A heatmap plot saved as a PDF file, visualizing the similarity matrix based on the provided genomic data.
    """
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
    plt.savefig(f'plots/similarity_matrix_expID_MaRe_min{min_value}_max_{max_value}.pdf', bbox_inches='tight')


def dataframe_reducing(dataframe):
    """
    Reduces the given DataFrame to include only a subset of columns. Specifically, it retains the 'ID' and 'POS' 
    columns along with a predefined set of sample columns. This function is useful for simplifying a larger DataFrame 
    by extracting only essential columns needed for further analysis or visualization.

    Args:
    dataframe (pandas.DataFrame): The DataFrame from which columns are to be selected.

    Returns:
    pandas.DataFrame: A new DataFrame containing only the selected columns ('ID', 'POS', and specific sample columns).
    """
    new_df_id_pos_samples = dataframe[['ID','POS', '201666260058_R06C01', '201666260058_R06C02', '201662330194_R02C01', '201662330194_R01C01', '201662330194_R01C02']]
    return new_df_id_pos_samples




def value_more_frequent(dataframe):
    """
    Determines the most frequent value in each row of a DataFrame, starting from the third column onwards. 
    The function calculates the mode for the numeric columns in each row and adds a new column named 
    'value_more_frequent' to the DataFrame that contains these most frequent values.

    Args:
    dataframe (pandas.DataFrame): A DataFrame with numeric data from which the most frequent value per row is calculated.

    Returns:
    pandas.DataFrame: The original DataFrame augmented with a new column 'value_more_frequent', which contains the most 
    frequent value for each row.
    """
    numeric_columns = dataframe.columns[2:]

    most_frequent_numbers = dataframe[numeric_columns].apply(lambda row: np.argmax(np.bincount(row)), axis=1)
    dataframe['value_more_frequent'] = most_frequent_numbers
    return dataframe



def frequence_majority_value(row, dataframe):
    """
    Calculates the frequency of the most frequent value in a specific row of a DataFrame. The function 
    counts how many times the value stored in the 'value_more_frequent' column of a row appears across the 
    numeric columns of that row.

    Args:
    row (pandas.Series): A row from a DataFrame, including the 'value_more_frequent' value.
    dataframe (pandas.DataFrame): The DataFrame from which the row is taken, used for column reference.

    Returns:
    int: The count of how frequently the most frequent value appears in the specified row.
    """
    counts = 0
    number_more_frequent = row['value_more_frequent']
    for col in dataframe.columns[2:-1]:
        if row[col] == number_more_frequent:
            counts += 1
    return counts



def missing_counts(row, dataframe):
    """
    Counts the number of missing values (coded as '5') in specified columns of a given row in a DataFrame. 
    The function focuses on a predefined set of columns identified by their labels.

    Args:
    row (pandas.Series): A single row from the DataFrame for which the missing value count is calculated.
    dataframe (pandas.DataFrame): The DataFrame from which the row is taken.

    Returns:
    int: The count of missing values in the specified columns of the row.
    """
    labels = ['201666260058_R06C01', '201666260058_R06C02', '201662330194_R02C01', '201662330194_R01C01', '201662330194_R01C02']
    num_samples = len(labels)+ 2
    counts = 0
    for col in dataframe.columns[2:num_samples]:
        if row[col] == 5:
            counts += 1
    return counts



def include_value(row, dataframe):
    """
    Determines whether a row from a DataFrame should be included based on specific criteria. The row is 
    included ('yes') if it has no missing calls ('number_no_calls' is 0), the most frequent value is not missing 
    ('value_more_frequent' is not 5), and the frequency of the majority value is exactly 5. Otherwise, it's excluded ('no').

    Args:
    row (pandas.Series): A row from the DataFrame being evaluated.
    dataframe (pandas.DataFrame): The DataFrame from which the row is taken.

    Returns:
    str: 'yes' if the row meets the inclusion criteria, otherwise 'no'.
    """
    if row['number_no_calls'] == 0 and row['value_more_frequent'] != 5 and row['frequence_majority_value'] == 5:
        return 'yes'
    return 'no'



def convert_back(dataframe):
    """
    Converts numeric values back to their original string representations in a DataFrame. The function applies 
    this conversion only to rows marked as 'yes' in the 'include' column. It uses a predefined dictionary for 
    mapping numeric values back to their original string format.

    Args:
    dataframe (pandas.DataFrame): The DataFrame containing numeric values and an 'include' column.

    Returns:
    pandas.DataFrame: A new DataFrame with converted values, including only rows marked as 'yes' in the 'include' column.
    """
    dictionary_opp ={0 :'0/0', 1: '0/1', 2:'1/1', 3:'1/2', 4:'2/2', 5:'./.'}
    labels = ['201666260058_R06C01', '201666260058_R06C02', '201662330194_R02C01', '201662330194_R01C01', '201662330194_R01C02']
    num_samples = len(labels)+ 3
    new_df_id_pos_included = dataframe[dataframe['include']== 'yes']
    columns_to_map = new_df_id_pos_included.columns[2 : num_samples] 
    new_df_id_pos_included[columns_to_map] = new_df_id_pos_included[columns_to_map].applymap(dictionary_opp.get)
    return new_df_id_pos_included

    
def create_gDNA_vector(dataframe):
    """
    Simplifies a DataFrame to include only the first three columns and renames the third column to 
    'gDNA_consensus'. This function is specifically designed to work with genomic DNA (gDNA) data, 
    transforming a broader DataFrame into a more focused one with a key gDNA consensus column.

    Args:
    dataframe (pandas.DataFrame): The DataFrame containing genomic DNA data.

    Returns:
    pandas.DataFrame: A simplified DataFrame with only the first three columns, where the third column is 
    renamed to 'gDNA_consensus'.
    """
    dataframe = dataframe.iloc[:,:3]
    dataframe.rename(columns={'201666260058_R06C01': 'gDNA_consensus'}, inplace=True)
    return dataframe


def create_consensus_gDNA_dataframe(dataframe):
    """
    Creates a consensus genomic DNA (gDNA) DataFrame. It filters a gDNA dataset based on a list of IDs, 
    retains only autosomal chromosomes, and simplifies the data. The function reads a merged gDNA VCF file, 
    filters it for autosomal chromosomes and IDs present in the provided DataFrame, and then retains only the 
    essential columns. The columns are then renamed and simplified for further analysis.

    Args:
    dataframe (pandas.DataFrame): A DataFrame containing IDs to filter the gDNA data.

    Returns:
    pandas.DataFrame: A DataFrame with consensus gDNA data, filtered and simplified for further analysis.
    """
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
            
    new_folders = ['reference_positions', 'plots', 'reference', 'map', 'logs']
    for directory in new_folders:
        if not os.path.exists(directory):
            os.makedirs(directory)

    # Setup logging
    logging.basicConfig(filename=os.path.join('logs', 'gDNA_processing.log'), 
                        level=logging.INFO, 
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        filemode='w')  # 'w' for overwrite mode, 'a' for append mode

    try:
        logging.info("Directory setup completed.")

        # Concatenate gDNA
        try:
            file_merged = concatenate_gDNA(file1='gDNA_1', file2='gDNA_2', file3='gDNA_3', file4='gDNA_4', file5='gDNA_5', path_data=data_dir)
            logging.info("Concatenation of gDNA files completed successfully.")
        except Exception as e:
            logging.error(f"Error in concatenating gDNA files: {e}")

        # Filter autosomes and clean data
        try:
            file_filtered = filter_autosomes_and_cleaning(dataframe=file_merged)
            logging.info("Autosomes filtering and cleaning completed successfully.")
        except Exception as e:
            logging.error(f"Error in filtering autosomes and cleaning: {e}")

        # Plot similarity matrix
        try:
            plot_similarity_matrix(data=file_filtered, labels=labels, min_value=0.986, max_value=1.00)
            logging.info("Similarity matrix plot generated successfully.")
        except Exception as e:
            logging.error(f"Error in generating similarity matrix plot: {e}")

        # Data reduction
        try:
            reduced_dataframe = dataframe_reducing(dataframe=file_filtered)
            logging.info("Data reduction completed successfully.")
        except Exception as e:
            logging.error(f"Error in data reduction: {e}")

        # Frequency calculations
        try:
            dataframe_freq = value_more_frequent(dataframe=reduced_dataframe)
            dataframe_freq['frequence_majority_value'] = dataframe_freq.apply(lambda row: frequence_majority_value(row, dataframe_freq), axis=1)
            dataframe_freq['number_no_calls'] = dataframe_freq.apply(lambda row: missing_counts(row, dataframe_freq), axis=1)
            dataframe_freq['include'] = dataframe_freq.apply(lambda row: include_value(row, dataframe_freq), axis=1)
            logging.info("Frequency calculations completed successfully.")
        except Exception as e:
            logging.error(f"Error in frequency calculations: {e}")

        # Counting 'yes' and 'no' in 'include' column
        try:
            count_no = dataframe_freq[dataframe_freq['include'] == 'no'].shape[0]
            count_yes = dataframe_freq[dataframe_freq['include'] == 'yes'].shape[0]
            logging.info(f"Count of 'no' in 'include' column: {count_no}")
            logging.info(f"Count of 'yes' in 'include' column: {count_yes}")
        except Exception as e:
            logging.error(f"Error in counting 'yes' and 'no': {e}")

        # Converting back and exporting data
        try:
            dataframe_back = convert_back(dataframe=dataframe_freq)
            dataframe_back.to_excel(os.path.join(results_dir, 'tables/gDNA_dataframe_consensus_raw.xlsx'), index=False)
            gDNA_vec = create_gDNA_vector(dataframe_back)
            gDNA_vec.to_excel(os.path.join(results_dir, 'tables/gDNA_vector_consensus.xlsx'), index=False)
            final_dataframe = create_consensus_gDNA_dataframe(dataframe=dataframe_back)
            logging.info("Data conversion and export completed successfully.")
        except Exception as e:
            logging.error(f"Error in data conversion and export: {e}")

    except Exception as e:
        logging.error(f"Unexpected error in the script: {e}")