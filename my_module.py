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
import pandas as pd
import random as rd
import subprocess
from sklearn.metrics import jaccard_score
import plotly.graph_objects as go
import numpy as np
import subprocess



def download_unzip_map(url, output_folder):
    """Downloading map files necessary to run the imputation

    Args:
        url (str): url from where you can get the data
        output_folder (str): name of the output folder
    """
    try:
        os.makedirs(output_folder, exist_ok=True)
        
        # download map files using wget
        wget_command = f"wget -c --timeout=120 --waitretry=60 --tries=10000 --retry-connrefused -P {output_folder} {url}"
        subprocess.run(wget_command, shell=True, check=True)

        # decompress with unzip
        zip_file_path = os.path.join(output_folder, os.path.basename(url))
        unzip_command = f"unzip {zip_file_path} -d {output_folder}"
        subprocess.run(unzip_command, shell=True, check=True)

        # remove zip file 
        os.remove(zip_file_path)
        print(f"Download and unzip successfully completed in: {output_folder}")

    except Exception as e:
        print(f"Error during downloading and unzipping: {str(e)}")



def download_chromosomes_vcf():
    """ Downloading beagle reference files for imputation
    """

    # check if the folder already exists
    if not os.path.exists("reference"):
        os.makedirs("reference")

    for chrom in range(1, 23):
        vcf_filename = f"chr{chrom}.1kg.phase3.v5a.vcf.gz"
        vcf_filepath = os.path.join("reference",
                                    vcf_filename
                                    )

        # check if the file is already inside the directory
        if os.path.exists(vcf_filepath):
            print(f"file {vcf_filename} is already present. Skipping the download...")
        else:
            url = f"https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/{vcf_filename}"
            wget_command = ["wget", "-c", "--timeout=120", "--waitretry=60", "--tries=10000", "--retry-connrefused", "-P", "reference", url]
            try:
                subprocess.run(wget_command, check=True)
                print(f"Download completato per {vcf_filename}")
            except subprocess.CalledProcessError:
                print(f"Error during the download of {vcf_filename}")



def get_snparray_positions(path_vcf):
    print('Running module: get_snparray_positions')
    """_summary_

    Args:
        path_vcf (_type_): _description_
    """
    output_positions_path = 'reference_positions/'
    
    input_vcf_file = pd.read_csv(path_vcf, header=0, sep='\t', dtype='object')

    chromosome_positions = {}

    for index, row in input_vcf_file.iterrows():
        chromosome = row['#CHROM']
        position = row['POS']
        
        # this way I can ignore X and Y
        if chromosome in ('X','Y'):
            continue

        if chromosome not in chromosome_positions:
            chromosome_positions[chromosome] = []
        chromosome_positions[chromosome].append(position)

    for chromosome, positions in chromosome_positions.items():
        output_file = os.path.join(output_positions_path, f"chromosome_position_{chromosome}_snparray.txt")
        with open(output_file, "w") as out_file:
            for position in positions:
                out_file.write(position + "\n")

            
def get_reference_positions():
    """_summary_
    """
    print('Running module: get_reference_positions')
    for chrom in range(1, 23):
        filename = f"reference/chr{chrom}.1kg.phase3.v5a.vcf.gz"
        print(f"Processing {filename}")

        cmd = f"zcat {filename} | grep -v '^#' | cut -f2"
        try:
            result = subprocess.run(cmd, shell=True, check=True, text=True, stdout=subprocess.PIPE)
            output = result.stdout

            output_filename = f"reference_positions/chromosome_position_{chrom}_reference.txt"

            # save the file
            with open(output_filename, "w") as output_file:
                output_file.write(output)

            print(f"Result saved in {output_filename}")
        except subprocess.CalledProcessError as e:
            print(f"Execution error in {filename}")
        

def difference_positions():
    """_summary_
    """
    print('Running module: difference_positions')

    for chromosome in list(range(1,23)):
        file_reference = f'reference_positions/chromosome_position_{chromosome}_reference.txt'
        file_snparray = f'reference_positions/chromosome_position_{chromosome}_snparray.txt'
        output_filename = f"reference_positions/chromosome_position_{chromosome}_difference.txt"

        # Leggi i numeri dai file di riferimento e array SNP
        with open(file_reference, 'r') as reference_file, open(file_snparray, 'r') as snparray_file:
            reference_numbers = set(line.strip() for line in reference_file)
            snparray_numbers = set(line.strip() for line in snparray_file)

        missing_positions = reference_numbers - snparray_numbers

        with open(output_filename, 'w') as output_file:
            for pos in missing_positions:
                # enter the right format chromosome:number necessary for beagle algorithm to exclude markers
                formatted_number = f"{chromosome}:{pos}"
                output_file.write(formatted_number + '\n')

        print(f"Differences saved in {output_filename}")


def concatenate_vcf(file1, file2, file3, file4, file5, kind, path_data):
    print('Running module: concatenate_vcf')
    """ This function creates a merged file both for single cell and bulk
    Args:
        file1 (str): gDNA_1 or SC_129
        file2 (str): gDNA_2 or SC_130
        file3 (str): gDNA_3 or SC_131
        file4 (str): gDNA_4 or SC_132
        file5 (str): gDNA_5 or SC_133
        kind (str): either 'gDNA' or 'SC' data
    """
    # load vcfs
    vcf_1 = pd.read_csv(os.path.join(path_data, f'raw/GM12878_{file1}.vcf.gz'), header=30, sep='\t', dtype='object')
    vcf_2 = pd.read_csv(os.path.join(path_data, f'raw/GM12878_{file2}.vcf.gz'), header=30, sep='\t', dtype='object')
    vcf_3 = pd.read_csv(os.path.join(path_data, f'raw/GM12878_{file3}.vcf.gz'), header=30, sep='\t', dtype='object')
    vcf_4 = pd.read_csv(os.path.join(path_data, f'raw/GM12878_{file4}.vcf.gz'), header=30, sep='\t', dtype='object')
    vcf_5 = pd.read_csv(os.path.join(path_data, f'raw/GM12878_{file5}.vcf.gz'), header=30, sep='\t', dtype='object')
    # merge vcf files
    file_merged = pd.concat([vcf_1, vcf_2.iloc[:,-1], vcf_3.iloc[:,-1], vcf_4.iloc[:,-1], vcf_5.iloc[:,-1]], axis=1)
    # save to processed directory
    file_merged.to_csv(os.path.join(path_data, f'processed/GM12878_merged_{kind}.vcf.gz'), sep='\t', index=False)


def consensus_gDNA_dataframe(path_processed):
    print('Running module: consensus_gDNA_dataframe')
    """_summary_

    Args:
        path_processed (_type_): _description_
    """
    df = pd.read_csv(os.path.join(path_processed, f'GM12878_merged_gDNA.vcf.gz'), header=0, sep='\t', dtype='object')
    df.loc[:, df.columns[9:]] = df.loc[:, df.columns[9:]].apply(lambda x: x.str.split(':').str.get(0))
    
    chromosomes_int = list(range(1, 23))
    chromosomes_str = [str(chrom) for chrom in chromosomes_int]
    new_df = df[df['#CHROM'].isin(chromosomes_str)]
    
    dictionary = {'0/0': 0, '0/1': 1, '1/1': 2, '1/2': 3, '2/2': 4, './.': 5}
    columns_to_map = new_df.columns[9:]  
    new_df[columns_to_map] = new_df[columns_to_map].applymap(dictionary.get)
    
    filtered_df = new_df[new_df.iloc[:, 9:].apply(lambda row: row.nunique() == 1, axis=1)]
    pos_list_proxy_gDNA = filtered_df['POS'].tolist()
    
    df_consensus = pd.read_csv(os.path.join(path_processed, f'GM12878_merged_gDNA.vcf.gz'), header=0, sep='\t', dtype='object')
    df_consensus_autosomes = df_consensus[df_consensus['#CHROM'].isin(chromosomes_str)]
    df_consensus_filtered = df_consensus_autosomes[df_consensus_autosomes['POS'].isin(pos_list_proxy_gDNA)]
    
    df_consensus_filtered.to_csv(os.path.join(path_processed, f'GM12878_consensus_gDNA.vcf.gz'), sep='\t', index=False)
    
    
    
def filtering_all_SC(path_processed):
    print('Running module: filtering_all_SC')
    """ This function filters the same positions of consensus gDNA on SC datfarames
    """
    chromosomes_str = [str(chrom) for chrom in list(range(1,23))]

    df_consensus_filtered = pd.read_csv(os.path.join(path_processed,f'GM12878_consensus_gDNA.vcf.gz'), header=0, sep='\t', dtype='object')
    df_sc = pd.read_csv(os.path.join(path_processed, f'GM12878_merged_SC.vcf.gz'), header=0, sep='\t', dtype='object')
    df_sc_autosomes = df_sc[df_sc['#CHROM'].isin(chromosomes_str)]
    print(df_consensus_filtered.shape[0], df_sc_autosomes.shape[0])
    
    positions_list = df_consensus_filtered['POS'].tolist()
    df_sc_filtered = df_sc_autosomes[df_sc_autosomes['POS'].isin(positions_list)]
    print(df_consensus_filtered.shape[0], df_sc_filtered.shape[0])
    df_sc_filtered.to_csv(os.path.join(path_processed, f'GM12878_merged_SC.vcf.gz'), sep='\t', index=False)
    

def splitting_chromosomes(kind, chrom, path_data, directory):
    print('Running module: splitting_chromosomes')
    """ This function split the dataframe for each chromosome and type of data (gDNA or SC)

    Args:
        kind (str): either consensus_gDNA or merged_SC
        chrom (str): list of chromosomes in string format
        directory (str): path
    """
    # load the dataset
    # vcf = pd.read_csv(os.path.join(data_dir, f'processed/GM12878_merged_{kind}.vcf.gz'), header=0, sep='\t', dtype='object')
    vcf = pd.read_csv(os.path.join(path_data, f'processed/GM12878_{kind}.vcf.gz'), header=0, sep='\t', dtype='object')
    # split for each chromosome
    vcf_chr = vcf[vcf['#CHROM'] == chrom]
    # rename the combinations and set the right format
    vcf_chr.columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'pairing_1', 'pairing_2', 'pairing_3', 'pairing_4', 'pairing_5']
    vcf_chr['FORMAT'] = 'GT'
    # keep only the genotypes
    #vcf_chr[:,9:] = vcf_chr[:,9:].apply(lambda x : x.str.split(':').str.get(0))
    vcf_chr.loc[:, vcf_chr.columns[9:]] = vcf_chr.loc[:, vcf_chr.columns[9:]].apply(lambda x: x.str.split(':').str.get(0))
    # save splitted chromosome for each dataset
    # vcf_chr.to_csv(os.path.join(directory, f'GM12878_merged_{kind}_chr{chrom}.vcf.gz'), sep='\t', index=False)
    vcf_chr.to_csv(os.path.join(directory, f'GM12878_{kind}_chr{chrom}.vcf.gz'), sep='\t', index=False)



def dataset_arrangement(chrom, seed, perc, path_result, path_data):
    print('Running module: dataset_arrangement')
    """This function  preparares the VCF dataframes to test and check

    Args:
        chrom (str): list of chromosomes
        seed (int): a list of integers
        perc (str): it is recoomended to use 10 and 20, if you mask more snps imputation is not garanteed
    """
    
    vcf_gDNA = pd.read_csv(os.path.join(f'data/processed/gDNA_chroms/GM12878_consensus_gDNA_chr{chrom}.vcf.gz'), header=0, sep='\t', dtype='object')
    vcf_SC = pd.read_csv(os.path.join(f'data/processed/SC_chroms/GM12878_merged_SC_chr{chrom}.vcf.gz'), header=0, sep='\t', dtype='object')
    # set the list of snps considering the package random, module seed 
    rd.seed(seed)
    pos = vcf_SC['POS'].to_list()
    pos_seed = rd.sample(pos, round((len(pos)*(int(perc)/100))))
    # save the snps considered in a text file
    pos_masked = os.path.join(path_result, f'lists/positions_chr{chrom}_seed{seed}_perc{perc}.txt')
    with open(pos_masked,'w') as output:
        output.write(str(pos_seed))
    vcf_SC_toimp = vcf_SC[~vcf_SC['POS'].isin(pos_seed)]
    vcf_gDNA_ctrl = vcf_gDNA[vcf_gDNA['POS'].isin(pos_seed)]
    vcf_SC_ctrl = vcf_SC[vcf_SC['POS'].isin(pos_seed)]
    # little check, shapes must be of the same lenght
    print(vcf_SC_ctrl.shape[0], vcf_gDNA_ctrl.shape[0])
    # save the new dataset to impute and test before imputation
    vcf_gDNA_ctrl = vcf_gDNA_ctrl.to_csv(os.path.join(path_data, f'processed/gDNA_control/GM12878_gDNA_merged_chr{chrom}_seed{seed}_perc{perc}_ctrl.vcf.gz'), sep='\t', index=False)
    vcf_SC_ctrl = vcf_SC_ctrl.to_csv(os.path.join(path_data, f'processed/SC_control/GM12878_SC_merged_chr{chrom}_seed{seed}_perc{perc}_ctrl.vcf.gz'), sep='\t', index=False)
    vcf_SC_toimp = vcf_SC_toimp.to_csv(os.path.join(path_data, f'processed/SC_toimp/GM12878_SC_merged_chr{chrom}_seed{seed}_perc{perc}_toimp.vcf.gz'), sep='\t', index=False)



def similarity_before(dictionary, seeds, percentages, path_processed, path_result):
    print('Running module: similarity_before')
    """This function calculate the similarity before the imputation

    Args:
        dictionary (dict): dictionary to map genotypes - should be like 1/1:1, 1/0:0...
    """

    # create an empty list to fill whit all jaccard scores
    j_list_pre = []
    pairings = ['pairing_1', 'pairing_2', 'pairing_3', 'pairing_4', 'pairing_5']
    for match in pairings:    
        dict_pre = []
        for chrom in [str(chrom) for chrom in list(range(1,23))]:
            for seed in seeds:
                for perc in percentages:
                    # load your dataframes splitted for each chromosome, seed and percentage
                    gDNA_merg_chr_ctrl = pd.read_csv(os.path.join(path_processed, f'gDNA_control/GM12878_gDNA_merged_chr{chrom}_seed{seed}_perc{perc}_ctrl.vcf.gz'), header=0, sep='\t', dtype='object')
                    SC_merg_chr_ctrl = pd.read_csv(os.path.join(path_processed, f'SC_control/GM12878_SC_merged_chr{chrom}_seed{seed}_perc{perc}_ctrl.vcf.gz'), header=0, sep='\t', dtype='object')
                    # define two vectors to compare
                    vec_gDNA = gDNA_merg_chr_ctrl[match].map(dictionary).tolist()
                    vec_SC = SC_merg_chr_ctrl[match].map(dictionary).tolist()
                    j_value = jaccard_score(vec_gDNA, vec_SC, average='micro')
                    j_list_pre.append(j_value)
                    dict_js = {'sample':match,
                               'chromosome':chrom,
                               'seed':seed,
                               'perc':perc,
                               'j_score':j_value,
                               'tot_snps': len(vec_SC)}
                    dict_pre.append(dict_js)
        # save the jaccard scores obtained for each pair/match
        pd.DataFrame(dict_pre).to_excel(os.path.join(path_result, f'tables/jaccard_scores_{match}_before.xlsx'), index=False)
     
        

def imputation(chrom, seed, perc, path_reference, path_map, path_data, path_output, path_positions):
    print('Running module: imputation')
    """This function run beagle imputation algorithm

    Args:
        chrom (str): list of chromosomes to impute in string format 
        seed (int): a list of integers
        perc (str): it is recoomended to use 10 and 20, or a number not greater than 20
    """
    # if Beagle does not work probably there is a problem with the paths, please check them out
    ref = os.path.join(path_reference, f'chr{chrom}.1kg.phase3.v5a.vcf.gz ')
    map = os.path.join(path_map, f'plink.chr{chrom}.GRCh37.map ')
    gt = os.path.join(path_data, f'processed/SC_toimp/GM12878_SC_merged_chr{chrom}_seed{seed}_perc{perc}_toimp.vcf.gz ')
    out = os.path.join(path_output, f'GM12878_SC_merged_chr{chrom}_seed{seed}_perc{perc}imputed')
    positions = os.path.join(path_positions, f'chromosome_position_{chrom}_difference.txt ')
    print(f'Imputation of single cell GM12878 chromosome:{chrom}, seed:{seed} and percentage:{perc}% is started, take a little break...')
    command = "java -jar beagle.22Jul22.46e.jar " + f"excludemarkers={positions}" + f"ref={ref}" + f"gt={gt}" + f"map={map}" + f"out={out}"
    result = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()



def similarity_after(dictionary_imp, dictionary_gen, seeds, percentages, path_output, path_processed, path_result):
    print('Running module: similarity_after')

    """This function calculates the jaccard score after filtering the dataframes.

    Args:
        dictionary_imp (dict): it contains genotypes associated with integers and separeted by '|'
        dictionary_gen (dict): it contains genotypes associated with integers and separeted by '/'
    """
    j_list_post = []
    
    pairings = ['pairing_1', 'pairing_2', 'pairing_3', 'pairing_4', 'pairing_5']

    for match in pairings:
        dict_post = []
        print(f'Analysing the pair: {match}')
        for chrom in list(range(1,23)):
            for seed in seeds:
                for perc in percentages:
                    # load your dataset
                    SC_merg_chr_imputed = pd.read_csv(os.path.join(path_output, f'GM12878_SC_merged_chr{chrom}_seed{seed}_perc{perc}imputed.vcf.gz'), header=8, sep='\t', dtype='object')
                    gDNA_merg_chr_ctrl = pd.read_csv(os.path.join(path_processed, f'gDNA_control/GM12878_gDNA_merged_chr{chrom}_seed{seed}_perc{perc}_ctrl.vcf.gz'), header=0, sep='\t', dtype='object')   
                    
                    # change the format
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
                    print(f'analising singlecell chromosome:{chrom} seed:{seed} percentage:{perc} imputed')         
                    
                    # define two vectors to compare and then calculate the jaccard score considering 'micro' option
                    vec_gDNA_final = gDNA_merg_chr_ctrl_nomiss_filtered[match].map(dictionary_gen).tolist()
                    vec_SC_final = SC_merg_chr_imputed_filtered[match].map(dictionary_imp).tolist()
                    j_value = jaccard_score(vec_gDNA_final, vec_SC_final, average='micro')
                    j_list_post.append(j_value)
                    dict_js = {'sample':match,
                               'chromosome':chrom,
                               'seed':seed,
                               'perc':perc,
                               'j_score':j_value,
                               'tot_snps':len(vec_SC_final)}
                    dict_post.append(dict_js)       
        # collect all the results inside a dataframe            
        pd.DataFrame(dict_post).to_excel(os.path.join(path_result, f'tables/jaccard_scores_{match}_after.xlsx'), index=False)



def concat_tables(path_result):
    print('Running module: concat_tables')
    """ This function concatenate tables of jaccard scores before and after imputation

    Args:
        results_dir (str): result path to directory
        pairings (str): names of couple created to match gDNA and SC
    """
    pairings = ['pairing_1', 'pairing_2', 'pairing_3', 'pairing_4', 'pairing_5']
    # tables before imputation
    file_names_pre = [os.path.join(path_result,f'tables/jaccard_scores_{match}_before.xlsx') for match in pairings]
    dfs_pre = [pd.read_excel(file) for file in file_names_pre]
    result_pre = pd.concat(dfs_pre)
    # tables after imputation
    file_names_post = [os.path.join(path_result,f'tables/jaccard_scores_{match}_after.xlsx') for match in pairings]
    dfs_post = [pd.read_excel(file) for file in file_names_post] 
    result_post = pd.concat(dfs_post)
    # prepare a concatenate table
    result_merged = result_post.copy()
    result_merged['j_score_pre'] = result_pre['j_score']
    result_merged['tot_snps_pre'] = result_pre['tot_snps']
    result_merged.rename(columns={'j_score': 'j_score_post', 'tot_snps':'tot_snps_post'}, inplace=True)
    result_merged = result_merged[['sample','chromosome','seed','perc','j_score_pre','j_score_post', 'tot_snps_pre', 'tot_snps_post']]
    result_merged.to_excel(os.path.join(path_result, f'tables/jaccard_scores_concatenated.xlsx'), index=False)
  
  
    
def create_violin_plot(results_dir, plots_dir, percentage):
    print('Running module: create_violin_plot')
    """ This function creates violin plot to show differences between before and after imputation 

    Args:
        table (dataframe): excel table containing jaccard scores before and after imputation
        plot_dir (str): path
        percentage (int): either 10 or 20
    """
    # start an empty figure
    fig = go.Figure()
    # load the table with the concatenated results
    table = pd.read_excel(os.path.join(results_dir, 'tables/jaccard_scores_concatenated.xlsx'))
    # create the violin plots
    fig.add_trace(go.Violin(x=table[table['perc']==percentage]['chromosome'], y=table[table['perc']==percentage]['j_score_pre'],
                            legendgroup='PRE', scalegroup='PRE', name='PRE IMPUTATION', side='negative', line_color='red'))
    fig.add_trace(go.Violin(x=table[table['perc']==percentage]['chromosome'], y=table[table['perc']==percentage]['j_score_post'],
                            legendgroup='POST', scalegroup='POST', name='POST IMPUTATION', side='positive', line_color='green'))
    fig.update_traces(meanline_visible=True)
    fig.update_layout(violingap=0, violinmode='overlay', 
                      title=f'Distribution of jaccard scores before and after imputation by chromosomes with {percentage} masked SNPs',
                      xaxis=dict(tickmode='array', tickvals=np.arange(1, 23)),
                      yaxis=dict(range=[0.4, 1]),
                      xaxis_title='Chromosomes',
                      yaxis_title='Jaccard scores',
                      width=1200,  # Larghezza della figura
                      height=400,
                      legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1))
    # save your fig
    pdf_file_path = os.path.join(plots_dir, f'figure_violinplot_samples_perc{percentage}.pdf')
    fig.write_image(pdf_file_path)
    

def calculate_frequencies(percentages, path_result):
    print('Running module: calculate_frequencies')
    """_summary_

    Args:
        percentages (_type_): _description_
    """
    table = pd.read_excel(os.path.join(path_result, 'tables/jaccard_scores_concatenated.xlsx'))
    # create relative frequencies
    table['ratio'] = table['tot_snps_post'] / table['tot_snps_pre']
    table['snps_lost'] = 1 - table['ratio']
    for perc in percentages:
        filtered_table = table[table['perc'] == perc]  
        filtered_table = filtered_table.groupby('chromosome')[['snps_lost']].mean()
        filtered_table.to_excel(os.path.join(path_result, 'statistics', f'statistics_frequencies_{perc}.xlsx'))



def calculate_statistics(group, path_result):
    print('Running module: calculate_statistics')
    """ This function creates to calculate same basics statistics like mean, median, max, min and standard deviation

    Args:
        group (str): either sample or chromosome
    """
    table = pd.read_excel(os.path.join(path_result, 'tables/jaccard_scores_concatenated.xlsx'))
    percentages = ['10','20']
    for perc in [int(perc) for perc in percentages]:
        filtered_table = table[table['perc'] == perc]
        # calculate the mean
        mean_scores = filtered_table.groupby(group)[['j_score_pre', 'j_score_post']].mean()
        mean_scores = mean_scores.rename(columns={'j_score_post': 'mean_scores_post', 'j_score_pre': 'mean_scores_pre'})
        # calculate the median
        median_scores = filtered_table.groupby(group)[['j_score_pre', 'j_score_post']].median()
        median_scores = median_scores.rename(columns={'j_score_post': 'median_scores_post', 'j_score_pre': 'median_scores_pre'})
        # calculate the maximum
        max_scores = filtered_table.groupby(group)[['j_score_pre', 'j_score_post']].max()
        max_scores = max_scores.rename(columns={'j_score_post': 'max_scores_post', 'j_score_pre': 'max_scores_pre'})
        # calculate the mainimum
        min_scores = filtered_table.groupby(group)[['j_score_pre', 'j_score_post']].min()
        min_scores = min_scores.rename(columns={'j_score_post': 'min_scores_post', 'j_score_pre': 'min_scores_pre'})
        # calculate the standard deviation
        std_scores = filtered_table.groupby(group)[['j_score_pre', 'j_score_post']].agg(lambda x: round(np.std(x, ddof=0), 3))
        std_scores = std_scores.rename(columns={'j_score_post': 'std_j_score_post', 'j_score_pre': 'std_j_score_pre'})
    
        mean_scores.to_excel(os.path.join(path_result, 'statistics', f'mean_{group}_perc{perc}.xlsx'))#index=False
        median_scores.to_excel(os.path.join(path_result, 'statistics', f'median_{group}_perc{perc}.xlsx'))
        max_scores.to_excel(os.path.join(path_result, 'statistics', f'max_{group}_perc{perc}.xlsx'))
        min_scores.to_excel(os.path.join(path_result, 'statistics', f'min_{group}_perc{perc}.xlsx'))
        std_scores.to_excel(os.path.join(path_result, 'statistics', f'std_{group}_perc{perc}.xlsx'))