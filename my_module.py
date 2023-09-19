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
import numpy as np
from sklearn.metrics import jaccard_score
from sklearn.metrics import recall_score
#import plotly.graph_objects as go
import plotly.subplots as sp
import plotly.graph_objs as go



pd.options.mode.chained_assignment = None


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
    """ This function creates txt files containing positions for each SNP within SNP array data.

    Args:
        path_vcf (_type_): path to reach SNP array data
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



def concatenate_SC(files, path_data, output_filename):
    """_summary_

    Args:
        files (_type_): _description_
        path_data (_type_): _description_
        output_filename (_type_): _description_
    """
    

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
    """_summary_

    Args:
        path_processed (_type_): _description_
    """
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
    """_summary_

    Args:
        kind (_type_): _description_
        chrom (_type_): _description_
        path_processed (_type_): _description_
    """

    vcf = pd.read_csv(os.path.join(path_processed, f'{kind}_consensus_dataframe.vcf.gz'), header=0, sep='\t', dtype='object')
    vcf_chr = vcf[vcf['#CHROM'] == chrom]
    vcf_chr['FORMAT'] = 'GT'
    vcf_chr.loc[:, vcf_chr.columns[9:]] = vcf_chr.loc[:, vcf_chr.columns[9:]].apply(lambda x: x.str.split(':').str.get(0))
    vcf_chr.to_csv(os.path.join(path_processed, f'{kind}_chroms/{kind}_consensus_chr{chrom}.vcf.gz'), sep='\t', index=False)



def dataset_arrangement(chrom, seed, path, path_result, perc='10'):
    """_summary_

    Args:
        chrom (_type_): _description_
        seed (_type_): _description_
        path (_type_): _description_
        path_result (_type_): _description_
        perc (str, optional): _description_. Defaults to '10'.
    """

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
    """_summary_

    Args:
        dictionary (_type_): _description_
        seeds (_type_): _description_
        path_processed (_type_): _description_
        path_result (_type_): _description_
        sc_samples (_type_): _description_
    """
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
     


def imputation(chrom, seed, path_reference, path_map, path_data, path_output, path_positions):
    """_summary_

    Args:
        chrom (_type_): _description_
        seed (_type_): _description_
        path_reference (_type_): _description_
        path_map (_type_): _description_
        path_data (_type_): _description_
        path_output (_type_): _description_
        path_positions (_type_): _description_
    """
    # if Beagle does not work probably there is a problem with the paths, please check them out
    ref = os.path.join(path_reference, f'chr{chrom}.1kg.phase3.v5a.vcf.gz ')
    map = os.path.join(path_map, f'plink.chr{chrom}.GRCh37.map ')
    gt = os.path.join(path_data, f'processed/SC_toimp/SC_cons_chr{chrom}_seed{seed}_toimp.vcf.gz ')
    out = os.path.join(path_output, f'SC_cons_chr{chrom}_seed{seed}_imputed')
    positions = os.path.join(path_positions, f'chromosome_position_{chrom}_difference.txt ')
    print(f'Imputation of single cell GM12878 chromosome:{chrom} and seed:{seed} is started, take a little break...')
    command = "java -jar beagle.22Jul22.46e.jar " + f"excludemarkers={positions}" + f"ref={ref}" + f"gt={gt}" + f"map={map}" + f"out={out}"
    result = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()




def similarity_recall_after(dictionary_imp, dictionary_gen, seeds, path_output, path_processed, path_result, sc_samples):
    """_summary_

    Args:
        dictionary_imp (_type_): _description_
        dictionary_gen (_type_): _description_
        seeds (_type_): _description_
        path_output (_type_): _description_
        path_processed (_type_): _description_
        path_result (_type_): _description_
        sc_samples (_type_): _description_
    """

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
    """_summary_

    Args:
        path_result (_type_): _description_
        sc_samples (_type_): _description_
    """
    
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



def delta_scores(name_delta, variable1, variable2, path_result):
    """Calculate and save delta scores based on specified variables.

    Args:
        name_delta (str): Name for the delta score column (e.g., jaccard_delta, recall_delta, snp_delta).
        variable1 (str): Name of the first variable to subtract from (e.g., j_score_post, recall_post, tot_snps_pre).
        variable2 (str): Name of the second variable to subtract (e.g., j_score_pre, recall_pre, tot_snps_post).
        path_result (str): Path to the directory where the result will be saved.
    """
    
    table = pd.read_excel(os.path.join(path_result, 'tables/similarity_recall_concatenated.xlsx'))
    table[name_delta] = table[variable1] - table[variable2]
    table.to_excel(os.path.join(path_result, f'tables/similarity_recall_delta.xlsx'), index=False)


    
def violin_plot(table_delta, variable, variable_title, dir_plots, dir_results, y_max=1, y_min=0.2):
    """_summary_

    Args:
        table_delta (_type_): _description_
        variable (_str_): either chromosome or sample_name
        variable_title (_str_): either chromosomes or samples
    """
    
    table_delta = pd.read_excel(os.path.join(dir_results, 'tables/similarity_recall_delta.xlsx'))
    
    fig = sp.make_subplots(rows=2, cols=1, shared_xaxes=True,
                       subplot_titles=['Jaccard Scores', 'Recall Scores'],
                       vertical_spacing=0.07)
    
    # plot 1: Jaccard Scores 
    fig.add_trace(go.Violin(x=table_delta[variable],y=table_delta['j_score_pre'],
                        box_visible=False,line_color='salmon',
                        meanline_visible=True,opacity=1,
                        name='PRE IMPUTATION'),
                        row=1,col=1)
    
    fig.add_trace(go.Violin(x=table_delta[variable],y=table_delta['j_score_post'],
                        box_visible=False,line_color='mediumseagreen',
                        meanline_visible=True,fillcolor='rgba(0, 128, 0, 0.3)',
                        opacity=1,name='POST IMPUTATION'),
                        row=1,col=1)
    
    # plot 2: Recall Scores 
    fig.add_trace(go.Violin(x=table_delta[variable],y=table_delta['recall_pre'],
                        box_visible=False,line_color='salmon',
                        meanline_visible=True,
                        opacity=1,name='PRE IMPUTATION'),
                        row=2,col=1)
    
    fig.add_trace(go.Violin(x=table_delta[variable],y=table_delta['recall_post'],
                        box_visible=False,line_color='mediumseagreen',
                        meanline_visible=True,fillcolor='rgba(0, 128, 0, 0.3)',
                        opacity=1,name='POST IMPUTATION'),
                        row=2,col=1)

    fig.update_layout(height=600,
                      width=900,
                      title=f'Comparison of Jaccard scores and Recall scores pre and post imputation in {variable_title}',
                      legend=dict(orientation="h",
                                  yanchor="bottom",
                                  y=1.04,
                                  xanchor="right",
                                  x=1)
                      )
    
    fig.update_yaxes(range=[y_min, y_max], row=1, col=1)
    fig.update_yaxes(range=[y_min, y_max], row=2, col=1)
    
    fig.write_image(os.path.join(dir_plots,f'violinplot_{variable_title}.pdf'))
    

#def calculate_frequencies(table, percentages, path_result):

    #table = pd.read_excel(os.path.join(path_result, 'tables/jaccard_scores_concatenated.xlsx'))
    # create relative frequencies
    #table['ratio'] = table['tot_snps_post'] / table['tot_snps_pre']
    #table['snps_lost'] = 1 - table['ratio']
    #for perc in percentages:
        #filtered_table = table[table['perc'] == perc]  
        #filtered_table = filtered_table.groupby('chromosome')[['snps_lost']].mean()
        #iltered_table.to_excel(os.path.join(path_result, 'statistics', f'statistics_frequencies_{perc}.xlsx'))


def calculate_statistics(table, group, path_result):
    """ This function creates to calculate same basics statistics like mean, median, max, min and standard deviation

    Args:
        group (str): either sample or chromosome
    """
    table = pd.read_excel(os.path.join(path_result, 'tables/jaccard_scores_concatenated.xlsx'))

    # calculate the mean
    mean_scores = table.groupby(group)[['j_score_pre', 'j_score_post']].mean()
    mean_scores = mean_scores.rename(columns={'j_score_post': 'mean_scores_post', 'j_score_pre': 'mean_scores_pre'})
    # calculate the median
    median_scores = table.groupby(group)[['j_score_pre', 'j_score_post']].median()
    median_scores = median_scores.rename(columns={'j_score_post': 'median_scores_post', 'j_score_pre': 'median_scores_pre'})
    # calculate the maximum
    max_scores = table.groupby(group)[['j_score_pre', 'j_score_post']].max()
    max_scores = max_scores.rename(columns={'j_score_post': 'max_scores_post', 'j_score_pre': 'max_scores_pre'})
    # calculate the mainimum
    min_scores = table.groupby(group)[['j_score_pre', 'j_score_post']].min()
    min_scores = min_scores.rename(columns={'j_score_post': 'min_scores_post', 'j_score_pre': 'min_scores_pre'})
    # calculate the standard deviation
    std_scores = table.groupby(group)[['j_score_pre', 'j_score_post']].agg(lambda x: round(np.std(x, ddof=0), 3))
    std_scores = std_scores.rename(columns={'j_score_post': 'std_j_score_post', 'j_score_pre': 'std_j_score_pre'})
    
    mean_scores.to_excel(os.path.join(path_result, 'statistics', f'mean_{group}.xlsx'))
    median_scores.to_excel(os.path.join(path_result, 'statistics', f'median_{group}.xlsx'))
    max_scores.to_excel(os.path.join(path_result, 'statistics', f'max_{group}.xlsx'))
    min_scores.to_excel(os.path.join(path_result, 'statistics', f'min_{group}.xlsx'))
    std_scores.to_excel(os.path.join(path_result, 'statistics', f'std_{group}.xlsx'))