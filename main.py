#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Module for testing Beagle imputation algorithm on single cell SNP array data

__author__ = Marco Reverenna
__copyright__ = Copyright 2022-2023
__version__ = 1.0
__maintainer__ = Marco Reverenna
__email__ = marcoreverenna@gmail.com
__status__ = Dev
"""

# load the packages
import os
import gzip
import pandas as pd
import random as rd
import subprocess
from sklearn.metrics import jaccard_score
import plotly.graph_objects as go
import statistics as stat
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
        vcf_filepath = os.path.join("reference", vcf_filename)

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

     
            
def get_reference_positions():
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
        
                
                
def get_snparray_positions(path_vcf):
    """_summary_

    Args:
        path_vcf (_type_): _description_
    """
    
    output_folder = 'reference_positions/'
    
    input_vcf_file = pd.read_csv(path_vcf, header=30, sep='\t', dtype='object')

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
        output_file = os.path.join(output_folder, f"chromosome_position_{chromosome}_snparray.txt")
        with open(output_file, "w") as out_file:
            for position in positions:
                out_file.write(position + "\n")



def difference_positions():
    for chromosome in chromosomes_int:
        file_reference = f'reference_positions/chromosome_position_{chromosome}.txt'
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



def concatenate_vcf(file1, file2, file3, file4, file5, kind):
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
    vcf_1 = pd.read_csv(os.path.join('data/raw', f'GM12878_{file1}.vcf.gz'), header=30, sep='\t', dtype='object')
    vcf_2 = pd.read_csv(os.path.join('data/raw', f'GM12878_{file2}.vcf.gz'), header=30, sep='\t', dtype='object')
    vcf_3 = pd.read_csv(os.path.join('data/raw', f'GM12878_{file3}.vcf.gz'), header=30, sep='\t', dtype='object')
    vcf_4 = pd.read_csv(os.path.join('data/raw', f'GM12878_{file4}.vcf.gz'), header=30, sep='\t', dtype='object')
    vcf_5 = pd.read_csv(os.path.join('data/raw', f'GM12878_{file5}.vcf.gz'), header=30, sep='\t', dtype='object')

    # merge vcf files
    file_merged = pd.concat([vcf_1, vcf_2.iloc[:,-1], vcf_3.iloc[:,-1], vcf_4.iloc[:,-1], vcf_5.iloc[:,-1]], axis=1)

    # save to processed directory
    file_merged.to_csv(os.path.join(data_dir, f'processed/GM12878_merged_{kind}.vcf.gz'), sep='\t', index=False)



def consensus_gDNA_dataframe():
    df = pd.read_csv(f'../code_paper_validation/data/processed/GM12878_merged_gDNA.vcf.gz', header=0, sep='\t', dtype='object')
    df.loc[:, df.columns[9:]] = df.loc[:, df.columns[9:]].apply(lambda x: x.str.split(':').str.get(0))
    
    chromosomes_int = list(range(1, 23))
    chromosomes_str = [str(chrom) for chrom in chromosomes_int]
    new_df = df[df['#CHROM'].isin(chromosomes_str)]
    
    dictionary = {'0/0': 0, '0/1': 1, '1/1': 2, '1/2': 3, '2/2': 4, './.': 5}
    columns_to_map = new_df.columns[9:]  
    new_df[columns_to_map] = new_df[columns_to_map].applymap(dictionary.get)
    
    filtered_df = new_df[new_df.iloc[:, 9:].apply(lambda row: row.nunique() == 1, axis=1)]
    pos_list_proxy_gDNA = filtered_df['POS'].tolist()
    
    df_consensus = pd.read_csv(f'../code_paper_validation/data/processed/GM12878_merged_gDNA.vcf.gz', header=0, sep='\t', dtype='object')
    df_consensus_autosomes = df_consensus[df_consensus['#CHROM'].isin(chromosomes_str)]
    df_consensus_filtered = df_consensus_autosomes[df_consensus_autosomes['POS'].isin(pos_list_proxy_gDNA)]
    
    df_consensus_filtered.to_csv(f'../code_paper_validation/data/processed/GM12878_consensus_gDNA.vcf.gz', sep='\t', index=False)
    


def splitting_chromosomes(kind, chrom, directory):
    """ This function split the dataframe for each chromosome and type of data (gDNA or SC)

    Args:
        kind (str): either consensus_gDNA or merged_SC
        chrom (str): list of chromosomes in string format
        directory (str): path
    """
    # load the dataset
    # vcf = pd.read_csv(os.path.join(data_dir, f'processed/GM12878_merged_{kind}.vcf.gz'), header=0, sep='\t', dtype='object')
    vcf = pd.read_csv(os.path.join(data_dir, f'processed/GM12878_{kind}.vcf.gz'), header=0, sep='\t', dtype='object')


    # split for each chromosome
    vcf_chr = vcf[vcf['#CHROM'] == chrom]

    # rename the combinations and set the right format
    vcf_chr.columns = column_names
    vcf_chr['FORMAT'] = 'GT'

    # keep only the genotypes
    #vcf_chr[:,9:] = vcf_chr[:,9:].apply(lambda x : x.str.split(':').str.get(0))
    vcf_chr.loc[:, vcf_chr.columns[9:]] = vcf_chr.loc[:, vcf_chr.columns[9:]].apply(lambda x: x.str.split(':').str.get(0))

    # save splitted chromosome for each dataset
    # vcf_chr.to_csv(os.path.join(directory, f'GM12878_merged_{kind}_chr{chrom}.vcf.gz'), sep='\t', index=False)
    vcf_chr.to_csv(os.path.join(directory, f'GM12878_{kind}_chr{chrom}.vcf.gz'), sep='\t', index=False)



def dataset_arrangement(chrom, seed, perc):
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
    pos_masked = os.path.join(results_dir, f'lists/positions_chr{chrom}_seed{seed}_perc{perc}.txt')
    with open(pos_masked,'w') as output:
        output.write(str(pos_seed))

    vcf_SC_toimp = vcf_SC[~vcf_SC['POS'].isin(pos_seed)]
    vcf_gDNA_ctrl = vcf_gDNA[vcf_gDNA['POS'].isin(pos_seed)]
    vcf_SC_ctrl = vcf_SC[vcf_SC['POS'].isin(pos_seed)]

    # little check, shapes must be of the same lenght
    print(vcf_SC_ctrl.shape[0], vcf_gDNA_ctrl.shape[0])

    # save the new dataset to impute and test before imputation
    vcf_gDNA_ctrl = vcf_gDNA_ctrl.to_csv(os.path.join(data_dir, f'processed/gDNA_control/GM12878_gDNA_merged_chr{chrom}_seed{seed}_perc{perc}_ctrl.vcf.gz'), sep='\t', index=False)
    vcf_SC_ctrl = vcf_SC_ctrl.to_csv(os.path.join(data_dir, f'processed/SC_control/GM12878_SC_merged_chr{chrom}_seed{seed}_perc{perc}_ctrl.vcf.gz'), sep='\t', index=False)
    vcf_SC_toimp = vcf_SC_toimp.to_csv(os.path.join(data_dir, f'processed/SC_toimp/GM12878_SC_merged_chr{chrom}_seed{seed}_perc{perc}_toimp.vcf.gz'), sep='\t', index=False)



def similarity_before(dictionary):
    """This function calculate the similarity before the imputation

    Args:
        dictionary (dict): dictionary to map genotypes - should be like 1/1:1, 1/0:0...
    """

    # create an empty list to fill whit all jaccard scores
    j_list_pre = []
    
    for match in pairings:    
        dict_pre = []

        for chrom in chromosomes_str:
            for seed in seeds:
                for perc in percentages:
                    
                    # load your dataframes splitted for each chromosome, seed and percentage
                    gDNA_merg_chr_ctrl = pd.read_csv(os.path.join(processed_dir, f'gDNA_control/GM12878_gDNA_merged_chr{chrom}_seed{seed}_perc{perc}_ctrl.vcf.gz'), header=0, sep='\t', dtype='object')
                    SC_merg_chr_ctrl = pd.read_csv(os.path.join(processed_dir, f'SC_control/GM12878_SC_merged_chr{chrom}_seed{seed}_perc{perc}_ctrl.vcf.gz'), header=0, sep='\t', dtype='object')
                    
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
        pd.DataFrame(dict_pre).to_excel(os.path.join(results_dir, f'tables/jaccard_scores_{match}_before.xlsx'), index=False)
     
        

def imputation(chrom, seed, perc):
    """This function run beagle imputation algorithm

    Args:
        chrom (str): list of chromosomes to impute in string format 
        seed (int): a list of integers
        perc (str): it is recoomended to use 10 and 20, or a number not greater than 20
    """
    # define all the paths
    # if Beagle does not work probably there is a problem with the paths, please check them out
    ref = os.path.join(reference_dir, f'chr{chrom}.1kg.phase3.v5a.vcf.gz ')
    map = os.path.join(map_dir, f'plink.chr{chrom}.GRCh37.map ')
    gt = os.path.join(data_dir, f'processed/SC_toimp/GM12878_SC_merged_chr{chrom}_seed{seed}_perc{perc}_toimp.vcf.gz ')
    out = os.path.join(output_dir, f'GM12878_SC_merged_chr{chrom}_seed{seed}_perc{perc}imputed')
    positions = f'../code_paper_validation/reference_positions/file_diff_{chrom}.txt '
    
    # print a statement
    print(f'Imputation of single cell GM12878 chromosome:{chrom}, seed:{seed} and percentage:{perc}% is started, take a little break...')
    command = "java -jar ../code_paper_validation/beagle.22Jul22.46e.jar " + f"excludemarkers={positions}" + f"ref={ref}" + f"gt={gt}" + f"map={map}" + f"out={out}"
    result = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()



def similarity_after(dictionary_imp, dictionary_gen):
    """This function calculates the jaccard score after filtering the dataframes.

    Args:
        dictionary_imp (dict): it contains genotypes associated with integers and separeted by '|'
        dictionary_gen (dict): it contains genotypes associated with integers and separeted by '/'
    """
    j_list_post = []
    
    for match in pairings:
        dict_post = []
        print(f'Analysing the pair: {match}')
        for chrom in chromosomes_int:
            for seed in seeds:
                for perc in percentages:
                    # load your dataset
                    SC_merg_chr_imputed = pd.read_csv(os.path.join(output_dir, f'GM12878_SC_merged_chr{chrom}_seed{seed}_perc{perc}imputed.vcf.gz'), header=8, sep='\t', dtype='object')
                    gDNA_merg_chr_ctrl = pd.read_csv(os.path.join(processed_dir, f'gDNA_control/GM12878_gDNA_merged_chr{chrom}_seed{seed}_perc{perc}_ctrl.vcf.gz'), header=0, sep='\t', dtype='object')
                    
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
        pd.DataFrame(dict_post).to_excel(os.path.join(results_dir, f'tables/jaccard_scores_{match}_after.xlsx'), index=False)



def concat_tables():
    """ This function concatenate tables of jaccard scores before and after imputation

    Args:
        results_dir (str): result path to directory
        pairings (str): names of couple created to match gDNA and SC
    """
    # tables before imputation
    file_names_pre = [os.path.join(results_dir,f'tables/jaccard_scores_{match}_before.xlsx') for match in pairings]
    dfs_pre = [pd.read_excel(file) for file in file_names_pre]
    result_pre = pd.concat(dfs_pre)
    
    # tables after imputation
    file_names_post = [os.path.join(results_dir,f'tables/jaccard_scores_{match}_after.xlsx') for match in pairings]
    dfs_post = [pd.read_excel(file) for file in file_names_post] 
    result_post = pd.concat(dfs_post)
    
    # prepare a concatenate table
    result_merged = result_post.copy()
    result_merged['j_score_pre'] = result_pre['j_score']
    result_merged['tot_snps_pre'] = result_pre['tot_snps']
    result_merged.rename(columns={'j_score': 'j_score_post', 'tot_snps':'tot_snps_post'}, inplace=True)
    result_merged = result_merged[['sample','chromosome','seed','perc','j_score_pre','j_score_post', 'tot_snps_pre', 'tot_snps_post']]
    result_merged.to_excel(os.path.join(results_dir, f'tables/jaccard_scores_concatenated.xlsx'), index=False)
  
  
    
def create_violin_plot(results_dir, plots_dir, percentage):
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
    

#def calculate_jaccard_minmax():
#    """_summary_
#    """
#
#   table = pd.read_excel(os.path.join(results_dir, 'tables/jaccard_scores_concatenated.xlsx'))
#    mean_j_score_post = table.groupby('chromosome')['j_score_post'].mean().reset_index()
#    mean_j_score_pre = table.groupby('chromosome')['j_score_pre'].mean().reset_index()

#    result_table = pd.merge(mean_j_score_post, mean_j_score_pre, on='chromosome')

#    result_table = result_table.rename(columns={'j_score_post': 'mean_j_score_post', 'j_score_pre': 'mean_j_score_pre'})

#    result_table['difference'] = result_table['mean_j_score_post'] - result_table['mean_j_score_pre']

     #maximum and minimum difference values and corresponding chromosomes
#    max_difference = result_table['difference'].max()
#    max_difference_index = result_table['difference'].idxmax()
#    max_chromosome = result_table.loc[max_difference_index, 'chromosome']

#    min_difference = result_table['difference'].min()
#    min_difference_index = result_table['difference'].idxmin()
#    min_chromosome = result_table.loc[min_difference_index, 'chromosome']

#    result_table = result_table.apply(lambda x: round(x, 3) if x.name != 'chromosome' else x)

#    result_table.to_excel(os.path.join(results_dir, 'statistics', 'table_similarity_scores.xlsx'), index=False)

     # save the print statements to a text file
#    with open(os.path.join(results_dir, 'statistics', 'statistics_dispersion.txt'), 'w') as f:
#        f.write(f'Considering all the pairings, seeds, and percentages, the greater increase in Jaccard score is {max_difference:.3f} in chromosome {max_chromosome}\n')
#        f.write(f'Considering all the pairings, seeds, and percentages, the lower increase in Jaccard score is {min_difference:.3f} in chromosome {min_chromosome}\n')


#def calculate_std_table():
#    """_summary_
#    """
#    table = pd.read_excel(os.path.join(results_dir, 'tables/jaccard_scores_concatenated.xlsx'))
#    std_by_chromosome_10 = table[table['perc'] == 10].groupby('chromosome')[['j_score_pre', 'j_score_post']].agg(lambda x: round(stat.stdev(x), 3))
#    std_by_chromosome_20 = table[table['perc'] == 20].groupby('chromosome')[['j_score_pre', 'j_score_post']].agg(lambda x: round(stat.stdev(x), 3))

#    std_by_chromosome_10 = std_by_chromosome_10.rename(columns={'j_score_pre': 'dev_standard_pre_10', 'j_score_post': 'dev_standard_post_10'})
#    std_by_chromosome_20 = std_by_chromosome_20.rename(columns={'j_score_pre': 'dev_standard_pre_20', 'j_score_post': 'dev_standard_post_20'})

#    chromosomes = table['chromosome'].unique()
#    all_chromosomes = pd.DataFrame({'chromosome': chromosomes})

#    result_table = pd.merge(all_chromosomes, std_by_chromosome_10, on='chromosome', how='left')
#    result_table = pd.merge(result_table, std_by_chromosome_20, on='chromosome', how='left')

#    result_table = result_table.fillna(0)

#    result_table['dev_standard_pre_10'] = result_table['dev_standard_pre_10'].apply(lambda x: f'{x:.3f}')
#    result_table['dev_standard_post_10'] = result_table['dev_standard_post_10'].apply(lambda x: f'{x:.3f}')
#    result_table['dev_standard_pre_20'] = result_table['dev_standard_pre_20'].apply(lambda x: f'{x:.3f}')
#    result_table['dev_standard_post_20'] = result_table['dev_standard_post_20'].apply(lambda x: f'{x:.3f}')
#    result_table.to_excel(os.path.join(results_dir, 'statistics', 'stdev_all_pairings.xlsx'), index=False)
"""
def calculate_frequencies():
    
    table = pd.read_excel(os.path.join(results_dir, 'tables/jaccard_scores_concatenated.xlsx'))
    table['ratio'] = table['tot_snps_post'] / table['tot_snps_pre']
    
    for perc in [int(perc) for perc in percentages]:
        
        max_ratio_10 = table[table['perc'] == 10]['ratio'].max()
        max_ratio_20 = table[table['perc'] == 20]['ratio'].max()
        min_ratio_10 = table[table['perc'] == 10]['ratio'].min()
        min_ratio_20 = table[table['perc'] == 20]['ratio'].min()

    grouped_table = table.groupby(['chromosome', 'perc'])['ratio'].agg(['max', 'min'])
    grouped_table = grouped_table.rename(columns={'max': 'max_ratio', 'min': 'min_ratio'})
    grouped_table = grouped_table.reset_index()
    
    pivoted_table = grouped_table.pivot(index='chromosome', columns='perc', values=['max_ratio', 'min_ratio'])
    pivoted_table.columns = ['max_ratio_10', 'min_ratio_10', 'max_ratio_20', 'min_ratio_20']
    
    pivoted_table = pivoted_table.apply(lambda x: round(x, 3) if x.name != 'chromosome' else x)
    
    pivoted_table.to_excel(os.path.join(results_dir, 'statistics', 'table_relative_frequencies.xlsx'), index=False)

    with open(os.path.join(results_dir, 'statistics', 'statistics_frequencies.txt'), 'w') as f:
        f.write(f'Maximum value for "ratio" for "perc" 10: {max_ratio_10:.3f}\n')
        f.write(f'Maximum value for "ratio" for "perc" 20: {max_ratio_20:.3f}\n')
        f.write(f'Minimum value for "ratio" for "perc" 10: {min_ratio_10:.3f}\n')
        f.write(f'Minimum value for "ratio" for "perc" 20: {min_ratio_20:.3f}\n')


def calculate_frequencies_2(percentages):

    table = pd.read_excel(os.path.join(results_dir, 'tables/jaccard_scores_concatenated.xlsx'))
    table['ratio'] = table['tot_snps_post'] / table['tot_snps_pre']
    
    results = {}  # un dizionario per memorizzare i risultati per ogni percentuale
    
    for perc in percentages:
        max_ratio = table[table['perc'] == perc]['ratio'].max()
        min_ratio = table[table['perc'] == perc]['ratio'].min()
        
        results[f'Relative frequency maximum value for perc {perc}'] = max_ratio
        results[f'Relative frequency minimum value for perc {perc}'] = min_ratio

    grouped_table = table.groupby(['chromosome', 'perc'])['ratio'].agg(['max', 'min'])
    grouped_table = grouped_table.rename(columns={'max': 'max_ratio', 'min': 'min_ratio'})
    grouped_table = grouped_table.reset_index()
    
    pivoted_table = grouped_table.pivot(index='chromosome', columns='perc', values=['max_ratio', 'min_ratio'])
    pivoted_table.columns = ['max_ratio_10', 'min_ratio_10', 'max_ratio_20', 'min_ratio_20']
    
    pivoted_table = pivoted_table.apply(lambda x: round(x, 3) if x.name != 'chromosome' else x)
    
    pivoted_table.to_excel(os.path.join(results_dir, 'statistics', 'table_relative_frequencies.xlsx'), index=False)

    #with open(os.path.join(results_dir, 'statistics', 'statistics_frequencies.txt'), 'w') as f:
    #    for key, value in results.items():
    #        f.write(f'{key}: {value:.3f}\n')
"""

def calculate_frequencies(percentages):
    """_summary_

    Args:
        percentages (_type_): _description_
    """

    table = pd.read_excel(os.path.join(results_dir, 'tables/jaccard_scores_concatenated.xlsx'))
    # create relative frequencies
    table['ratio'] = table['tot_snps_post'] / table['tot_snps_pre']
    table['snps_lost'] = 1 - table['ratio']
    
    for perc in percentages:
        filtered_table = table[table['perc'] == perc]  
        filtered_table = filtered_table.groupby('chromosome')[['snps_lost']].mean()
        filtered_table.to_excel(os.path.join(results_dir, 'statistics', f'statistics_frequencies_{perc}.xlsx'))



def calculate_statistics(group):
    """ This function creates to calculate same basics statistics like mean, median, max, min and standard deviation

    Args:
        group (str): either sample or chromosome
    """
    table = pd.read_excel(os.path.join(results_dir, 'tables/jaccard_scores_concatenated.xlsx'))
    
    for perc in [int(perc) for perc in percentages]:
        
        filtered_table = table[table['perc'] == perc]
        
        mean_scores = filtered_table.groupby(group)[['j_score_pre', 'j_score_post']].mean()
        mean_scores = mean_scores.rename(columns={'j_score_post': 'mean_scores_post', 'j_score_pre': 'mean_scores_pre'})
        
        median_scores = filtered_table.groupby(group)[['j_score_pre', 'j_score_post']].median()
        median_scores = median_scores.rename(columns={'j_score_post': 'median_scores_post', 'j_score_pre': 'median_scores_pre'})

        max_scores = filtered_table.groupby(group)[['j_score_pre', 'j_score_post']].max()
        max_scores = max_scores.rename(columns={'j_score_post': 'mex_scores_post', 'j_score_pre': 'max_scores_pre'})
        
        min_scores = filtered_table.groupby(group)[['j_score_pre', 'j_score_post']].min()
        min_scores = min_scores.rename(columns={'j_score_post': 'min_scores_post', 'j_score_pre': 'min_scores_pre'})
        
        std_scores = filtered_table.groupby(group)[['j_score_pre', 'j_score_post']].agg(lambda x: round(np.std(x, ddof=0), 3))
        std_scores = std_scores.rename(columns={'j_score_post': 'std_j_score_post', 'j_score_pre': 'std_j_score_pre'})
    
        mean_scores.to_excel(os.path.join(results_dir, 'statistics', f'mean_{group}_perc{perc}.xlsx'))#index=False
        median_scores.to_excel(os.path.join(results_dir, 'statistics', f'median_{group}_perc{perc}.xlsx'))
        max_scores.to_excel(os.path.join(results_dir, 'statistics', f'max_{group}_perc{perc}.xlsx'))
        min_scores.to_excel(os.path.join(results_dir, 'statistics', f'min_{group}_perc{perc}.xlsx'))
        std_scores.to_excel(os.path.join(results_dir, 'statistics', f'std_{group}_perc{perc}.xlsx'))



if __name__=='__main__':
    # avoid warning
    pd.options.mode.chained_assignment = None
    
    # define the paths
    data_dir = '../code_paper_validation/data/'
    map_dir = '../code_paper_validation/map/'
    plots_dir = '../code_paper_validation/plots/'
    results_dir = '../code_paper_validation/results/'
    output_dir = '../code_paper_validation/output/'
    reference_dir = '../code_paper_validation/reference/'
    ref_pos_dir = '../code_paper_validation/reference_positions/'
    processed_dir = '../code_paper_validation/data/processed/'
    
    # define global variables
    chromosomes_int = list(range(1,23))
    chromosomes_str = [str(chrom) for chrom in chromosomes_int]
    column_names = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'pairing_1', 'pairing_2', 'pairing_3', 'pairing_4', 'pairing_5']
    pairings = ['pairing_1', 'pairing_2', 'pairing_3', 'pairing_4', 'pairing_5']
    seeds = range(1, 3)
    percentages = ['10','20']
            
    # create the output and reference_positions directories
    new_folders = ['reference_positions', 'plots', 'reference','map']
    for directory in new_folders:
        if not os.path.exists(directory):
            os.makedirs(directory)
    
    # create sub-folder inside 
    directories_results = ['lists', 'tables', 'statistics']
    for directory in directories_results:
        dir_path = os.path.join(results_dir, directory)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
    
    # create the folders inside data directory
    # create the directories to store new processed data and final outcomes
    directories_data = ['processed', 'output']
    for directory in directories_data:
        dir_path = os.path.join(data_dir, directory)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
       
    # create other step-folders inside processed-folder
    directories_processed = ['gDNA_chroms', 'SC_chroms', 'gDNA_control', 'SC_toimp', 'SC_control']
    for directory in directories_processed:
        dir_path = os.path.join('data/processed/', directory)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
        
    
    # run the functions
    """
    Download material for imputation
    """
    download_unzip_map("https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip", "map")
    download_chromosomes_vcf()
    
    for chrom in range(1,23):
        # enter the url
        url = f"https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr{chrom}.1kg.phase3.v5a.vcf.gz"
        # wget to download the files
        command = ["wget", "-c", "--timeout=120", "--waitretry=60", "--tries=10000", "--retry-connrefused", "-P", "reference", url]
        subprocess.run(command)

    
    """
    Dataframes arrangements
    """
    concatenate_vcf('gDNA_1', 'gDNA_2', 'gDNA_3', 'gDNA_4', 'gDNA_5', 'gDNA')
    concatenate_vcf('SC_129', 'SC_130', 'SC_131', 'SC_132', 'SC_133', 'SC')
    
    consensus_gDNA_dataframe()

    for chromosome in chromosomes_str:
        splitting_chromosomes(kind='consensus_gDNA', chrom = chromosome, directory= 'data/processed/gDNA_chroms')
        splitting_chromosomes(kind='merged_SC', chrom = chromosome, directory= 'data/processed/SC_chroms')

    for chromosome in chromosomes_str:
        for percentage in percentages:
            for sd in seeds:
                dataset_arrangement(chrom=chromosome, seed=sd, perc=percentage)
       
    """
    Similarity scores and imputation
    """           
    similarity_before(dictionary={'0/0': 0, '0/1': 1, '1/0': 1, '1/1': 2, '1/2': 3, '2/2': 4, './.': 5})
    
    for chromosome in chromosomes_str:
        for sd in seeds:
            for percentage in percentages:
                imputation(chrom=chromosome, seed=sd, perc=percentage)        
    
    similarity_after(dictionary_imp={'0|0': 0, '0|1': 1, '1|0': 1, '1|1': 2, '1|2': 3, '2|1': 3, '2|2': 4,'2|0': 5, '0|2': 5},
                     dictionary_gen={'0/0': 0, '0/1': 1, '1/0': 1, '1/1': 2, '1/2': 3, '2/1': 3, '2/2': 4,'2/0': 5, '0/2': 5})
    
    """
    Results
    """
    concat_tables()
    
    create_violin_plot(results_dir, plots_dir, 20)
    create_violin_plot(results_dir, plots_dir, 10)
    
    # calculate statistics
    calculate_statistics(group='chromosome')
    calculate_statistics(group='sample')
    calculate_frequencies(percentages=[10, 20])
    
    """
    Preparing positions files
    """
    #get_reference_positions()
    #get_snparray_positions(path_vcf='data/raw/GM12878_gDNA_1.vcf.gz')
    #difference_positions()
       
    #readme_content = """
    # Reference position informations.
    #In this folder are present all the files created considering the SNP positions that we exclude from each chromosome.
    #This way the imputation analysis save several time considering only SNPs of the SNP array and so we can avoid to filter positions subsequently.
    #"""
    #folder_path = 'reference_positions' 
    #readme_file_path = os.path.join(folder_path, 'README.txt')
    #with open(readme_file_path, 'w') as readme_file:
    #    readme_file.write(readme_content)