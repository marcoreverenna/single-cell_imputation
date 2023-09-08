#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Testing Beagle imputation algorithm on single cell SNP array data

__author__ = Marco Reverenna
__copyright__ = Copyright 2022-2023
__version__ = 1.0
__maintainer__ = Marco Reverenna
__email__ = marcoreverenna@gmail.com
__status__ = Dev
"""

import my_module as module


# avoid warning
module.pd.options.mode.chained_assignment = None
    
# define the paths
data_dir = '../test1/data/'
map_dir = '../test1/map/'
plots_dir = '../test1/plots/'
results_dir = '../test1/results/'
output_dir = '../test1/data/output/'
reference_dir = '../test1/reference/'
ref_pos_dir = '../test1/reference_positions/'
processed_dir = '../test1/data/processed/'
    
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
    if not module.os.path.exists(directory):
        module.os.makedirs(directory)
    
# create sub-folder inside 
directories_results = ['lists', 'tables', 'statistics']
for directory in directories_results:
    dir_path = module.os.path.join(results_dir,
                                   directory
                                   )
    if not module.os.path.exists(dir_path):
        module.os.makedirs(dir_path)
    
# create the folders inside data directory
# create the directories to store new processed data and final outcomes
directories_data = ['processed', 'output']
for directory in directories_data:
    dir_path = module.os.path.join(data_dir,
                                   directory
                                   )
    if not module.os.path.exists(dir_path):
        module.os.makedirs(dir_path)
       
# create other step-folders inside processed-folder
directories_processed = ['gDNA_chroms', 'SC_chroms', 'gDNA_control', 'SC_toimp', 'SC_control']
for directory in directories_processed:
    dir_path = module.os.path.join('data/processed/',
                                   directory
                                   )
    if not module.os.path.exists(dir_path):
        module.os.makedirs(dir_path)
        
    
# run the functions from module 

module.download_unzip_map("https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip", "map")

module.download_chromosomes_vcf()

module.concatenate_vcf(file1= 'gDNA_1', file2='gDNA_2', file3='gDNA_3', file4='gDNA_4', file5='gDNA_5', kind='gDNA', path_data= data_dir)

module.concatenate_vcf(file1='SC_129', file2='SC_130', file3='SC_131', file4='SC_132', file5='SC_133', kind='SC', path_data= data_dir)
 
module.consensus_gDNA_dataframe(path_processed= processed_dir)

module.filtering_all_SC(path_processed = processed_dir)

for chromosome in chromosomes_str:
    module.splitting_chromosomes(kind='consensus_gDNA',
                                 chrom = chromosome,
                                 path_data = data_dir,
                                 directory= 'data/processed/gDNA_chroms'
                                 )
    module.splitting_chromosomes(kind='merged_SC',
                                 chrom = chromosome,
                                 path_data = data_dir, 
                                 directory= 'data/processed/SC_chroms'
                                 )

for chromosome in chromosomes_str:
    for percentage in percentages:
        for sd in seeds:
            module.dataset_arrangement(chrom=chromosome,
                                       seed=sd,
                                       perc=percentage,
                                       path_result = results_dir,
                                       path_data = data_dir
                                       )
       
module.similarity_before(dictionary={'0/0': 0, '0/1': 1, '1/0': 1, '1/1': 2, '1/2': 3, '2/2': 4, './.': 5},
                         seeds= range(1, 3),
                         percentages= ['10','20'],
                         path_processed= processed_dir,
                         path_result= results_dir
                         )


module.get_snparray_positions(path_vcf='data/processed/GM12878_consensus_gDNA.vcf.gz')

module.get_reference_positions()

module.difference_positions()
    
for chromosome in chromosomes_str:
    for sd in seeds:
        for percentage in percentages:
            module.imputation(chrom=chromosome,
                              seed=sd,
                              perc=percentage,
                              path_reference = reference_dir,
                              path_map= map_dir,
                              path_data= data_dir,
                              path_output = output_dir,
                              path_positions = ref_pos_dir
                              )        
    
module.similarity_after(dictionary_imp= {'0|0': 0, '0|1': 1, '1|0': 1, '1|1': 2, '1|2': 3, '2|1': 3, '2|2': 4,'2|0': 5, '0|2': 5},
                        dictionary_gen= {'0/0': 0, '0/1': 1, '1/0': 1, '1/1': 2, '1/2': 3, '2/1': 3, '2/2': 4,'2/0': 5, '0/2': 5},
                        seeds= range(1,3),
                        percentages= ['10','20'],
                        path_output= output_dir,
                        path_processed= processed_dir,
                        path_result =results_dir
                        )
    

module.concat_tables(path_result=results_dir)
    
module.create_violin_plot(results_dir,
                          plots_dir,
                          20
                          )
module.create_violin_plot(results_dir,
                          plots_dir,
                          10
                          )

module.calculate_statistics(group='chromosome',
                            path_result=results_dir
                            )
module.calculate_statistics(group='sample',
                            path_result=results_dir
                            )
module.calculate_frequencies(percentages=[10, 20],
                             path_result=results_dir
                             )
    
      
#readme_content =
# Reference position informations.
#In this folder are present all the files created considering the SNP positions that we exclude from each chromosome.
#This way the imputation analysis save several time considering only SNPs of the SNP array and so we can avoid to filter positions subsequently.
#"""
#folder_path = 'reference_positions' 
#readme_file_path = os.path.join(folder_path, 'README.txt')
#with open(readme_file_path, 'w') as readme_file:
#    readme_file.write(readme_content)

#for chrom in range(1,23):
#    # enter the url
#    url = f"https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr{chrom}.1kg.phase3.v5a.vcf.gz"
#    # wget to download the files
#    command = ["wget", "-c", "--timeout=120", "--waitretry=60", "--tries=10000", "--retry-connrefused", "-P", "reference", url]
#    module.subprocess.run(command)