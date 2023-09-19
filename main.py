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

sc_samples = ['201666260058_R01C01','201285550024_R06C01',
              '201285550024_R04C01', '201285550024_R02C01',
              '201666260058_R02C01', '201285550024_R06C02',
              '201285550024_R04C02', '201285550024_R02C02',
              '201666260058_R03C01', '201662330032_R06C01',
              '201662330032_R04C01', '201662330032_R02C01',
              '201666260058_R04C01', '201662330032_R06C02',
              '201662330032_R04C02', '201662330032_R02C02',
              '201666260058_R05C01', '201662330098_R06C01',
              '201662330098_R04C01', '201662330098_R02C01']


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
print('Running module: download_unzip_map')
module.download_unzip_map("https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip", "map")

print('Running module: download_chromosomes_vcf')
module.download_chromosomes_vcf()

print('Running module: concatenate_SC')
module.concatenate_SC(files = ["SC_129", "SC_131", "SC_133", "SC_135", "SC_137",
                               "SC_139", "SC_141", "SC_143", "SC_145", "SC_147",
                               "SC_149", "SC_151", "SC_153", "SC_155", "SC_157",
                               "SC_159", "SC_161", "SC_163", "SC_165", "SC_167"],
                               path_data = data_dir,
                               output_filename = "GM12878_SC_merged.vcf.gz"
                               )

print('Running module: filtering_all_SC')
module.filtering_all_SC(path_processed = processed_dir)

print('Running module: splitting_chromosomes SC')
for chromosome in chromosomes_str:
    module.splitting_chromosomes(kind='SC',
                                 chrom = chromosome,
                                 path_processed = processed_dir
                                 )    
                                 
print('Running module: splitting_chromosomes gDNA')    
for chromosome in chromosomes_str:
    module.splitting_chromosomes(kind='gDNA',
                                 chrom = chromosome,
                                 path_processed = processed_dir
                                 )  
  
print('Running module: dataset_arrangement')            
for chromosome in chromosomes_str:
    print(f'chromosome: {chromosome}')
    for percentage in percentages:
        for sd in seeds:
            module.dataset_arrangement(chrom=chromosome,
                                       seed=sd,
                                       path = processed_dir,
                                       path_result = results_dir,
                                       perc='10'
                                       )

print('Running module: similarity_recall_before')
module.similarity_recall_before(dictionary={'0/0': 0, '0/1': 1, '1/0': 1, '1/1': 2, '1/2': 3, '2/2': 4, './.': 5},
                         seeds= range(1, 3),
                         path_processed= processed_dir,
                         path_result= results_dir,
                         sc_samples= sc_samples
                         )


print('Running module: get_snparray_positions')
module.get_snparray_positions(path_vcf='data/processed/GM12878_consensus_gDNA.vcf.gz')

print('Running module: get_reference_positions')
module.get_reference_positions()

print('Running module: difference_positions')
module.difference_positions()

print('Running module: imputation')
for chromosome in chromosomes_str:
    for sd in seeds:
        for percentage in percentages:
            module.imputation(chrom=chromosome,
                              seed=sd,
                              path_reference = reference_dir,
                              path_map= map_dir,
                              path_data= data_dir,
                              path_output = output_dir,
                              path_positions = ref_pos_dir
                              )        
print('Running module: similarity_recall_after')    
module.similarity_recall_after(dictionary_imp= {'0|0': 0, '0|1': 1, '1|0': 1, '1|1': 2, '1|2': 3, '2|1': 3, '2|2': 4,'2|0': 5, '0|2': 5},
                        dictionary_gen= {'0/0': 0, '0/1': 1, '1/0': 1, '1/1': 2, '1/2': 3, '2/1': 3, '2/2': 4,'2/0': 5, '0/2': 5},
                        seeds= range(1, 3),
                        path_output= output_dir,
                        path_processed= processed_dir,
                        path_result =results_dir,
                        sc_samples=sc_samples
                        )
   
print('Running module: concat_tables')    
module.concat_tables(path_result=results_dir, sc_samples = sc_samples)

print('Running module: violin_plot')    
module.violin_plot(table_delta = module.pd.read_excel(module.os.path.join(results_dir, 'tables/similarity_recall_delta.xlsx')),
            variable = 'chromosome',
            variable_title = 'chromosomes',
            y_max=1,
            y_min=0.2,
            dir_results = results_dir,
            dir_plots= plots_dir
            )

module.violin_plot(table_delta = module.pd.read_excel(module.os.path.join(results_dir, 'tables/similarity_recall_delta.xlsx')),
            variable = 'sample_name',
            variable_title = 'samples',
            y_max=1,
            y_min=0.2,
            dir_results = results_dir,
            dir_plots= plots_dir
            )


print('Running module: violin_plot')    
module.calculate_statistics(table = module.pd.read_excel(module.os.path.join(results_dir, 'tables/similarity_recall_delta.xlsx')),
                            group = 'chromosome',
                            path_result = results_dir
                            )

module.calculate_statistics(table = module.pd.read_excel(module.os.path.join(results_dir, 'tables/similarity_recall_delta.xlsx')),
                            group = 'sample_name',
                            path_result = results_dir
                            )
      
      
#module.calculate_frequencies(table = module.pd.read_excel(module.os.path.join(results_dir, 'tables/similarity_recall_delta.xlsx')),
#                             percentages=[10, 20],
#                             path_result=results_dir)
#readme_content =
# Reference position informations.
#In this folder are present all the files created considering the SNP positions that we exclude from each chromosome.
#This way the imputation analysis save several time considering only SNPs of the SNP array and so we can avoid to filter positions subsequently.
"""
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
"""