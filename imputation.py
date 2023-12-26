#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Function for running Beagle algorithm.

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
import subprocess

def imputation(chrom, seed, path_reference, path_map, path_data, path_output, path_positions):
    """
    Performs genetic imputation for a specific chromosome and seed using Beagle. This function constructs and executes a command to run Beagle, a genetic imputation tool. 
    It expects specific paths for reference files, genetic map, genetic data, output, and chromosome positions.
    If Beagle fails to run, it's likely due to issues with the provided paths, so they should be verified.

    Args:
        chrom (int/str): The chromosome number for which imputation is to be performed.
        seed (int): The seed value used for randomization in the imputation process.
        path_reference (str): The file path to the reference VCF file.
        path_map (str): The file path to the genetic map.
        path_data (str): The file path to the genetic data.
        path_output (str): The directory path where the output will be stored.
        path_positions (str): The file path to the chromosome positions file.
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

if __name__ == "__main__":
    try:
        logging.info("Running function: imputation")
        for chromosome in [str(chrom) for chrom in list(range(1,23))]:
            for sd in range(1, 3):
                for percentage in ['10']:
                    imputation(chrom=chromosome,
                               seed=sd,
                               path_reference = 'reference/',
                               path_map= 'map/',
                               path_data= 'data/',
                               path_output = 'data/output/',
                               path_positions = 'reference_positions/'
                               )
    except Exception as e:
        logging.error(f"Error in imputation: {e}")