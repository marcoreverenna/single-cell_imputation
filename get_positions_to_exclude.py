#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Fucntions for getting the lists of positions to start the imputation and make imputation faster.

__author__ = Marco Reverenna
__copyright__ = Copyright 2022-2023
__version__ = 1.0
__maintainer__ = Marco Reverenna
__email__ = marcoreverenna@gmail.com
__status__ = Dev
"""

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

pd.options.mode.chained_assignment = None


# load the packages
import os
import logging
import pandas as pd
import subprocess


def get_snparray_positions(path_vcf):
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


if __name__ == '__main__':
    processed_dir = 'data/processed/'
    try:
        logging.info("Running function: get_snparray_positions")
        get_snparray_positions(path_vcf= os.path.join(processed_dir, 'gDNA_consensus_dataframe.vcf.gz'))
    except Exception as e:
        logging.error(f"Error in function get_snparray_positions: {e}")
    
    try:
        logging.info("Running function: get_reference_positions")
        get_reference_positions()
    except Exception as e:
        logging.error(f"Error in function get_reference_positions: {e}")
        
    try:
        logging.info("Running function: difference_positions")
        difference_positions()
    except Exception as e:
        logging.error(f"Error in function difference_positions: {e}")