#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Fucntions for downloading reference VCF and map file for imputation on single-cell SNP array data.

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
import subprocess


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

pd.options.mode.chained_assignment = None


def download_unzip_map(url, output_folder):
    if os.path.exists(output_folder) and os.listdir(output_folder):
        logging.info(f"Folder '{output_folder}' already exists and contains data. Skip the download.")
    else:
        logging.info(f"Start the download from {url}")
        
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



def download_references_vcf():
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


if __name__ == '__main__':
    # run the functions from module
    try:
        download_unzip_map("https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip", "map")
        logging.info("Running module: download_unzip_map")
    except Exception as e:
        logging.error(f"Error in download_unzip_map: {e}")

    try:
        download_references_vcf()
        logging.info("Running module: download_references_vcf")
    except Exception as e:
        logging.error(f"Error in module download_references_vcf{e}:")