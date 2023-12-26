# Imputation on single-cell SNP array data
This repository has designed to perform imputation on SNP-array data to fill gaps in single-cell genomic data. This is particularly crucial in addressing the challenges posed by Whole Genome Amplification (WGA), a common technique used in single-cell genomics that can introduce significant background noise and result in missing genetic information. 


### Table of contents
- [Description](#description) - OAverview of the project's purpose and goals
- [Getting started](#getting-started) - Instructions on how to begin with this project
- [Bioinformatic parameters](#bioinformatic-parameters) - Explanation and details of the bioinformatic parameters used throughout the pipeline
- [Repository structure](#repository-structure) - A layout of the repository's architecture, describing the purpose of each file or directory
- [References](#references) - Tools used in the project
- [Authors](#authors) - List of contributors to the project
- [Acknowledgments](#acknowledgement) - Credits and thanks to those who helped with the project


### Description <a name = "description"></a>
The repository shows some statistics related on SNP-array data in order to understand if imputation is able to
reintegrate the loss information presente in single-cell data. The analysis begins with the creation of a bulk reference considering
five different bulk data, followed by un pre processing dei dati. Subsequently, l'analisi procede con il calcolo di coefficienti di similarita
e recall per comparare le situazione che precede e segue l'imputazione. Finally, vengono compiute delle statistiche descrittive e creati
dei plots per mostrare i risultati ottenuti.

### Getting started <a name = "getting-started"></a>
1. Set up the conda environment so install required libraries (specified in `requirements.txt` file).
2. Use the functions from `get_gdna_consensus.py` to manipulate and analyze genomic DNA (gDNA) data: they perform various operations ranging from data concatenation, filtering, cleaning and analysis to visualization and data transformation.
3. Use the functions from `get_references_map.py` for downloading large genomic data files: they automate the process of downloading, unzipping, and organizing genomic data files into specified directories.
4. Use the functions from `data_processing_pre_imputation.py` for processing, filtering, and analyzing genomic data, particularly focused on single-cell (SC) genomics and consensus genomic DNA (gDNA) data.
5. Use the functions from `get_positions_to_exclude.py` to .
6. Use the functions from `imputation.py` to performs genetic imputation for each chromosome.
7. Use the functions from `data_processing_post_imputation.py` to .
8. Use the functions from `creating_statistics.py` to .
9. Use the functions from `creating_plots.py` to .

### Bioinformatic parameters <a name = "bioinformatic-parameters"></a>

### Repository structure <a name = "repository-structure"></a>
|File |Description|
|:---:|-----------|
|[data/](data/)|This folder must contain another folder called "raw" in which there should be your personal input data included single-cell and bulk VCF files|
|[requeriments.txt](requeriments.txt)|File with names and versions of packages installed in the virtual environment to run the imputation|
|[beagle.22Jul22.46e.jar](beagle.22Jul22.46e.jar)|Beagle imputation tool to perform the imputation|


## References <a name = "references"></a>
- [Beagle5.4](http://faculty.washington.edu/browning/beagle/beagle.html)

## Authors <a name = "authors"></a>
Contact me at [marcor@dtu.dk](https://github.com/marcoreverenna) for more detail or explanations.

## Acknowledgements <a name = "acknowledgement"></a>
I would like to extend my heartfelt gratitude to [KU](https://www.ku.dk/english/) and [CCS](https://ccs.ku.dk)(Center for Chromosome Stability) for providing the essential resources and support that have been fundamental in the development and success of [Eva Hoffmann group](https://icmm.ku.dk/english/research-groups/hoffmann-group/) projects.



