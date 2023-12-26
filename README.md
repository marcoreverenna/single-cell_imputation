# Statistical analysis on single-cell SNP array data and imputation
Scripts to implement imputation on SNP-array data to fill the gaps on single-cell data and recover genetic information due to background noise caused by WGA 
(whole-genome-amplification).


### Table of contents
- [Project](#project) - OAverview of the project's purpose and goals
- [Getting started](#getting-started) - Instructions on how to begin with this project
- [Bioinformatic para	Bmeters](#bioinformatic-parameters) - Explanation and details of the bioinformatic parameters used throughout the pipeline
- [Repository structure](#repository-structure) - A layout of the repository's architecture, describing the purpose of each file or directory
- [References](#references) - Tools used in the project
- [Authors](#authors) - List of contributors to the project
- [Acknowledgments](#acknowledgement) - Credits and thanks to those who helped with the project


### Project <a name = "project"></a>
The repository shows some statistics related on SNP-array data in order to understand if imputation is able to
reintegrate the loss information presente in single-cell data. The analysis begins with the creation of a bulk reference considering
five different bulk data, followed by un pre processing dei dati. Subsequently, l'analisi procede con il calcolo di coefficienti di similarita
e recall per comparare le situazione che precede e segue l'imputazione. Finally, vengono compiute delle statistiche descrittive e creati
dei plots per mostrare i risultati ottenuti.

### Getting started
1. Set up the conda environment so install required libraries (specified in `requirements.txt` file).
2. Use the functions from `get_gdna_consensus.py` to .
3. Use the functions from `get_references_map.py` to .
4. Use the functions from `data_processing_pre_imputation.py` to .
5. Use the functions from `get_positions_to_exclude.py` to .
6. Use the functions from `imputation.py` to .
7. Use the functions from `data_processing_post_imputation.py` to .
8. Use the functions from `creating_statistics.py` to .
9. Use the functions from `creating_plots.py` to .

### Repository structure
|File |Description|
|:---:|-----------|
|[data/](data/)|This folder must contain another folder called "raw" in which there should be your personal input data included single-cell and bulk VCF files|
|[requeriments.txt](requeriments.txt)|File with names and versions of packages installed in the virtual environment to run the imputation|
|[beagle.22Jul22.46e.jar](beagle.22Jul22.46e.jar)|Beagle imputation tool to perform the imputation|


## References <a name = references"></a>
- [Beagle5.4](https://www.beagle.com)

## Authors <a name = "authors"></a
Contact me at marcor@dtu.dk for more detail or explanations.
- [marcor@dtu.dk](https://github.com/marcoreverenna)

## Acknowledgements <a name = "acknowledgement"></a>
I would like to extend my heartfelt gratitude to [KU](https://www.ku.dk/) and CCS (Center for Chromosome Stability) for providing the essential resources
and support that have been fundamental in the development and success of [Eva Hoffmann group](https://www.hoffmann-group-lab.dk) projects.



