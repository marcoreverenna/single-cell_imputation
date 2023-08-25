# Validation of single-cell using Beagle5.4 algorithm

<p align="left"> 
    <br> Validation of imputation algorithm in singlecell SNP array data
</p>

---

## 📝 Table of Contents

- [About](#about)
- [Getting Started](#getting_started)
- [Prerequisites and installing](#prerequisites_and_installing)
- [Deployment](#deployment)
- [Usage](#usage)
- [Built Using](#built_using)
- [Authors](#authors)
- [Acknowledgments](#acknowledgement)

## 🧐 About <a name = "about"></a>
The objective of this project is to restore unsatisfactory genotypes using an imputation algorithm. Analyzing single-cell SNP array data poses challenges due to the amplification of DNA, which is necessary because of the limited genetic material available. After quality control, unsuitable genotypes are detected and eliminated. Subsequently, they are imputed again to enhance the accuracy of downstream analysis.

## 🏁 Getting Started <a name = "getting_started"></a>
These instructions will enable you to have a copy of the project up and running on your genomedk profile for development, testing and playback of the workflow for all Leuven PGD036 family.

### 🔧 Prerequisites and installing <a name = "prerequisites_and_installing"></a>
This workflow is set up to be executed on genomedk cluster. Therefore, the only prerequisite is to have a genomedk login and to be a member of the meiomap project group on genomedk. Genomedk uses package manager called [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html "Conda") and [Git](https://github.com/git-guides/install-git "Git") to download the entire project. Date are confidential; reference panel and map are in the folder on the server but downloadable from [Beagle5.4](https://faculty.washington.edu/browning/beagle/beagle.html).


1. Install dependencies into isolated environment
```
conda env create --name imp_proj --file imp_env.yaml
```
2. Activate environment
```
source activate imp_proj
```
3. Download data for imputation
```
sh get_data.sh
```
4. Process parquet file
```
python processing_parquet.py
```
5. Execute workflow, i.e.
```
snakemake --snakefile Snakefile --cores 10
```

## 🚀 Deployment <a name = "deployment"></a>
Add notes...

## 🎈 Usage <a name="usage"></a>
Add notes about how to use the system.

Add additional notes about how to deploy this on a live system.
## ⛏️ Built Using <a name = "built_using"></a>
- [GenomeDK](https://genome.au.dk/) - Server Environment
- [Beagle5.4](https://faculty.washington.edu/browning/beagle/beagle.html) - Imputation algorithm
## ✍️ Authors <a name = "authors"></a>
- [@marcoreverenna](https://github.com/marcoreverenna) -
- [@ivanvogel](https://github.com/puko818) -
## 🎉 Acknowledgements <a name = "acknowledgement"></a>
Add notes...
