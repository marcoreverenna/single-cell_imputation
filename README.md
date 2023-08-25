# Validation of single-cell using Beagle5.4 algorithm

<p align="left"> 
    <br> Validation of imputation algorithm in single-cell SNP array data
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
The objective of this project is to assess the performance of the Beagle 5.4 imputation algorithm on single-cell SNP array data. Single-cell SNP array data is less reliable compared to that obtained from bulk samples (gDNA) due to the altered signal received from the array technology, which is often not properly recognized. This results in a higher prevalence of missing values and values that do not match the more reliable ones found in gDNA. Therefore, an initial validation approach of the algorithm involves considering all autosomes, masking different percentages of the SNPs, considering 25 different lists of positions present in the array data, re-imputing them, and comparing the condition before and after imputation using similarity coefficients.

## 🏁 Getting Started <a name = "getting_started"></a>
These instructions provide the means to establish a functional copy of the project on either your local device or your genomedk profile. Given the computationally intensive nature of the imputation algorithm in the validation phase, it is highly recommended to leverage a High-Performance Computing (HPC) system, such as [Genomedk](https://genome.au.dk/), to optimize and expedite the workflow as much as possible.

### 🔧 Prerequisites and installing <a name = "prerequisites_and_installing"></a>
As previously mentioned, it is advisable to execute the script on the genomedk cluster. For this project, a package manager called [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html "Conda") is employed, containing all the necessary packages to derive results solely from the shareable raw data stored in the repository. [Git](https://github.com/git-guides/install-git "Git") is utilized to download the entire project. For further details regarding the reference materials and map files, please refer to the algorithm's page on [Beagle5.4](https://faculty.washington.edu/browning/beagle/beagle.html).

1. Install dependencies into isolated environment
```
conda env create --name imp_proj --file environment.yml
```
2. Request resources
```
srun --job-name=imputation_beagle --account your-account --partition normal --mem-per-cpu 10G --cpus-per-task 10 --time 1:00:00 --pty /bin/bash
```
3. Activate environment
```
source activate imp_proj
```
4. Run main.py
```
python main.py
```
## 🚀 Deployment <a name = "deployment"></a>
Add notes...

## 🎈 Usage <a name="usage"></a>
Add notes about how to use the system.

Add additional notes about how to deploy this on a live system.
## ⛏️ Built Using <a name = "built_using"></a>
- [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html "Conda") - Conda environment
- [Git](https://github.com/git-guides/install-git "Git") - Git system
- [GenomeDK](https://genome.au.dk/) - Server Environment
- [Beagle5.4](https://faculty.washington.edu/browning/beagle/beagle.html) - Imputation algorithm
## ✍️ Authors <a name = "authors"></a>
- [@marcoreverenna](https://github.com/marcoreverenna) -
- [@ivanvogel](https://github.com/puko818) -
## 🎉 Acknowledgements <a name = "acknowledgement"></a>
Add notes...
