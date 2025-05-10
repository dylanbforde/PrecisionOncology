# PrecisionOncology

## Project Setup

This repository contains tools for precision oncology data analysis focused on colorectal cancer datasets from TCGA.

### Prerequisites

- Git
- Conda package manager

# Project Goals
Explore relevant TCGA data from patients that can be used for unsupervised clustering. Next a gene expression analysis for the population, then look at gene expression patterns on the discovered subtypes.
Then link patient survival data and see if there are any subtypes with noticeably different health outcomes.
If there are groups with high/low survival you can look at genes in those subgroups and determine if there are possible treatment options that are not currently considered.
For example if HER2 is expressed in subgroups that did not survive, it could be targeted by an oncologist as an area that is not being treated in the cancer.
If there are genes that do not currently have treatments, and their function is unknown, we can suggest those for further study.

# To-Do
- [x] Decide on what data points to use for unsupervised learning
**Decided to use data_cna, data_mrna_seq_v2_rsem.txt, data_genetic_ancestry.txt, and data_mutations.txt.**
- [x] Prepare data for unsupervised learning
- [x] Do unsupervised learning
- [x] Plot tsne of feature embeddings
- [x] Feature Comparison for subgroups, what explains subgroups
- [x] Visualistaion for something like mutations in subgroups
- [ ] How can we test for this subgroup
- [ ] Suggest treatment for that subgroup

### Installation

1. **Clone the repository**
   ```
   git clone https://github.com/dylanbforde/PrecisionOncology.git
   cd PrecisionOncology
   ```

2. **Download the required dataset**
   - Go to: https://www.cbioportal.org/study/summary?id=coadread_tcga_pan_can_atlas_2018
   - Download the dataset
   - Save the zip file to the `data` folder
   - Extract the contents into `data/tcga_data/`

3. **Create the conda environment**
   ```
   conda env create -f conda_environment.yml
   conda activate PrecisionOncology
   ```

   **For Continued Use**
   ```
   conda activate PrecisionOncology in cmd
   Ctrl + Shift + P -> Python:SelectInterpreter -> The one with Precision Oncology
   ```

4. **Verify installation**
   ```
   python3 src/test_setup.py
   ```
   If the test fails, please contact with the error message.

### Project Structure

```
PrecisionOncology/
├── data/
│   ├── adapted_data/
│   │   ├── ancestry_matrix.csv
│   │   ├── cna_matrix.csv
│   │   ├── merged_data.csv
│   │   ├── mrna_matrix.csv
│   │   └── mutation_matrix.csv
│   └── tcga_data/
│       └── coadread_tcga_pan_can_atlas_2018/
├── scripts/
|   └── run_data_preparation.sh
├── src/
│   ├── test_setup.py
│   └── data_preparation/
│       ├── merge_data.py
│       ├── prepare_all_data.py
│       ├── prepare_ancestry.py
│       ├── prepare_cna.py
│       ├── prepare_mrna.py
│       └── prepare_mutations.py
├── .gitignore
├── environment.yml
└── README.md
```