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
- [ ] Decide on what data points to use for unsupervised learning
- [ ] Prepare data for unsupervised learning
- [ ] Do unsupervised learning
- [ ] Plot tsne of feature embeddings
- [ ] Feature Comparison for subgroups, what explains subgroups
- [ ] Visualistaion for something like mutations in subgroups
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
│   └── tcga_data/
│       └── coadread_tcga_pan_can_atlas_2018/
├── scripts/
├── src/
│   └── test_setup.py
├── .gitignore
├── environment.yml
└── README.md
```