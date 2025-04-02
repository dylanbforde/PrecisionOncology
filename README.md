# PrecisionOncology

## Project Setup

This repository contains tools for precision oncology data analysis focused on colorectal cancer datasets from TCGA.

### Prerequisites

- Git
- Conda package manager

### Installation

1. **Clone the repository**
   ```
   git clone https://github.com/username/PrecisionOncology.git
   cd PrecisionOncology
   ```

2. **Download the required dataset**
   - Go to: https://www.cbioportal.org/study/summary?id=coadread_tcga_pan_can_atlas_2018
   - Download the dataset
   - Save the zip file to the `data` folder
   - Extract the contents into `data/tcga_data/coadread_tcga_pan_can_atlas_2018`

3. **Create the conda environment**
   ```
   conda env create -f environment.yml
   conda activate PrecisionOncology
   ```

4. **Verify installation**
   ```
   python src/test_setup.py
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

### Troubleshooting

If you encounter any issues during setup, please create an issue in the repository with details about the error.