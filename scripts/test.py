''' 
This is a script to test you correctly imported data
'''

import pandas as pd

try:
    pd.read_csv('./data/tcga_data/coadread_tcga_pan_can_atlas_2018/data_cna.txt', delimiter=' ')
except FileNotFoundError:
    print("Data not copied correctly")
except ModuleNotFoundError:
    print("Please use the conda environment")

print("worked correctly")
