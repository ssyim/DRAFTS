#!/usr/bin/env python

# Code to 1) compute abundances of DNA and RNA barcode counts and 2) transcription levels

import os, sys

import pandas as pd
import numpy as np
from numpy import mean

from natsort import natsorted

### helper functions ###

### user input ###
try:
    ref_csv = str(sys.argv[1]) # csv file of library metadata
    DNA_dir = str(sys.argv[2]) # directory to read processed DNA barcode counts data from
    RNA_dir = str(sys.argv[3]) # directory to read processed RNA barcode counts data from
    out_dir = str(sys.argv[4]) # directory to write to

except IndexError:
    print ("Usage:")
    print ("./03_DRAFTS_compute_tx.py [ref_csv] [dna_bc_directory] [rna_bc_directory] [out_directory]")
    exit(1)

# read library metadata 
lib_key = pd.read_csv(ref_csv, index_col = "Barcode") # barcode sequences in metadata is reverse complemented

# list up DNA-seq data and RNA-seq data
DNA_file_list = []
for subdir, dirs, files in os.walk(DNA_dir):
    for file in files:
        if file.endswith("_bccounts.csv"):
            DNA_file_list.append(os.path.join(subdir, file))
DNA_file_list = natsorted(DNA_file_list)

RNA_file_list = []
for subdir, dirs, files in os.walk(RNA_dir):
    for file in files:
        if file.endswith("_bccounts.csv"):
            RNA_file_list.append(os.path.join(subdir, file))
RNA_file_list = natsorted(RNA_file_list)

