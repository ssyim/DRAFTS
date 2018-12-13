#!/usr/bin/env python

# Code to 1) compute abundances of DNA and RNA barcode counts and 2) transcription levels

import os, sys
from multiprocessing import Pool as ThreadPool

import pandas as pd
import numpy as np

from natsort import natsorted

### helper functions ###

def filter_lib(file, library):
	read = pd.read_csv(file,index_col="Barcode")
	filter_key = read.loc[library.index] #filtering library barcodes	
	return filter_key

def filter_tx_log10(file, library):
	read = pd.read_csv(file,index_col="Barcode")
	tx_log10 = read[read.columns[-1]]
	filter_key = tx_log10.loc[library.index] #filtering library barcodes	
	return filter_key

### user input ###
try:
    ref_csv = str(sys.argv[1]) # csv file of library metadata
    DNA_dir = str(sys.argv[2]) # directory to read processed DNA barcode counts data from
    RNA_dir = str(sys.argv[3]) # directory to read processed RNA barcode counts data from
    out_dir = str(sys.argv[4]) # directory to write to

    # [important]
    # out_dir should have a folder named 01_tx

except IndexError:
    print ("Usage:")
    print ("./03_DRAFTS_compute_tx.py [ref_csv] [dna_bc_directory] [rna_bc_directory] [out_dir]")
    exit(1)

# read reference barcodes of library
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

### main script body ###

# compute tx for each sample
def compute_tx(file_DNA, file_RNA):

	counts = []
	counts.append(filter_lib(file_DNA, lib_key))
	counts.append(filter_lib(file_RNA, lib_key))
	combined_counts = pd.concat(counts, axis=1)
	
	colnames = []
	DNA_col_name = file_DNA.split('/')[-1].split('_bccounts.csv')[0]
	RNA_col_name = file_RNA.split('/')[-1].split('_bccounts.csv')[0]
	colnames.append(DNA_col_name)
	colnames.append(RNA_col_name)
	combined_counts.columns = colnames

	n_tx = '_'.join([combined_counts.columns[0],combined_counts.columns[1],'tx'])
	n_tx_log10 = '_'.join([combined_counts.columns[0],combined_counts.columns[1],'tx','log10'])

	# compute tx values by normalizing relative abundance of RNA by relative abundance of DNA
	relative_DNA = combined_counts[combined_counts.columns[0]]/combined_counts[combined_counts.columns[0]].sum()
	relative_RNA = combined_counts[combined_counts.columns[1]]/combined_counts[combined_counts.columns[1]].sum()
	combined_counts[n_tx] = relative_RNA/relative_DNA

	# filter based on DNA and RNA counts
	combined_counts[n_tx] = combined_counts[n_tx].mask(np.isnan(combined_counts[combined_counts.columns[1]]), "no_RNA_counts")
	combined_counts[n_tx] = combined_counts[n_tx].mask(np.isnan(combined_counts[combined_counts.columns[0]]), "no_DNA_counts")
	combined_counts[n_tx] = combined_counts[n_tx].mask(combined_counts[combined_counts.columns[0]] < 10, "low_DNA_counts")
	
	# log10 transformation
	combined_counts[n_tx_log10] = np.log10(combined_counts[n_tx].apply(pd.to_numeric, errors='coerce'))

	combined_counts.index.rename('Barcode', inplace=True)
	combined_counts.to_csv(out_dir + "/01_tx/" + n_tx + ".csv")

pool = ThreadPool(36)
pairwise_results = pool.starmap(compute_tx, zip(DNA_file_list,RNA_file_list)) # assuming the DNA-RNA pairs are in the same order in each filelist
pool.close()
pool.join()

# report summary
tx_file_list = []
for subdir, dirs, files in os.walk(out_dir):
    for file in files:
        if file.endswith("_tx.csv"):
            tx_file_list.append(os.path.join(subdir, file))
tx_file_list = natsorted(tx_file_list)

tx_log10 = []
for file in tx_file_list:
    tx_log10.append(filter_tx_log10(file, lib_key))

# combining all data
combined_tx_log10 = pd.concat(tx_log10, axis=1)
combined_tx_log10.to_csv(out_dir + "/01_tx/" + "/summary_log10_tx.csv")

print ("Done")