# DRAFTS

DNA Regulatory element Analysis by cell-Free Transcription and Sequencing

<p>Code and materials from paper "<b>Multiplex transcriptional characterizations across diverse and hybrid bacterial cell-free expression systems</b>" Yim SS<sup>*</sup>, Johns NI<sup>*</sup>, Park J, Gomes ALC, McBee RM, Richardson M, Ronda C, Chen SP, Garenne D, Noireaux V, Wang HH. JOURNAL TO BE UPDATED (YEAR TO BE UPDATED). <sup>*</sup>denotes equal contribution</p>

<p>The full paper and supplementary information can be accessed <a href="http://wanglab.c2b2.columbia.edu/publications/">here</a>. [LINK TO BE UPDATED]</p>

<p>Raw sequencing data can be found at NCBI SRA under <a href="https://www.ncbi.nlm.nih.gov/bioproject/PRJNA509603">PRJNA509603</a>.</p>

## dependencies
The following must be installed prior to executing the code in this repository. For Python packages, it may be convenient to obtain these through a distribution such as Anaconda. Installation should only take a few minutes.
<UL>
<LI>python 3.6.X, ipython/jupyter
<UL>
<LI>biopython
<LI>pandas
<LI>numpy
<LI>scipy
<LI>matplotlib
<LI>seaborn
</UL>
<LI>bbmerge
</UL>

## 1. processing of raw sequencing data
`01_DRAFTS_process_raw.sh`
- expects nextseq/miseq raw data folder, where each folder has 2 files of R1 and R2 (paired-end reads) and files sequenced from different lanes of flowcell are separated in four different folders labeled with _L00n
- assumes foldernames are Samplename_L001_*, Samplename_L002_*, Samplename_L003_*, or Samplename_L004_*
Samplename here is SampleID for each sample in the sample sheet for illumna sequencing run

run `01_DRAFTS_process_raw.sh` 1) to find and combine raw nextseq data in search_dir, 2) unzip them to the out_dir, then 3) assemble paired-end reads

```
sudo bash 01_DRAFTS_process_raw.sh [search_dir] [out_dir (optional)]
```
after running `01_DRAFTS_process_raw.sh`, group DNA-seq and RNA-seq reads in seperate folders for further analysis

## 2. error filtering and barcode counting
`02_DRAFTS_extract_data.py`
- out_dir should contain a folder named 01_bccounts with 2 empty folders insde named [01_dna_bccounts, 02_rna_bccounts],
- and a folder named 02_log with 10 empty folders inside named [01_bccounts, 02_lowq, 03_missingadapter, 04_badbc, 05_goodbc_badalign, 06_frag, 07_goodbc_perfectalign, 08_goodbc_goodalign, 09_goodbc_perfectalign_bccounts, 10_goodbc_goodalign_bccounts, 11_log_files]

run `02_DRAFTS_extract_data.py` to 1) filter errors in oligo library synthesis or sequencing, 2) extract barcode counts and 3) other info for qc and additional analysis

```
sudo python 02_DRAFTS_extract_data.py [ref_csv] [dna_directory] [rna_directory] [out_dir]
```

## 3. calculation of transcription levels
`03_DRAFTS_compute_tx.py`
- out_dir should contain a folder named 01_tx

run `03_DRAFTS_compute_tx.py` to 1) compute abundances of DNA and RNA barcode counts and 2) transcription levels

```
sudo python 03_DRAFTS_compute_tx.py [ref_csv] [dna_bc_directory] [rna_bc_directory] [out_dir]
```
