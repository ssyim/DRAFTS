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
run 01_DRAFTS_process_raw.sh
1) to find and combine raw nextseq data in search_dir, 2) unzip them to the out_dir, 3) then assemble paired-end reads


## 2. error filtering and barcode counting

## 3. calculation of transcription levels
