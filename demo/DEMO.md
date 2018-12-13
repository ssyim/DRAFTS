## demo analysis

1. run `01_DRAFTS_process_raw.sh`

```
sudo bash 01_DRAFTS_process_raw.sh path_to/01_demo_raw_dir path_to/02_demo_out_dir/01_parsed
```

2. move processed reads of DNA and RNA into 02_dna and 02_rna in 02_demo_out_dir

3. run `02_DRAFTS_extract_data.py`

```
sudo python 02_DRAFTS_extract_data.py path_to/demo_lib_ref.csv path_to/02_demo_out_dir/02_dna path_to/02_demo_out_dir/02_rna path_to/02_demo_out_dir/03_processed
```

4. run `03_DRAFTS_compute_tx.py`

```
sudo python 03_DRAFTS_compute_tx.py path_to/demo_lib_ref.csv path_to/02_demo_out_dir/03_processed/01_bccounts/01_dna_bccounts path_to/demo/02_demo_out_dir/03_processed/01_bccounts/02_rna_bccounts path_to/demo/02_demo_out_dir/04_data
```
