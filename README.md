
![GCMP Logo](GCMP_Logo_FJP_V2_JZ_rightward_arrows_r2_FJP-01.png)
[![NSF-1942647](https://img.shields.io/badge/NSF-1942647-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1942647).
# GCMP Global Disease
This project hold the workflow for the Global Coral Microbiome Project's global comparison of coral disease susceptibility and microbiome structure.

## *Project Organization*:
 
 
### Analysis Folders
The analysis folder will be organized into subdirectories by analytical step. Each analytical step is required to have 3 subfolders: input output and procedure


- analysis/
    - organelle_removal/
        - input/
        - output/
        - procedure/
        
 
**input folder:** Typically, necessary data files cannot be checked in to GitHub. Therefore, scripts at each analytical step should download the files they need as a first step (e.g. with wget or curl). 

**output folder:** The output folder holds all results of the analysis. *All files in output should be able to be regenerated (deleted and then recreated) from the data in input using the scripts in procedure*. This is essential to ensure outputs stay up-to-date if upstream QC steps change (which is common).

**procedure:**. The procedure should be written using Jupyter Notebooks when possible, but R markdown or bash scripts may be more convenient for some steps. If there is more than one procedure file, please create an index file that explains what order to run them in and what each is doing.


#### Workflow

Each step is labelled with **output** that is used in later steps, and **products** that are used in Figures or Supp. Figures in the paper
To recreate the analysis we will run the procedures in the following analysis folders in order:

1. coral_disease_data
      **input**: raw data from the FRRP, HICORDIS and Lamb et al, disease datasets; the Huang & Roy coral phylogeny
      **output**: merged_disease_table.tsv
      **products**: disease prevalence graphs
2. metadata
      **input**:  GCMP_EMP_map_r28_no_empty_samples.txt
      **output**: one_hot_encoding_metadata.tsv
3. organelle_removal
      **input**: 
         GCMP_EMP_map_r28_no_empty_samples.txt -- the sample metadata
         all.seqs.fa -- the sequences from QIITA
         all.biom -- the feature table from QIITA
      **output**:
         silva_metaxa2_reference_taxonomy.qza -- taxonomic annotations that better detect mitochondrial sequences
         effects_of_rarefaction/
             feature_table_silva_metaxa_2_all.qza 
         
7. phylogeny_insertion
8. decontamination
9. PICRUSt2
10. core_analysis
11. coral_disease_vs_adiv 
12. GDM
