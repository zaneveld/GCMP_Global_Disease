
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


