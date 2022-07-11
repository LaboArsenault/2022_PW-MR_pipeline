# PW-MR pipeline for cis-MR analyses

PW-MR pipeline for cis-MR analyses <br>
Version : 1.0 <br>
Created by : Jerome Bourgault <br>
Contact : jerome.bourgault@criucpq.ulaval.ca <br>

## Download
- cd your/folder/
- git clone ...
- Rscript PWAS.R [options]

## Prerequisites
### R packages :
- data.table
    - TwoSampleMR
    - MRInstruments
    - coloc
    - MRPRESSO
    - MendelianRandomization
    - hyprcoloc
    - ggplot2
    - ggrepel
    - dplyr
    - magrittr
    - TissueEnrich
    - R2BGLiMS
    - ieugwasr
### Software :
    - gcta64 (in PATH)

## Important informations
#### *Make sure all included '.R' and '.sh' scripts are executable and all included folders are writable.*

### To run this script properly (not mandatory, but easier to run with default options) :
#### Data folder should include (in './data/') :
    - proteins information file (should contain at least the follow columns with identical names : id, trait, prot_name, file; All proteins included in the file will be used for the analysis.)
    - genes reference file (should contain the following columns in the same order (Chromosome/scaffold name must have the same name) : Gene stable ID version, Gene name, Chromosome/scaffold name, Gene start (bp), Gene end (bp), Gene stable ID)
    - chromosome lengths file (should contain at least the following columns in the same order with the same names : Chr, Total_length_db)
    - LD reference file for MR clumping (PLINK format)
    - (optional) rsids reference file (should contain the following columns in the same order, NO HEADER : rsid, pos, chr)
#### Cojo folder should include (in './cojo/ldref/') :
    - Your LD reference files in PLINK format, splitted by chromosomes if possible (to reduce computing time)
           LD reference file in ./data/ can also be used for CoJo, if samplesize is sufficient.

## The sequence to run a full analysis goes like this :
- 1.PWAS.R
- 2.results_table.R (only works with results in MRres folder, not with sensitivity analyses results)
- 3.PWAS_sensitivity.R

## References and methods used by this script
- TwoSampleMR : Hemani et al., 2018 (doi: 10.7554/eLife.34408)
- JAM : Newcombe et al., 2016 (doi: 10.1002/gepi.21953)
- PCA for PW-MR : Burgess et al., 2017 (doi: 10.1002/gepi.22077)
- CoJo : Yang et al., 2012 (doi: 10.1038/ng.2213)
- coloc : Wallace, 2021 (doi: 10.1101/2021.02.23.432421)
- MR-PRESSO : Verbanck et al., 2018 (doi: 10.1038/s41588-018-0099-7)
- HyPrColoc : Foley et al., 2020 (doi: 10.5281/zenodo.4293559)
- TissueEnrich : Jain et al., 2018 (doi: 10.1093/bioinformatics/bty890)
