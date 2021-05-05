# HDA TDS Group 4
## About this project 📽️

This project is prepared for the Translational Data Science module as part of MSc Health Data Analytics and Machine Learning. Group 4 has been tasked to utilise UK Biobank (UKB) and implement one of the most cutting-edge techniques Mendelian Randomisation in Computational Epidemiology to assess the causal relationship between Inflammatory Bowel Disease (IBD), its subtypes Crohn's Disease (CD) and Ulcerative Colitis (UC), and Colon Cancer. Our main population comes from the European cohort of Inflammatory Bowel Disease Genetics Consortium (IBDGC), and the validation was performed in the UKB cohort. Hospital Episode Statistics (HES) was provided as part of the UKB data.

For the detailed description of the data sources, data processing and analysis please see the full report.

## Project team 🧑‍🤝‍🧑

* CID1: 01999787
* CID2: 01912853
* CID3: 01906946
* CID4: 01993406
* CID5: 00600463

## Getting started ✈️

### Prerequisites 💻

Below are the list of packages and the version numbers required to reproduce the project.

  ```sh
  npm install npm@latest -g
  ```
### Key files 📂

The following files identifies UKB participants by eid:
| File name           | Description                           |
|---------------------|---------------------------------------|
| ukb26390.csv        | Main UKB data                         |
| ukb27725.csv        | HES data                              |
| GWAS_PCs.rds        | First 10 PCs of the UKB genotype data |
| hesin_diag.txt      | HES medical diagnosis                 |
| w19266_20200204.csv | List of withdrawn participants        |

The following files come from IIBDGC:
| File name                          | Description                                           |
|------------------------------------|-------------------------------------------------------|
| UR.CD.gwas_info03_filtered.assoc   | European IBD 1000 Genome imputed GWAS summarised data |
| EUR.IBD.gwas_info03_filtered.assoc | European CD 1000 Genome imputed GWAS summarised data  |
| EUR.UC.gwas_info03_filtered.assoc  | European UC 1000 Genome imputed GWAS summarised data  |


### Key scripts 📜

The following R scripts should be run sequentially by the respective .sh scripts:
| R script                            | Corresponding .sh Script |
|-------------------------------------|--------------------------|
| 1_SNP_List_Creation.R               |                          |
| 2_Genetic_Data_Creation.R           |                          |
| 3_Outcome_Covariate_Data_Creation.R |                          |
| 4_Final_Data_Creation.R             |                          |
| 5_Analysis_and_Visualisation.R      |                          |
| 5plus_Analysis_and_Visualisation.R  |                          |
| 6.1_Sensitivity_Analysis.R          |                          |
| 6.2_Sensitivity_Analysis.R          |                          |

### Installation 🖱️

1. Clone the repo.
   ```sh
   git clone https://github.com/xinzhao13/HDA_TDS.git
   ```
2. Install dependencies listed above 👆.
3. Locate the UKB and HES data.
4. Download the summarised data from IBD Genetics Consortium from [here](https://www.ibdgenetics.org/downloads.html)
5. Run the pipeline.
    ```sh
    ````


### HPC folder structure 🗄️
```bash
.
├── 1_SNP_List_Creation.R
├── 2_Genetic_Data_Creation.R
├── 3_Outcome_Covariate_Data_Creation.R
├── 4_Final_Data_Creation.R
├── 5_Analysis_and_Visualisation.R
├── 5plus_Analysis_and_Visualisation.R
├── 6.1_Sensitivity_Analysis.R
├── 6.2_Sensitivity_Analysis.R
├── analysis
├── data
├── Job_Submission_DataCreation.sh
├── Job_Submission_Step1Only.sh
├── Job_Submission_Step2Only.sh
├── Job_Submission_Step3Only.sh
├── Job_Submission_Step4Only.sh
├── Job_Submission_Step5Only.sh
├── Job_Submission_Step5plusOnly.sh
├── Job_Submission_Step6Only.sh
├── Job_Submission_Step6plusOnly.sh
├── original
└── README.md
```

## Relevant resources 🏖️

* [UK Biobank](https://www.ukbiobank.ac.uk/)
* [UK Biobank - List of Cancer Cases](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=40006)
* [Understanding Hospital Episode Statistics (HES) data](https://biobank.ndph.ox.ac.uk/showcase/showcase/docs/HospitalEpisodeStatistics.pdf)
* [International Inflammatory Bowel Disease Genetics Consortium (IIBDGC)](https://www.ibdgenetics.org/)
* [MendelianRandomisation: An R Package for Performing Mendelian Randomization Analyses Using Summarized Data](https://cran.r-project.org/web/packages/MendelianRandomization/vignettes/Vignette_MR.pdf)
* [Getting Started with HPC at Imperial](https://www.imperial.ac.uk/admin-services/ict/self-service/research-support/rcs/support/getting-started/)
* [PhenoScannerV2](http://www.phenoscanner.medschl.cam.ac.uk/)