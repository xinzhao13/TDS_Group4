# HDA TDS Group 4
## About this project 📽️

This project is prepared for the Translational Data Science module as part of MSc Health Data Analytics and Machine Learning. Group 4 has been tasked to utilise UK Biobank (UKB) and implement one of the most cutting-edge techniques Mendelian Randomisation (MR) in epidemiology to assess the causal relationship between Inflammatory Bowel Disease (IBD), its subtypes Crohn's Disease (CD) and Ulcerative Colitis (UC), and Colon Cancer. 

We adopted the Polygenic Risk Scores (PRS) as instrumental variables to infer the causal relationship between IBD (and UC and CD respectively) and Colon Cancer.

We used the European cohort of Inflammatory Bowel Disease Genetics Consortium (IIBDGC) to select the IBD-associated SNPs and to generate the weights of those SNPS, the causal inference between IBD and colon cancer was performed in the UKB cohort. We used the UKB genotype data, cancer diagnosis, lifestyle factors, medication, dietary habits, biomarkers. Hospital Episode Statistics (HES) was linked to the UKB data by the unique eid identifier of participants, from which we used the non-cancer disease outcomes.

In addition to the observational analysis and the main MR analysis, we have studied potential confounding variables using LASSO, sPLS and Random Forest models, used the selected/important covariates to attenuate the observational associations between IBD and colon cancer. We have queried the selected SNPS and checked for potential confounders using PhenoScanner, performed logistic regressions between the individual IBD-associated SNPs and colon cancer, followed by a set of sensitivity analyses (WM, IVW and MR-Egger etc) to assess the assumptions of the MR analysis.

For the detailed description of the data sources, data processing and analysis, please see the full report.

&nbsp;


## Project team 🧑‍🤝‍🧑

* CID1: 01999787
* CID2: 01912853
* CID3: 01906946
* CID4: 01993406
* CID5: 00600463 

&nbsp;


## Key files 📂

The following files identifies UKB participants by eid:
| File name           | Description                           |
|---------------------|---------------------------------------|
| ukb26390.csv        | Main UKB data                         |
| ukb27725.csv        | UKB biomarker data                    |
| GWAS_PCs.rds        | First 10 PCs of the UKB genotype data |
| hesin_diag.txt      | HES medical diagnosis                 |
| w19266_20200204.csv | List of withdrawn participants        |

&nbsp;

The following files come from IIBDGC:
| File name                          | Description                                           |
|------------------------------------|-------------------------------------------------------|
| UR.CD.gwas_info03_filtered.assoc   | European IBD 1000 Genome imputed GWAS summarised data |
| EUR.IBD.gwas_info03_filtered.assoc | European CD 1000 Genome imputed GWAS summarised data  |
| EUR.UC.gwas_info03_filtered.assoc  | European UC 1000 Genome imputed GWAS summarised data  |

&nbsp;


## Key scripts 📜

The following R scripts will be run sequentially by the respective bash scripts:
| R script                            | Description                                                                                                                                                                                                                                     |
|-------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 1_SNP_List_Creation.R               | Selects the IBD-associated SNPs (and also for UC, CD separately) from IIBDGC, this includes clumping and the removal of rare variants.                                                                                                          |
| 2_Genetic_Data_Creation.R           | Matches effect alleles and major/minor alleles between UKB and IIBDGC, generates genetic data on selected SNPs from UKB genotype data, and calculates the PRS scores.                                                                           |
| 3_Outcome_Covariate_Data_Creation.R | Extracts disease status including colon cancer, other digestive cancer types, non-cancer diseases and covariates from the UKB data.                                                                                                             |
| 4_Final_Data_Creation.R             | Joins all disease status and covariates with PRS scores.                                                                                                                                                                                        |
| 5_Analysis_and_Visualisation.R      | Produces table 1, performs the main MR analysis, performs mean/mode imputation of any remaining missing values, performs sPLS and Random Forest for confounder analysis, visualises the PRS scores and the results of the confounder analysis. |
| 5plus_Analysis_and_Visualisation.R  | Performs LASSO and attenuation analysis, visualises the results.                                                                                                                                                                                |
| 6.1_Sensitivity_Analysis.R          | Prepares for multivariable logistic regressions at SNP level, queries PhennoScanner.                                                                                                                                                            |
| 6.2_Sensitivity_Analysis.R          | Performs the multivariable logistic regressions at SNP level, consolidates the results, visualises the 2-sample MR and sensitivity analysis.                                                                                                    |


&nbsp;

The following bash scripts calls the R scripts above:

| Bash script                     | Description |
|---------------------------------|-------------|
| 
| Job_Submission_Step3Only.sh     |             |
| Job_Submission_Step4Only.sh     |             |
| Job_Submission_Step5Only.sh     |             |
| Job_Submission_Step5plusOnly.sh |             |
| Job_Submission_Step6Only.sh     |             |


&nbsp;

## Getting started ✈️

### Prerequisites 💻

Below are the list of packages and the version numbers required to reproduce the project.

  ```sh
  npm install npm@latest -g
  ```

### Installation 🖱️

1. Clone the repo.
   ```sh
   git clone https://github.com/xinzhao13/TDS_Group4.git
   cd TDS_Group4
   ```
2. Install dependencies listed above 👆.
3. Create "data" folder under the directory.
4. Locate the UKB and HES data and save in the data folder.
5. Download the summarised data from IBD Genetics Consortium from [here](https://www.ibdgenetics.org/downloads.html) and save in the data folder.
6. Change
7. Run the pipeline in the following order.
    ```sh
    qsub Job_Submission_DataCreation.sh
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
├── Job_Submission_DataCreation.sh.e3221321
├── Job_Submission_DataCreation.sh.o3221321
├── Job_Submission_SensitivityAnalysis1.sh
├── Job_Submission_SensitivityAnalysis2.sh
├── Job_Submission_SensitivityAnalysis3.sh
├── Job_Submission_Step5Only.sh
├── Job_Submission_Step5plusOnly.sh
├── Job_Submission_Step5plusOnly.sh.e3483701
├── Job_Submission_Step5plusOnly.sh.o3483701
├── original
├── RandomForestTest.R
├── RandomForestTest.sh
├── RandomForestTest.sh.e3484204
├── RandomForestTest.sh.e3484228
├── RandomForestTest.sh.o3484204
├── RandomForestTest.sh.o3484228
├── README.md
└── Ss

4 directories, 25 files
bash-4.2$ tree -L 2
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
│   ├── 1-2-phenoscanner_output.csv
│   ├── appendix;\ forest.jpg
│   ├── appendix;\ scatter.jpg
│   ├── MR_logistic_regressions_results.txt
│   ├── MR_plot.jpg
│   ├── PRSvalidation_logistic_regressions_results.txt
│   ├── PRS_validation_plot.jpg
│   ├── results;\ forest+scatter.jpg
│   └── Vennspls.pdf
├── data
├── original
│   ├── example_extraction
│   ├── generating_toy_dataset.R
│   ├── getting_started.R
│   ├── README.txt
│   ├── recoding_disease.R
│   └── Ss
├── RandomForestTest.R
├── RandomForestTest.sh
├── RandomForestTest.sh.e3484204
├── RandomForestTest.sh.e3484228
├── RandomForestTest.sh.o3484204
├── RandomForestTest.sh.o3484228
├── README.md
└── Ss
    ├── Job_Submission_Step1Only.sh
    ├── Job_Submission_Step2Only.sh
    ├── Job_Submission_Step3Only.sh
    ├── Job_Submission_Step4Only.sh
    └── Job_Submission_Step6plusOnly-copy.sh
```

&nbsp;


## Relevant resources 🏖️

* [UK Biobank](https://www.ukbiobank.ac.uk/)
* [UK Biobank - List of Cancer Cases](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=40006)
* [Understanding Hospital Episode Statistics (HES) data](https://biobank.ndph.ox.ac.uk/showcase/showcase/docs/HospitalEpisodeStatistics.pdf)
* [International Inflammatory Bowel Disease Genetics Consortium (IIBDGC)](https://www.ibdgenetics.org/)
* [MendelianRandomisation: An R Package for Performing Mendelian Randomization Analyses Using Summarized Data](https://cran.r-project.org/web/packages/MendelianRandomization/vignettes/Vignette_MR.pdf)
* [Getting Started with HPC at Imperial](https://www.imperial.ac.uk/admin-services/ict/self-service/research-support/rcs/support/getting-started/)
* [PhenoScannerV2](http://www.phenoscanner.medschl.cam.ac.uk/)
