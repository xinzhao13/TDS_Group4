# HDA TDS Group 4
## About this project 📽️

This project is prepared for the Translational Data Science module as part of MSc Health Data Analytics and Machine Learning. Group 4 has been tasked to utilise UK Biobank and implement one of the most cutting-edge techniques Mendelian Randomisation in Computational Epidemiology to assess the causal relationship between Inflammatory Bowel Disease (IBD), its subtypes Crohn's Disease (CD) and Ulcerative Colitis (UC), and Colon Cancer.

## Project team 🧑‍🤝‍🧑

* CID1
* CID2
* CID3
* CID4
* CID5

## Getting started ✈️

### Prerequisites 💻

Below are the list of packages and the version numbers required to reproduce the project.

  ```sh
  npm install npm@latest -g
  ```

### Installation 🖱️

1. Clone the repo.
   ```sh
   git clone https://github.com/xinzhao13/HDA_TDS.git
   ```
2. Install dependencies listed above 👆.
3. Run the pipeline.
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
### Key files 📂

## Relevant resources 🏖️

* [UK Biobank](https://www.ukbiobank.ac.uk/)
* [UK Biobank - List of Cancer Cases](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=40006)
* [Understanding Hospital Episode Statistics (HES) data](https://biobank.ndph.ox.ac.uk/showcase/showcase/docs/HospitalEpisodeStatistics.pdf)
* [International Inflammatory Bowel Disease Genetics Consortium (IIBDGC)](https://www.ibdgenetics.org/)
* [MendelianRandomisation: An R Package for Performing Mendelian Randomization Analyses Using Summarized Data](https://cran.r-project.org/web/packages/MendelianRandomization/vignettes/Vignette_MR.pdf)
* [Getting Started with HPC at Imperial](https://www.imperial.ac.uk/admin-services/ict/self-service/research-support/rcs/support/getting-started/)
