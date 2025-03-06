
This pipeline processes FASTQ sequencing data generated with the SEALigHTS technique to detect splicing anomalies. It consists of 3 scripts that must be executed sequentially.
Each script will begin with require you to select the folder of the run to analyse with an easygui to allow easy folder selection of the run you want to analyse.

## Prerequisites
Ensure you have the following dependencies installed:
- Python 3.x
- pandas
- numpy
- matplotlib
- easygui
- openpyxl
- xlwings

Please, download the all repositories tpo get the following architecture :
- \FASTQ\Run_name
- \Probes_Design
- \out\Run_name (created by the script if it doesn't already exist)

## Execution Order


### 1. Demultiplexing Step
**Script:** `01_SEALigHTS_Demultiplexing_V4.py`

**Description:**
- Takes raw FASTQ sequencing data.
- Demultiplexes reads based on probe sequences.
- Counts junction occurrences.
- Generates output files in `out/<Run_name>` directory.

**Input Requirements:**
- Place FASTQ files in a directory named after the run.
- Ensure the `Probes_Design` folder is in the same directory.
- Ensure the `list_genes` corresponds to the genes in the probe mix and in the Probes_Design folder.

**Output:**
- Demultiplexed read counts stored by patient and gene in `out/<Run_name>`.

### 2. Analysis Step
**Script:** `02_SEALigHTS_Analysis_V4_UMI.py`

**Description:**
- Processes demultiplexed data.
- Normalizes junction counts.
- Identifies potential splicing anomalies.
- Generates quality control plots.

Critieria used for the search of splicing anomalies in V4:
- highly expressed junction definition ("Jonction fortes") : linear canonical splicing junction with a median count > 50UMI
- Exonskipping :
                - Normalised count of the junction >q75+1.5IQR  
                - Normalised count of the junction > 1/25 of the q50 of the canonical junctions inside the skipping (1/50 for PMS2)
                - fold change > 2
                - more than 10% of canonical junctions inside the jump are < q25
                - junction excluded if count < 5 UMI 
- Dup_orcircRNA:
                - Normalised count of the junction >q75+1.5IQR 
                - Normalised count of the junction > 1/3 of the q50 of the canonical junctions inside the circRNA or duplication
                - if q50 =0 : the anomaly "Type" is named Dup+++_or_circRNA
- low junction:
                - Normalised count of the junction < q25 - 1.5IQR
                - Normalised count of the junction < 75% of the q50 of this junction
                - the juncton must be categorised as a "Jonction_forte"
                - Score 11 = ratio of the junction compared to the upstream and downstream junctions is < q25-1.5IQR
                - Score 10 = only the ratio with the upstream junction is < q25
                - Score 01 = only the ratio with the downstream junction is < q25
- Low_gene_expression: 
                Calculation of the ratio of each junction of a gene between the patient and the median for the run 
                - the median of each ratio must be < 2/3 (66%)
                - more than 80% of the junctions must be below < q25
       

**Input Requirements:**
- The demultiplexed data from Step 1.
- The `Probes_Design` folder.

**Output:**
- Normalized junction matrices.
- Anomaly reports in `out/<Run_name>/Results`.

### 3. Result Formatting
**Script:** `03_Result_formatting.py`

**Description:**
- Formats anomaly reports into an Excel file.
- Adds hyperlinks for easy access to figures and matrices.

**Input Requirements:**
- Anomaly reports from Step 2 (`out/<Run_name>/Results`).

**Output:**
- A formatted Excel file saved as `Results_<Run_name>_formatted.xlsx`.

## Notes
- Ensure all scripts are executed in order.
- The scripts prompt user interaction for selecting directories.
- The formatted results file contains hyperlinks to key figures and matrices.

For further details, refer to the inline comments within each script.
