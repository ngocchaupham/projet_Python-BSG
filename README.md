# Python-project : Comparison of Nucleotide and Protein Sequences Between Eukaryotes and Prokaryotes

---

### Overview
In this project, my objective is to compare nucleotide and protein sequences between a eukaryotic organism and a prokaryotic organism. To achieve this, I utilized the Python programming language and bash commands to prepare and download the necessary data, then extract the sequences of genes of interest, analyze their characteristics, and compare the results between the two organisms.
This folder contains a collection of bash and Python scripts to extract data from FASTA and GTF files of an organism, analyze sequences, and compare certain aspects between two organisms. To do this:

1. Open the Terminal.
2. Run the [`Load_Data.sh`](/Load_Data.sh) script: This will download the necessary data from the internet.
3. Then, run the [`run.sh`](/run.sh) script: This will execute the Python scripts. Please wait a few moments for the process to complete.
4. Check the results in the [`data`](/data) folder: You will find the newly generated FASTA files, the comparison report, and figures.
   
---

# Data Source Descriptions

The organisms chosen for this project were Homo sapiens (chromosome 1) as the eukaryotic organism and E. coli as the prokaryotic organism. Two types of files in FASTA and GTF formats were downloaded for each organism from the Ensembl website:
- The FASTA files contain the complete genomic sequence of each organism.
- The GTF files provide annotations for each genomic sequence, including the locations of genes and exons.

The definition of the working directory and the downloading of source files were carried out via a bash script named `Load_data.sh`. All data are stored in a folder named `data`.

# Software and Python Libraries Used

For this project, I used the Pycharm code editor to write the Python program. Several Python libraries were also used to perform various project steps, including: os, sys, random, statistics, pandas, and matplotlib.
Comments are included within each script to explain the code and the methodology I used to extract the sequences of the genes of interest, analyze their characteristics, and compare the results.


