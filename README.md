# DEAPR Overview
Dr. Sue Rathe created a novel approach to prioritize the gene and pathway
changes in controlled experiments (DEAPR).

DEAPR compares gene expression from small groups of experimental samples and
their associated controls to identify genes that are differentially expressed.
It then applies a weighting strategy based on both fold change and FPKM difference
to rank the genes from most differentially expressed to least.

This repository contains two simple Python utilities that implement
analysis tools for DEAPR and for pathways.

The code was written by Jeremy White under direction by Sue Rathe.

# Running the programs

The programs are simple python utilities.  Anyone familiar with a command line
environment should find them simple to run.  Each utility has full documentation of it's usage. 

On Windows, these directions may be helpful:

## First time setup:
 - Visit python.org and download Python
 - Install python - make sure to click the 'add Python to the path' checkbox
 - Open a command prompt, and type:
      `pip install numpy`
   This installs the Python numerical analysis library

## Each time you run:
Then, each time you want to run it, do the following:
 - Put the raw data file, proteins data file, and the program into a directory somewhere.  For my testing, I put them in c:\users\throw\deapr, but it can be anywhere.
 - Open a command prompt
 - Change to the directory you made; in my case that was:
      `cd c:\users\throw\deapr`
  - Run the program.  You will need to pass the parameters as described in the deapr.py program. Here is one example:
      `python deapr.py Raw_data.txt Protein_coding_genes.txt S13,S14,S15 "First group" S16,S17,S18 "Second group" --output result.csv`
    That put the resulting output into a file called result.csv.
