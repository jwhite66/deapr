#!/usr/bin/env python3 -B
#----------------------------------------------------------------------------
# deapr
#  Perform DEAPR process on selected samples of raw data.
# This process was described by Sue Rathe with initial code by Jeremy White.
#
# Parameters:
#   rawdata         A tab delimited text file that includes a header row, and
#                    an 'Ensembl ID', a 'Gene Name', and then a variable number
#                    of samples, each one labeled in the header row.
#   proteins        A tab delimited text file that indicates which ensembles have
#                    a protein coding.  There is no header.  Column 1 is ensemble id,
#                    column 2 is 'protein coding' if the ensemble has a protein coding.
#   first samples   A comma delimited list of samples to use as the first sample set
#                   The sample names must match the headers given in the rawdata.
#   second samples  A comma delimited list of samples to use as the second sample set
#
# Options:
#   --debug         Produce debug statements and intermediate data files.
#                   Debug output will be written to stderr.
# Output:
#  This will write a csv file suitable for import into Excel to stdout.
#
# Approach:
#Select samples,In this case we are comparing S13-S15 to S16-S18
#,Assign the gene type by referencing the PC Genes worksheet (using VLOOKUP)
#Set min value,Sort Select samples by gene type
#,Copy protein coding genes into set min value worksheet
#,Calculate maximum value for each gene and delete MAX < 1
#,Change all values < 0.01 to 0.01
#Apply DELV-SRMM logic,Cull out genes with low values and low fold changes
#,   Delete any genes where both group averages are < 1
#,   Delete any genes with fold change between -2 and +2
#,Determine if they are DELV of SRMM candidates
#,For DELV:
#,   Look for cases in each group where the MAX/MIN if less than 2 (low variability)
#,   We want to designate as KEEP (2 DELVs?) when both groups are DELVs
#,For SRMM:
#,   We want to designate SRMM genes as keep if:
#,   The SRMM calculation is > 1.5 OR < -1.5 AND the minimum value in the higher group is > 1
#Apply weights,Copy over the samples from the Apply DELV-SRMM logic 
#,   where the 2 DELVs? is 'KEEP' OR the SRMM? Is 'KEEP'
#,Calculate the absolute fold change and the absolute difference of the values
#,Sort the absolute value of the fold change (descending) and apply a rank (starting with 1)
#,Sort the absolute value of the FPKM difference (descending) and apply a rank (starting with 1)
#,Combine 90% of the fold change rank with 10% of the difference rank to get the weighted rank
#,   (We should probably make the % modifiable by the user)
#Final report,For publication
##----------------------------------------------------------------------------
## deapr
##  Perform DEAPR process on selected samples of raw data.
## This process was described by Sue Rathe with initial code by Jeremy White.
##
## Parameters:
##   rawdata         A tab delimited text file that includes a header row, and
##                    an 'Ensembl ID', a 'Gene Name', and then a variable number
##                    of samples, each one labeled in the header row.
##   proteins        A tab delimited text file that indicates which ensembles have
##                    a protein coding.  There is no header.  Column 1 is ensemble id,
##                    column 2 is 'protein coding' if the ensemble has a protein coding.
##   first samples   A comma delimited list of samples to use as the first sample set
##                   The sample names must match the headers given in the rawdata.
##   second samples  A comma delimited list of samples to use as the second sample set
##
## Options:
##   --debug         Produce debug statements and intermediate data files.
##                   Debug output will be written to stderr.
## Output:
##  This will write a csv file suitable for import into Excel to stdout.
##
## Approach:
##  
#
#Select samples,In this case we are comparing S13-S15 to S16-S18
#,Assign the gene type by referencing the PC Genes worksheet (using VLOOKUP)
#Set min value,Sort Select samples by gene type
#,Copy protein coding genes into set min value worksheet
#,Calculate maximum value for each gene and delete MAX < 1
#,Change all values < 0.01 to 0.01
#Apply DELV-SRMM logic,Cull out genes with low values and low fold changes
#,   Delete any genes where both group averages are < 1
#,   Delete any genes with fold change between -2 and +2
#,Determine if they are DELV of SRMM candidates
#,For DELV:
#,   Look for cases in each group where the MAX/MIN if less than 2 (low variability)
#,   We want to designate as KEEP (2 DELVs?) when both groups are DELVs
#,For SRMM:
#,   We want to designate SRMM genes as keep if:
#,   The SRMM calculation is > 1.5 OR < -1.5 AND the minimum value in the higher group is > 1
#Apply weights,Copy over the samples from the Apply DELV-SRMM logic 
#,   where the 2 DELVs? is 'KEEP' OR the SRMM? Is 'KEEP'
#,Calculate the absolute fold change and the absolute difference of the values
#,Sort the absolute value of the fold change (descending) and apply a rank (starting with 1)
#,Sort the absolute value of the FPKM difference (descending) and apply a rank (starting with 1)
#,Combine 90% of the fold change rank with 10% of the difference rank to get the weighted rank
#,   (We should probably make the % modifiable by the user)
#Final report,For publication
#----------------------------------------------------------------------------

