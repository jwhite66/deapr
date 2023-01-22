#!/usr/bin/env python3
"""
#----------------------------------------------------------------------------
# deapr
#  Perform DEAPR process on selected samples of raw data.
# This process was described by Sue Rathe with initial code by Jeremy White.
#
# Parameters:
#   raw data        A tab delimited text file that includes a header row, and
#                    an 'Ensembl ID', a 'Gene Name', and then a variable number
#                    of samples, each one labeled in the header row.
#   proteins        A tab delimited text file that indicates which ensembles have
#                    a protein coding.  There is no header.  Column 1 is ensemble id,
#                    column 2 is 'protein coding' if the ensemble has a protein coding.
#   first group     A comma delimited list of samples to use as the first group
#                    The sample names must match the headers given in the rawdata.
#   first name      A name for the first group.
#   second group    A comma delimited list of samples to use as the second sample set
#                    The sample names must match the headers given in the rawdata.
#   second name     A nmae for the second group
#
# Options:
#   --debug             Produce debug statements and intermediate data files.
#                       Debug output will be written to stderr.
#   --fold-weight=x     Adjust the default percent fold change weight, default 90
#   --minimum-fpkm=x    Establish the default FPKM value, default 0.01
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
"""

import sys
import csv
import locale
import argparse

class Data:
    """ Hold the raw data provided by the files """
    def __init__(self):
        """ Construct the object """
        self.raw = []
        self.proteins = []
    def append_raw(self, row):
        """ Add a row """
        self.raw.append(row)
    def append_proteins(self, row):
        """ Add a row """
        self.proteins.append(row)

def read_data(args):
    """ Load the expected data from the raw data file and
        the protein encodings file and return an object
        that includes all of the raw data """
    data = Data()
    try:
        with open(args.rawdata, encoding=locale.getpreferredencoding()) as csvfile:
            reader = csv.DictReader(csvfile, dialect='excel-tab')
            for row in reader:
                data.append_raw(row)
    except FileNotFoundError:
        print(f"Error: could not read {args.rawdata}", file=sys.stderr)
        return None
    csvfile.close()

    try:
        with open(args.proteins, encoding=locale.getpreferredencoding()) as csvfile:
            reader = csv.DictReader(csvfile, fieldnames=['Ensembl ID', 'protein'],
                                    dialect='excel-tab')
            for row in reader:
                data.append_proteins(row)
    except FileNotFoundError:
        print(f"Error: could not read {args.proteins}", file=sys.stderr)
        return None
    csvfile.close()

    return data

def validate_args(args, data):
    """ Make sure groups are specified correctly """
    return True


def parse_args():
    """ Parse our command line arguments """
    parser = argparse.ArgumentParser(prog = 'deapr',
                        description = 'What the program does')

    parser.add_argument('rawdata')
    parser.add_argument('proteins')
    parser.add_argument('group1')
    parser.add_argument('group1name')
    parser.add_argument('group2')
    parser.add_argument('group2name')
    parser.add_argument('--fold-weight', action='store', default=90)
    parser.add_argument('--minimum-fpkm', action='store', default=0.01)
    parser.add_argument('-d', '--debug')

    args = parser.parse_args()

    return args

def main(args):
    """ The mainline of execution for this module """
    data = read_data(args)
    if not data:
        sys.exit(1)

    if not validate_args(args, data):
        sys.exit(2)

    print(data.proteins[0])

    print("Hi mom")
    print(args)
    print(args.rawdata)
    print(args.fold_weight)

if __name__ == "__main__":
    main(parse_args())
