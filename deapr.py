#!/usr/bin/env python3
"""
#----------------------------------------------------------------------------
# deapr
#  Perform DEAPR process on selected samples of raw data.
#
# This process was described by Sue Rathe and initial code was written
# by Jeremy White.
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
#   first name      A name for the first group. To be displayed in the final report.
#   second group    A comma delimited list of samples to use as the second sample set
#                    The sample names must match the headers given in the rawdata.
#   second name     A name for the second group. To be displayed in the final report.
#
# Options:
#   --debug=lvl         Produce debug statements and intermediate data files.
#                       Debug output will be written to stderr.
#                       Level is verbosity, starting at 1, higher means more
#   --fold-weight=x     Adjust the default percent fold change weight, default 90
#   --minimum-fpkm=x    Establish the default FPKM value, default 0.01
# Output:
#  This will write a csv file suitable for import into Excel to stdout.
#
# Approach:
#  Select only rows that have a protein coding
#  Select the two groups specified
#  Apply the minimum fpkm to all groups
#  If the maximum for a given pair of groups is below a threshold,
#   do not include it
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
#
#
#----------------------------------------------------------------------------
# Copyright (c) 2023 Jeremy P White
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#----------------------------------------------------------------------------
"""

import sys
import csv
import locale
import argparse
import numpy as np

# We have a few constants we use that are not parameterized
#  We define them here

# If neither selected group has a value greater or equal to this amount,
#  then we do not include the ensemble
SAMPLE_CUTOFF = 1.0

class Ensemble:
    """ Hold the data for a specific Ensemble ID """
    def __init__(self, row, group1, group2):
        """ Make a record for a specific data row """
        self.ensemble_id = row['Ensembl ID']
        self.gene_name = row['Gene Name']
        self.group1 = np.array(list(float(row[sel]) for sel in group1))
        self.group2 = np.array(list(float(row[sel]) for sel in group2))
        self.max_value = np.amax(np.concatenate([self.group1, self.group2]))

class Data:
    """ Hold the raw data provided by the files """
    def __init__(self):
        """ Construct the object """
        self.raw = []
        self.proteins = []
        self.selected = []

    def select(self, args):
        """ Create a subset of the raw data that includes only proteins,
            and the data requested by the groups.  Prune and adjust it a bit as well """
        for row in self.raw:
            if row['Ensembl ID'] in self.proteins:
                ensemble = Ensemble(row, args.group1.split(","), args.group2.split(","))
                for i, value in enumerate(ensemble.group1):
                    if value < args.minimum_fpkm:
                        ensemble.group1[i] = args.minimum_fpkm
                for i, value in enumerate(ensemble.group2):
                    if value < args.minimum_fpkm:
                        ensemble.group2[i] = args.minimum_fpkm

                if ensemble.max_value >= SAMPLE_CUTOFF:
                    self.selected.append(ensemble)
                elif args.debug > 2:
                    print("Pruning " + row['Ensembl ID'] + "; max is " + str(ensemble.max_value))
            elif args.debug > 2:
                print("Not a protein " + row['Ensembl ID'])

    def write_selected(self, args, fname):
        """ For debug purposes, write the selected data out. """
        with open(fname, "w", encoding=locale.getpreferredencoding()) as outfile:
            outfile.write("Ensembl ID,Gene Name,")
            for sample in args.group1.split(","):
                outfile.write(sample + ",")
            for sample in args.group2.split(","):
                outfile.write(sample + ",")
            outfile.write("Max\n")

            for ensemble in self.selected:
                outfile.write(f"{ensemble.ensemble_id},{ensemble.gene_name},")
                for sample in ensemble.group1:
                    outfile.write(f"{sample},")
                for sample in ensemble.group2:
                    outfile.write(f"{sample},")
                outfile.write(str(ensemble.max_value) + "\n")
            outfile.close()

def read_data(args):
    """ Load the expected data from the raw data file and
        the protein encodings file and return an object
        that includes all of the raw data """
    data = Data()
    try:
        with open(args.rawdata, encoding=locale.getpreferredencoding()) as csvfile:
            reader = csv.DictReader(csvfile, dialect='excel-tab')
            for row in reader:
                data.raw.append(row)
    except FileNotFoundError:
        print(f"Error: could not read {args.rawdata}", file=sys.stderr)
        return None
    csvfile.close()

    try:
        with open(args.proteins, encoding=locale.getpreferredencoding()) as csvfile:
            reader = csv.DictReader(csvfile, fieldnames=['Ensembl ID', 'protein'],
                                    dialect='excel-tab')
            for row in reader:
                if row['protein'] == 'protein coding':
                    data.proteins.append(row['Ensembl ID'])
    except FileNotFoundError:
        print(f"Error: could not read {args.proteins}", file=sys.stderr)
        return None
    csvfile.close()

    return data

def validate_args(args, data):
    """ Make sure groups are specified correctly """
    valid_names = []
    for name in data.raw[0]:
        valid_names.append(name)
    if not (valid_names[0] == 'Ensembl ID' and valid_names[1] == 'Gene Name'):
        print(f"Error: {args.rawdata} does not have expected column headers", file=sys.stderr)
        return False

    valid_names = valid_names[2:]
    for sample in args.group1.split(","):
        if not sample in valid_names:
            print(f"Error: {sample} is not a valid sample name", file=sys.stderr)
            return False

    for sample in args.group2.split(","):
        if not sample in valid_names:
            print(f"Error: {sample} is not a valid sample name", file=sys.stderr)
            return False
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
    parser.add_argument('--minimum-fpkm', action='store', type=float, default=0.01)
    parser.add_argument('-d', '--debug', action='store', type=int, default=0)

    args = parser.parse_args()

    return args

def main(args):
    """ The mainline of execution for this module """
    data = read_data(args)
    if not data:
        sys.exit(1)

    if not validate_args(args, data):
        sys.exit(2)

    data.select(args)

    if args.debug > 1:
        data.write_selected(args, "selected.csv")

    #print(data.selected[18].ensemble_id)
    #print(data.selected[18].gene_name)
    #print(data.selected[18].group1)
    #print(data.selected[18].group2)

if __name__ == "__main__":
    main(parse_args())
