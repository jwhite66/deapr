#!/usr/bin/env python3
"""
#----------------------------------------------------------------------------
#   deapr - Perform DEAPR process on selected samples of raw data.
#
# DEAPR compares gene expression from small groups of experimental samples and
# their associated controls to identify genes that are differentially expressed.
# It then applies a weighting strategy based on both fold change and FPKM difference
# to rank the genes from most differentially expressed to least.
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
#  Apply DELV-SRMM logic, cull out genes with low values and low fold changes
#    Delete any genes where both group averages are < 1
#    Delete any genes with fold change between -2 and +2
#  Determine if they are DELV of SRMM candidates
#  For DELV:
#    Look for cases in each group where the MAX/MIN if less than 2 (low variability)
#    We want to designate as KEEP (2 DELVs?) when both groups are DELVs
#  For SRMM:
#    We want to designate SRMM genes as keep if:
#    The SRMM calculation is > 1.5 OR < -1.5 AND the minimum value in the higher group is > 1
#  Apply weights
#    First, keep only samples with SRMM to keep or 2 DELVs
#    Calculate the absolute fold change and the absolute difference of the values
#    Sort the absolute value of the fold change (descending) and apply a rank (starting with 1)
#    Sort the absolute value of the FPKM difference (descending) and apply a rank (starting with 1)
#    Combine 90% of the fold change rank with 10% of the difference rank to get the weighted rank
#
#  Write out a final report
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
#  then we do not include the gene
SAMPLE_THRESHOLD = 1.0

# If the averages of neither sample pass muster we do not include
#  the row in the second pass results
MAX_AVG_CUTOFF = 1.0

# If absolute value of the fold change is not more than this cutoff,
#  we do not include it in the second pass results.
FOLD_CHANGE_CUTOFF = 2.0

# If there is a lot of variability in the sample, we may not want
#  it.  This is tracked as 'DELV'
DELV_CUTOFF = 2.0

# To keep the SRMM, we require that the samples included in it
#  have a minimum above 1
SRMM_MIN_CUTOFF = 1.0

# And we want the SRMM itself to be 'interesting'
SRMM_THRESHOLD = 1.5

#------------------------------------------------------------------------------
#  The Ensemble and Data class hold the bulk of the 'business logic'
#   They are the meat of this application.
#------------------------------------------------------------------------------
class Ensemble:
    """ Hold the data for a specific Ensemble ID aka Gene aka row in the table """
    def __init__(self, row, group1, group2, minimum_fpkm):
        """ Make a record for a specific data row """
        self.eid = row['Ensembl ID']
        self.gene_name = row['Gene Name']

        self.group1 = np.array(list(float(row[sel]) for sel in group1))
        self.group1[self.group1 < minimum_fpkm] = minimum_fpkm
        self.group1_min = np.amin(self.group1)
        self.group1_max = np.amax(self.group1)
        self.group1_avg = np.average(self.group1)

        self.group2 = np.array(list(float(row[sel]) for sel in group2))
        self.group2[self.group2 < minimum_fpkm] = minimum_fpkm
        self.group2_min = np.amin(self.group2)
        self.group2_max = np.amax(self.group2)
        self.group2_avg = np.average(self.group2)

        self.max_avg = max(self.group1_avg, self.group2_avg)
        self.max_value = max(self.group1_max, self.group2_max)
        self.avg_diff = abs(self.group1_avg - self.group2_avg)

        self.keep_srmm = True
        self.srmm = 0.0

    def calculate_srmm(self):
        """ Calculate the srmm and indicate if we should keep it """
        self.keep_srmm = True
        if self.group2_min > self.group1_max:
            self.srmm = self.group2_min / self.group1_max
            if self.group2_min < SRMM_MIN_CUTOFF:
                self.keep_srmm = False
        elif self.group1_min > self.group2_max:
            self.srmm = -1 * self.group1_min / self.group2_max
            if self.group1_min < SRMM_MIN_CUTOFF:
                self.keep_srmm = False
        else:
            self.srmm = 1
        if abs(self.srmm) <= SRMM_THRESHOLD:
            self.keep_srmm = False


def sort_fold(ensemble):
    """ Support sorting by fold, ideally descending """
    return abs(ensemble.fold_change)

def sort_diff(ensemble):
    """ Support sorting by the delta in the averages, ideally descending """
    return ensemble.avg_diff

def sort_weighted(ensemble):
    """ Sort by the weighted rank """
    return ensemble.weighted_rank

class Data:
    """ Hold the raw data provided by the files """
    def __init__(self):
        """ Construct the object """
        self.raw = []
        self.proteins = []
        self.pass1 = []
        self.pass2 = []
        self.pass3 = []

    def run_pass1(self, args):
        """ Create a subset of the raw data that includes only proteins,
            and the data requested by the groups.  Prune and adjust it a bit as well """
        for row in self.raw:
            if row['Ensembl ID'] in self.proteins:
                ensemble = Ensemble(row, args.group1.split(","), args.group2.split(","),
                                    args.minimum_fpkm)

                if ensemble.max_value >= SAMPLE_THRESHOLD:
                    self.pass1.append(ensemble)
                elif args.debug > 2:
                    print("Pruning " + row['Ensembl ID'] +
                          "; max is " + str(ensemble.max_value), file=sys.stderr)
            elif args.debug > 2:
                print("Not a protein " + row['Ensembl ID'], file=sys.stderr)

    def run_pass2(self, args):
        """ Create a subset of the selected set that applies DELV and SRMM logic. """
        for ensemble in self.pass1:
            if ensemble.max_avg <= MAX_AVG_CUTOFF:
                if args.debug > 2:
                    print("Pruning " + ensemble.eid +
                          "; max_avg is " + str(ensemble.max_avg), file=sys.stderr)
                continue

            if ensemble.group2_avg > ensemble.group1_avg:
                ensemble.fold_change = ensemble.group2_avg / ensemble.group1_avg
            else:
                ensemble.fold_change = ensemble.group1_avg / ensemble.group2_avg * -1.0

            if abs(ensemble.fold_change) <= FOLD_CHANGE_CUTOFF:
                if args.debug > 2:
                    print("Pruning " + ensemble.eid +
                          "; fold change is " + str(ensemble.fold_change), file=sys.stderr)
                continue

            ensemble.calculate_srmm()

            ensemble.delv1 = False
            if ensemble.group1_max / ensemble.group1_min < DELV_CUTOFF:
                ensemble.delv1 = True
            ensemble.delv2 = False
            if ensemble.group2_max / ensemble.group2_min < DELV_CUTOFF:
                ensemble.delv2 = True

            self.pass2.append(ensemble)

    def run_pass3(self, args):
        """ Create a subset that has either both DELVs or SRMM
            and then establish a ranking """
        for ensemble in self.pass2:
            if not ensemble.keep_srmm and not (ensemble.delv1 and ensemble.delv2):
                if args.debug > 2:
                    print(f"Pruning {ensemble.eid}; keep_srmm {ensemble.keep_srmm}; " +
                          f"delve {ensemble.delv1}/{ensemble.delv2}",
                          file=sys.stderr)
                continue


            self.pass3.append(ensemble)

        self.pass3.sort(key=sort_fold, reverse=True)
        for i, ensemble in enumerate(self.pass3):
            ensemble.fold_rank = i + 1

        self.pass3.sort(key=sort_diff, reverse=True)
        for i, ensemble in enumerate(self.pass3):
            ensemble.diff_rank = i + 1
            ensemble.weighted_rank = (args.fold_weight_float * ensemble.fold_rank) + \
                              (1.0 - args.fold_weight_float) * ensemble.diff_rank

        self.pass3.sort(key=sort_weighted)

def write_report(args, outfile, ensembles):
    """ Produce the final report """
    outfile.write("\n,,,FPKMs\n")
    outfile.write(f",,,{args.group1name},,,{args.group2name},,,,,,,,,,,,\n")
    outfile.write(",Selected by,,,,,,,,,,,,,,Rankings\n")

    outfile.write("Gene Name,DELV?,SRMM?,")
    for sample in args.group1.split(","):
        outfile.write(sample + ",")
    for sample in args.group2.split(","):
        outfile.write(sample + ",")
    outfile.write("Avg Grp 1,Avg Grp 2,")
    outfile.write("Fold Chg,Abs Fold Chg,Abs Diff,")
    outfile.write("Abs Fold Chg,Abs Diff,Weighted\n")

    for ensemble in ensembles:
        outfile.write(f"{ensemble.gene_name},")
        delv = ""
        if ensemble.delv1 and ensemble.delv2:
            delv = "Yes"
        srmm = ""
        if ensemble.keep_srmm:
            srmm = "Yes"
        outfile.write(f"{delv},{srmm},")

        for sample in ensemble.group1:
            outfile.write(f"{sample:.3f},")
        for sample in ensemble.group2:
            outfile.write(f"{sample:.3f},")

        outfile.write(f"{ensemble.group1_avg:.3f},{ensemble.group2_avg:.3f},")
        outfile.write(f"{ensemble.fold_change:.3f},{abs(ensemble.fold_change):.3f},")
        outfile.write(f"{abs(ensemble.avg_diff):.3f},")
        outfile.write(f"{ensemble.fold_rank},")
        outfile.write(f"{ensemble.diff_rank},")
        outfile.write(f"{abs(ensemble.weighted_rank):.1f}")
        outfile.write("\n")



#------------------------------------------------------------------------------
#  The following functions are just for debug purpose; they let us create
#   .csv files which show the intermediate steps.  That is particularly
#   useful for comparing to Sue's original spreadsheets.
#------------------------------------------------------------------------------
def write_selected(list_to_use, args, fname, use_pass2, use_pass3):
    """ For debug purposes, write the selected data out. """
    with open(fname, "w", encoding=locale.getpreferredencoding()) as outfile:
        outfile.write("Ensembl ID,Gene Name,")
        for sample in args.group1.split(","):
            outfile.write(sample + ",")
        for sample in args.group2.split(","):
            outfile.write(sample + ",")

        if not use_pass2:
            outfile.write("Max")
        else:
            outfile.write("Avg Grp 1,Avg Grp 2,Fold chg (2+/-),Max Avg > 1," +
                          "Min,Max,Max/Min < 2,DELV?,Min,Max ,Max/Min < 2," +
                          "DELV?,2 DELVs?,SRMM,SRMM?")
        if use_pass3:
            outfile.write(",Abs fold,Abs diff,Rank fold,Rank diff,Wgt rank")
        outfile.write("\n")

        for ensemble in list_to_use:
            write_pass1(outfile, ensemble, not use_pass2)
            if use_pass2:
                write_pass2(outfile, ensemble)
            if use_pass3:
                write_pass3(outfile, ensemble)
            outfile.write("\n")
        outfile.close()

def write_pass1(outfile, ensemble, write_max):
    """ Write out basic information common to all passes """
    outfile.write(f"{ensemble.eid},{ensemble.gene_name},")
    for sample in ensemble.group1:
        outfile.write(f"{sample},")
    for sample in ensemble.group2:
        outfile.write(f"{sample},")
    if write_max:
        outfile.write(f"{ensemble.max_value},")

def write_pass2(outfile, ensemble):
    """ Debug function to write out the values after pass2 """
    outfile.write(f"{ensemble.group1_avg},")
    outfile.write(f"{ensemble.group2_avg},")
    outfile.write(f"{ensemble.fold_change},")
    outfile.write(f"{ensemble.max_avg},")
    outfile.write(f"{ensemble.group1_min},")
    outfile.write(f"{ensemble.group1_max},")
    mcalc = ensemble.group1_max / ensemble.group1_min
    outfile.write(f"{mcalc},")
    if ensemble.delv1:
        outfile.write("KEEP,")
    else:
        outfile.write(",")

    outfile.write(f"{ensemble.group2_min},")
    outfile.write(f"{ensemble.group2_max},")
    mcalc = ensemble.group2_max / ensemble.group2_min
    outfile.write(f"{mcalc},")
    if ensemble.delv2:
        outfile.write("KEEP,")
    else:
        outfile.write(",")
    if ensemble.delv1 and ensemble.delv2:
        outfile.write("KEEP,")
    else:
        outfile.write(",")
    outfile.write(f"{ensemble.srmm},")
    if ensemble.keep_srmm:
        outfile.write("KEEP")

def write_pass3(outfile, ensemble):
    """ Write data out generated in pass 3 """
    outfile.write(",")
    outfile.write(f"{abs(ensemble.fold_change)},")
    outfile.write(f"{abs(ensemble.avg_diff)},")
    outfile.write(f"{ensemble.fold_rank},")
    outfile.write(f"{ensemble.diff_rank},")
    outfile.write(f"{abs(ensemble.weighted_rank)}")


#------------------------------------------------------------------------------
#  Mainline and utility functions.
#   Below are the functions for parsing arguments, reading data, and so on,
#   including the program mainline.
#------------------------------------------------------------------------------
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
    except OSError:
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
    except OSError:
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
                        description = 'Apply a process to compare gene expression samples.')

    parser.add_argument('rawdata')
    parser.add_argument('proteins')
    parser.add_argument('group1')
    parser.add_argument('group1name')
    parser.add_argument('group2')
    parser.add_argument('group2name')
    parser.add_argument('--output', action='store', default='-')
    parser.add_argument('--fold-weight', action='store', default=90)
    parser.add_argument('--minimum-fpkm', action='store', type=float, default=0.01)
    parser.add_argument('-d', '--debug', action='store', type=int, default=0)

    args = parser.parse_args()

    args.fold_weight_float = float(args.fold_weight) / 100.0

    return args

def main(args):
    """ The mainline of execution for this module """
    data = read_data(args)
    if not data:
        sys.exit(1)

    if not validate_args(args, data):
        sys.exit(2)

    if args.output == "-":
        outfile = sys.stdout
    else:
        try:
            outfile = open(args.output, "w", encoding=locale.getpreferredencoding())
        except OSError:
            print(f"Error: unable to open {args.output}", file=sys.stderr)
            sys.exit(3)

    data.run_pass1(args)

    if args.debug > 1:
        write_selected(data.pass1, args, "pass1.csv", False, False)

    data.run_pass2(args)
    if args.debug > 1:
        write_selected(data.pass2, args, "pass2.csv", True, False)

    data.run_pass3(args)
    if args.debug > 1:
        write_selected(data.pass3, args, "pass3.csv", True, True)

    write_report(args, outfile, data.pass3)


if __name__ == "__main__":
    main(parse_args())
