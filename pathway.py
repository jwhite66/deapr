#!/usr/bin/env python3
"""
#----------------------------------------------------------------------------
#   pathway - Pathway ranking for DEAPR proceesed genes.
#
# This process was described by Sue Rathe and initial code was written
# by Jeremy White.
#
# Parameters:
#   deapr output    A tab delimted file that was produced by the deapr python program.
#   GA output       A tab delimited file that contains a Score, a SuperPath name,
#                    a superpath count of genes, a superpath count of matched genes,
#                    and then a quote enlcosed comma delimited list of matched genes
#
# Options:
#   --debug=lvl         Produce debug statements and intermediate data files.
#                       Debug output will be written to stderr.
#                       Level is verbosity, starting at 1, higher means more
#   --output=filename   Use a file other than stdout for the report
#
# Output:
#  This will write a csv file suitable for import into Excel to stdout.
#  It will return 0 on success, > 0 for a failure, with an error written to stderr
#
# Process:
#  For each DEAPR result, we invert the rank to get a weight.  Higher is better.
#  For each pathway, we add up the weights of each involved gene to get a total gene weight.
#  We then compute an average weight across all pathways.  Finally, we use the Score
#  given by the GA process and multiply it by the ratio of a pathways weight to the average
#  weight to get a weighted score.
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

# Maximum number of genes taken from the deapr output
MAX_GENES = 400

#------------------------------------------------------------------------------
#  The Data class holds the information we read and synthesize
#------------------------------------------------------------------------------
class Data:
    """ Hold the raw data provided by the files """
    def __init__(self):
        """ Construct the object """
        self.pathways = []
        self.deapr = []
        self.rank = {}
        self.weight = {}
        self.total_weight = 0.0
        self.total_count = 0.0
        self.average_ratio = 0.0

    def build_map(self):
        """ Build a map with gene weights """
        all_genes = {}
        for pathway in self.pathways:
            for gene in pathway['genes']:
                all_genes[gene.upper().strip()] = True

        max_weight = float(self.deapr[-1]['Weighted'])

        for gene in self.deapr:
            weighted = float(gene['Weighted'])
            if gene['Gene Name'] in all_genes:
                self.rank[gene['Gene Name']] = weighted
                self.weight[gene['Gene Name']] = max_weight + 1 - weighted

    def run_pass1(self):
        """ Deconstruct the data """
        for pathway in self.pathways:
            weight = 0.0
            for gene in pathway['genes']:
                if gene in self.weight:
                    weight += self.weight[gene]
                else:
                    print("Missing gene: " + gene, file=sys.stderr)
            pathway['total_weight'] = weight
            self.total_count += float(pathway['# SuperPath Matched Genes'])
            self.total_weight += weight

        self.average_ratio = self.total_weight / self.total_count

        for pathway in self.pathways:
            wgt_score = pathway['total_weight'] / float(pathway['# SuperPath Matched Genes'])
            wgt_score /= self.average_ratio
            wgt_score *= float(pathway['Score'])
            pathway['wgt_score'] = wgt_score

        self.pathways.sort(key=sort_weighted, reverse=True)


    def print_report(self, outfile):
        """ Print the desired report """
        print("SuperPath Name,Cnt,Ttl Gene Wgt,GA Score,Path Ttl,Wgt Score", file=outfile)
        for pathway in self.pathways:
            path = pathway['SuperPath Name']
            if path.find(",") >= 0:
                path = '"' + path + '"'
            print("{},{},{},{},{},{}".format(path,
               pathway['# SuperPath Matched Genes'],
               pathway['total_weight'],
               pathway['Score'],
               pathway['# SuperPath Total Genes'],
               pathway['wgt_score']), file=outfile)
        print("TOTAL,{},{}".format(self.total_count, self.total_weight), file=outfile)
        print("AVERAGE,,{}".format(self.average_ratio), file=outfile)


    def print_pathways(self):
        """ Debug method to print the pathways """
        print("SuperPath Name, Gene")
        for pathway in self.pathways:
            for gene in pathway['genes']:
                print("{},{}".format(pathway['SuperPath Name'], gene))

    def print_assignments(self):
        """ Debug method to print the assignments """
        print("Gene,Rank,SuperPath Name,Gene Wgt (Max rank-Rank+1)")
        for pathway in self.pathways:
            for gene in pathway['genes']:
                if gene in self.rank:
                    path = pathway['SuperPath Name']
                    if path.find(",") >= 0:
                        path = '"' + path + '"'
                    print("{},{},{},{}".format(gene, self.rank[gene], path, self.weight[gene]))

def sort_weighted(pathway):
    """ Sort by the weighted score """
    return pathway['wgt_score']

#------------------------------------------------------------------------------
#  Mainline and utility functions.
#   Below are the functions for parsing arguments, reading data, and so on,
#   including the program mainline.
#------------------------------------------------------------------------------
def read_data(args):
    """ Load the expected data from the raw data file """
    data = Data()
    try:
        with open(args.pathway, encoding=locale.getpreferredencoding()) as csvfile:
            reader = csv.DictReader(csvfile, dialect='excel-tab')
            for row in reader:
                row['genes'] = row['Matched Genes (Symbols)'].split(",")
                for i, _ in enumerate(row['genes']):
                    row['genes'][i] = row['genes'][i].strip()
                data.pathways.append(row)
    except OSError:
        print(f"Error: could not read {args.pathway}", file=sys.stderr)
        return None
    csvfile.close()

    try:
        with open(args.deapr, encoding=locale.getpreferredencoding()) as csvfile:
            fields = None
            for i in range(0, 5):
                line = csvfile.readline()
                if line.startswith("Gene Name"):
                    fields = line.strip().split(",")
                    break
            if fields is None:
                print(f"Error: {args.deapr} does not appear to have gene data", file=sys.stderr)
                sys.exit(1)
            reader = csv.DictReader(csvfile, fields)
            for row in reader:
                data.deapr.append(row)
            data.deapr = data.deapr[:MAX_GENES]
    except OSError:
        print(f"Error: could not read {args.deapr}", file=sys.stderr)
        return None
    csvfile.close()

    return data

def parse_args():
    """ Parse our command line arguments """
    parser = argparse.ArgumentParser(prog = 'pathway',
                        description = 'Perform pathway analysys in a DEAPR set.')

    parser.add_argument('deapr')
    parser.add_argument('pathway')
    parser.add_argument('--output', action='store', default='-')
    parser.add_argument('-d', '--debug', action='store', type=int, default=0)

    args = parser.parse_args()

    return args

def main(args):
    """ The mainline of execution for this module """
    data = read_data(args)
    if not data:
        sys.exit(1)

    if args.output == "-":
        outfile = sys.stdout
    else:
        try:
            outfile = open(args.output, "w", encoding=locale.getpreferredencoding())
        except OSError:
            print(f"Error: unable to open {args.output}", file=sys.stderr)
            sys.exit(3)

    data.build_map()
    #data.print_assignments()
    data.run_pass1()
    data.print_report(outfile)

if __name__ == "__main__":
    main(parse_args())
