#!/usr/bin/env python

## Usage: python printEcoliIDs.py -i CP001396.gbk > MC4100_IDs.csv

import argparse
import re
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(description='print csv file of IDs in Genbank file.')

    parser.add_argument("-i", "--input", help='Input genbank file',required=True)
    args = parser.parse_args()

    genome = next(SeqIO.parse(args.input, "genbank"))
    print(','.join(['Gene','locus_tag','blattner','gene_length','product']))
    for feat in genome.features:
        ## just print out protein-coding genes, skip the rest.
        if feat.type != 'CDS':
            continue
        length = feat.location.end - feat.location.start
        locus_tag = feat.qualifiers['locus_tag'].pop()
        try:
            gene = feat.qualifiers['gene'].pop()
        except KeyError:
            gene = locus_tag
        try:
            product = feat.qualifiers['product'].pop()
        except:
            product = 'NA'
        try:
            note = feat.qualifiers['note'].pop()
            blattner = note.split('MG1655 equivalent: ')[-1].split()[-1]
            blattner = re.sub('[,;()]', '', blattner)
        except:
            blattner = 'NA'

        ## strip all punctuation that could cause parsing problems.
        product = re.sub('[,;()]', '', product)
        ## I should not have to do this next-- fix these corner cases later.
        blattner = re.sub('[,;()]', '', blattner)
        
        print(','.join([str(x) for x in (gene,locus_tag,blattner,length,product)]))

main()
