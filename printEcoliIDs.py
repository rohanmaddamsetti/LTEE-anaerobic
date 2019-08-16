#!/usr/bin/env python

## Usage: python printEcoliIDs.py -i CP001396.gbk > MC4100_IDs.csv

import argparse
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(description='print csv file of IDs in Genbank file.')

    parser.add_argument("-i", "--input", help='Input genbank file',required=True)
    args = parser.parse_args()

    genome = next(SeqIO.parse(args.input, "genbank"))
    print(','.join(['gene','locus_tag','blattner','length','product']))
    for feat in genome.features:
        ## just print out protein-coding genes, skip the rest.
        if feat.type != 'CDS':
            continue
        try:
            gene = feat.qualifiers['gene'].pop()
        except KeyError:
            gene = 'NA'
        locus_tag = feat.qualifiers['locus_tag'].pop()
        product = feat.qualifiers['product'].pop()
        try:
            note = feat.qualifiers['note'].pop()
            blattner = note.split('MG1655 equivalent: ')[-1]
        except:
            blattner = 'NA'
        length = feat.location.end - feat.location.start
        
            
        print(','.join([str(x) for x in (gene,locus_tag,blattner,length,product)]))

main()
