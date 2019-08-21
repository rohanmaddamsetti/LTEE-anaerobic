#!/usr/bin/env python

''' 
arcA_motif_fasta.py by Rohan Maddamsetti.

Usage: python arcA_motif_to_fasta.py > ../results/Park2013_arcA_motifs.fasta

## Then upload Park2013_arcA_motifs.fasta to http://meme-suite.org/tools/meme.



'''

def main():

    arcA_motifs_f = "../data/Park2013-arcA-binding-sites.csv"
    arcA_motifs_fh = open(arcA_motifs_f,'r')
    for i,l in enumerate(arcA_motifs_fh):
        if i == 0: ## skip header
            continue
        l = l.strip()
        fields = l.split(',')
        gene, strand, coordinate, sequence, ri_bits = fields
        ##header = '>' + '|'.join(['Gene='+gene,'Strand='+strand,'Coordinate='+coordinate,'Ri_bits='+ri_bits])
        header = '>' + 'Gene='+ gene+'|'+'Coordinate='+coordinate
        print(header)
        print(sequence)
    
main()
