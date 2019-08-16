#!/usr/bin/env python

''' 
make_gene_to_blattner_csv.py by Rohan Maddamsetti.

this script makes a table of gene names to blattner numbers using the annotation
downloaded from the EcoCyc database in the file called
All_instances_of_Genes_in_Escherichia_coli_K-12_substr._MG1655.txt.

Usage: python make_gene_to_blattner_csv.py > Ecocyc_gene_to_blattner.csv

'''

import re

def main():

    gene_to_blattner = {}
    
    f = "All_instances_of_Genes_in_Escherichia_coli_K-12_substr._MG1655.txt"
    fh = open(f,'r')
    for i, line in enumerate(fh):
        if i == 0: ## skip header.
            continue
        line = line.strip()
        ## capture all words in double-quotes.
        words = re.findall(r'"(.*?)"', line)
        
        if len(words) < 2: ## no blattner? then skip.
            continue
        ## blattner is second to last word in the line.
        my_blattner = words[-2]
        if not my_blattner.startswith('b'): ## skip if clearly not a blattner.
            continue
        for w in words:
            gene_to_blattner[w] = my_blattner

    print('gene,blattner')
    for k,v in gene_to_blattner.items():
        print(','.join([k,v]))

main()
