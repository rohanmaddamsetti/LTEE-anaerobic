#!/usr/bin/env python

'''  align_K12_arcA_to_REL606.py by Rohan Maddamsetti.

1) Automatically make blast dbs.

2) Find and print REL606 coordinates for each arcA binding region reported by
   Park et al. 2013.

Usage: python align_K12_arcA_to_REL606.py > ../results/K12_arcA_motifs_in_REL606.csv

'''

from os.path import join, basename, exists
from os import listdir, makedirs, chdir, getcwd
import sys
from pprint import pprint
from copy import deepcopy
import subprocess
import argparse

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Blast.Applications import NcbiblastnCommandline as BLAST
from Bio.Blast import NCBIXML
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import MafftCommandline as MAFFT
from io import StringIO
from Bio import AlignIO

def setup_blastdb():

    projdir = "/Users/Rohandinho/BoxSync/active-projects/LTEE-anaerobic/"
    resultsdir = join(projdir,"results/")

    REL606fasta = join(projdir,"data/REL606.7.fasta")
    title = 'REL606'
    db = join(resultsdir,title)
    
    chdir(resultsdir)
    ## if blastdbs haven't been made, make blast dbs.
    if not exists("REL606.fasta.nhr"):
        REL606_args = ['makeblastdb','-in', REL606fasta,
                          '-dbtype','nucl', '-out', db, '-title', title]
        subprocess.run(REL606_args)

def printREL606_coordinates():
    
    ## The BLAST queries come from the Park et al. 2013 paper.
    projdir = "/Users/Rohandinho/BoxSync/active-projects/LTEE-anaerobic/"
    resultsdir = join(projdir,"results/")
    query_file = join(projdir, 'results/Park2013_arcA_motifs.fasta')

    ## parameters need to be optimized for short query sequences.
    ## set task to blastn-short.
    
    the_db = join(projdir,"results/REL606")
    blast_out = join(projdir,"results/arcA_blast_results.xml")
    this_blast = BLAST(cmd='blastn',
                       task='blastn-short',
                       query=query_file,
                       db=the_db,
                       out=blast_out,
                       outfmt=5,
                       max_target_seqs=1)
    ## run BLASTN from the results dir.
    chdir(resultsdir)
    this_blast()
    ## Now parse the output to finish populating the fields in alignment_dict.
    result_handle = open(blast_out,'r')
    blast_records = NCBIXML.parse(result_handle)
    ## print a header.
    print('K12_arcA_motif,REL606_start,REL606_end')
    for rec in blast_records:
        ## if no hits, write na.
        if len(rec.alignments) == 0:
            print(','.join([rec.query,'NA','NA']))
            continue
        ## we only want the first hit for each query.
        tophit = rec.alignments[0].hsps[0]
        tophitname = rec.alignments[0].title
        query_start = tophit.sbjct_start
        ## watch the indexing! have to subtract one.
        query_end = tophit.sbjct_start + tophit.align_length - 1
        ## print the K12 query, and start and end coordinates in REL606.
        print(','.join([rec.query,str(query_start),str(query_end)]))

def main():

    projdir = "/Users/Rohandinho/BoxSync/active-projects/LTEE-anaerobic/"
    ##setup_blastdb()
    printREL606_coordinates()

main()
