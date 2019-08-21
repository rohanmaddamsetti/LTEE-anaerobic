#!/usr/bin python

'''
print_REL606_noncoding_seqs.py by Rohan Maddamsetti.

Usage: python print_REL606_noncoding_seqs.py > ../results/REL606_noncoding_seqs.fasta
'''

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

def main():

    genome_rec = SeqIO.read('../data/REL606.7.gbk','genbank')
    genome_dna = genome_rec.seq
    starting_len = len(genome_dna)
    genome_dna_list = [x for x in genome_rec.seq]
    ## remember python conventions for slices:
    ## start at 0, and end-index is NOT included
    ## (half-open interval).
    intergenic_start = 0
    for feat in genome_rec.features:
        if feat.type == 'source': ## don't mask the whole genome!
            continue
        start = int(feat.location.start)
        end = int(feat.location.end)

        ## mask the feature.
        for i in range(start,end):
            genome_dna_list[i] = 'N'

    all_intergenic_dna = ''.join([x for x in genome_dna_list if x!= 'N'])
    genome_dna = ''.join([x for x in genome_dna_list])

    ## print out intergenic regions.
    start = 0
    end = 0
    in_intergenic = True
    seq_buffer = []
    for i, c in enumerate(genome_dna_list):
        if in_intergenic:
            if c == 'N': ## flush buffer.
                in_intergenic = False
                end = i
                header = ''.join(['>REL606:',str(start),'-',str(end)])
                my_seq = ''.join(seq_buffer)
                
                print(header)
                print(my_seq)
                seq_buffer = []
            else:
                seq_buffer.append(c)
        else: ## in a gene or some other feature
            if c != 'N': ## then now in intergenic region.
                in_intergenic = True
                start = i
                assert seq_buffer == []
                seq_buffer.append(c)
    if in_intergenic: ## flush the buffer now the loop is done.
        assert(len(seq_buffer)>0)
        end = len(genome_dna)
        header = ''.join(['>REL606:',str(start),'-',str(end)])
        my_seq = ''.join(seq_buffer)
        print(header)
        print(my_seq)
  
                  
        ##intergenic_end = feat.location.start
        ## print Genbank coordinates (starts at 1), and end-index is included--
        ## therefore, add 1 to python start index but not to python end index.
        ##header = '>REL606:' + str(intergenic_start+1)+ '-' + str(intergenic_end)
        ##print(header)
        ##intergenic_dna = genome_dna[intergenic_start:intergenic_end]
        ##print(intergenic_dna)
        ##intergenic_start = feat.location.end
        
        
        ##mask = 'N'* len(feat)
        ##genome_dna = genome_dna[:int(feat.location.start)] + mask + genome_dna[int(feat.location.end):]
    assert len(genome_dna) == starting_len
main()
