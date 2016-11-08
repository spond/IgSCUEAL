from __future__ import division, print_function


from Bio import SeqIO, AlignIO
import Bio.SeqRecord
import Bio.Alphabet
import Bio.Seq
import io
import subprocess
import os
import argparse
import sys
import json
import csv

from BioExt import ( BLOSUM62,PAM250)

from hy454 import align_to_refseq, to_positional

        
################################################################################

def run_group_alignment (seqrecords, refseqs):
    
    
    sm = BLOSUM62.load()

    seq_scores = {'names':{}}
    
    sequence_names = {}
    for seq in seqrecords:
        sequence_names[seq.id] = len (sequence_names)
        
    seqrecords.insert (0, 0)
    
    
    for type in refseqs:
        seq_scores [type] = {}
        seq_scores['names'][type] = []
    
        for seq in refseqs[type]:
            seqrecords[0] = seq
            seq_scores['names'][type].append (seq.id)
            
            msa, discarded = align_to_refseq(
                seq,
                seqrecords,
                score_matrix=sm,
                codon=True,
                expected_identity=0.4,
                keep_insertions=True,
                quiet=True
            )
            
            if (discarded):
                print (discarded)    
                
            max_score = msa[0].annotations['_pbpscore']
            max_len   = msa[0].annotations['_nbpidentical']
            seq_scores[type][seq.id] = [None for idx in range (len (sequence_names))]
            msa.pop (0)
            for successful in msa:
                idx = sequence_names [successful.id]
                seq_scores[type][seq.id][idx] = (successful.annotations['_pbpscore']/max_score+successful.annotations['_nbpidentical']/max_len)*0.5
            
        
    return seq_scores


def main ():
    argument_parser = argparse.ArgumentParser (description='Process clonotypes from a binned JSON.')
    argument_parser.add_argument('-i', '--input',       help = 'The input FASTA file',   required = True, type=argparse.FileType('r'))
    argument_parser.add_argument('-b', '--bnab',        help = 'The FASTA file with target bNab sequences',  required = True, type=argparse.FileType('r'))
    argument_parser.add_argument('-g', '--germline',    help = 'The FASTA file with reference germline sequences',  required = True, type=argparse.FileType('r'))
    argument_parser.add_argument('-o', '--output',      help = 'The output JSON file',  required = True, type=argparse.FileType('w'))
    
    cli_args = argument_parser.parse_args()
    
    query_seqs = []
    for seq in Bio.SeqIO.parse (cli_args.input, "fasta"):
        query_seqs.append(seq)
        
    ref_seqs = {'bnab':[],'germline': []}
    for seq in Bio.SeqIO.parse (cli_args.bnab, "fasta"):
        ref_seqs['bnab'].append(seq)
        
    for seq in Bio.SeqIO.parse (cli_args.germline, "fasta"):
        ref_seqs['germline'].append(seq)

    json.dump (run_group_alignment (query_seqs,  ref_seqs), cli_args.output)

main ()
