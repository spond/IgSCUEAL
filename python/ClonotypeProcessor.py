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

from  BioExt.scorematrices import ( BLOSUM62,PAM250)

from BioExt.uds import align_to_refseq
from BioExt.misc import gapless

################################################################################

def define_clonotypes (seqs, distances):


    if len (seqs) == 1:
        return   [{'size': 1, 'members': [seqs[0].id], 'centroid': seqs[0].format('fasta')}]

    edges              = {}
    nuc_distances      = {}

    for line in distances:
        if len(line) < 3: continue
        dist = float (line[2])

        ids = [0,0]
        for k in (0,1):
            ids[k] = int(line[k])
            if (ids[k] not in edges):
                edges[ids[k]] = []
                nuc_distances[ids[k]] = []


        edges[ids[1]].append(ids[0])
        edges[ids[0]].append(ids[1])
        nuc_distances[ids[0]].append (dist)
        nuc_distances[ids[1]].append (dist)



    id_to_cluster_assignment = {}
    clusters                 = []

    for id in range (len (seqs)):

        if (id not in edges):
            clusters.append([id])
            continue

        if id not in id_to_cluster_assignment:
            considered_clusters = {}

            for idx, neighbor in enumerate(edges[id]):
                if neighbor < id:
                    try:
                        cluster_id = id_to_cluster_assignment[neighbor]
                        if cluster_id not in considered_clusters:
                            considered_clusters [cluster_id] = 1
                            for node in clusters[cluster_id]:
                                nuc_distances[id][edges[id].index (node)]

                        id_to_cluster_assignment[id] = cluster_id
                        clusters[cluster_id].append(id)
                        break

                    except (KeyError, ValueError):
                        pass

            if id not in id_to_cluster_assignment:
                id_to_cluster_assignment [id] = len (clusters)
                clusters.append ([id])

    return_clusters = []
    for c in clusters:
        lengths = [len(seqs[k].seq) for k in c]
        return_clusters.append({'size': len (c), 'members': [seqs[k].id for k in c], 'centroid': seqs[lengths.index(max(lengths))].format('fasta')})

    return return_clusters

################################################################################

def make_file (string):
    try:
        fh = open (string, 'r+')
    except:
        fh = open (string, 'w+')

    return fh

################################################################################

def run_TN93 (sequences, dist_thresh = 0.02, overlap = 100):

    if len (sequences['seqs']) == 1:
       return ['ID1,ID2,Distance']

    with open ("/dev/null", "w") as devnull:
        tn93_process = subprocess.Popen(['/usr/local/bin/tn93', '-q', '-t', str (dist_thresh), '-a', 'resolve', '-f', 'csvn',  '-l',  str (overlap)], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=devnull, universal_newlines = True)
        stdout, stderr = tn93_process.communicate (sequences['alignment'])
        if stderr is not None:
            raise Exception (stderr)
        ig_reader = csv.reader (stdout.split('\n'), delimiter = ',')
        column_headings = next (ig_reader)
        return [row for row in ig_reader]

################################################################################

def run_group_alignment (sequence_group):

    print ("%d sequences with matching JUNCTION regions" % (len (sequence_group)  - 1))
    seqrecords = []

    for seq_id in sequence_group:
        #print ("Step 1\n%s" % sequence_group[seq_id])
        massaged_string = sequence_group[seq_id].replace ('NNN','').replace ('---','').replace ('-','N')
        #print ("Step 2\n%s" % massaged_string)
        if len (massaged_string) % 3:
            massaged_string = massaged_string [:len (massaged_string) - len (massaged_string) % 3]
            #print ("Step 3\n%s" % massaged_string)

        seqrecords.append(gapless(Bio.SeqRecord.SeqRecord (Bio.Seq.Seq(massaged_string), id = seq_id, name = seq_id, description = '')))

    if len (seqrecords) == 1:
        refseq = seqrecords[0].format ('fasta')
        return {'ref': refseq, 'alignment': refseq, 'seqs': seqrecords}

    # find the longest sequence
    seq_lengths = [len(record.seq) for record in seqrecords]
    refseq_id = seq_lengths.index(max(seq_lengths))
    refseq = seqrecords.pop(refseq_id)


    #print (len (seqrecords))

    if len(refseq.seq) % 3:
        seqrecords = [s for s in seqrecords]
        print (">ref\n%s" % str(refseq.seq))
        print ('\n'.join ([">%s\n%s" % (str(k.id), str(k.seq)) for k in seqrecords]))

    sm = BLOSUM62.load()

    msa, discarded = align_to_refseq(
        refseq,
        seqrecords,
        score_matrix=sm,
        do_codon=True,
        reverse_complement=False,
        #expected_identity=0.6,
        keep_insertions=False,
    )

    if len (discarded):
        print (">ref\n%s" % str(refseq.seq))
        print ('\n'.join ([">%s\n%s" % (str(k.id), str(k.seq)) for k in seqrecords]))
        print (discarded)
        raise Exception ("Non-empty discarded")
        sys.exit (1)

    string_buffer = io.StringIO ()
    Bio.SeqIO.write (msa, string_buffer, "fasta")
    all_lines = string_buffer.getvalue()
    string_buffer.close()
    return {'ref': refseq.format ('fasta'), 'alignment': all_lines, 'seqs': seqrecords}


def main ():
    argument_parser = argparse.ArgumentParser (description='Process clonotypes from a binned JSON.')
    argument_parser.add_argument('-i', '--input',    help = 'The input JSON file',  nargs = '?', required = True, type=argparse.FileType('r'), default = sys.stdin)
    argument_parser.add_argument('-o', '--output',   help = 'The output JSON file', nargs = '?', required = True, type=make_file)

    cli_args = argument_parser.parse_args()

    json_input     =  json.load (cli_args.input)
    all_clonotypes = []
    try:
        cli_args.output.seek(0)
        all_clonotypes = json.load (cli_args.output)
        print ("[Loaded %d already processed results]" % len (all_clonotypes))
    except:
        pass

    unique_cdr3 = sorted(json_input.keys ())

    for idx, cdr3 in enumerate(unique_cdr3):
        if idx < len (all_clonotypes):
            continue
        #if len (json_input[cdr3]) > 5000:
        #    continue
        seqs = run_group_alignment(json_input[cdr3])
        distances = run_TN93 (seqs)
        this_batch = define_clonotypes(seqs['seqs'],distances)
        all_clonotypes.append (this_batch)
        if idx and idx % 100 == 0 or len (json_input[cdr3]) >= 500:
            cli_args.output.seek(0)
            json.dump (all_clonotypes, cli_args.output, indent=4)


        print ("%d --> %d (%d/%d)" % (len (seqs['seqs']), len (this_batch), len (all_clonotypes), len (json_input)))

    cli_args.output.seek(0)
    json.dump (all_clonotypes, cli_args.output, indent=4)

main ()
