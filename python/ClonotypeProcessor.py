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
import tempfile
import random

from  BioExt.scorematrices import ( BLOSUM62,PAM250)

from BioExt.uds import align_to_refseq
from BioExt.misc import gapless

################################################################################

def define_clonotypes (sequences, dist_thresh = 0.02, overlap = 100, file = None):

    seqs = sequences['seqs']
    equivalence = sequences['equivalence']
    id_to_seq = {}

    if 'mapping' in sequences:
    	reverse_id_map = {v:k for (k,v) in sequences['mapping'].items()}
    else:
    	reverse_id_map = {s.id : s.id for s in seqs}

    id_to_seq = {s.id : s for s in seqs}

    if len (seqs) == 1:
        members = [seqs[0].id]
        members.extend (equivalence [seqs[0].id])
        return   [{'size': 1 + len(equivalence [seqs[0].id]), 'members': members , 'centroid': seqs[0].format('fasta')}]

    with open ("/dev/null", "w") as devnull:
        if file is not None:
            tn93_process = subprocess.Popen(['/usr/local/bin/tn93-cluster', '-q', '-t', str (dist_thresh), '-a', 'resolve', '-m', 'json',  '-l',  str (overlap), file], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=devnull, universal_newlines = True)
        else:
            tn93_process = subprocess.Popen(['/usr/local/bin/tn93-cluster', '-q', '-t', str (dist_thresh), '-a', 'resolve', '-m', 'json',  '-l',  str (overlap)], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=devnull, universal_newlines = True)
        stdout, stderr = tn93_process.communicate (sequences['alignment'] if file is None else '')

        if stderr is not None:
            raise Exception (stderr)

        clusters = [c['members'] for c in json.loads (stdout)]

    return_clusters = []


    for c in clusters:
        lengths = [len(id_to_seq[reverse_id_map[k]]) for k in c]
        members = []
        for k in c:
            members.append (reverse_id_map[k])
            members.extend (equivalence[reverse_id_map[k]])
        return_clusters.append({'size': len (members), 'members': members, 'centroid': id_to_seq[reverse_id_map[c[lengths.index(max(lengths))]]].format('fasta')})

    return return_clusters

################################################################################

def run_bealign (ref, sequences):
    spooled_fasta = tempfile.NamedTemporaryFile(dir = os.getcwd(), mode = 'w', delete = False)

    used_ids = set ()
    digits = range (10)

    def random_id ():
        while True:
            test_id = ''.join([str (k) for k in random.sample(range(1000000), 5)])
            if test_id not in used_ids:
                used_ids.add (test_id)
                return test_id

    id_map = {ref.id : random_id()}
    print (">%s\n%s\n" % (id_map[ref.id], ref.seq), file = spooled_fasta)
    for s in sequences:
        id_map [s.id] = random_id()
        print (">%s\n%s\n" % (id_map [s.id], s.seq), file = spooled_fasta)

    bamfile = spooled_fasta.name + ".bam"
    msafile = spooled_fasta.name + ".msa"
    spooled_fasta.close()

    subprocess.call(['/opt/python-3.4.3/bin/bealign', '-K', '-m', 'BLOSUM62', spooled_fasta.name, bamfile])
    subprocess.call(['/opt/python-3.4.3/bin/bam2msa', bamfile, msafile])

    os.remove (spooled_fasta.name)
    os.remove (bamfile)
    os.remove (bamfile + ".bai")


    return msafile, id_map


################################################################################

def run_TN93 (sequences, dist_thresh = 0.02, overlap = 100, file = None):

    if len (sequences['seqs']) == 1:
       return ['ID1,ID2,Distance']

    with open ("/dev/null", "w") as devnull:
        if file is not None:
            tn93_process = subprocess.Popen(['/usr/local/bin/tn93', '-q', '-t', str (dist_thresh), '-a', 'resolve', '-f', 'csvn',  '-l',  str (overlap), '-d', "'~'", file], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=devnull, universal_newlines = True)
        else:
            tn93_process = subprocess.Popen(['/usr/local/bin/tn93', '-q', '-t', str (dist_thresh), '-a', 'resolve', '-f', 'csvn',  '-l',  str (overlap), '-d', "'~'"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=devnull, universal_newlines = True)
        stdout, stderr = tn93_process.communicate (sequences['alignment'] if file is None else '')
        if stderr is not None:
            raise Exception (stderr)
        ig_reader = csv.reader (stdout.split('\n'), delimiter = ',')
        column_headings = next (ig_reader)
        return [row for row in ig_reader]

    #if file is not None:
    #    os.remove (file)

################################################################################

def run_group_alignment (sequence_group):

    print ("%d sequences with matching JUNCTION regions" % (len (sequence_group) ))

    unique_sequences = {}
    seqrecords       = []
    id_eqivalence    = {}

    for seq_id in sequence_group:
        #print ("Step 1\n%s" % sequence_group[seq_id])
        massaged_string = sequence_group[seq_id].replace ('NNN','').replace ('---','').replace ('-','N')
        #print ("Step 2\n%s" % massaged_string)
        if len (massaged_string) % 3:
            massaged_string = massaged_string [:len (massaged_string) - len (massaged_string) % 3]
            #print ("Step 3\n%s" % massaged_string)


        #seqrecords.append(gapless(Bio.SeqRecord.SeqRecord (Bio.Seq.Seq(massaged_string), id = seq_id, name = seq_id, description = '')))

        if not massaged_string in unique_sequences:
            unique_sequences [massaged_string] = []

        unique_sequences [massaged_string].append (seq_id)

    stable_keys = sorted (unique_sequences.keys())

    for sequence in stable_keys:
        ids = unique_sequences[sequence]
        id_eqivalence [ids[0]] = ids[1:]
        seqrecords.append(gapless(Bio.SeqRecord.SeqRecord (Bio.Seq.Seq(sequence), id = ids[0], name = ids[0], description = '')))


    if len (seqrecords) == 1:
        refseq = seqrecords[0].format ('fasta')
        return {'ref': refseq, 'alignment': refseq, 'seqs': seqrecords, 'equivalence' : id_eqivalence}

    # find the longest sequence
    seq_lengths = [len(record.seq) for record in seqrecords]
    refseq_id = seq_lengths.index(max(seq_lengths))
    refseq = seqrecords.pop(refseq_id)


    #print (len (seqrecords))

    if len(refseq.seq) % 3: ## error condition (not divisible by 3)
        seqrecords = [s for s in seqrecords]
        print (">ref\n%s" % str(refseq.seq))
        print ('\n'.join ([">%s\n%s" % (str(k.id), str(k.seq)) for k in seqrecords]))

    sm = BLOSUM62.load()

    if len (sequence_group) >= 2048:
        # handle via an external call
        file, mapping = run_bealign (refseq, seqrecords)
        return {'ref': refseq.format ('fasta'), 'file': file, 'mapping': mapping, 'seqs': seqrecords, 'equivalence' : id_eqivalence}
    else:
        if len (sequence_group) < 32:
            os.environ['NCPU'] = '1'
        elif len (sequence_group) < 128:
            os.environ['NCPU'] = '4'
        elif len (sequence_group) < 512:
            os.environ['NCPU'] = '8'
        else:
            os.environ['NCPU'] = '16'



        #os.environ['NCPU'] = '1' if  len (sequence_group) < 32 else '-1'

        #print (os.environ.get ('NCPU', -1))

        msa, discarded = align_to_refseq(
            refseq,
            seqrecords,
            score_matrix=sm,
            do_codon=True,
            reverse_complement=False,
            #expected_identity=0.6,
            keep_insertions=False,
        )

        print ("Done aligning %d sequences with matching JUNCTION regions" % (len (sequence_group) ))

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

        return {'ref': refseq.format ('fasta'), 'alignment': all_lines, 'seqs': seqrecords, 'equivalence' : id_eqivalence}


def main ():
    argument_parser = argparse.ArgumentParser (description='Process clonotypes from a binned JSON.')
    argument_parser.add_argument('-i', '--input',    help = 'The input JSON file',  nargs = '?', required = True, type=argparse.FileType('r'), default = sys.stdin)
    argument_parser.add_argument('-o', '--output',   help = 'The output JSON file', nargs = '?', required = True, type=str)
    argument_parser.add_argument('-c', '--collapse',   help = 'Collapse reads that match the rearrangement and everything except N in the CDR3 region into a single clustrer', action = "store_true")

    cli_args = argument_parser.parse_args()

    if not os.path.isfile (cli_args.output):
        open (cli_args.output, 'w').close()

    json_input     =  json.load (cli_args.input)
    all_clonotypes = []
    with open (cli_args.output, 'r') as fh:
        try:
            all_clonotypes = json.load (fh)
            print ("[Loaded %d already processed results]" % len (all_clonotypes))
        except Exception as e:
            print ("Error loading cache: %s" % str (e), file = sys.stderr)

    unique_cdr3 = sorted(json_input.keys ())

    if cli_args.collapse:
        print ("[Collapsing clusters matching up to N's in CDR3]")
        collapsed_json = {}
        # first generate the list of all resolved keys, i.e. those that do not have an N in the CDR3 region

        ambig_map = str.maketrans ("ACGT", "ACGT", "N-")
        ambig_keys = set ()
        resolved_keys = {}

        for cluster_tag in unique_cdr3:
            rearrangement , cdr3 = cluster_tag.split ('|')
            if cdr3 != cdr3.translate (ambig_map):
                ambig_keys.add (cluster_tag)
            else:
                if rearrangement not in resolved_keys:
                    resolved_keys [rearrangement] = {}

                if len (cdr3) not in resolved_keys [rearrangement]:
                    resolved_keys [rearrangement] [len (cdr3)] = {}

                resolved_keys [rearrangement][len(cdr3)][cdr3] = [cluster_tag, ];

        if len (ambig_keys):

            def compare_with_ambigs (s1, s2):
                if len (s1) == len (s2):
                    for i in range (len (s1)):
                        if s1 [i] != s2 [i]:
                            if s1 [i] == 'N' or s1 [i] == '-' or s2 [i] == 'N' or s2 [i] == '-':
                                continue
                            return False
                    return True
                return False

            unmatched_ambigs = set ()
            for ambig_tag in ambig_keys:
                rearrangement , cdr3 = ambig_tag.split ('|')
                try:
                    candidate_set = resolved_keys [rearrangement][len(cdr3)]
                    try:
                        match = next (filter (lambda key : compare_with_ambigs (cdr3, key), resolved_keys [rearrangement][len(cdr3)].keys ()))
                        resolved_keys [rearrangement][len(cdr3)][match].append (ambig_tag)

                    except StopIteration:
                        unmatched_ambigs.add (ambig_tag)

                except KeyError:
                    unmatched_ambigs.add (ambig_tag)


            reaggregated = {}

            for r,rr in resolved_keys.items():
                for l,ll in rr.items ():
                    for c, cc in ll.items():
                        reaggregated [cc[0]] = {}
                        for tag in cc:
                            reaggregated [cc[0]].update (json_input [tag])

            for tag in unmatched_ambigs:
                reaggregated[tag] = json_input[tag]

            old_l = len (unique_cdr3)

            json_input = reaggregated
            unique_cdr3 = sorted(json_input.keys ())

            print ("[REDUCED FROM %d CLUSTERS TO %d CLUSTERS]" % (old_l, len (unique_cdr3)))


    for idx, cdr3 in enumerate(unique_cdr3):
        if idx < len (all_clonotypes):
            continue
        if len (json_input[cdr3]) > 2000:
            print ("••••••••••••• Large (%d) CDR3 cluster for %s" % (len (json_input[cdr3]), cdr3))
        seqs = run_group_alignment(json_input[cdr3])
        #distances = run_TN93 (seqs, file = seqs['file'] if 'file' in seqs else None)
        this_batch = define_clonotypes(seqs,file = seqs['file'] if 'file' in seqs else None)
        for k in this_batch:
            k['tag'] = cdr3

        all_clonotypes.append (this_batch)
        print ("%d --> %d (%d/%d)" % (len(json_input[cdr3]), len (this_batch), len (all_clonotypes), len (json_input)))

        if 'file' in seqs:
            os.unlink (seqs['file'])

        if (idx > 0 and idx % 1000 == 0) or len (json_input[cdr3]) >= 1000:
            with open (cli_args.output, 'w') as fh:
                print ("... Writing JSON to file" )
                json.dump (all_clonotypes, fh, indent=4)



    with open (cli_args.output, 'w') as fh:
        print ("... Writing JSON to file")
        json.dump (all_clonotypes, fh, indent=4)

main ()
