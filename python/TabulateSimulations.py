import json, sys, csv, argparse, operator, re

def echo (calls, region_names, total, latex):
    if latex:
        for k in range (3):
            calls[k][0] += calls[k][1]
            calls[k][2] += calls[k][3]

        calls = [["%.2f" % (k/total*100) for k in i] for i in calls]
        
        print (' & '.join (["%s (%s)" % (calls[0][0],calls[0][1]), "%s (%s)" % (calls[1][0],calls[1][1]), "%s (%s)" % (calls[2][0],calls[2][1]),
                            "%s (%s)" % (calls[0][2],calls[0][3]), "%s (%s)" % (calls[1][2],calls[1][3]), "%s (%s)" % (calls[2][2],calls[2][3]),
                            calls[0][4], calls[1][4], calls[2][4]]))
    else:
        for k, r in enumerate (calls):
            print ("Results for %s" % region_names[k])
            print ("\tCorrect             %d (%.3g %%)" % (r[0],100*r[0]/total))
            print ("\tAlternative         %d (%.3g %%)" % (r[1],100*r[1]/total))
            print ("\tMismatch (allele)   %d (%.3g %%)" % (r[2],100*r[2]/total))
            print ("\tMismatch (gene)     %d (%.3g %%)" % (r[3],100*r[3]/total))
            print ("\tNo result           %d (%.3g %%)" % (r[4],100*r[4]/total))

allele_id = re.compile ('\*[0-9]+')

parser = argparse.ArgumentParser(
    description='Parse IgSCUEAL'
)

parser.add_argument (
    '-i', '--igscueal',
    type=argparse.FileType('r'),
    help='IgSCUEAL main TSV file',
    required = True,
)

parser.add_argument (
    '-r', '--rearrangement',
    type=argparse.FileType('r'),
    help='IgSCUEAL rearrangement TSV file',
    required = True,
)

parser.add_argument (
    '-w', '--weak',
    action = 'store_true',
    help='Report weak matches',
    required = False,
    default = False
)

parser.add_argument (
    '-e', '--errors',
    action = 'store_true',
    help='Report mis-matches',
    required = False,
    default = False
)

parser.add_argument (
    '-l', '--latex',
    action = 'store_true',
    help='Output results in a LaTeX friendly format',
    required = False,
    default = False
)

    

args = None
retcode = -1
args = parser.parse_args()

main_tsv = csv.reader (args.igscueal, delimiter = '\t')
r_tsv    = csv.reader (args.rearrangement, delimiter = '\t')

next (main_tsv)
next (r_tsv)

total = 0

call_cloud = {}
for line in r_tsv:
    if len (line) == 3:
        if line[0] not in call_cloud:
            call_cloud [line[0]] = []
        call_cloud [line[0]].append (line[1].split(','))

#print (call_cloud)
calls = [[0,0,0,0,0] for k in range (3)]

#    calls [region] = [0,0,0] # correct-best, correct-in-set, incorrect

region_names = ["V","D","J"]
equvalence_map = [{},{},{}]
equvalence_map[0] = {"IGHV3-66*04" : "IGHV3-66*01",
                     "IGHV4-28*03" : "IGHV4-28*01",
                     "IGHV1-69*01" : "IGHV1-69D*01",
                     "IGHV3-30*04" : "IGHV3-30-3*03"
                     }
                  

for line in main_tsv:
    correct  = line[1].split ('|')
    if len (correct) > 3:
        name = correct.pop (0)
    else:
        name = correct[0]
    try:
        for k in range (3):
            if correct[k] in equvalence_map[k]:
                correct[k] = equvalence_map[k][correct[k]]
                
        #correct[0] = correct[0].replace ("/OR","-OR")
                
        inferred = line[2].split (',')
        d_index = 19 # if float (line[22]) > float (line[26]) else 23
        if (len (inferred) > 1):
            inferred.insert (1,line[d_index])
        
        
        for region,assignment in enumerate (correct):
             if region >= 3: continue
             index = assignment.find (inferred[region])
             if index >= 0 and index + len (inferred[region]) == len (assignment):
                calls[region][0] += 1
             else:
                partial_match = False
                if line[1] in call_cloud:
                    if region != 1: #not D
                        cloud_index = 0 if region == 0 else 1
                        for other_calls in call_cloud[line[1]]:
                            if len (other_calls) > cloud_index:
                                if  assignment.find (other_calls[cloud_index]) >= 0:
                                    partial_match = True
                                    if args.weak:
                                        print ("%s: WEAK match in region %s: %s -> %s" % (name, region_names[region], correct[region], inferred[region]))
                                    break
                                
                    else:
                        alternatives = inferred[region].split ('|')
                        if len (alternatives) > 1:
                            for da in alternatives:
                                if assignment.find (da) >= 0:
                                    partial_match = True
                            
                if partial_match:
                    calls[region][1] += 1
                    continue
               
                
                # try stripping *0X and match genes
                
                no_allele_correct = allele_id.sub ("",assignment)
                no_allele_assigned = allele_id.sub ("",inferred[region])
                index = no_allele_correct.find (no_allele_assigned)
                if index >= 0 and index + len (no_allele_assigned) == len (no_allele_correct):
                    calls[region][2] += 1
                else:
                    calls[region][3] += 1
                if args.errors and region != 1:
                    print ("%s: ERROR in region %s: %s -> %s" % (name, region_names[region], correct[region], inferred[region]))
        
    except:
        #raise
        pass
        
    total += 1

echo (calls, region_names, total, args.latex)
       



