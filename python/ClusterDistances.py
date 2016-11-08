import csv, argparse, sys, re, json, math

clonotype_weight = re.compile ('_([0-9]+)$')

####################################

def numeric_range (value, lower = 0, upper = 1) :
    value = float (value)
    if (value < lower or value > upper):
        raise argparse.ArgumentTypeError('Expected a value in [%f, %f]. Was given %f' % (lower, upper, value))
    return value

####################################

def connected_component (seed, links, extra_links):
    while True:
        current_set = set (seed)
        for k in seed:
            link_to = [l for l in links if l[0] == k or l[1] == k]
            link_to.extend ([l for l in extra_links if l[0] == k or l[1] == k])
            current_set.update ([l[0] if l[1] == k else l[1] for l in link_to])
        if len (current_set) == len (seed):
            break
            
        seed = current_set
    return seed
    

####################################

def float01 (value):
    return numeric_range (value)
    
####################################

def make_graphviz_rep (clusters, full_clusters, FH, special = None, filter = None, id_to_sequence_map = None):
    print ('graph G { outputorder = egdesfirst;', file = FH)
    for idx, c in enumerate (clusters['clusters']):
        if (filter and idx not in filter): continue
        extra = ', style = "filled", fillcolor = gray'
        if special is not None:
            for k in special:
                if k [0] is not None:
                    member_nodes = []
                    for n in k[0]:
                        if n in full_clusters[idx]:
                           member_nodes.append (n)
                    if len (member_nodes):
                         extra = ', style = "filled", fillcolor = %s, label = "%s", fontcolor = "white", labelfontsize = 20' % (k[1]," ,".join(sorted([id_to_sequence_map[n] for n in member_nodes])))
            
        print ('%d [label = "", width = %g, height = %g%s];' % (idx, math.sqrt (c), math.sqrt (c), extra), file = FH)
    
    for link in clusters['links']:
        if (filter and link[0] not in filter and link[1] not in filter): continue
        print ("%d -- %d [penwidth = %g];" % (link[0], link[1], 1 + math.log(link[2],10)), file = FH)

    for link in clusters['external']:
        if (filter and link[0] not in filter and link[1] not in filter or link[0] == link[1]): continue
        print ('%d -- %d [penwidth = %d, color = "red"];' % (link[0], link[1], max(1,link[2])), file = FH)
        
    print ("}", file = FH)
    
    
####################################

def process_a_node (node_id,sequence_id_map,patterns,id_to_sequence_map,edges,nuc_distances):
    if node_id not in sequence_id_map:
        for idx, p in enumerate(patterns):
            if p.search(node_id):
                match_sets[idx].add (len(sequence_id_map))
                       
        sequence_id_map[node_id] = len (sequence_id_map)
        id_to_sequence_map.append (node_id)
        edges.append ([])
        nuc_distances.append ([])
    return sequence_id_map[node_id]
    
####################################

####################################

argument_parser = argparse.ArgumentParser (description='Make a representation of Ig reads.')
argument_parser.add_argument('-i', '--input',   help = 'The comma-separated file produced by TN93. Must have exactly 3 columns', required = True, type=argparse.FileType('r'))
argument_parser.add_argument('-f', '--fasta',   help = 'The aligned file', required = True, type=argparse.FileType('r'))
argument_parser.add_argument('-c', '--cluster',   help = 'Maximal distance allowed to cluster sequences into a blob [default 0.02]', type=float01, default = 0.02)
argument_parser.add_argument('-l', '--link',   help = 'Maximal distance allowed to link clusters [default 0.05]', type=float01, default = 0.05)
argument_parser.add_argument('-o', '--output',   help = 'Write the resulting dot file here', type=argparse.FileType('w'), default = sys.stdout)
argument_parser.add_argument('-p', '--pattern',   help = ' sequence ids', nargs = '*', type = str)
argument_parser.add_argument('-a', '--additional',   help = ' Json files listing additional paths', nargs = '*', type = argparse.FileType('r'))

####################################

cli_args = argument_parser.parse_args()

ig_reader = csv.reader (cli_args.input, delimiter = ',')
column_headings = next (ig_reader)

if len (column_headings) != 3:
    raise BaseException ('The input file has an incorrect number of columns')
    
sequence_id_map    = {}
id_to_sequence_map = [];
edges              = []
nuc_distances      = []

idx = 0

patterns   = []
match_sets = []
colors     = []

if cli_args.pattern:
    for k in range(0,len(cli_args.pattern),2):
        patterns.append (re.compile(cli_args.pattern[k]))
        colors.append (cli_args.pattern[k+1])
        match_sets.append (set())
        
        

for line in ig_reader:
    dist = float (line[2])
    
    if (dist > cli_args.link): continue
    
    ids = [0,0]
    for k in (0,1):
        ids[k] = process_a_node(line[k],sequence_id_map,patterns,id_to_sequence_map,edges,nuc_distances)
        
    edges[ids[1]].append(ids[0])
    edges[ids[0]].append(ids[1])
    nuc_distances[ids[0]].append (dist)
    nuc_distances[ids[1]].append (dist)
    
    idx += 1
        
    if (idx % 100000 == 0):
        print (idx, file = sys.stderr)
            
id_to_cluster_assignment = {}
clusters                 = []

additional_links         = []
additional_links_weights = []

if cli_args.additional:
    for a_path_set in cli_args.additional:
        path_sets = json.load (a_path_set)
        for a_path in path_sets.items():
            last_node = None
            for node in a_path[1]:
                node_id = process_a_node (node,sequence_id_map,patterns,id_to_sequence_map,edges,nuc_distances)
                if last_node:
                    this_pair = sorted([last_node, node_id])
                    if this_pair not in additional_links:
                        additional_links.append (this_pair)
                        additional_links_weights.append(1)
                    else:
                        additional_links_weights[additional_links.index(this_pair)] += 1
                        
                last_node = node_id    
                
for k in range(len(additional_links)):
    additional_links[k].append (additional_links_weights[k])
 
for id in range (len (id_to_sequence_map)):

    if id not in id_to_cluster_assignment:
        considered_clusters = {}
    
        for idx, neighbor in enumerate(edges[id]):
            if neighbor < id and nuc_distances[id][idx] < cli_args.cluster: # could cluster with this node
                try:
                    cluster_id = id_to_cluster_assignment[neighbor]
                    if cluster_id not in considered_clusters:
                        considered_clusters [cluster_id] = 1
                        for node in clusters[cluster_id]:
                            if nuc_distances[id][edges[id].index (node)] > cli_args.cluster:
                                raise KeyError
                    
                    id_to_cluster_assignment[id] = cluster_id
                    clusters[cluster_id].append(id)
                    break
                    
                except (KeyError, ValueError):
                    pass
            
        if id not in id_to_cluster_assignment:    
            id_to_cluster_assignment [id] = len (clusters)
            clusters.append ([id])
            
        
link_clusters = {}

summary_dict = {'clusters' : [], 'links': [], 'external': [[id_to_cluster_assignment[links[k]] if k < 2 else links[k] for k in range (len(links))] for links in additional_links]}

#print (additional_links, "\n", summary_dict ['external'], file = sys.stderr)

for c in clusters:
    summary_dict ['clusters'].append (len(c))
    size = 0
    for seq in c:
        seq_name = id_to_sequence_map[seq]
        try:
            size += int(clonotype_weight.match (seq_name).group(1))
        except:
            size += 1

             
for index,eset in enumerate(edges):
    edges[index] = set (eset)

for c1 in range (len(clusters)):
    print ("%d/%d" % (c1, len(clusters)), file = sys.stderr)
    for c2 in range (c1+1, len(clusters)): 
        total_links = 0       
        for n1 in clusters[c1]:
            for n2 in clusters[c2]:
                try:
                    if n2 in edges[n1]:
                        total_links += 1
                        
                except (KeyError):
                    pass
                    
        if (total_links > 0):
            summary_dict ['links'].append ([c1,c2,total_links])
            #print ("%d -- %d (%d)" % (c1, c2, total_links), file = sys.stderr)
            


#print (summary_dict, file = sys.stderr)
        
make_graphviz_rep (summary_dict, clusters, cli_args.output, [(match_sets[k], colors[k]) for k in range (len(colors))],id_to_sequence_map)                                     
      
if len (match_sets):       
    seed = set ()
    for i,c in enumerate(clusters):
        for s in match_sets:
            if len(s.intersection (set (c))):
                seed.add (i)
                break
        
    #print (summary_dict ['external'], seed, file = sys.stderr)
    make_graphviz_rep (summary_dict, clusters, sys.stdout, [(match_sets[k], colors[k]) for k in range (len(colors))],connected_component (seed, summary_dict['links'], summary_dict['external']),id_to_sequence_map)                                                  
