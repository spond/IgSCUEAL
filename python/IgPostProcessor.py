import csv, argparse, sys, re, operator, json, random


########### START GLOBALS #############

required_headers = ['Index', 'Name', 'Best Rearrangement', 'Support', 'Sequence']
column_mapper    = {}
cdr3_id          = 0
is_light_chain   = False
has_constant     = False

V_processor      = re.compile ('^V([^\-\*\,]+)([^\*\,]+)?([^\,]+)?')
J_processor      = re.compile ('\,J([^\*\,]+)([^\,]+)?')
D_processor      = re.compile ('^D([^\-\*\,]+)([^\*\,]+)?([^\,]+)?')
CH_processor     = re.compile ('\,C([^\*\,]+)([^\,]+)?')

########### END GLOBALS #############

class matchNStr (str):

    def __new__ (self, string):
        return super(matchNStr, self).__new__(self, string)

    def __hash__ (self):
        return len(self).__hash__()

    def __lt__ (self,other):
        if len (self) != len (other):
            return super(matchNStr,self).__lt__(other)

        for char in range (len (self)):
            if self [char] == 'N' or other [char] == 'N':
                continue
            if self[char] < other[char]:
                return True
            if self[char] > other[char]:
                return False

        return False

    def __eq__  (self,other):
        if len (self) != len (other):
            return False

        for char in range (len (self)):
            if self [char] == 'N' or other [char] == 'N':
                continue
            if self[char] != other[char]:
                return False

        return True


'''class filter_csv (csv.reader):

    def __init__ (line_validator = None, *args):
        super.__init__ (args)
        self.validator = line_validator
    def next ():
        if self.line_validator is None:
            return super.next()
        while True:
            next_line = super.next()
            if self.line_validator (next_line):
                return next_line

        return None'''

####################################

class regExpClass:
    regexp    = None
    column_id = 0
    next_filter = None

    def __init__ (self, double, column_mapper, nf = None):
        self.column_id = column_mapper.index(double[0])
        self.regexp    = re.compile (double[1])
        self.next_filter = nf

    def check_line (self,line):
        if self.regexp.search (line[self.column_id]) is not None:
            if (self.next_filter is not None): return self.next_filter.check_line (line)
            return True
        return False


####################################
class idFilterClass:
    filter_set    = None
    column_id = 0
    next_filter = None

    def __init__ (self, id_tag, the_set, column_mapper, nf = None):
        self.column_id = column_mapper.index(id_tag)
        self.filter_set = the_set
        self.next_filter = nf

    def check_line (self,line):
        if line[self.column_id] in self.filter_set:
            if (self.next_filter is not None): return self.next_filter.check_line (line)
            return True
        return False


####################################

class colLengthClass:
    length        = None
    column_id     = 0
    comp_function = None

    def __init__ (self, triple, column_mapper, nf = None):
        self.column_id = column_mapper.index(triple[0])
        if triple[1] == 'l':
            self.comp_function = operator.lt
        elif triple[1] == 'g':
            self.comp_function = operator.gt
        else:
            self.comp_function = operator.eq

        self.length = int (triple[2])
        self.next_filter = nf

    def check_line (self,line):
        if len (line) > self.column_id and self.comp_function (len(line[self.column_id]), self.length):
            if (self.next_filter is not None): return self.next_filter.check_line (line)
            return True
        return False

####################################

def add_random_suffix (name, sep = '_', length = 10):
    return name + sep + "".join (random.sample ("0123456789ABDCDE", length))

####################################

def make_sequence_id (line):
    return "%s:%s" % (line[column_mapper['Name']], line[column_mapper['Index']])

####################################

def reduceRearrangementFromString (rearrangement):
    v  = ['','','']
    j  = ['','']
    ch = ['','']

    v_match  = V_processor.match (rearrangement)
    for i in range (3):
        v[i] = v_match.group (i+1)
        if v[i] is None:
            v[i] = '?'

    v[0] = 'V' + v[0]

    j_match = J_processor.search (rearrangement)
    #print (rearrangement, j_match)
    if j_match is not None and j_match.group(1) is not None:
        for i in range (2):
            j[i] = j_match.group (i+1)
            if j[i] is None:
                j[i] = '?'
        j[0] = 'J' + j[0]

    ch_match = CH_processor.search (rearrangement)
    #print (rearrangement, j_match)
    if ch_match is not None and ch_match.group(1) is not None:
        for i in range (2):
            ch[i] = ch_match.group (i+1)
            if ch[i] is None:
                ch[i] = '?'
        ch[0] = 'C' + ch[0]


    return v, j, ch


def reduceRearrangementToV (line):
    v, j, ch = reduceRearrangementFromString (line[column_mapper['Best Rearrangement']])
    if len (v) > 1:
        return v[0] + v[1]

    return v[0]


def reduceRearrangement (line):
    v, j, ch = reduceRearrangementFromString (line[column_mapper['Best Rearrangement']])

    if not is_light_chain:
        d_alleles = line[column_mapper['D_ALLELE']].split ('|')
        if (d_alleles == '0'):
            d_alleles = ''
        else:
            if cli_args.inverse_d:
                if line [column_mapper['D/INV_SCORE']] > line [column_mapper['D_SCORE']]:
                    d_alleles = [ k + "/I" for k in line[column_mapper['D/INV_ALLELE']].split ('|')]

        assignment = []
        for d in d_alleles:
            d_match = D_processor.search (d)
            if d_match:
                assignment.append ([v, ['D' + d_match.group (1), d_match.group (2) if d_match.group (2) else '?', d_match.group (3) if d_match.group (3) else '?'],j, ch])
            else:
               assignment.append ([v, ['D?', '?','?'],j, ch])


    else:
        assignment = ','.join (['|'.join (v),'|'.join (j)])

    return [['|'.join([','.join (b) for b in rearr]), 1/len (assignment)] for rearr in assignment]

####################################

def cdr3 (line):
    return len (line[cdr3_id])

####################################

def increment_dict_key (dict, key, value):
    if not key in dict:
        dict [key] = value
    else:
        dict [key] += value

####################################

def compute_diversity (array):
    group_weights = 0
    total = 0
    for v in array:
        group_weights += v * (v - 1)
        total += v

    return 1 - group_weights/total/(total-1)

####################################

def outputUsage (line_filter, support, support_col, mappingFunction):
    usage_report = {}

    if type (mappingFunction) == int:
        column_index = mappingFunction
        mappingFunction = lambda line : line [column_index]

    for line in ig_reader:
        if len ([k for k in line if len (k) >0]) > 10: # not a failed alignment
            if float (line[support_col]) >= support:
                if line_filter is not None:
                    if not line_filter.check_line (line):
                         continue

                try:
                    tag = mappingFunction (line)
                except:
                    continue

                increment_dict_key (usage_report, tag, 1.)

    to_display = sorted ([[key, value] for key, value in usage_report.items()], key = lambda x : x[1])
    return '\n'.join (['\t'.join ([str(value) for value in row]) for row in to_display]) + '\nGini-Simpson diversity %g' % compute_diversity ([k[1] for k in to_display])

####################################

def outputSummary (line_filter, support, support_col, rearrangement_col, mappingFunction = None, histogram = None, read_alternatives = None, coverages = None, cluster_sizes = None):


    binByRearrangement = {}
    total        = 0
    pass_support = 0
    mapped       = 0
    filtered_in  = 0
    alternatives = {}
    coverage_mapping = {}
    read_lengths  = {'Low support' : {},
                     'Passed filter' : {},
                     'Failed filter' : {} }



    read_column_id = column_mapper ["Sequence"]

    def add_to_counts (line, tag):
        length = len(line[read_column_id])
        if length not in read_lengths [tag]:
            read_lengths[tag][length] = 0
        read_lengths[tag][length] += 1

    if coverages is not None:
        coverages = coverages.split (',')
        for key in coverages:
            coverage_mapping [key] = 0

    if read_alternatives:
        alt_reader = csv.reader (read_alternatives, delimiter = '\t')
        next (alt_reader)
        for l in alt_reader:
            read_name = l[0]

            if not read_name in alternatives:
                alternatives[read_name] = []

            try:
                alternatives[read_name].append ([[','.join (k) for k in reduceRearrangementFromString (l[1])], float (l[2])])
            except:
                continue

    for line in ig_reader:
        total += 1
        if len ([k for k in line if len (k) >0]) > 10: # not a failed alignment
            mapped += 1
            if float (line[support_col]) >= support:
                pass_support += 1

                if line_filter is not None:
                    if not line_filter.check_line (line):
                        add_to_counts (line, 'Failed filter')
                        continue

                try:
                    if mappingFunction:
                        alleles = mappingFunction (line)
                    else:
                        alleles = [[line[rearrangement_col],1]]
                except:
                    add_to_counts (line, 'Failed filter')
                    continue

                filtered_in += 1

                if read_alternatives and line[column_mapper['Name']] in alternatives:
                    new_alleles = []
                    for d in alleles:
                        d_a = d[0].split ('|')[1]
                        for vj in alternatives[line[column_mapper['Name']]]:
                            wt = d[1] * vj[1]
                            new_alleles.append (["|".join ([vj[0][0], d_a, vj[0][1], vj[0][2]]), wt])
                    rescale = sum ([k[1] for k in new_alleles])
                    new_alleles = [[k[0], k[1]/rescale] for k in new_alleles]

                    alleles = new_alleles

                add_to_counts (line, 'Passed filter')


                if coverages is not None:
                    for key in coverages:
                        if len (line[column_mapper[key]]):
                            coverage_mapping[key] += 1


                if (histogram):
                    prop   = histogram (line)

                    for a in alleles:
                        if a[0] not in binByRearrangement:
                            binByRearrangement [a[0]] = {'Count' : a[1]}
                        else:
                            binByRearrangement [a[0]]['Count'] += a[1]

                        if prop not in binByRearrangement [a[0]]:
                            binByRearrangement [a[0]][prop] = a[1]
                        else:
                            binByRearrangement [a[0]][prop] += a[1]

                else:
                    for a in alleles:
                        if a[0] not in binByRearrangement:
                            binByRearrangement [a[0]] = a[1]
                        else:
                            binByRearrangement [a[0]] += a[1]
            else:
                add_to_counts (line, 'Low support')


    binByRearrangement ['Summary']     = [['Total', total], ['Mapped', mapped], ['Passed support', pass_support], ['Passed filter', filtered_in], ['Support', support], ['Read Lengths', read_lengths]]
    if coverages is not None:
        binByRearrangement ['Summary'].append ( ['Coverage', coverage_mapping])
    if cluster_sizes is not None:
        binByRearrangement ['Summary'].append ( ['Cluster sizes', cluster_sizes])


    return json.dumps (binByRearrangement, indent = 2)

####################################

def ouputLengths (line_filter, support, support_col, binning_column):
    binByLength = {}
    for line in ig_reader:
        if len ([k for k in line if len (k) >0]) > 10:
            if float (line[support_col]) >= support:
                if line_filter is not None:
                    if not line_filter.check_line (line): continue

                my_len = len(line[binning_column])
                if my_len not in binByLength:
                    binByLength [my_len] = 1
                else:
                    binByLength [my_len] += 1

    return json.dumps (binByLength)

####################################

def ouputSequences (line_filter, support, support_col, name_col, sequence_col):

    for line in ig_reader:
        if len ([k for k in line if len (k) >0]) > 10:
            if float (line[support_col]) >= support:
                if line_filter is not None:
                    if not line_filter.check_line (line): continue
                if type(name_col) == int:
                    name = line[name_col]
                else:
                    name = name_col(line)

                print (">%s\n%s\n" % (name, line[sequence_col]))


####################################

def outputJSON (line_filter, support, support_col, group_by = None):

    if group_by is None:
       output = []
       for line in ig_reader:
            if len ([k for k in line if len (k) >0]) > 10:
                if float (line[support_col]) >= support:
                    if line_filter is not None:
                        if not line_filter.check_line (line): continue
                    output.append (line)
    else:
        output = {}
        for line in ig_reader:
            if len ([k for k in line if len (k) >0]) > 10:
                if float (line[support_col]) >= support:
                    if line_filter is not None:
                        try:
                            if not line_filter.check_line (line): continue
                        except IndexError:
                            continue

                    #bin_by = matchNStr(line[group_by[1]])
                    bin_by = "|".join([line[k] for k in group_by[1]])
                    #print (group_by, file = sys.stderr)
                    if bin_by not in output:
                        output [bin_by] = {}

                    name_col = group_by[0]
                    if type(name_col) == int:
                        name = line[name_col]
                    else:
                        name = name_col(line)


                    output[bin_by][name] = line[group_by[2]]

    return json.dumps (output, sort_keys=True, indent=4)


####################################

def numeric_range (value, lower = 0, upper = 1) :
    value = float (value)
    if (value < lower or value > upper):
        raise argparse.ArgumentTypeError('Expected a value in [%f, %f]. Was given %f' % (lower, upper, value))
    return value

####################################

def float01 (value):
    return numeric_range (value)

####################################

argument_parser = argparse.ArgumentParser (description='Post-process IgSCUEAL output files.')
argument_parser.add_argument('-i', '--input',        help = 'The tab-separated file produced by IgSCUEAL. Must have at least 10 columns', nargs = '?', required = True, type=argparse.FileType('r'), default = sys.stdin)
argument_parser.add_argument('-s', '--support',      help = 'Minimum assignment support required to process a read (default is 0.9)', type=float01, default = 0.9)
argument_parser.add_argument('-p', '--pattern',      help = 'A regular expression pattern for pulling out reads that match a particular pattern in a particular column, e.g. "Rearrangement V3-21,J3"', type=str, nargs = 2, action = 'append')
argument_parser.add_argument('-c', '--column',       help = 'A column [l,g,e] length argument to select lines based on the length of the string in a given column, e.g. JUNCTION_AA g 10 (all lines with JUNCTIONS_AA with more than 10 chars)', nargs = 3, action = 'append')
argument_parser.add_argument('-x', '--cluster',      help = 'Only keep sequences which represent unique clonotypes', required = False, type=argparse.FileType('r'))
argument_parser.add_argument('-o', '--output',       help = 'Output mode (default fasta): either summarize filtered reads by column value (supply the column name), print FASTA for the reads matching selection criteria (fasta), or output JSON (json) for matching reads. '+
                                                          'Finally specify "json output bin [bin2...]" to output a json of "ID": "output" for each unique value of "bin+bin2" (group by)', nargs = '+', default = ['fasta'], type = str)
argument_parser.add_argument('-l', '--length',       help = 'Summarize the lengths of sequences in a given column', type = str, default = '')
argument_parser.add_argument('-r', '--rearrangement',help = 'Summarize all reads by [H (family, gene, allele) : D (gene-allele) : J (gene, allele)]', action = 'store_const', const = True)
argument_parser.add_argument('-g', '--coverage',       help = 'Count how many reads passing the filter contain a mapped region; the columns to be counted are provided as a comma separated list; works in conjunction with --rearrgement', default = None)
argument_parser.add_argument('-t', '--light_chain', help='If set, assume that the reads come from a light Ig chain (no d-region)', action = 'store_const', const = True, default = False);
argument_parser.add_argument('-a', '--alternatives',  help='If set, read alternative rearrangements from this file and use those for rearranegment binning', required = False, type=argparse.FileType('r'));
argument_parser.add_argument('-v', '--inverse_d',   help='If set (for heavy chains only), produce the D-region assignment by also considering inverted regions', action = 'store_const', const = True, default = False);


cli_args = argument_parser.parse_args()

ig_reader = csv.reader (cli_args.input, delimiter = '\t')
column_headings = next (ig_reader)

if len (column_headings) < 10:
    raise BaseException ('The input file had too few columns')

for col in column_headings:
    column_mapper [col] = column_headings.index(col)

for req in required_headers:
    try:
        column_mapper [req]
    except KeyError as e:
        print ("Missing required column '%s'" % req)
        sys.exit (1)

line_filter             = None

is_light_chain = cli_args.light_chain

if cli_args.pattern is not None:
    for col_double in cli_args.pattern:
        line_filter = regExpClass (col_double, column_headings, line_filter)

if cli_args.column is not None:
    for col_triple in cli_args.column:
        line_filter = colLengthClass (col_triple, column_headings, line_filter)

cluster_sizes = None

if cli_args.cluster is not None:
    filter_set = set ()

    cluster_sizes = {}

    for group in json.load (cli_args.cluster):
        for cluster in group:
            filter_set.add (cluster['centroid'].split('\n')[0][1:])
            if not cluster['size'] in cluster_sizes:
                cluster_sizes[cluster['size']] = 0
            cluster_sizes[cluster['size']] += 1

    line_filter = idFilterClass ("Name", filter_set, column_headings, line_filter)


if cli_args.rearrangement is not None:
    cdr3_id = column_mapper["CDR3_AA"]
    print (outputSummary (line_filter, cli_args.support, column_mapper['Support'], None,reduceRearrangement,cdr3, cli_args.alternatives, coverages = cli_args.coverage, cluster_sizes = cluster_sizes))

else:
    if cli_args.length != '':
        if cli_args.length not in column_headings:
            raise BaseException ("%s is not a valid column name (must be one of %s)" % (cli_args.length, str(column_headings)))

        print (ouputLengths (line_filter, cli_args.support, column_mapper['Support'], column_headings.index (cli_args.length)))

    else:

        if cli_args.output[0] == 'fasta':
            ouputSequences (line_filter, cli_args.support, column_mapper['Support'], make_sequence_id, column_mapper['Mapped Read'])
        elif cli_args.output[0] == 'json':
            if len (cli_args.output) > 1:
                print (outputJSON (line_filter, cli_args.support, column_mapper['Support'],[make_sequence_id, [column_mapper[k] for k in cli_args.output[2:]], column_mapper[cli_args.output[1]]]))
            else:
                print (outputJSON (line_filter, cli_args.support, column_mapper['Support']))
        elif cli_args.output[0] == 'usage':
            if len (cli_args.output) > 1:
                print (outputUsage (line_filter, cli_args.support, column_mapper['Support'], column_mapper[cli_args.output[1]] ))
            else:
                print (outputUsage (line_filter, cli_args.support, column_mapper['Support'], reduceRearrangementToV ))
        else:
            print (outputSummary (line_filter, cli_args.support, column_mapper['Support'], column_mapper[cli_args.output[0]]))




