import csv, argparse, sys, re, operator, json


########### START GLOBALS #############

required_headers = ['Name', 'Best Rearrangement', 'Support', 'Sequence']
column_mapper    = {}
cdr3_id          = 0
is_light_chain   = False

V_processor      = re.compile ('^V([0-9]+|\-ancestral)\-?([^\,]+)?')
J_processor      = re.compile ('\,J(\-ancestral|[0-9]+)')

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
        if self.comp_function (len(line[self.column_id]), self.length):
            if (self.next_filter is not None): return self.next_filter.check_line (line)
            return True
        return False        


####################################

def reduceRearrangement (line):
    rearrangement  = line[column_mapper['Best Rearrangement']]
    
    if is_light_chain:
        assignment = ['','','']    
        j_id = 2
    else:
        assignment = ['','','','']
        j_id = 3
    
    v_match  = V_processor.match (rearrangement)
    assignment[0] = 'V' + v_match.group(1).replace('-','').replace('ancestral','?')
    
    if v_match.group (2) is not None:
        assignment[1] = assignment[0] + '-' + v_match.group (2).replace('-','')
    else:
        assignment[1] = assignment[0] + '-?'
    
    j_match = J_processor.search (rearrangement)
    if j_match is not None and j_match.group(1) is not None:
        assignment[j_id] = 'J' + j_match.group(1).replace('-','').replace('ancestral','?');
        
    if not is_light_chain:
        assignment[2] = line[column_mapper['D_ALLELE']]
        if (assignment[2] == '0'):
            assignment[2] = ''
        else:
            assignment[2] = 'D' + assignment[2][1:].replace('_','-')
    
    return '|'.join(assignment)

####################################

def cdr3 (line):
    return len (line[cdr3_id])

####################################
        
def ouputSummary (line_filter, support, support_col, rearrangement_col, mappingFunction = None, histogram = None):
    binByRearrangement = {}
    total        = 0
    low_support  = 0
    mapped       = 0
    filtered_in  = 0
    
    for line in ig_reader:
        total += 1
        if len (line) > 10: # not a failed alignment
            mapped += 1
            if float (line[support_col]) >= support:
                if line_filter is not None:
                    if not line_filter.check_line (line): 
                        continue
                filtered_in += 1
                if mappingFunction:                    
                    alleles = mappingFunction (line)
                else:
                    alleles = line[rearrangement_col]
                    
                if (histogram):                   
                    if alleles not in binByRearrangement:
                        binByRearrangement [alleles] = {'Count' : 1}
                    else:
                        binByRearrangement [alleles]['Count'] += 1
                    
                    prop = histogram (line)
                    if prop not in binByRearrangement [alleles]:
                        binByRearrangement [alleles][prop] = 1
                    else:
                        binByRearrangement [alleles][prop] += 1
                                                
                else:
                    if alleles not in binByRearrangement:
                        binByRearrangement [alleles] = 1
                    else:
                        binByRearrangement [alleles] += 1
            else:
                low_support += 1

    binByRearrangement ['Summary']     = [['Total', total], ['Mapped', mapped], ['Low support', low_support], ['Passed filter', filtered_in]] 
                                        
    return json.dumps (binByRearrangement, indent = 2)                                    

####################################
        
def ouputLengths (line_filter, support, support_col, binning_column):
    binByLength = {}
    for line in ig_reader:
        if len (line) > 10: # not a failed alignment
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
        if len (line) > 10: # not a failed alignment
            if float (line[support_col]) >= support:
                if line_filter is not None:
                    if not line_filter.check_line (line): continue
                print (">%s\n%s\n" % (line[name_col], line[sequence_col]))       
                          

####################################

def outputJSON (line_filter, support, support_col, group_by = None):
    
    if group_by is None:
       output = []
       for line in ig_reader:
            if len (line) > 10: # not a failed alignment
                if float (line[support_col]) >= support:
                    if line_filter is not None:
                        if not line_filter.check_line (line): continue
                    output.append (line)
    else:
        output = {}
        for line in ig_reader:
            if len (line) > 10: # not a failed alignment
                if float (line[support_col]) >= support:
                    if line_filter is not None:
                        if not line_filter.check_line (line): continue
                        
                    #bin_by = matchNStr(line[group_by[1]])
                    bin_by = "|".join([line[k] for k in group_by[1]])
                    #print (bin_by, file = sys.stderr)
                    if bin_by not in output:
                        output [bin_by] = {}
                    output[bin_by][line[group_by[0]]] = line[group_by[2]]
                
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

argument_parser = argparse.ArgumentParser (description='Read IgSCUEAL output files.')
argument_parser.add_argument('-i', '--input',        help = 'The tab-separated file produced by IgSCUEAL. Must have at least 10 columns', nargs = '?', required = True, type=argparse.FileType('r'), default = sys.stdin)
argument_parser.add_argument('-s', '--support',      help = 'Minimum assignment support required to process a read (default is 0.9)', type=float01, default = 0.9)
argument_parser.add_argument('-p', '--pattern',      help = 'A regular expression pattern for pulling out reads that match a particular pattern in a particular column, e.g. "Rearrangement V3-21,J3"', type=str, nargs = 2, action = 'append')
argument_parser.add_argument('-c', '--column',       help = 'A column [l,g,e] length argument to select lines based on the length of the string in a given column, e.g. JUNCTION_AA g 10 (all lines with JUNCTIONS_AA with more than 10 chars)', nargs = 3, action = 'append')
argument_parser.add_argument('-x', '--cluster',      help = 'Only keep sequences which represent unique clonotypes', required = False, type=argparse.FileType('r'))
argument_parser.add_argument('-o', '--output',       help = 'Output mode (default fasta): either summarize filtered reads by column value (supply the column name), print FASTA for the reads matching selection criteria (fasta), or output JSON (json) for matching reads. '+ 
                                                          'Finally specify "json output bin [bin2...]" to output a json of "ID": "output" for each unique value of "bin+bin2" (group by)', nargs = '+', default = ['fasta'], type = str)
argument_parser.add_argument('-l', '--length',       help = 'Summarize the lengths of sequences in a given column', type = str, default = '')
argument_parser.add_argument('-r', '--rearrangement',help = 'Summarize all reads by [H-family,H-allele,D-allele,J-allele]', action = 'store_const', const = True)
argument_parser.add_argument('-t', '--light_chain', help='If set, assume that the reads come from a light Ig chain (no d-region)', action = 'store_const', const = True, default = False);


cli_args = argument_parser.parse_args()

ig_reader = csv.reader (cli_args.input, delimiter = '\t')
column_headings = next (ig_reader)

if len (column_headings) < 10:
    raise BaseException ('The input file had too few columns')
    
for col in column_headings:
    column_mapper [col] = column_headings.index(col)
    
line_filter             = None

is_light_chain = cli_args.light_chain

if cli_args.pattern is not None:
    for col_double in cli_args.pattern:
        line_filter = regExpClass (col_double, column_headings, line_filter)

if cli_args.column is not None:
    for col_triple in cli_args.column:
        line_filter = colLengthClass (col_triple, column_headings, line_filter)
    
if cli_args.cluster is not None:
    filter_set = set ()
    for group in json.load (cli_args.cluster):
        for cluster in group:
           filter_set.add (cluster['centroid'].split('\n')[0][1:])
        
    line_filter = idFilterClass ("Name", filter_set, column_headings, line_filter)

if cli_args.rearrangement is not None:
    cdr3_id = column_mapper["JUNCTION_AA"]
    print (ouputSummary (line_filter, cli_args.support, column_mapper['Support'], None,reduceRearrangement,cdr3))
    
else:    
    if cli_args.length != '':
        if cli_args.length not in column_headings:
            raise BaseException ("%s is not a valid column name (must be one of %s)" % (cli_args.length, str(column_headings)))
        
        print (ouputLengths (line_filter, cli_args.support, column_mapper['Support'], column_headings.index (cli_args.length)))
        
    else:
            
        if cli_args.output[0] == 'fasta':
            ouputSequences (line_filter, cli_args.support, column_mapper['Support'], column_mapper['Name'], column_mapper['Mapped Read'])
        elif cli_args.output[0] == 'json':
            if len (cli_args.output) > 1:
                print (outputJSON (line_filter, cli_args.support, column_mapper['Support'],[column_mapper["Name"], [column_mapper[k] for k in cli_args.output[2:]], column_mapper[cli_args.output[1]]]))                
            else:
                print (outputJSON (line_filter, cli_args.support, column_mapper['Support']))
        else:
            print (ouputSummary (line_filter, cli_args.support, column_mapper['Support'], column_mapper[cli_args.output]))
    
    
	

