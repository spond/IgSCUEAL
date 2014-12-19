LoadFunctionLibrary ("GrabBag");

filesIn = {{"FR1","CDR1","FR2","CDR2","FR3","CDR3","J"}};
fprintf (stdout, "The directory containing V-region amino-acid alignments:");
fscanf (stdin, "String", base_directory);
extension = ".fas.aligned";

segments = Columns (filesIn);

for (k = 0; k < segments; k+=1) {
    filesIn [k] = base_directory + "/" + filesIn [k] + extension;
}

fprintf (stdout, "Loading ", filesIn[0], "\n");

DataSet 			first 		= ReadDataFile (filesIn[0]);
GetString			(masterOrdering, first, -1);

segment_lengths = {segments, 1};
segment_lengths [0] = first.sites -1;

for (fileID = 1; fileID < segments; fileID += 1)
{
    fprintf (stdout, "Loading ", filesIn[fileID], "\n");
	DataSet 		current = ReadDataFile (filesIn[fileID]);
	
	if (fileID < segments-1) {
        GetString		(currentOrdering, current, -1);
        mapping	= mapSets 		(masterOrdering, currentOrdering);
        DataSetFilter	dsf		= CreateFilter (current, 1, "", Join (",", mapping));
        Export 			(filterStr, dsf);
        DataSet 		current = ReadFromString (filterStr);
        DataSet			first	= Concatenate 	 (purge, first, current);
        segment_lengths [fileID] = first.sites;
    } else {
        segment_lengths [fileID] = segment_lengths [fileID-1] + current.sites;
    }
}

DataSetFilter dsf = CreateFilter (first, 1);
SetDialogPrompt ("Write merged protein alignment to");
fprintf (PROMPT_FOR_FILE,CLEAR_FILE, dsf); 
aa_config = LAST_FILE_PATH + ".partition_info";

fprintf (aa_config, CLEAR_FILE, "
IG_MRCA_sequence_length = ",segment_lengths[segments-1]*3,";
_IG_nuc_annotation = {1,IG_MRCA_sequence_length};
aa_breaks = ", segment_lengths, ";
jMotifRegExp       = \"[FW]G.G\";\n");
