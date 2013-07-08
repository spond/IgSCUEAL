LoadFunctionLibrary ("GrabBag");
LoadFunctionLibrary ("ReadDelimitedFiles");

SetDialogPrompt ("IMGT V-region protein alignment");
fscanf (PROMPT_FOR_FILE, "Lines", imgt_aa);

fprintf (stdout, "Regular expression to filter annotation with (e.g. F for functional, F|ORF for function or ORF):");
fscanf  (stdin, "String", regExpFilter);

sequences = {};

for (k = 0; k < Columns (imgt_aa); k+=2) {
    annotation = splitStringByTab (imgt_aa[k]);
    match = annotation[5] $ regExpFilter;
    if (match[0] == 0 && match[1] == Abs (annotation[5])-1) {
        sequence = splitOnRegExp(imgt_aa[k+1]^{{"\\.",""}}, "\ +");
        if (Abs (sequence) == 11) {
            sequences [imgt_aa[k]] =  {"FR1": sequence[0] + sequence[1],
                       "CDR1": sequence[2], 
                       "FR2": sequence[3]+sequence[4],
                       "CDR2": sequence[5] + sequence[6],
                       "FR3": sequence[6] + sequence[7] + sequence[8] + sequence [9],
                       "CDR3": sequence[10]};
                       
            if (Abs (sequences) == 1) {
                regions = Rows (sequences [imgt_aa[k]]);
            }
        }
    }
}

assert (Abs (sequences) > 0, "No 11-segment sequences matched the functional filter");

fprintf (stdout, "Save partitioned alignments to this directory (no trailing slash necessary):");
fscanf  (stdin, "String", outDirPath);



/*----------------------------------------------------------------------------------------------------------*/

