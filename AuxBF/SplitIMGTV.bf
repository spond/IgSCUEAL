LoadFunctionLibrary ("GrabBag");
LoadFunctionLibrary ("ReadDelimitedFiles");

SetDialogPrompt ("IMGT V-region protein alignment");
fscanf (PROMPT_FOR_FILE, "Lines", imgt_aa);

fprintf (stdout, "Regular expression to filter annotation with (e.g. F for functional, F|ORF for function or ORF):");
fscanf  (stdin, "String", regExpFilter);

sequences = {};

for (k = 0; k < Columns (imgt_aa); k+=2) {
    annotation = splitStringByTab (imgt_aa[k]);

    if (Abs (annotation) < 6) {
        fprintf (stdout, "[WARNING] Sequence description line ", imgt_aa[k], " does not contain 6 or more tab-separated fields and will be skipped\n");
        continue;
    }
    match = annotation[5] $ regExpFilter;
    if (match[0] == 0 && match[1] == Abs (annotation[5])-1) {
        v_region = imgt_aa[k+1]^{{"\\.",""}};
        sequence = splitOnRegExp(v_region, " +");
        if (Abs (sequence) == 11) {
            has_stops = (v_region $ "\\*");
            if ( has_stops [0] >= 0) {
                if (has_stops[0] == Abs (v_region)-1) {
                    fprintf (stdout, "[WARNING] Automatically truncated the trailing '*' in sequence ", 
                                      imgt_aa[k], "\n");
                    sequence [10] = sequence[10] ^ {{"\\*",""}};
                } else {               
                    fprintf (stdout, "[WARNING] Discarding sequence ", imgt_aa[k], 
                                " because the translation contained '*' characters\n");
                    continue;
                }
            }
            sequences [imgt_aa[k]] =  {"FR1": sequence[0] + sequence[1],
                       "CDR1": sequence[2], 
                       "FR2": sequence[3]+sequence[4],
                       "CDR2": sequence[5],
                       "FR3": sequence[6] + sequence[7] + sequence[8] + sequence [9],
                       "CDR3": sequence[10]};
                       
            if (Abs (sequences) == 1) {
                regions = Rows (sequences [imgt_aa[k]]);
            }
        } else {
            fprintf (stdout, "Sequence ", imgt_aa[k+1], " contained ", Abs (sequences), " space-separated blocks (expected 11) and has been skipped\n");     
        }
    } else {
        fprintf (stdout, "Sequence ", imgt_aa[k], " did not match the functional filter and has been skipped\n");
    }
}

assert (Abs (sequences) > 0, "No 11-segment sequences matched the functional filter");

fprintf (stdout, "Save partitioned alignments to this directory (no trailing slash necessary):");
fscanf  (stdin, "String", outDirPath);

for (rtw = 0; rtw < Columns (regions); rtw += 1) {
    region_to_write = regions[rtw];
    file_path = outDirPath + "/" + region_to_write + ".fas";
    fprintf (file_path, CLEAR_FILE);
    sequences ["segment_writer"][""];
}

function segment_writer (key, value) {
    fprintf (file_path, ">", key, "\n", value[region_to_write], "\n");
    return 0;
}


/*----------------------------------------------------------------------------------------------------------*/

