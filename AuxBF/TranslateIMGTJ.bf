LoadFunctionLibrary ("GrabBag");
LoadFunctionLibrary ("ReadDelimitedFiles");
LoadFunctionLibrary ("chooseGeneticCode", {"0": "Universal"});

SetDialogPrompt ("IMGT J-region nucleotide sequences");
DataSet j_region = ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter j_region_filter = CreateFilter (j_region,1);

sequences     = {};
nuc_sequences = {};

codon_mapping = defineCodonToAA ();

for (k = 0; k < j_region_filter.species; k+=1) {
    GetString (seqName, j_region_filter, k);
    seq_annotation = splitOnRegExp (seqName, "\\|");
    if (Abs (seq_annotation) != 16) {
        fprintf (stdout, "[WARNING] Sequence `seqName` is not fully annotated and will be skipped\n");
        continue;
    }
    frame = Max (0, (-1) + seq_annotation [7]);
    GetDataInfo (nuc_sequence, j_region_filter, k);
    coding_region = nuc_sequence[frame][frame+(Abs(nuc_sequence)$3*3)-1];
    aa_sequence   = translateCodonToAA (coding_region, codon_mapping, 0);
    if (Abs (aa_sequence) == 0 || (aa_sequence $ "\\?")[0] >= 0) {
        fprintf (stdout, "[WARNING] Sequence `seqName` appears to have premature stop codons or is empty and will be skipped\n");
        continue;
    }
    sequences [seqName] = aa_sequence;
    nuc_sequences[seqName] = coding_region;
}

assert (Abs (sequences) > 0, "No fully annotated sequences found");

SetDialogPrompt ("Save translated J amino-acid sequences (FASTA format) to:");
fprintf (PROMPT_FOR_FILE, CLEAR_FILE);
sequences["segment_writer"][""];

SetDialogPrompt ("Save clipped J nucleotide sequences (FASTA format) to:");
fprintf (PROMPT_FOR_FILE, CLEAR_FILE);
nuc_sequences["segment_writer"][""];

function segment_writer (key, value) {
    fprintf (LAST_FILE_PATH, ">", key, "\n", value, "\n");
    return 0;
}


/*----------------------------------------------------------------------------------------------------------*/

