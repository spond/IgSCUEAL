LoadFunctionLibrary ("chooseGeneticCode", {"0": "Universal"});
LoadFunctionLibrary ("GrabBag");
LoadFunctionLibrary ("ReadDelimitedFiles");


_nucAccessionExtractRE = "^([0-9A-Z]+)\\|([^\\|]+)";

SetDialogPrompt                     ("Amino-acid alignment:");
DataSet  protein                    = ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter  proteinF             = CreateFilter (protein,1);
SetDialogPrompt                     ("Nucleotide sequences (unaligned)");
DataSet  nucleotides                = ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter  nucleotidesF         = CreateFilter (nucleotides,1);

nucleotideAccessions    = {};
allelename              = {};
nuc_offset              = {};

for (seq_id = 0; seq_id < nucleotides.species; seq_id += 1) {
    GetString (seqName, nucleotides, seq_id);
    if (seq_id == 0) {
        is_j = (seqName $ "J-REGION")[0]>=0;
    }
    match = seqName$_nucAccessionExtractRE;
    if (match[0] >= 0) {
        bits  = splitOnRegExp (seqName, "\\|");
        nucleotideAccessions + ((seqName[match[4]][match[5]] + "|" + seqName[match[2]][match[3]])&&6);
        allelename +  seqName[match[4]][match[5]];
        if (is_j) {
            nuc_offset + 0;
        } else {
            nuc_offset + (-1+(bits[7]));
        }
    }
}

used_names = {};

SetDialogPrompt ("Save the V-region nucleotide alignment (as FASTA) to:");
fprintf (PROMPT_FOR_FILE, CLEAR_FILE, KEEP_OPEN);

for (seq_id = 0; seq_id < protein.species; seq_id += 1) {
    GetString   (seqName, protein, seq_id);
    
    if (is_j) {
        split_name = splitOnRegExp(seqName, "\\|");
         seqName = split_name[1] + "|" + split_name[0];
    } else {
        split_name = splitOnRegExp(seqName, "[\\ \t]");
        seqName = split_name[2] + "|" + split_name[3];
    }
    idx = matchStringToSetOfPatterns (seqName, nucleotideAccessions);

    
    assert (idx >= 0, "No nucleotide sequence for `seqName`");
    
    GetDataInfo (nucsequence,nucleotidesF,idx);
    GetDataInfo (protsequence,proteinF,seq_id);
    nucsequence = nucsequence[nuc_offset[idx]][Abs(nucsequence)-1];
   
    normalized_name = normalizeSequenceID (allelename[idx] + "_CRF_1", "used_names");
   
    nucsequence = nucsequence^{{"[^A-Z]",""}};
    
    fprintf (LAST_FILE_PATH, ">", normalized_name, "\n", mapCodonsToAA (nucsequence,protsequence),"\n");
}

fprintf (LAST_FILE_PATH, CLOSE_FILE);


/*----------------------------------------------------------------------------------------------------------*/

