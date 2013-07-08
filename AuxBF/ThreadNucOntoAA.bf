LoadFunctionLibrary ("chooseGeneticCode", {"0": "Universal"});
LoadFunctionLibrary ("GrabBag");
LoadFunctionLibrary ("ReadDelimitedFiles");

_nucAccessionExtractRE = "^([0-9A-Z]+)\\|([^\\|]+)";

DataSet  protein                    = ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter  proteinF             = CreateFilter (protein,1);
DataSet  nucleotides                = ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter  nucleotidesF         = CreateFilter (nucleotides,1);

nucleotideAccessions    = {};
allelename              = {};
nuc_offset              = {};

for (seq_id = 0; seq_id < nucleotides.species; seq_id += 1) {
    GetString (seqName, nucleotides, seq_id);
    match = seqName$_nucAccessionExtractRE;
    if (match[0] >= 0) {
        bits  = splitOnRegExp (seqName, "\\|");
        nucleotideAccessions + ((seqName[match[4]][match[5]] + "|" + seqName[match[2]][match[3]])&&6);
        allelename +  seqName[match[4]][match[5]];
        nuc_offset + (-1+(bits[7]));
    }
}

used_names = {};

for (seq_id = 0; seq_id < protein.species; seq_id += 1) {
    GetString   (seqName, protein, seq_id);
    
    split_name = splitStringByTab(seqName);
    seqName = split_name[2] + "|" + split_name[3];
    idx = matchStringToSetOfPatterns (seqName, nucleotideAccessions);

    
    assert (idx >= 0, "No nucleotide sequence for `seqName`");
    
    GetDataInfo (nucsequence,nucleotidesF,idx);
    GetDataInfo (protsequence,proteinF,seq_id);
    nucsequence = nucsequence[nuc_offset[idx]][Abs(nucsequence)-1];
   
    normalized_name = normalizeSequenceID (allelename[idx] + "_CRF_1", "used_names");
   
    nucsequence = nucsequence^{{"[^A-Z]",""}};
    
    //fprintf (stdout, "\n****\n", nucsequence, "\n");
    
    fprintf (stdout, ">", normalized_name, "\n", mapCodonsToAA (nucsequence,protsequence),"\n");
}


/*----------------------------------------------------------------------------------------------------------*/

