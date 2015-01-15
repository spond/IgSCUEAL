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
reference_alleles       = {};
index_reference         = {};
valid_nuc_sequences     = {};

for (seq_id = 0; seq_id < nucleotides.species; seq_id += 1) {
    GetString (seqName, nucleotides, seq_id);
    if (seq_id == 0) {
        is_j = (seqName $ "J-REGION")[0]>=0;
    }
    match = seqName$_nucAccessionExtractRE;
    if (match[0] >= 0) {
        bits  = splitOnRegExp (seqName, "\\|");
        allele = seqName[match[4]][match[5]];
        accession = seqName[match[2]][match[3]];
        is_01 = (allele $ "\\*01$")[0];
        
        if (is_01 >= 0) {
            reference_alleles[allele[0][is_01-1]] = Abs (nucleotideAccessions);
            index_reference [Abs (nucleotideAccessions)] = allele[0][is_01-1];
        }
        
        
        nucleotideAccessions + ((allele + "|" + accession)&&6);
        allelename +  seqName[match[4]][match[5]];
        valid_nuc_sequences + seq_id;
        
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

matched_accessions = {};

function pad_to_length (nl, pl) {
    ps = "";
    for (nli = nl; nli < 3*pl; nli += 1) {
        ps += "-";
    }
    return ps;
}

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
        
    GetDataInfo (nucsequence,nucleotidesF,valid_nuc_sequences[idx]);
    GetDataInfo (protsequence,proteinF,seq_id);
    nucsequence = nucsequence[nuc_offset[idx]][Abs(nucsequence)-1];
    
    matched_accessions [nucleotideAccessions[idx]] = protsequence;
   
    normalized_name = normalizeSequenceID (allelename[idx] + "_CRF_1", "used_names");
   
    nucsequence = mapCodonsToAA (nucsequence^{{"[^A-Z]",""}},protsequence);
    if (Abs (index_reference[idx])) {
        reference_alleles[index_reference[idx]] = protsequence;
    }
    
    fprintf (LAST_FILE_PATH, ">", normalized_name, "\n", nucsequence, pad_to_length (Abs (nucsequence), Abs(protsequence)), "\n");
}


if (Abs (matched_accessions) < Abs (nucleotideAccessions)) {
    for (si = 0; si < Abs (nucleotideAccessions); si += 1) {   
        if (si >= 0) {
            if (Abs (matched_accessions[nucleotideAccessions[si]]) == 0) {
                idx = nucleotideAccessions[si] $ "^[^\\]+";
                idx = (nucleotideAccessions[si])[0][idx[1]];
                simple = reference_alleles[idx];
                if (Type (simple) == "String") {
                    GetDataInfo (nucsequence,nucleotidesF,valid_nuc_sequences[si]);
                    nucsequence = nucsequence^{{"[^A-Z]",""}};
                    nucsequence_mapped = mapCodonsToAAFuzzy (nucsequence,simple,10);
                    if (None != nucsequence_mapped) {
                        normalized_name = normalizeSequenceID (allelename[si] + "_CRF_1", "used_names");
                        fprintf (LAST_FILE_PATH, ">", normalized_name, "\n", nucsequence_mapped, pad_to_length (Abs (nucsequence_mapped), Abs(protsequence)), "\n");
                    } else {
                        fprintf (stdout, "Failed to map ", allelename[si], "\n");
                    }
                }
            }
        }
    }
}

fprintf (LAST_FILE_PATH, CLOSE_FILE);


/*----------------------------------------------------------------------------------------------------------*/

