
_nucleotide_rc = {};
_nucleotide_rc["A"] = "T";
_nucleotide_rc["C"] = "G";
_nucleotide_rc["G"] = "C";
_nucleotide_rc["T"] = "A";
_nucleotide_rc["M"] = "K";
_nucleotide_rc["R"] = "Y";
_nucleotide_rc["W"] = "W";
_nucleotide_rc["S"] = "S";
_nucleotide_rc["Y"] = "R";
_nucleotide_rc["K"] = "M";
_nucleotide_rc["B"] = "V";  /* not A */
_nucleotide_rc["D"] = "H";  /* not C */
_nucleotide_rc["H"] = "D";  /* not G */
_nucleotide_rc["V"] = "B";  /* not T */
_nucleotide_rc["N"] = "N";

alignOptions = {};

alignOptions ["SEQ_ALIGN_CHARACTER_MAP"]="ARNDCQEGHILKMFPSTWYV";


if (useBlosum62) {
    scoreMatrix = 
    {
    {                 6,                -3,                -4,                -4,                -2,                -2,                -2,                -1,                -3,                -3,                -3,                -2,                -2,                -4,                -2,                 0,                -1,                -5,                -3,                -1,                -4,                -2,                -2,                -7}
    {                -3,                 8,                -2,                -4,                -6,                 0,                -2,                -5,                -2,                -6,                -4,                 1,                -3,                -5,                -4,                -2,                -3,                -5,                -3,                -5,                -3,                -1,                -2,                -7}
    {                -4,                -2,                 8,                 0,                -5,                -1,                -2,                -2,                 0,                -6,                -6,                -1,                -4,                -5,                -4,                 0,                -1,                -7,                -4,                -5,                 6,                -2,                -2,                -7}
    {                -4,                -4,                 0,                 8,                -6,                -2,                 0,                -3,                -3,                -5,                -6,                -2,                -6,                -6,                -3,                -1,                -3,                -7,                -6,                -6,                 6,                 0,                -3,                -7}
    {                -2,                -6,                -5,                -6,                10,                -5,                -7,                -5,                -5,                -3,                -3,                -6,                -3,                -5,                -5,                -2,                -2,                -4,                -4,                -2,                -5,                -6,                -4,                -7}
    {                -2,                 0,                -1,                -2,                -5,                 8,                 1,                -4,                 0,                -6,                -4,                 0,                -1,                -6,                -3,                -1,                -2,                -3,                -3,                -4,                -1,                 6,                -2,                -7}
    {                -2,                -2,                -2,                 0,                -7,                 1,                 7,                -4,                -1,                -6,                -5,                 0,                -4,                -6,                -3,                -1,                -2,                -5,                -4,                -4,                 0,                 6,                -2,                -7}
    {                -1,                -5,                -2,                -3,                -5,                -4,                -4,                 7,                -4,                -7,                -6,                -3,                -5,                -5,                -4,                -2,                -4,                -4,                -5,                -6,                -2,                -4,                -4,                -7}
    {                -3,                -2,                 0,                -3,                -5,                 0,                -1,                -4,                10,                -6,                -5,                -2,                -3,                -3,                -4,                -2,                -4,                -5,                 0,                -6,                -1,                -1,                -3,                -7}
    {                -3,                -6,                -6,                -5,                -3,                -6,                -6,                -7,                -6,                 6,                 0,                -5,                 0,                -1,                -5,                -5,                -2,                -5,                -3,                 2,                -5,                -6,                -2,                -7}
    {                -3,                -4,                -6,                -6,                -3,                -4,                -5,                -6,                -5,                 0,                 6,                -5,                 1,                -1,                -5,                -5,                -3,                -3,                -3,                 0,                -6,                -5,                -2,                -7}
    {                -2,                 1,                -1,                -2,                -6,                 0,                 0,                -3,                -2,                -5,                -5,                 7,                -3,                -6,                -2,                -1,                -2,                -5,                -3,                -4,                -2,                 0,                -2,                -7}
    {                -2,                -3,                -4,                -6,                -3,                -1,                -4,                -5,                -3,                 0,                 1,                -3,                 9,                -1,                -5,                -3,                -2,                -3,                -3,                 0,                -5,                -2,                -1,                -7}
    {                -4,                -5,                -5,                -6,                -5,                -6,                -6,                -5,                -3,                -1,                -1,                -6,                -1,                 8,                -6,                -4,                -4,                 0,                 1,                -3,                -6,                -6,                -3,                -7}
    {                -2,                -4,                -4,                -3,                -5,                -3,                -3,                -4,                -4,                -5,                -5,                -2,                -5,                -6,                 9,                -2,                -3,                -6,                -5,                -4,                -4,                -3,                -4,                -7}
    {                 0,                -2,                 0,                -1,                -2,                -1,                -1,                -2,                -2,                -5,                -5,                -1,                -3,                -4,                -2,                 7,                 0,                -5,                -3,                -4,                -1,                -1,                -2,                -7}
    {                -1,                -3,                -1,                -3,                -2,                -2,                -2,                -4,                -4,                -2,                -3,                -2,                -2,                -4,                -3,                 0,                 7,                -4,                -3,                -1,                -2,                -2,                -2,                -7}
    {                -5,                -5,                -7,                -7,                -4,                -3,                -5,                -4,                -5,                -5,                -3,                -5,                -3,                 0,                -6,                -5,                -4,                12,                 0,                -6,                -7,                -4,                -4,                -7}
    {                -3,                -3,                -4,                -6,                -4,                -3,                -4,                -5,                 0,                -3,                -3,                -3,                -3,                 1,                -5,                -3,                -3,                 0,                 9,                -3,                -5,                -3,                -2,                -7}
    {                -1,                -5,                -5,                -6,                -2,                -4,                -4,                -6,                -6,                 2,                 0,                -4,                 0,                -3,                -4,                -4,                -1,                -6,                -3,                 6,                -6,                -4,                -2,                -7}
    {                -4,                -3,                 6,                 6,                -5,                -1,                 0,                -2,                -1,                -5,                -6,                -2,                -5,                -6,                -4,                -1,                -2,                -7,                -5,                -6,                 7,                -1,                -3,                -7}
    {                -2,                -1,                -2,                 0,                -6,                 6,                 6,                -4,                -1,                -6,                -5,                 0,                -2,                -6,                -3,                -1,                -2,                -4,                -3,                -4,                -1,                 7,                -2,                -7}
    {                -2,                -2,                -2,                -3,                -4,                -2,                -2,                -4,                -3,                -2,                -2,                -2,                -1,                -3,                -4,                -2,                -2,                -4,                -2,                -2,                -3,                -2,                -2,                -7}
    {                -7,                -7,                -7,                -7,                -7,                -7,                -7,                -7,                -7,                -7,                -7,                -7,                -7,                -7,                -7,                -7,                -7,                -7,                -7,                -7,                -7,                -7,                -7,                 1}
    };
    
    _protBaseFrequencies = {{         0.074}
                            {         0.025}
                            {         0.054}
                            {         0.054}
                            {         0.047}
                            {         0.074}
                            {         0.026}
                            {         0.068}
                            {         0.058}
                            {         0.099}
                            {         0.025}
                            {         0.045}
                            {         0.039}
                            {         0.034}
                            {         0.052}
                            {         0.057}
                            {         0.051}
                            {         0.073}
                            {         0.013}
                            {         0.032}
                            };
}
else
{
    scoreMatrix = 
    {
    {                 7,                -7,                -7,                -4,               -10,               -11,                -4,                -3,               -10,                -6,                -9,                -9,                -7,               -13,                -3,                -2,                 1,               -16,               -15,                 0,                -5,                -5,                -3,               -17}
    {                -7,                 7,                -5,               -11,                -8,                -2,                -7,                -2,                 0,                -6,                -6,                 2,                -3,               -12,                -4,                -2,                -2,                -5,                -9,               -10,                -7,                -3,                -3,               -17}
    {                -7,                -5,                 8,                 2,                -9,                -6,                -6,                -7,                 0,                -6,               -12,                 0,               -10,               -12,                -9,                 1,                 0,               -17,                -3,               -10,                 6,                -6,                -3,               -17}
    {                -4,               -11,                 2,                 8,               -14,               -10,                 0,                -2,                -3,               -11,               -15,                -7,               -13,               -15,               -13,                -5,                -6,               -16,                -6,                -5,                 7,                 0,                -3,               -17}
    {               -10,                -8,                -9,               -14,                11,               -16,               -15,                -5,                -7,               -11,                -9,               -13,               -14,                 0,               -12,                -1,                -6,                -2,                 0,                -8,               -10,               -16,                -5,               -17}
    {               -11,                -2,                -6,               -10,               -16,                 8,                -2,               -10,                 0,               -12,                -4,                 0,                -8,               -12,                -1,                -9,                -8,               -14,                -9,               -13,                -7,                 6,                -4,               -17}
    {                -4,                -7,                -6,                 0,               -15,                -2,                 7,                -1,                -9,               -12,               -15,                -1,               -10,               -17,               -13,               -11,                -8,               -15,               -12,                -5,                 0,                 6,                -4,               -17}
    {                -3,                -2,                -7,                -2,                -5,               -10,                -1,                 7,               -10,               -11,               -14,                -6,               -12,                -9,               -11,                -1,                -7,                -5,               -14,                -5,                -4,                -3,                -4,               -17}
    {               -10,                 0,                 0,                -3,                -7,                 0,                -9,               -10,                10,               -10,                -4,                -5,               -10,                -6,                -3,                -6,                -6,               -11,                 2,               -14,                -1,                -2,                -3,               -17}
    {                -6,                -6,                -6,               -11,               -11,               -12,               -12,               -11,               -10,                 7,                 0,                -7,                 0,                -2,               -10,                -4,                 0,               -14,                -9,                 2,                -7,               -12,                -2,               -17}
    {                -9,                -6,               -12,               -15,                -9,                -4,               -15,               -14,                -4,                 0,                 6,               -10,                 0,                 0,                -3,                -5,                -8,                -6,                -8,                -4,               -13,                -6,                -4,               -17}
    {                -9,                 2,                 0,                -7,               -13,                 0,                -1,                -6,                -5,                -7,               -10,                 7,                -4,               -14,                -9,                -5,                -1,               -12,               -13,                -9,                -1,                -1,                -2,               -17}
    {                -7,                -3,               -10,               -13,               -14,                -8,               -10,               -12,               -10,                 0,                 0,                -4,                10,                -7,               -11,                -9,                -1,               -11,               -15,                 0,               -11,                -9,                -3,               -17}
    {               -13,               -12,               -12,               -15,                 0,               -12,               -17,                -9,                -6,                -2,                 0,               -14,                -7,                10,               -11,                -5,               -10,                -5,                 1,                -5,               -13,               -14,                -3,               -17}
    {                -3,                -4,                -9,               -13,               -12,                -1,               -13,               -11,                -3,               -10,                -3,                -9,               -11,               -11,                 8,                -1,                -3,               -13,               -11,               -12,               -10,                -3,                -5,               -17}
    {                -2,                -2,                 1,                -5,                -1,                -9,               -11,                -1,                -6,                -4,                -5,                -5,                -9,                -5,                -1,                 8,                 0,               -12,                -6,                -9,                 0,               -10,                -3,               -17}
    {                 1,                -2,                 0,                -6,                -6,                -8,                -8,                -7,                -6,                 0,                -8,                -1,                -1,               -10,                -3,                 0,                 7,               -16,               -10,                -4,                -2,                -8,                -2,               -17}
    {               -16,                -5,               -17,               -16,                -2,               -14,               -15,                -5,               -11,               -14,                -6,               -12,               -11,                -5,               -13,               -12,               -16,                10,                -4,               -16,               -16,               -14,                -8,               -17}
    {               -15,                -9,                -3,                -6,                 0,                -9,               -12,               -14,                 2,                -9,                -8,               -13,               -15,                 1,               -11,                -6,               -10,                -4,                10,               -12,                -4,               -10,                -4,               -17}
    {                 0,               -10,               -10,                -5,                -8,               -13,                -5,                -5,               -14,                 2,                -4,                -9,                 0,                -5,               -12,                -9,                -4,               -16,               -12,                 7,                -7,                -7,                -3,               -17}
    {                -5,                -7,                 6,                 7,               -10,                -7,                 0,                -4,                -1,                -7,               -13,                -1,               -11,               -13,               -10,                 0,                -2,               -16,                -4,                -7,                 7,                -2,                -4,               -17}
    {                -5,                -3,                -6,                 0,               -16,                 6,                 6,                -3,                -2,               -12,                -6,                -1,                -9,               -14,                -3,               -10,                -8,               -14,               -10,                -7,                -2,                 6,                -4,               -17}
    {                -3,                -3,                -3,                -3,                -5,                -4,                -4,                -4,                -3,                -2,                -4,                -2,                -3,                -3,                -5,                -3,                -2,                -8,                -4,                -3,                -4,                -4,                -3,               -17}
    {               -17,               -17,               -17,               -17,               -17,               -17,               -17,               -17,               -17,               -17,               -17,               -17,               -17,               -17,               -17,               -17,               -17,               -17,               -17,               -17,               -17,               -17,               -17,                 1}
    };	
    _protBaseFrequencies = {
                            {0.060490222}
                            {0.020075899}
                            {0.042109048}
                            {0.071567447}
                            {0.028809447}
                            {0.072308239}
                            {0.022293943}
                            {0.069730629}
                            {0.056968211}
                            {0.098851122}
                            {0.019768318}
                            {0.044127815}
                            {0.046025282}
                            {0.053606488}
                            {0.066039665}
                            {0.050604330}
                            {0.053636813}
                            {0.061625237}
                            {0.033011601}
                            {0.028350243}
                            };
}

alignOptions ["SEQ_ALIGN_SCORE_MATRIX"] = 	scoreMatrix[{{0,0}}][{{19,19}}];
alignOptions ["SEQ_ALIGN_GAP_OPEN"]		= 	40;
alignOptions ["SEQ_ALIGN_GAP_OPEN2"]	= 	20;
alignOptions ["SEQ_ALIGN_GAP_EXTEND"]	= 	10;
alignOptions ["SEQ_ALIGN_GAP_EXTEND2"]	= 	5;
alignOptions ["SEQ_ALIGN_AFFINE"]		=   1;
alignOptions ["SEQ_ALIGN_NO_TP"]		=   1;

outputAlignment = "";

LoadFunctionLibrary  ("chooseGeneticCode");
 
protScoreMatrix = scoreMatrix;
scoreMatrix = {64,64};
_HY_NUC_CODON_HAVE_SCORE_MATRIX = 1;
ExecuteAFile (HYPHY_LIB_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"SeqAlignmentCodonShared.ibf");
protScoreMatrix = alignOptions["SEQ_ALIGN_SCORE_MATRIX"];
alignOptions = {};


protLetters = "ARNDCQEGHILKMFPSTWYV";
_cdnaln_cdnScoreMatrix = pSM2cSM(protScoreMatrix, protLetters);


alignOptions ["SEQ_ALIGN_SCORE_MATRIX"] = 	_cdnaln_cdnScoreMatrix;
maxScore = Max (protScoreMatrix,0);
minScore = Min (protScoreMatrix,0);

    
alignOptions ["SEQ_ALIGN_GAP_OPEN"]		= 	1.5*Max(maxScore,-minScore);
alignOptions ["SEQ_ALIGN_GAP_OPEN2"]	= 	1.5*Max(maxScore,-minScore);
alignOptions ["SEQ_ALIGN_GAP_EXTEND"]	= 	1;
alignOptions ["SEQ_ALIGN_GAP_EXTEND2"]	= 	1;
alignOptions ["SEQ_ALIGN_FRAMESHIFT"]	= 	2*Max(maxScore,-minScore);
alignOptions ["SEQ_ALIGN_CODON_ALIGN"]	= 	1;


alignOptions ["SEQ_ALIGN_CHARACTER_MAP"]=  "ACGT";
alignOptions ["SEQ_ALIGN_NO_TP"]		=   1;
alignOptions ["SEQ_ALIGN_AFFINE"]		=   1;

shift_penalty = computeExpectedPerBaseScore (.2,protScoreMatrix,_protBaseFrequencies);
//fprintf (stdout, "\n****", shift_penalty, "******\n");
_cdnaln_partialScoreMatrices = cSM2partialSMs(_cdnaln_cdnScoreMatrix, {{shift_penalty__*1.5,shift_penalty__,shift_penalty__,shift_penalty*1.5}});

alignOptions ["SEQ_ALIGN_PARTIAL_3x1_SCORES"] = _cdnaln_partialScoreMatrices["3x1"];
alignOptions ["SEQ_ALIGN_PARTIAL_3x2_SCORES"] = _cdnaln_partialScoreMatrices["3x2"];
alignOptions ["SEQ_ALIGN_PARTIAL_3x4_SCORES"] = _cdnaln_partialScoreMatrices["3x4"];
alignOptions ["SEQ_ALIGN_PARTIAL_3x5_SCORES"] = _cdnaln_partialScoreMatrices["3x5"];	

if (_RUN_ALIGNER_AS_A_LIBRARY) {
    return 0;
}

SetDialogPrompt 			("Reference alignment:");
DataSet ref_ds 			  =  ReadDataFile (PROMPT_FOR_FILE);
referenceAlignmentPath	  = LAST_FILE_PATH;
if (verboseFlag) {
	fprintf 				   (stdout, "Read ", ref_ds.species-1, " reference sequences.\n");
}

refTree	= DATAFILE_TREE;

fscanf (stdin,"String",sequenceName);
fscanf (stdin,"String",sequenceData);

inputAlignment = ">"+sequenceName+"\n"+sequenceData;

DataSet ds_to_align = ReadFromString (inputAlignment);

seqCount			= ds_to_align.species;

if (verboseFlag) {
	fprintf (stdout, "Read ", seqCount, " sequences.\n");
}

DataSetFilter	ref_filter		= CreateFilter (ref_ds,1);
DataSetFilter	qry_filter		= CreateFilter (ds_to_align,1);

GetInformation 	(refSeqs, ref_filter);
GetInformation 	(qrySeqs, qry_filter);

refSequence	     = refSeqs[ref_ds.species-1];
_qry_sequence_id = 0;
qrySequence = qrySeqs[_qry_sequence_id];
additional_sequences = {};

function handle_alignment_templates (file_path) {
    afn = "../data/" + file_path;
    if (!afn) {
        LoadFunctionLibrary         ("ReadDelimitedFiles");
        DataSet                     ancestors = ReadDataFile (afn);
        DataSetFilter               ancestral_filter = CreateFilter (ancestors, 1);   
        GetInformation 	            (ancSeqs, ancestral_filter);
        
 
        _subtypeAssignmentByNode    ["NODE1"] = "MRCA"; 
        
        
        GetString (seq_names, ancestral_filter, -1);
        
        for (pattern = 0; pattern < Abs (_alignmentTemplates); pattern += 1) {
            _additional_references_to_consider = {"0" : {}, "1" : {}};
            for (s = 0; s < ancestral_filter.species; s += 1) {
                
                seq_type = _subtypeAssignmentByNode[seq_names[s] && 1];
                if (Type (seq_type) == "String") {
                    pattern_id = matchStringToSetOfPatterns (seq_type, (_alignmentTemplates[pattern])["Components"]);
                    if (pattern_id >= 0) {
                         (_additional_references_to_consider[pattern_id])[s] = 1;
                    }   
                }
            }
            
            bp = (_alignmentTemplates[pattern])["BP"];
            
            seq_ids  = Rows (_additional_references_to_consider[0]);
            seq_ids2 = Rows (_additional_references_to_consider[1]);
            
            for (i1 = 0; i1 <  Abs (_additional_references_to_consider[0]); i1+=1) {
                vs = ancSeqs[0 + seq_ids[i1]];
                for (i2 = 0; i2 < Abs (_additional_references_to_consider[1]); i2+=1) {
                    js = ancSeqs[0 + seq_ids2[i2]];
                    igg = ( vs[0][bp-1] + js [bp][Abs(js)-1]) ;
                    igg_stripped = (igg ^ {{"---", ""}}) ^ {{"-","N"}};
                    additional_sequences[igg_stripped] = {2,1};
                    (additional_sequences[igg_stripped])[0] = _subtypeAssignmentByNode[seq_names[0 + seq_ids[i1]]] + "-" + _subtypeAssignmentByNode[seq_names[0 + seq_ids2[i2]]];
                    (additional_sequences[igg_stripped])[1] = igg;
                }
            }
        }
    }
}

if (Type (_alignmentTemplates) == "AssociativeList") {
    handle_alignment_templates (ancestralAlignmentFileName);
    handle_alignment_templates (referenceAlignmentFileName);
}

//fprintf (stdout, Abs (additional_sequences, "\n");

GetString 	(qryName, ds_to_align, _qry_sequence_id);
toAlignDS	= ">REFERENCE\n" + refSequence + "\n>" + qryName + "\n"+qrySequence;
	
DataSet		   unal = ReadFromString (toAlignDS);
DataSetFilter  filteredData 	= CreateFilter	(unal,1);
GetInformation (CleanedSequences,filteredData);
	/* preprocess sequences */
	
for (seqCounter = 0; seqCounter < Columns(CleanedSequences); seqCounter += 1) {
    aSeq = CleanedSequences[seqCounter];
    CleanedSequences[seqCounter] = aSeq^{{"[^a-zA-Z]",""}};
    
    if (_allowNsInSequences == 0) {
        CleanedSequences[seqCounter] = CleanedSequences[seqCounter]^{{"^N+",""}};
        CleanedSequences[seqCounter] = CleanedSequences[seqCounter]^{{"N+$",""}};
    }
}

_reference_sequence = CleanedSequences[0];
_query_sequence     = CleanedSequences[1];
	    
_expected_score = Max(0,computeExpectedPerBaseScore (.4,protScoreMatrix,_protBaseFrequencies))*Abs(aSeq)$3;
bestScore = _expected_score;		
		
inStr = {{_reference_sequence,_query_sequence}};	
		
AlignSequences(aligned, inStr, alignOptions);
aligned = aligned[0];

if (aligned[0] >= bestScore) {
    bestScore   = aligned[0];
    bestAlignment = aligned;
}


additional_sequences ["try_alignment"][""];
function try_alignment (key, value) {
    inStr [0] = key;	
    //fprintf (stdout, key, ":", Abs (key), "\n\n");
    AlignSequences(aligned, inStr, alignOptions);
    aligned = aligned[0];
    if (aligned[0] >= bestScore) {
        bestScore   = aligned[0];
        
        rs = aligned[1];
        qs = aligned[2];
        
        //fprintf (stdout, "\n\n\n>r\n", aligned[1], "\n>q\n", aligned[2], "\n");
        
        aligned[1] = "";
        aligned[2] = "";
        aligned[1] * 400;
        aligned[2] * 400;
        
        letters = {};
        lc      = -1;
        seq_in_ref = value[1];
        
        for (cc = 0; cc < Abs (seq_in_ref); cc+=1) {
            if (seq_in_ref[cc] != "-") {
                lc += 1;
            } else {
                letters[lc] += 1;
            }
        }
        
        //fprintf (stdout, letters, "\n");
        
        lc = -1;
        rc = -1;
        
        for (padder = 0; padder < letters[lc]; padder += 1) {
            aligned[1] += "N";
            aligned[2] += "-";           
        }
        
        while (rc < Abs (rs)) {
            rc += 1;
            aligned[1] * rs[rc];
            aligned[2] * qs[rc];
            if (rs[rc] != "-") {
                lc += 1;
                for (padder = 0; padder < letters[lc]; padder += 1) {
                    aligned[1] * "N";
                    aligned[2] * "-";           
                }
            }
        }
        
        aligned[1] * 0;
        aligned[2] * 0;
        
        //fprintf (stdout, "\n>r\n", aligned[1], "\n>q\n", aligned[2], "\n");
        //fprintf (stdout, value[0], "\n");
        bestAlignment = aligned;
        
    }    
}

if (tryReverseComplement) {
    inStr [1] = nucleotideReverseComplement (_query_sequence);
    AlignSequences(aligned, inStr, alignOptions);
    aligned = aligned[0];
    if (aligned[0] >= bestScore) {
        bestScore   = aligned[0];
        bestAlignment = aligned;
    }
    additional_sequences ["try_alignment"][""];
}

if (bestScore == _expected_score) {
    fprintf (stdout, "Poor alignment scores : homology too low");
    return 1;
}
				

aligned_reference = bestAlignment[1];
aligned_query     = bestAlignment[2];
//fprintf (stdout, "\n\n", aligned_reference, "\n\n", aligned_query, "\n\n");
cleaned_up = correctReadUsingCodonAlignedData (aligned_reference, aligned_query);
gappedSeqN  = cleaned_up ["QRY"];
fullRefSeq  = cleaned_up ["REF"];
ref_shift    = cleaned_up ["OFFSET"];
//fprintf (stdout, "\n\n", fullRefSeq, "\n\n", gappedSeqN, "\n\n");

if (ref_shift) {
    aligned_reference    = aligned_reference[0][ref_shift-1] + fullRefSeq;
    aligned_query    = aligned2[0][ref_shift-1] + gappedSeqN;
    fullRefSeq = aligned_reference;
    gappedSeqN = aligned_query;
} else {
    aligned_reference    = fullRefSeq;
    aligned_query       = gappedSeqN;
}

	
if (Abs (_annotateSequenceByAlignment)) {
    _extraResult = Eval (_annotateSequenceByAlignment + "(aligned_reference,aligned_query)");
}
else {
    _extraResult = 0;
}
	
regExpS = gappedSeqN $ "^\\-+";
if (regExpS[0]>=0 && KEEP_ALL_GAPS_IN == 0) {
    startFrom 		= (regExpS		[1]+1);		
} else {
    startFrom = 0;
}

regExpE	= gappedSeqN $ "\\-+$";

if (regExpE[0]>=0 && KEEP_ALL_GAPS_IN == 0) {
    endAt		= regExpE[0];
} else {
    endAt		= Abs (gappedSeqN);
}	 

outSeqs = {};
for (k=0; k< ref_ds.species; k +=1) {
    outSeqs[k] = "";
    outSeqs[k] * 128;		
}
	
shift = 0;


//fprintf (stdout, "\n\n", fullRefSeq, "\n", gappedSeqN, "\n\n", KEEP_ALL_GAPS_IN, "\n\n", startFrom, ":", endAt, "\n\n");

gappedSeqN_Stripped = ""; gappedSeqN_Stripped * 128;

for (s=startFrom; s<endAt; s+=1) {
    //fprintf (stdout, fullRefSeq[s], gappedSeqN[s], ":", shift, "\n");
    if (fullRefSeq[s] == "-") {
        shift += 1;
    }
    else {
        gappedSeqN_Stripped * gappedSeqN[s];
        /* only insert reference characters if 
           the query DOES NOT have indels in 
           that position */
         
        if (gappedSeqN[s] != "-" || KEEP_ALL_GAPS_IN == 1) {
            idx = (s-shift);
            for (k=0; k< ref_ds.species; k+=1) {
                outSeqs[k] * (refSeqs[k])[idx];		
            }
        }
    }
}

for (k=0; k< ref_ds.species; k+=1) {
    outSeqs[k] * 0;		
}

gappedSeqN_Stripped * 0;
outputAlignment * 256;



for (k=0; k< ref_ds.species-1; k=k+1) {
    GetString (refName, ref_ds,k);
    outputAlignment *( ">"+ refName+ "\n"+ outSeqs[k]+ "\n");
}	

if (KEEP_ALL_GAPS_IN == 1) {
    alignedQuerySeq = gappedSeqN_Stripped;
    outputAlignment * (">" + qryName + "\n" + alignedQuerySeq[startFrom][endAt-1] + "\n" + refTree + "\n");	
}
else {		
    alignedQuerySeq = gappedSeqN_Stripped ^ {{"\\-"}{""}};
    outputAlignment * (">" + qryName + "\n" + alignedQuerySeq + "\n" + refTree + "\n");	
}

outputAlignment * 0;

//fprintf ("/Volumes/sergei-raid/Desktop/igscueal.fas", CLEAR_FILE, outputAlignment);
//assert (0);
//fprintf (stdout, outputAlignment, "\n");
	

/*---------------------------------------------
reverse complement a nucleotide string
---------------------------------------------*/


function nucleotideReverseComplement (seqIn)
{
	_seqOut = "";_seqOut*128;
	_seqL   = Abs(seqIn);
	for (_r = _seqL-1; _r >=0 ; _r = _r-1)
	{
		_seqOut *_nucleotide_rc[seqIn[_r]];
	}
	_seqOut*0;
	return _seqOut;
}

// -------------------------------------------------------------------------- //

function computeExpectedPerBaseScore( _expectedIdentity, _cdnaln_scorematrix, _cdnaln_base_freqs ) {
    meanScore = 0;

    for (_aa1 = 0; _aa1 < 20; _aa1 += 1) {
        for (_aa2 = 0; _aa2 < 20; _aa2 += 1) {
            if ( _aa1 != _aa2 ) {
                meanScore += ( 1 - _expectedIdentity ) * _cdnaln_scorematrix[_aa1][_aa2] * _cdnaln_base_freqs[_aa1] * _cdnaln_base_freqs[_aa2];
            } else {
                meanScore += _expectedIdentity * _cdnaln_scorematrix[_aa1][_aa1] * _cdnaln_base_freqs[_aa1];
            }
        }
    }

    return meanScore;
}


/*---------------------------------------------------------------------*/

function mapStrings (sourceStr,targetStr)
// source ID -> target ID (-1 means no correspondence)

{
	mapping 	  = {};
	targetIndexing = {};
	_d = Abs(targetStr);
	
	for (_i = 0; _i < _d; _i += 1)
	{
		targetIndexing [targetStr[_i]] = _i + 1;
	}
	_d = Abs (targetStr);
	for (_i = 0; _i < _d; _i += 1)
	{
		mapping [_i] = targetIndexing[sourceStr[_i]] - 1;
	}
	
	return mapping;
}

// -------------------------------------------------------------------------- //

function _igg_alignment_cleanup (reference, query, offset_nuc) {
    too_short = 0;
    too_long  = 0;
    span      = 0; // how many nucleotides in the reference were covered by non-gaps
    _seqL     = Abs (reference);
        
    ref_cleaned = ""; ref_cleaned * 128; 
    qry_cleaned = ""; qry_cleaned * 128;
    
    _codon_in_reference = 0;
    
	for ( _rcidx = 0; _rcidx < _seqL; _rcidx += 1 ) {
	    _del1 = reference [_rcidx] != (reference [_rcidx]&&1);
	    if (_del1) {
	        too_short += 1;
	        _codon_in_reference += 1;
	        ref_cleaned * (reference [_rcidx]&&1);
	        qry_cleaned * (query [_rcidx]&&1);
	    } else {
	        _del1 = query [_rcidx] != (query [_rcidx]&&1);
	        if (_del1) {
	            if (_seqL-_rcidx < 3 && _codon_in_reference % 3 == 0) {
	                break;
	            }
	            too_long += 1;
	        } else {
                ref_cleaned * (reference [_rcidx]&&1);
                qry_cleaned * (query [_rcidx]&&1);
                span += 1;	            
                _codon_in_reference +=1;
	        }
	    }
	}
	ref_cleaned * 0; qry_cleaned * 0;
	
	return {"REF": ref_cleaned, "QRY": qry_cleaned, "TOO_SHORT" : too_short, "TOO_LONG": too_long, "SPAN": span, "OFFSET_AA" :  offset_nuc$3 + (offset_nuc % 3 > 0),"OFFSET" :  offset_nuc, "AA" : translateToAA (qry_cleaned, (3-offset_nuc%3)%3), "AA_REF" : translateToAA (ref_cleaned, (3-offset_nuc%3)%3)};
}

/*-------------------------------------------------*/

function	computeCorrection (str)
{
	correctionFore	 = (str$"^\\-+")[1]+1;
	correctionAft	 = (str$"\\-+$")[0];
	if (correctionAft >= 0)
	{
		correctionAft = Abs(str)-correctionAft;
	}
	else
	{
		correctionAft = 0;
	}
	return {{correctionFore__,correctionAft__}};
}        

/*-------------------------------------------------*/

function makeAAMap ()
{
	codonToAAMap = {};
	codeToAA 	 = "FLIMVSPTAYXHQNKDECWRG";
	
	nucChars = "ACGT";
	
	for (p1=0; p1<64; p1=p1+1)
	{
		codon 				= nucChars[p1$16]+nucChars[p1%16$4]+nucChars[p1%4];
		ccode 				= _Genetic_Code[p1];
		codonToAAMap[codon] = codeToAA[ccode];
	}
	return codonToAAMap;
}

/*-------------------------------------------------*/

function translateToAA (aSeq, offset)
{
	seqLen	= Abs (aSeq)-2;
	translString = "";
	translString * (seqLen/3+1);
	for (seqPos = offset; seqPos < seqLen; seqPos = seqPos+3)
	{
		codon = aSeq[seqPos][seqPos+2];
		prot = codonToAAMap[codon];
		if (Abs(prot))
		{
			translString * prot;
		}
		else
		{
			translString * "?";
		}
	} 
	translString * 0;
	translString = translString^{{"X$","?"}};
	return translString;
}

/*-------------------------------------------------*/

function correctReadUsingCodonAlignedData (aligned1, aligned2) {
      
    alL = computeCorrection(aligned1);
            
    /*alL is the starting,ending nucleotide on the reference relative to the read. if reference is longer than the read, then both are 0*/
    
    offsetFrom = (aligned2$"^\\-+")[1]+1;
    offsetTo   = (aligned2$"\\-+$")[0]-1;
    
    /* the $ looks for the regular expression in bestAl[2] and returns a 2x1 array with the starting and ending 0-based positions of the regular expression. in this case multiple indels, -. returns -1 for both if the regular expression is not found. 
        i.e. 0-based index leading indels start at (bestAl[2]$"^\\-+")[0] and end at (bestAl[2]$"^\\-+")[1]; trailing indels start at (bestAl[2]$"\\-+$")[0] and end at (bestAl[2]$"\\-+$")[0];		
        
        so offSetFrom to offSetTo will return the reference sequence co-ordinates overlapping with the read.
    */
    
    
    if (offsetTo < 0) {
        offsetTo = Abs(aligned2)-1; /*if no trailing indels then to end of read*/
    }
    
    seqOffset  = offsetFrom;          /*set the offset of the read relative to the reference. ie the number of indels needed on the read to align to the reference */
    offsetFrom +=  alL[0];           /*if the read starts before the reference then shift to start of reference ie. by alL[0] */
    offsetTo    =  offsetTo	- alL[1];           /*if the read extends beyond the reference then shift to end of reference ie. by alL[1] */
    
    theSeq     = aligned2;
    theSeq	   = theSeq[alL[0]][Abs(theSeq)-alL[1]-1]; /*the nucleotide sequence of the read that overlaps with the reference sequence */
    
    
    nucSeq	   = (aligned2)[offsetFrom][offsetTo]; /*read sequence pruned to exactly overlapping region*/
    nucSeqRef  = (aligned1)[offsetFrom][offsetTo]; /*reference sequence pruned to exactly overlapping region*/
    
    LoadFunctionLibrary ("chooseGeneticCode", {"0":"Universal"});
    codonToAAMap = makeAAMap();
    
    return _igg_alignment_cleanup (nucSeqRef, nucSeq,seqOffset);
       
}
