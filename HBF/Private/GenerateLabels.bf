respect_allele = 1;

ExecuteAFile ("../utils.ibf");
LoadFunctionLibrary ("ReadDelimitedFiles");
SetDialogPrompt 				  ("Existing labels");
DataSet 		ds = ReadDataFile (PROMPT_FOR_FILE);
labels  = LAST_FILE_PATH + ".labels";
ExecuteAFile (labels);
sequenceLabels = _subtypeAssignmentByNode;

Tree exportT = DATAFILE_TREE;
//Tree exportT = "(((IGLV9_49_01_CRF_1:0.2994294106132673,(IGLV4_3_01_CRF_1:0.1542971064991972,(IGLV4_60_01_CRF_1:0.09032441876115867,IGLV4_69_01_CRF_1:0.03846677780948428)Node7:0.05779590223201311)Node5:0.07981574288480114)Node3:0.2045417061467103,(IGLV5_52_01_CRF_1:0.1088980942123382,(IGLV5_37_01_CRF_1:0.05765744990257021,(IGLV5_39_01_CRF_1:0.04539309336044683,IGLV5_45_01_CRF_1:0.01904218435059701)Node14:0.03838555863321481)Node12:0.05757716183100315)Node10:0.1998290982961171)Node2:0.139745793539219,(IGLV8_61_01_CRF_1:0.1906384181779668,(IGLV7_43_01_CRF_1:0.03172826581346387,IGLV7_46_01_CRF_1:0.05643764219208338)Node19:0.2050631902570158)Node17:0.2019409991532234,((IGLJ4_01_CRF_1:0.2398090318059009,(IGLJ7_02_CRF_1:0.02734042376614373,IGLJ7_01_CRF_1:1.153489754209784e-09)Node25:0.1464224299274031)Node23:0.08942682201190405,IGLJ2_01_CRF_1:0,(IGLJ3_02_CRF_1:0,((IGLJ6_01_CRF_1:0.09657655585406927,IGLJ1_01_CRF_1:0.1349577323835641)Node32:0.175029466943911,(IGLJ5_02_CRF_1:0.02887722652149212,IGLJ5_01_CRF_1:0)Node35:0.1509311660687902)Node31:0)Node29:0.08595430966004458)Node22:0.0001786117019070815,((IGLV10_54_01_CRF_1:0.2735019294027625,((IGLV1_40_01_CRF_1:0.05297251010162618,(IGLV1_51_01_CRF_1:0.09286102694421611,(IGLV1_47_01_CRF_1:0.009952886800500453,(IGLV1_36_01_CRF_1:0.06434270668278068,IGLV1_44_01_CRF_1:0)Node48:0.01449671342503711)Node46:0.06293971155153423)Node44:0.02418833867509778)Node42:0.06842461575079496,(IGLV6_57_01_CRF_1:0.2821250642911635,(IGLV2_11_01_CRF_1:0.01260309837309609,(IGLV2_8_01_CRF_1:0.01058317980084184,(IGLV2_18_01_CRF_1:0.02947893894305976,(IGLV2_14_01_CRF_1:0.004631231746300492,IGLV2_23_01_CRF_1:0.04918726065758569)Node59:0.02531163641918093)Node57:0.01649850203869457)Node55:0.008906051462981788)Node53:0.1560809209481158)Node51:0.03549888055060379)Node41:0.04901262905871992)Node39:0.03351561904409638,(IGLV3_19_01_CRF_1:0.165979404977528,((IGLV3_9_01_CRF_1:0.04188750406494332,(IGLV3_12_01_CRF_1:0.06033837636226477,IGLV3_21_01_CRF_1:0.03068320288521998)Node67:0.02245430006835932)Node65:0.06669968391698507,(IGLV3_1_01_CRF_1:0.07100792659364495,(IGLV3_22_01_CRF_1:0.141348829824659,(IGLV3_27_01_CRF_1:0.06985451753027111,(IGLV3_10_01_CRF_1:0.07997679455101768,(IGLV3_16_01_CRF_1:0.03207404562849189,IGLV3_25_01_CRF_1:0.0181895656078146)Node78:0.03105255238657365)Node76:0.04503788839217135)Node74:0.02953381874984161)Node72:0.04841997827621977)Node70:0.04357322937528683)Node64:0.05993984265275056)Node62:0.08670210616533451)Node38:0.09640814131946361)";

fprintf (stdout, 		"\nAuto-generating internal node labels\n"); 

tc = TipCount (exportT);
for (k=0; k < tc; k += 1) {
	nodeName 	  = TipName (exportT,k);
	subexp = extractAllExpressions (nodeName^ {{"^IG[A-Z]"}{""}}, "[^_]+", "");
	sequenceLabels [nodeName&&1] = subexp[0];
	
    for (k2 = 1; k2 < Abs (subexp) - 3; k2 += 1) {
        sequenceLabels [nodeName&&1] += "-" + subexp [k2];
    }
    
    if (respect_allele) {
        sequenceLabels [nodeName&&1] += "*" + subexp[Abs (subexp) - 3];
    }
	    
	fprintf (stdout, nodeName, "->", sequenceLabels [nodeName&&1], "\n");
}

_doLabelGeneration (respect_allele);

fprintf (stdout, sequenceLabels, "\n");
fprintf (labels, CLEAR_FILE, "_subtypeAssignmentByNode = \n", sequenceLabels, ";\n");

