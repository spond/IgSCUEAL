RequireVersion ("2.22");

if (MPI_NODE_COUNT <= 1) {
	fprintf (stdout, "[ERROR] This script requires MPI \n");
	return 0;
}

/* sequence indices being processed */
MPI_NODE_STATUS = {MPI_NODE_COUNT-1,1}; 

ExecuteAFile 			("../Configs/settings.ibf"):
	
SetDialogPrompt 		("A sequence file to screen:");
DataSet ds_in 			= ReadDataFile (PROMPT_FOR_FILE);
DataSetFilter ds_fil 	= CreateFilter (ds_in,1);
GetString				(sequenceNames, ds_fil, -1);


fprintf 						(stdout, "\nRead ", ds_in.species, " sequences\n");
SetDialogPrompt					("Write a tsv file to:");

/* check if the output file already exists */
fprintf 						(PROMPT_FOR_FILE,CLEAR_FILE,KEEP_OPEN);
resultsFile				= 		LAST_FILE_PATH;

SetDialogPrompt					("Write a tsv file with all significantly supported rearrangements for each sequence:");
fprintf 						(PROMPT_FOR_FILE,CLEAR_FILE,KEEP_OPEN);
rearrangementFile	   = 		LAST_FILE_PATH;

SetDialogPrompt					("Write a tsv file with all branch support values for all screened sequences:");
fprintf 						(PROMPT_FOR_FILE,CLEAR_FILE,KEEP_OPEN);
treeAssignmentFile	   = 		LAST_FILE_PATH;

headerString			= 		"Index\tName\tBest Rearrangement\tSupport\tSequence";

for (k = 0; k < Abs(_extraOutputColumns); k+=1) {	
	headerString += "\t"+_extraOutputColumns[k];
}
headerString += "\tV-length\tJ-length";


fprintf 						(resultsFile,headerString);
fprintf                         (rearrangementFile,  KEEP_OPEN, "Name\tRearrangement\tSupport")l
fprintf                         (treeAssignmentFile, KEEP_OPEN, "Name\tBranch\tSupport");

fprintf (stdout, "\n");

jobsFinished = 0;

for (seqID = 0; seqID < ds_in.species; seqID += 1) {
	SendAJob (seqID);
}

/* clean up MPI jobs */

howManyPending = 0;
for (mpiNode = 0; mpiNode < MPI_NODE_COUNT-1; mpiNode += 1){
	if (MPI_NODE_STATUS[mpiNode]) {
		howManyPending += 1;
	}
}

for (; howManyPending; howManyPending = howManyPending-1) {
	ReceiveAJob (0);
}

fprintf                         (rearrangementFile, CLOSE_FILE);
fprintf                         (resultsFile, CLOSE_FILE);
fprintf                         (treeAssignmentFile, CLOSE_FILE);

/*------------------------------------------------------------------------*/

function SendAJob (sequenceID)
{
	inOptions 	   = {};
	inOptions["0"] = sequenceNames[sequenceID];
	GetDataInfo (theSeq, ds_fil, sequenceID);
	inOptions["1"] = theSeq;
	inOptions["2"] = ""+alignmentType;
	
	for (mpiNode = 0; mpiNode < MPI_NODE_COUNT-1; mpiNode = mpiNode+1)
	{
		if (MPI_NODE_STATUS[mpiNode] == 0) /* free node */
		{
			break;
		}
	}
	if (mpiNode == MPI_NODE_COUNT-1) /* all busy */
	{
		mpiNode = ReceiveAJob (0);
	}
	fprintf (stdout, "[SEND] Sequence ", inOptions["0"], " to MPI node ", mpiNode + 1, "\n");
	MPI_NODE_STATUS [mpiNode] = sequenceID+1;
	MPISend (mpiNode+1, "../HBF/MPI_Wrapper.bf", inOptions);
	return 0;
}

/*------------------------------------------------------------------------*/

function ReceiveAJob (dummy)
{
	MPIReceive 		(-1, whichNode, returnValue);
	whichNode	  = whichNode-1;
	processedID   =	MPI_NODE_STATUS [whichNode]-1; 
	processedName = sequenceNames[processedID];
	
	MPI_NODE_STATUS [whichNode] = 0;
	ExecuteCommands	("returnAVL = " + returnValue);
	jobsFinished    = jobsFinished + 1;
	fprintf (stdout, "[RECEIVE] Sequence ", processedName, " from node ", whichNode + 1, " (", (ds_in.species-jobsFinished), " alignments remaining)");
	rearr_found = returnAVL["BEST_REARRANGEMENT"];
	if (Abs(rearr_found) == 0) /* error */
	{
		fprintf (stdout, ": Error/ alignment failed\n");
		fprintf (resultsFile, "\n", processedID+1, "\t", processedName, "\tError: alignment failed\t\t");
		if (Abs (_extraOutputColumns)) {
            for (k = 0; k < Abs(_extraOutputColumns); k+=1){	
                    fprintf (resultsFile, "\t");
            }
		}
	}
	else
	{
		fprintf (stdout, ": ", rearr_found, "\n");
		
		fprintf (resultsFile,   "\n", processedID+1, 
		                        "\t", processedName, 
		                        "\t", rearr_found, 
		                        "\t", returnAVL["SUPPORT"], 
		                        "\t", returnAVL["SEQUENCE"]);
		                        
	    (returnAVL["REARRANGEMENTS"])["_write_rearrangements"][""];
	    (returnAVL["BRANCH_SUPPORT"])["_write_branch_support"][""];
		
		extra			= returnAVL["EXTRA"];
		
		if (Abs(extra)) {
			for (k = 0; k < Abs(_extraOutputColumns); k+=1) {	
				fprintf (resultsFile, "\t", extra[_extraOutputColumns[k]]);
			}
		} 
		


		_bl_dict = returnAVL["BRANCH_LENGTHS"];
		_bl_vect = {{0,0}};
		if (Type (_bl_dict) == "AssociativeList") {
		    hasJ = Rows (_bl_dict["TO_SISTER"])* Columns (_bl_dict["TO_SISTER"]) > 1;
		    _bl_vect[0] = (_bl_dict["TO_SISTER"])[0] + (_bl_dict["TO_PARENT"])[0] + (_bl_dict["TO_QUERY"])[0];
		    if (hasJ) {
		        _bl_vect[1] = (_bl_dict["TO_SISTER"])[1] + (_bl_dict["TO_PARENT"])[1] + (_bl_dict["TO_QUERY"])[1];
		    }
		}
		fprintf (resultsFile, "\t", _bl_vect[0], "\t", _bl_vect[1]);
		
	}
	return 			whichNode;
}

function _write_rearrangements (key, value) {
    fprintf (rearrangementFile, "\n", processedName, "\t", key, "\t", value);
}

function _write_branch_support (key, value) {
    fprintf (treeAssignmentFile, "\n", processedName, "\t", key, "\t", value);
}

