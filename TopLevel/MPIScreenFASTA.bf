LoadFunctionLibrary ("libv3/convenience/regexp");

function            PadWithZeros (value) {
    if (value < 10) {
        return "0" + value;
    }
    return "" + value;
}

function            FormatTime (seconds) {
    return PadWithZeros (seconds $ 3600) + ":" + PadWithZeros (seconds % 3600 $ 60) + ":" + PadWithZeros (seconds % 60);
}

function            CheckAndLoadTSV (file_name, id_index, min_fields, header_line, cross_check) {
    cached_ids = {};

    if (!file_name) { // file exists
        fscanf (file_name, REWIND, "Lines", cached_lines);
        fprintf(file_name, CLEAR_FILE, KEEP_OPEN);
        if (Abs (header_line)) {
            fprintf(file_name, cached_lines[0], "\n");
        }
        for (line_id = (Abs (header_line) > 0); line_id < Columns (cached_lines); line_id += 1) {
            fields = regexp.Split (cached_lines[line_id], "\t");
            if (Abs (fields) >= min_fields) {
                if (Type (cross_check) == "AssociativeList") {
                    if (Abs (cross_check [fields[id_index]]) == 0) {
                        continue;
                    }
                }
                cached_ids [fields[id_index]] = 1;
                fprintf(file_name, cached_lines[line_id], "\n");
            }
        }
        DeleteObject (cached_lines);
        fprintf (stdout, "\nLoaded ", Abs (cached_ids), " cached records\n");

    } else {
        fprintf (stdout, "\nCreated file ", file_name, "\n");
        fprintf(file_name, CLEAR_FILE, KEEP_OPEN);
        if (Abs (header_line)) {
            fprintf(file_name, header_line, "\n");
        }
    }
    return cached_ids;
}

RequireVersion ("2.26");

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


headerString			= 		"Index\tName\tBest Rearrangement\tSupport\tSequence";

for (k = 0; k < Abs(_extraOutputColumns); k+=1) {
	headerString += "\t"+_extraOutputColumns[k];
}
headerString += "\tV-length\tJ-length";

fprintf                         (stdout, "Write the TSV main file to: ");
fscanf                          (stdin, "String", resultsFile);
skip_these = CheckAndLoadTSV (resultsFile, 1, 5 + Abs (_extraOutputColumns), headerString, None);

fprintf                         (stdout, "Write a tsv file with all significantly supported rearrangements for each sequence to:");
fscanf                          (stdin, "String", rearrangementFile);
CheckAndLoadTSV (rearrangementFile, 0, 3, "Name\tRearrangement\tSupport", skip_these);

fprintf                         (stdout, "Write a tsv file with all branch support values for all screened sequences: ");
fscanf                          (stdin, "String", treeAssignmentFile);
CheckAndLoadTSV (treeAssignmentFile, 0, 3, "Name\tBranch\tSupport", skip_these);

fprintf (stdout, "\n");

jobsFinished = 0;

start_time = Time (1);

for (seqID = 0; seqID < ds_in.species; seqID += 1) {
    if (skip_these[sequenceNames[seqID]] == 0) {
        SendAJob (seqID);
    } else {
        jobsFinished += 1;
    }
}

return 0;

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
	fprintf (stdout, "[RECEIVE] Sequence ", processedName, " from node ", whichNode + 1, " (", (ds_in.species-jobsFinished), " alignments remaining, ",
	        FormatTime ((Time(1) - start_time) * ((ds_in.species-jobsFinished)/jobsFinished)), " time remaining)");
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

