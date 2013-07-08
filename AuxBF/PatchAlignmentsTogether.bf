LoadFunctionLibrary ("TreeTools");
LoadFunctionLibrary  ("GrabBag");
SetDialogPrompt     ("List of files");
fscanf              (PROMPT_FOR_FILE, "Lines", fileList);

fileCount           = Columns (fileList);
seqLengths          = {fileCount,1};

SetDialogPrompt ("Write to this file:");
fprintf (PROMPT_FOR_FILE, CLEAR_FILE, KEEP_OPEN);
saveTo = LAST_FILE_PATH;

for (fileID = 0; fileID < fileCount; fileID += 1) {
    ExecuteCommands ("DataSet component_" + fileID + " = ReadDataFile( " + fileList [fileID] + ");");
    //fprintf (stdout, "Read ", Eval ("component_" + fileID +".species"), " sequences from file ", fileList [fileID], "\n");
    seqLengths [fileID] = Eval ("component_" + fileID +".sites");
}

jointTreeString = ""; jointTreeString * 128;

for (fileID = 0; fileID < fileCount; fileID += 1) {
    ExecuteCommands ("DataSetFilter thisFilter = CreateFilter (component_" + fileID + ", 1);");
    if (fileID > 0) {
        prefixPad = returnPadder (+seqLengths[{{0,0}}][{{fileID-1,0}}]);
    } else {
        prefixPad = "";
    }
    if (fileID < fileCount - 1) {
        suffixPad = returnPadder (+seqLengths[{{fileID+1,0}}][{{fileCount-1,0}}]);
    } else {
        suffixPad = "";
    }
    
    for (seqID = 0; seqID < thisFilter.species; seqID += 1) {
        GetString (seqName, thisFilter, seqID);
        GetDataInfo (seqData, thisFilter, seqID);
        if (thisFilter.species < 4) {
            jointTreeString * ("," + seqName);
        }
        padded = ">`seqName`\n`prefixPad``seqData``suffixPad`\n";
        fprintf (saveTo, padded);
    }
    if (thisFilter.species >= 4) {
        pathInfo = splitFilePath (fileList [fileID]);
        pathInfo = pathInfo["DIRECTORY"] + DIRECTORY_SEPARATOR + pathInfo["FILENAME"] + ".tree";
        fscanf (pathInfo, "String", treeString);
        jointTreeString * ("," + treeString);
    }
    
 }

jointTreeString * 0;
fprintf (saveTo, "\n(", jointTreeString[1][Abs(jointTreeString)-1], ")", CLOSE_FILE);


function returnPadder (length) {
    padder = ""; padder * 128;
    
    for (gapc = 0; gapc < length; gapc += 1) {
        padder * "-";
    }
    
    padder * 0;
    return padder;
}

