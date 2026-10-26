//@File(label = "Input directory", style = "directory") inputDir
//@File(label = "Output directory", style = "directory") outputDir
//@String (label = "File suffix", value = ".nd2") fileSuffix
//@String (label = "File name contains:", value = "MaxIP") containString
//@Integer (label = "Channel for vessels/CD31:", value = 1, style = "spinner") vesselChannel
//@Integer (label = "Channel for measuring Claudin5:", value = 2, style = "spinner") measureChannel

// PFC_BBB_damage_batch.ijm
// ImageJ/Fiji script to process a batch of images to quantify BBB damage
// Segments vessels and measures intensity of another channel within the vessel area
// Tested on Ch1 = green = CD31, 2 = red = Cld5

// Theresa Swayne, 2026
//  -------- Suggested text for acknowledgement -----------
//   "These studies used the services of the Confocal and Specialized Microscopy Shared Resource 
//   of the Herbert Irving Comprehensive Cancer Center at Columbia University, 
//   funded in part through the NIH/NCI Cancer Center Support Grant P30CA013696."

// TO USE: Place all input images in the input folder.
// 	Create a folder for the output files. 
//  Run the script in Fiji. 
//	Limitations -- cannot have >1 dots in the filename


// TODO: 1) CD31 (vessel) mean intensity, 2) Cld5 +Area and 3) colocalization of Cl5/CD31 (Cl5 +fraction of CD31 +vascular area). 
// need to define what is cld5+
// All measurements are for the CD31+ area except the Cld5_Area and Whole Image modes.
// TODO: Fix headers and measrements so they make sense 

// ---- Setup ----

while (nImages>0) { // clean up open images
	selectImage(nImages);
	close();
}
run("Clear Results");
print("\\Clear"); // clear Log window
roiManager("reset");
//setBatchMode(true); // faster performance
run("Bio-Formats Macro Extensions"); // support native microscope files
run("Set Measurements...", "area mean modal integrated display redirect=None decimal=3");
run("Input/Output...", "jpeg=85 gif=-1 file=.csv copy_row save_column save_row");

run("Options...", "iterations=1 count=1 black"); // ensure that black is interpreted as background


// get date and time for timestamped results
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
startTime = getTime();
month = month+1;
timeString = "" + year + "-" + month + "-" + dayOfMonth + "-" + hour + "-" + minute; // must start with an empty string

// get time for runtime measurement
startTime = getTime();

// ---- Run ----

print("Processing images in",inputDir);
print("Detecting vessels in channel",vesselChannel,"and additionally measuring channel",measureChannel);

// Set up a results file
resultsName = timeString + "_results.csv";
resultsFile = outputDir + File.separator + resultsName;
resultsHeader = "Cld5_Label,Vessel_Area,Cld5_Mean,Cld5_Mode,Cld5_IntDen,Cld5_RawIntDen,Cld5_WholeImageMode,CD31_Label,CD31_Area, CD31_Mean,CD31_Mode,CD31_IntDen,CD31_RawIntDen,CD31_WholeImageMode,Cld5_Area,Cld5_AreaWithinCD31";
if (File.exists(resultsFile)==false) { // start the file with headers
	File.append(resultsHeader, resultsFile);	
	print("Created results file");
    }
    
// Call the processFolder function
filenum = processFolder(inputDir, outputDir, fileSuffix, containString, vesselChannel, measureChannel, resultsFile);

// Clean up images and get out of batch mode
while (nImages > 0) { // clean up open images
	selectImage(nImages);
	close(); 
}
roiManager("reset");
setBatchMode(false);
run("Clear Results");
elapsedTime = (getTime() - startTime)/1000;
print("Finished",filenum,"files in",elapsedTime,"seconds"); 

// save Log
selectWindow("Log");
saveAs("text", outputDir + File.separator + "Analysis_Log.txt");


// ---- Functions ----

function processFolder(input, output, suffix, contain, vesselChannel, measureChannel, resultsFile) {

	// this function searches folders for files matching the criteria and sends them to the processFile function
	
	filenum = 0;
	print("Processing folder", input);
	testString = ".*"+contain+".*";
	
	// scan folder tree to find files with correct names
	
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i])) {
			processFolder(input + File.separator + list[i], output, suffix, contain, vesselChannel, measureChannel, resultsFile); // handles nested folders
		}
		if(endsWith(list[i], suffix)) { // check for suffix
			if(matches(list[i], testString)) { // check for another string in the filename
				//print("File", list[i], "is a match!");
				filenum = filenum + 1;
				processFile(input, output, list[i], filenum, vesselChannel, measureChannel, resultsFile); // passes the filename and parameters to the processFile function
			}
		}
	}
	
	return filenum;
} // end of processFolder function

function processFile(inputFolder, outputFolder, fileName, fileNumber, vesselChannel, measureChannel, resultsFile) {
	
	// this function processes a single image
	
	path = inputFolder + File.separator + fileName;
	print("Processing file number",fileNumber," at path" ,path);	

	// determine the name of the file without extension
	dotIndex = lastIndexOf(fileName, ".");
	basename = substring(fileName, 0, dotIndex); 
	extension = substring(fileName, dotIndex);
	
	print("File name is",fileName,"and basename is",basename);
	
	// open the file
	run("Bio-Formats", "open=&path");
	run("Split Channels");
	//print("After opening image",fileNumber," we have", nResults, "results");

	// ---- Define vessel area using the selected channel ----
	
	vesselName = "C"+vesselChannel+"-"+fileName;
	vesCopyName = vesselName + "_copy";
	selectWindow(vesselName);
	run("Duplicate...", "title=["+vesCopyName+"]");
	selectWindow(vesCopyName);
	run("Gaussian Blur...", "sigma=2"); // smooth the image to get a cleaner segmentation
	run("Auto Threshold", "method=Huang white"); 
	run("Analyze Particles...", "size=20.00-Infinity show=Nothing clear add"); // excludes objects < 20 µm2 and adds to ROI manager

	// combine all the vessel parts into one ROI
	numRois = roiManager("count");
	print("After detecting the vessels we have",numRois,"ROIs");
	roiIndexes = Array.getSequence(numRois); // 0 to numROIs-1
	roiManager("deselect");
	roiManager("select", roiIndexes);
	roiManager("combine");
	roiManager("add"); // add the combined ROI to the manager
	roiManager("deselect");
	roiManager("select", numRois); // this is the last ROI
	rName = Roi.getName(); 
	print("The name of the last ROI is", rName);
	roiManager("rename", "vessels");
	roiManager("deselect");
	roiManager("select", numRois); // this is the last ROI
	rName = Roi.getName(); 
	print("The name of the last ROI is", rName);
	
	// delete the smaller ROIs
	numRois = roiManager("count");
	roisToDelete = Array.getSequence(numRois);
	roisToDelete = Array.slice(roisToDelete,0,roisToDelete.length-1); // don't select the last ROI which is the vessels
	min = roisToDelete[0];
	max = roisToDelete[roisToDelete.length - 1];
	print("We now have",numRois,"ROIs, and will be deleting from",min,"to", max);

	roiManager("select", roisToDelete);
	roiManager("delete"); // delete the individual vessel "particles" and keep the combined one
	Overlay.remove; // clear out the overlay from the other ROIs
	
	numRois = roiManager("count");
	print("After deleting we now have",numRois,"ROIs");
	roiManager("select", numRois-1); // this is the last ROI
	rName = Roi.getName(); 
	print("The name of the last ROI is", rName);

	// ---- Measure intensity in Cld5 channel ----
	
	measureName = "C"+measureChannel+"-"+fileName;
	selectWindow(measureName);
	
	// measure whole image; we will approximate the background value as the Mode 
	roiManager("deselect");
	run("Select All");
	run("Measure");

	cld5mode = getResult("Mode", nResults-1); // we'll add this to the results table later
	print("Mode of the Cld5 image is",cld5mode);

	// measure the vessel area; use for measuring junction intensity
	roiManager("select", 0); // indices start at 0 and we should have only 1 ROI now
	run("Measure");
		
	// add the image mode to the results
	setResult("Cld5_WholeImageMode", nResults-1, cld5mode);
	
	// get the results as comma-separated values
	headings = split(String.getResultsHeadings);
	row = nResults - 1;
	label = getResultString("Label",row); // required for non-numeric columns
	//print("Label is",label);
	resultString = label;
	for (a=1; a<lengthOf(headings); a++) {
		head = headings[a];
		val = getResult(head, row);
		//print("Retrieving heading",head,"with result", val);
	    resultString = resultString + "," + val;
	    //print("Result string is",resultString);
	}
	// ---- Measure intensity in CD31 channel ----
	
	//run("Clear Results");
	selectWindow(vesselName);
	
	// measure whole image; we will approximate the background value as the Mode 
	roiManager("deselect");
	run("Select All");
	run("Measure");

	cd31mode = getResult("Mode", nResults-1); // we'll add this to the results table later
	print("Mode of the CD31 image is",cd31mode);

	// measure the vessels
	selectWindow(vesselName);
	roiManager("select", 0); // indices start at 0 and we should have only 1 ROI now
	run("Measure");
		
	// add the image mode to the results
	setResult("CD31_WholeImageMode", nResults-1, cd31mode);
	
	// get the results as comma-separated values
	headings = split(String.getResultsHeadings);
	row = nResults - 1;
	label = getResultString("Label",row); // required for non-numeric columns
	//print("Label is",label);
	resultString = resultString + "," + label;
	for (a=1; a<lengthOf(headings); a++) {
		head = headings[a];
		val = getResult(head, row);
		//print("Retrieving heading",head,"with result", val);
	    resultString = resultString + "," + val;
	    //print("Result string is",resultString);
	
	}
	
	// ---- Measure CLD5+ area ----
	
	// Threshold CLD5

	measCopyName = measureName + "_copy";
	selectWindow(measureName);
	run("Duplicate...", "title=["+measCopyName+"]");
	selectWindow(measCopyName);
	run("Gaussian Blur...", "sigma=2"); // smooth the image to get a cleaner segmentation
	run("Auto Threshold", "method=Huang white"); 
	run("Analyze Particles...", "size=20.00-Infinity show=Nothing add"); // excludes objects < 20 µm2 and adds to ROI manager

	// combine all the CLD5+ parts into one ROI
	numRois = roiManager("count");
	print("After detecting Cld5 we now have",numRois,"ROIs");
	roisToCombine = Array.getSequence(numRois); // 0 to numROIs-1
	roisToCombine = Array.slice(roisToCombine,1,numRois-1); // don't select the 1st ROI which is the vessels
	roiManager("deselect");
	roiManager("select", roisToCombine);
	roiManager("combine");
	roiManager("add"); // add the combined ROI to the manager
	roiManager("deselect");
	roiManager("select", numRois); // this is the last ROI
	rName = Roi.getName(); 
	print("The name of the last ROI is", rName);
	roiManager("rename", "Cld5");
	roiManager("select", numRois); // this is the last ROI
	rName = Roi.getName(); 
	print("The name of the last ROI is", rName);

	numRois = roiManager("count"); // after adding the combined ROI for Cld5
	roisToDelete = Array.getSequence(numRois);
	roisToDelete = Array.slice(roisToDelete,1,numRois-1); // don't select the 1st ROI which is the vessels
	min = roisToDelete[0];
	max = roisToDelete[roisToDelete.length - 1];
	print("We now have",numRois,"ROIs, and will be deleting from",min,"to", max);
	roiManager("select", roisToDelete);
	roiManager("delete"); // delete the individual Cld5 "particles" and keep the combined one
	Overlay.remove; // clear out the overlay from the other ROIs
	
	//roiManager("select", 1); // indices start at 0 and we should have only 2 ROIs now


	// ---- Measure overlapped area (fraction of CD31 that is Cld5+) ----
	
	// TODO:
	// Do AND with ROIs
	// Measure area (getStatistics is fine)
	// Append area to results string
	
	// ---- Save data ----
	
	// append to the results file
	File.append(resultString, resultsFile);	

	// save ROI and snapshot for checking segmentation ----	
	roiManager("deselect");
	run("Select None");
	roiName = basename + "_ROIs.zip";
	roiManager("save", outputFolder +  File.separator + roiName);
	snapName = basename + "_snapshot.png";
	selectWindow(vesselName);
	run("Enhance Contrast", "saturated=0.35"); // Same as Auto B&C -- for visualization
	roiManager("show none");
	roiManager("show all without labels");
	//roiManager("select", 0); // this is the last ROI
	run("Flatten"); // creates an RGB image including the ROI line
	flatName = "C"+vesselChannel+"-"+basename + "-1" + extension;
	selectWindow(flatName);
	saveAs("png", outputFolder + File.separator + snapName);

	// ---- Clean up ----
	while (nImages>0) { // clean up open images
		selectImage(nImages);
	close();
	}
	roiManager("reset");
	
	//print("After closing image",fileNumber," we have", nResults, "results");
	
} // end of processFile function

function thresholdImage(imageName, roiName) {
	// thresholds an image, analyzes particles, combines these into 1 ROI
	
	numRois = roiManager("count");
	print("Before detecting",roiName," we have",numRois,"ROIs");
	
	copyName = imageName + "_copy";
	selectWindow(imageName);
	run("Duplicate...", "title=[" + copyName + "]"); // brackets avoid issues with complicated filenames
	selectWindow(copyName);
	run("Gaussian Blur...", "sigma=2"); // smooth the image to get a cleaner segmentation
	run("Auto Threshold", "method=Huang white"); 
	run("Analyze Particles...", "size=20.00-Infinity show=Nothing add"); // excludes objects < 20 µm2 and adds to ROI manager

	// combine all the particles into one ROI
	numRois = roiManager("count");
	print("After detecting",roiName," we have",numRois,"ROIs");
	
	// find the index of the named ROI
	
	
	roiIndexes = Array.getSequence(numRois); // 0 to numROIs-1
	roiManager("deselect");
	roiManager("select", roiIndexes);
	roiManager("combine");
	roiManager("add"); // add the combined ROI to the manager
	roiManager("deselect");
	roiManager("select", numRois); // this is the last ROI
	rName = Roi.getName(); 
	print("The name of the last ROI is", rName);
	roiManager("rename", "vessels");
	roiManager("deselect");
	roiManager("select", numRois); // this is the last ROI
	rName = Roi.getName(); 
	print("The name of the last ROI is", rName);
	
	// delete the smaller ROIs
	numRois = roiManager("count");
	roisToDelete = Array.getSequence(numRois);
	roisToDelete = Array.slice(roisToDelete,0,roisToDelete.length-2); // don't select the last ROI which is the vessels
	min = roisToDelete[0];
	max = roisToDelete[roisToDelete.length - 1];
	print("We now have",numRois,"ROIs, and will be deleting from",min,"to", max);

	roiManager("select", roisToDelete);
	roiManager("delete"); // delete the individual vessel "particles" and keep the combined one
	Overlay.remove; // clear out the overlay from the other ROIs
	
	numRois = roiManager("count");
	print("After deleting we now have",numRois,"ROIs");
	roiManager("select", numRois-1); // this is the last ROI
	rName = Roi.getName(); 
	print("The name of the last ROI is", rName);
}

function findRoiWithName(roiName) { // thanks to Olivier Burri https://forum.image.sc/t/selecting-roi-based-on-name/3809/2
	
	nR = roiManager("Count"); 
 
	for (i=0; i<nR; i++) { 
		roiManager("Select", i); 
		rName = Roi.getName(); 
		if (matches(rName, roiName)) { 
			return i; 
		} 
	} 
	return -1; 
}