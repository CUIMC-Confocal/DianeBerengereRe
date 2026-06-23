//@File(label = "Input directory", style = "directory") inputDir
//@File(label = "Output directory", style = "directory") outputDir
//@String (label = "File suffix", value = ".nd2") fileSuffix
//@String (label = "File name contains:", value = "MaxIP") containString

// PFC_BBB_damage_batch.ijm
// ImageJ/Fiji script to process a batch of images to quantify BBB damage
// Segments vessels using channel 3, and measures intensity of channel 2 within the vessel area

// Theresa Swayne, 2026
//  -------- Suggested text for acknowledgement -----------
//   "These studies used the services of the Confocal and Specialized Microscopy Shared Resource 
//   of the Herbert Irving Comprehensive Cancer Center at Columbia University, 
//   funded in part through the NIH/NCI Cancer Center Support Grant P30CA013696."

// TO USE: Place all input images in the input folder.
// 	Create a folder for the output files. 
//  Run the script in Fiji. 
//	Limitations -- cannot have >1 dots in the filename
//		Assumes channels in order: 1) DAPI, 2) green (channel of interest), 3) red (vessel marker)


// ---- Setup ----

while (nImages>0) { // clean up open images
	selectImage(nImages);
	close();
}
run("Clear Results");
print("\\Clear"); // clear Log window
roiManager("reset");
setBatchMode(true); // faster performance
run("Bio-Formats Macro Extensions"); // support native microscope files
run("Set Measurements...", "area mean modal integrated display redirect=None decimal=3");
run("Input/Output...", "jpeg=85 gif=-1 file=.csv copy_row save_column save_row");

// get date and time for timestamped results
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
startTime = getTime();
month = month+1;
timeString = "" + year + "-" + month + "-" + dayOfMonth + "-" + hour + "-" + minute; // must start with an empty string

// get time for runtime measurement
startTime = getTime();

// ---- Run ----

print("Starting");

// Set up a results file
resultsName = timeString + "_results.csv";
resultsFile = outputDir + File.separator + resultsName;
resultsHeader = "Label,Area,Mean,Mode,IntDen,RawIntDen,WholeImageMode";
if (File.exists(resultsFile)==false) { // start the file with headers
	File.append(resultsHeader, resultsFile);	
	print("Created results file");
    }
    
// Call the processFolder function
processFolder(inputDir, outputDir, fileSuffix, containString, resultsFile);


// Clean up images and get out of batch mode
while (nImages > 0) { // clean up open images
	selectImage(nImages);
	close(); 
}
setBatchMode(false);
run("Clear Results");
elapsedTime = (getTime() - startTime)/1000;
print("Finished in",elapsedTime,"seconds"); 

// ---- Functions ----

function processFolder(input, output, suffix, contain, resultsFile) {

	// this function searches folders for files matching the criteria and sends them to the processFile function
	
	filenum = 0;
	print("Processing folder", input);
	testString = ".*"+contain+".*";
	
	// scan folder tree to find files with correct names
	
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i])) {
			processFolder(input + File.separator + list[i], output, suffix, contain, resultsFile); // handles nested folders
		}
		if(endsWith(list[i], suffix)) { // check for suffix
			if(matches(list[i], testString)) { // check for another string in the filename
				filenum = filenum + 1;
				processFile(input, output, list[i], filenum, resultsFile); // passes the filename and parameters to the processFile function
			}
		}
	}
	
} // end of processFolder function

function processFile(inputFolder, outputFolder, fileName, fileNumber, resultsFile) {
	
	// this function processes a single image
	
	path = inputFolder + File.separator + fileName;
	print("Processing file number",fileNumber," at path" ,path);	

	// determine the name of the file without extension
	dotIndex = lastIndexOf(fileName, ".");
	basename = substring(fileName, 0, dotIndex); 
	extension = substring(fileName, dotIndex);
	
	//print("File basename is",basename);
	
	// open the file
	run("Bio-Formats", "open=&path");
	run("Split Channels");
	//print("After opening image",fileNumber," we have", nResults, "results");

	// ---- Define vessel area using channel 3 ----
	
	selectWindow("C3-" + fileName);
	copyName = "C3-" + fileName + "_copy";
	run("Duplicate...", "title="+copyName);
	selectWindow(copyName);
	run("Gaussian Blur...", "sigma=2"); // smooth the image to get a cleaner segmentation
	run("Auto Threshold", "method=Huang white"); 
	run("Analyze Particles...", "size=20.00-Infinity show=Nothing clear add"); // excludes objects < 20 µm2 and adds to ROI manager

	// combine all the vessel parts into one ROI
	numRois = roiManager("count");
	roiIndexes = Array.getSequence(numRois); // 0 to numROIs-1
	roiManager("deselect");
	roiManager("select", roiIndexes);
	roiManager("combine");
	roiManager("add"); // add the combined ROI to the manager
	roiManager("deselect");
	roiManager("show none");
	roiManager("select", numRois); // this is the last ROI
	
	roisToDelete = Array.getSequence(numRois);
	roiManager("select", roisToDelete);
	roiManager("delete"); // delete the individual vessel "particles" and keep the combined one
	Overlay.remove; // clear out the overlay from the other ROIs

	// ---- Measure intensity in channel 2 ----
	
	selectWindow("C2-"+fileName);
	
	// measure whole image; we will approximate the background value as the Mode 
	roiManager("deselect");
	run("Select All");
	run("Measure");

	mode = getResult("Mode", nResults-1); // we'll add this to the results table later
	print("Mode of the image is",mode);

	// measure the vessel area; use for measuring junction intensity
	roiManager("select", 0); // indices start at 0 and we should have only 1 ROI now
	run("Measure");
		
	// add the image mode to the results
	setResult("WholeImageMode", nResults-1, mode);
	
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
		
	// ---- Save data ----
	
	// append to the results file
	File.append(resultString, resultsFile);	

	// save ROI and snapshot for checking segmentation ----	
	roiName = basename + "_vessel_ROI.roi";
	roiManager("save", outputFolder +  File.separator + roiName);
	snapName = basename + "_snapshot.png";
	selectWindow("C3-"+fileName);
	run("Enhance Contrast", "saturated=0.35"); // Same as Auto B&C -- for visualization
	roiManager("show none");
	roiManager("show all without labels");
	//roiManager("select", 0); // this is the last ROI
	run("Flatten"); // creates an RGB image including the ROI line
	selectWindow("C3-"+basename + "-1" + extension);
	saveAs("png", outputFolder + File.separator + snapName);

	// ---- Clean up ----
	while (nImages>0) { // clean up open images
		selectImage(nImages);
	close();
	}
	roiManager("reset");
	
	//print("After closing image",fileNumber," we have", nResults, "results");
	
} // end of processFile function


	