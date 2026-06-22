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

// 	

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
	
// get date and time for timestamped results
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
startTime = getTime();
month = month+1;
timeString = "" + year + "-" + month + "-" + dayOfMonth + "-" + hour + "-" + minute; // must start with an empty string

// get time for runtime measurement
startTime = getTime();

// ---- Run ----

print("Starting");

// Call the processFolder function
processFolder(inputDir, outputDir, fileSuffix, containString);

// Save results
resultsName = timeString + "_results.csv";
saveAs("Results", outputDir + File.separator + resultsName);

// Clean up images and get out of batch mode
while (nImages > 0) { // clean up open images
	selectImage(nImages);
	close(); 
}
setBatchMode(false);
elapsedTime = (getTime() - startTime)/1000;
print("Finished in",elapsedTime,"seconds"); 

// ---- Functions ----

function processFolder(input, output, suffix, contain) {

	// this function searches folders for files matching the criteria and sends them to the processFile function
	
	filenum = 0;
	print("Processing folder", input);
	testString = ".*"+contain+".*";
	
	// scan folder tree to find files with correct names
	
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i])) {
			processFolder(input + File.separator + list[i], output, suffix); // handles nested folders
		}
		if(endsWith(list[i], suffix)) { // check for suffix
			if(matches(list[i], testString)) { // check for another string in the filename
				filenum = filenum + 1;
				processFile(input, output, list[i], filenum); // passes the filename and parameters to the processFile function
			}
		}
	}
	
} // end of processFolder function

function processFile(inputFolder, outputFolder, fileName, fileNumber) {
	
	// this function processes a single image
	
	path = inputFolder + File.separator + fileName;
	print("Processing file",fileNumber," at path" ,path);	

	// determine the name of the file without extension
	dotIndex = lastIndexOf(fileName, ".");
	basename = substring(fileName, 0, dotIndex); 
	extension = substring(fileName, dotIndex);
	
	print("File basename is",basename);
	
	// open the file
	run("Bio-Formats", "open=&path");
	run("Split Channels");
	print("After opening image",fileNumber," we have", nResults, "results");

	// ---- Define vessel area using channel 3 ----
	
	selectWindow("C3-"+fileName);
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
	roiManager("select", numRois); // this is the last ROI
	
	roisToDelete = Array.getSequence(numRois);
	roiManager("select", roisToDelete);
	roiManager("delete"); // delete the individual vessel "particles" and keep the combined one

	// ---- Measure intensity in channel 2 ----
	
	selectWindow("C2-"+fileName);
	
	// measure whole image; we will approximate the background value as the Mode 
	roiManager("deselect");
	run("Select All");
	run("Measure");
	print("After measuring whole image",fileNumber," we have", nResults, "results");
	mode = getResult("Mode", nResults-1); // we'll add this to the results table later
	print("Mode of the image is",mode);
	//IJ.deleteRows( nResults-1, nResults-1 ); // delete the last row of the results table while preserving the rest
	//print("After deleting whole image data from image",fileNumber," we have", nResults, "results");
	
	// measure the vessel area; use for measuring junction intensity
	roiManager("select", 0); // indices start at 0 and we should have only 1 ROI now
	run("Measure");
	print("After measuring vessels in image",fileNumber,"we have", nResults, "results");
		
	// add the image mode to the results
	setResult("WholeImageMode", nResults-1, mode);
	
	// ---- Save ROI and snapshot for checking segmentation ----
	
	roiName = basename + "_vessel_ROI.roi";
	roiManager("save", outputFolder +  File.separator + roiName);
	snapName = basename + "_snapshot.png";
	selectWindow("C3-"+fileName);
	roiManager("show none");
	roiManager("select", 0); // this is the last ROI
	run("Flatten");
	selectWindow("C3-"+basename + "-1" + extension);
	saveAs("png", outputFolder + File.separator + snapName);

	// ---- Clean up ----
	while (nImages>0) { // clean up open images
		selectImage(nImages);
	close();
	}

	print("After closing image",fileNumber," we have", nResults, "results");
	
} // end of processFile function


	