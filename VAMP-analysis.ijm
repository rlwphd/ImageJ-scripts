function boutonCheck(boutonImage, thresh, boutonType, headers, pathSave, nameBase, picSave) { 
// Checking to make sure that there are boutons present for each threshold, if not creating a zero result
	selectWindow(boutonImage);
	run("Duplicate...", "title=check duplicate");
	run("8-bit");
	valueThreshold = thresh/4095*255;
	setThreshold(valueThreshold, 255);
	run("Convert to Mask", "method=Default background=Dark black");
	run("Analyze Particles...", "size=0-infinity circularity=0.00-1.00 show=Nothing clear stack");
	numOfparts = nResults;
	close("check");

	dataList = newArray("Volume (micron^3)", "Nb of obj. voxels", "Nb of surf. voxels", "IntDen", "Mean", "StdDev", "Median", "Min", "Max", "XM", "YM", "ZM", "BX", "BY", "BZ", "B-width", "B-height", "B-depth");
	// If there are no particles then create a Results Table of 0s
	if (numOfparts == 0) {
		Table.create("Mouse-Results");
		for (i = 0; i < dataList.length; i++) {
			Table.set(dataList[i], numOfparts, 0);
		}
		Table.set("Threshold", numOfparts, thresh);
		Table.set("Data Type", numOfparts, "Empty");
		Table.set("Category", numOfparts, headers[0]);
		Table.set("Mouse", numOfparts, headers[1]);
		Table.set("Synapse", numOfparts, boutonType[1]);
		Table.set("Vamp", numOfparts, boutonType[2]);
		
	} else {
		// Otherwise analyze the boutons
			//print("analyzing the boutons");
			// Running bouton analysis by finding boutons in FITC channel and measuring values in the FITC and Cy3 channels
			if (thresh == 1500) {
				initialize = true;
			} else {
				initialize = false;
			}
			boutonAnalysis(boutonImage, "ch2", thresh, boutonType[1], boutonType, headers, dataList, picSave, pathSave, nameBase, initialize);
			boutonAnalysis(boutonImage, "ch3", thresh, boutonType[2], boutonType, headers, dataList, picSave, pathSave, nameBase, false);

			if (initialize) {
				// Getting the background data for determining cutoff values
				selectWindow("ch2");
				run("Flip Horizontally", "stack");
				run("Flip Vertically", "stack");
				background = boutonType[1]+" Background";
				boutonAnalysis(boutonImage, "ch2", thresh, background, boutonType, headers, dataList, false, pathSave, nameBase, false);
				selectWindow("ch2");
				run("Flip Horizontally", "stack");
				run("Flip Vertically", "stack");
			
				selectWindow("ch3");
				run("Flip Horizontally", "stack");
				run("Flip Vertically", "stack");
				background = boutonType[2]+" Background";
				boutonAnalysis(boutonImage, "ch3", thresh, background, boutonType, headers, dataList, false, pathSave, nameBase, false);
				selectWindow("ch3");
				run("Flip Horizontally", "stack");
				run("Flip Vertically", "stack");
			}
	}
}


// bouton analysis function
// redirectTitle is the name of the channel to analyze, thresh is the threshold value for the 3D Objects Counter
// boutonType is a text for specifying the name of the bouton, imgBool is a T/F for if you want to save the overlayed image
// pathSave is a passthrough of savePath, nameBase is baseName passthrough
function boutonAnalysis(boutonImage, redirectTitle, thresh, dataType, boutonType, headers, colNames, imgBool, pathSave, nameBase, initialize) {
	selectWindow(boutonImage);
	run("Select None");
	run("3D OC Options", "volume nb_of_obj._voxels nb_of_surf._voxels integrated_density mean_gray_value std_dev_gray_value median_gray_value minimum_gray_value maximum_gray_value centre_of_mass bounding_box show_masked_image_(redirection_requiered) dots_size=5 font_size=10 store_results_within_a_table_named_after_the_image_(macro_friendly) redirect_to="+redirectTitle);
	run("3D Objects Counter", "threshold="+thresh+" slice=1 min.=10 max.=100663296 statistics");
	Table.deleteColumn("Surface (micron^2)");
	Table.deleteColumn("X");
	Table.deleteColumn("Y");
	Table.deleteColumn("Z");
	Table.deleteColumn("Mean dist. to surf. (micron)");
	Table.deleteColumn("SD dist. to surf. (micron)");
	Table.deleteColumn("Median dist. to surf. (micron)");
	Table.rename("Statistics for boutons redirect to "+redirectTitle, "New");

	// Adding in the parameters from the image file name to the Results Table
	thresholdArray = newArray(Table.size);
	Array.fill(thresholdArray, thresh);
	dataArray = newArray(Table.size);
	categoryArray = newArray(Table.size);
	mouseArray = newArray(Table.size);
	synapseArray = newArray(Table.size);
	vampArray = newArray(Table.size);
	for (i=0; i<Table.size; i++) {
		dataArray[i] = dataType;
		categoryArray[i] = headers[0];
		mouseArray[i] = headers[1];
		synapseArray[i] = boutonType[1];
		vampArray[i] = boutonType[2];
	}
	Table.setColumn("Threshold", thresholdArray);
	Table.setColumn("Data Type", dataArray);
	Table.setColumn("Category", categoryArray);
	Table.setColumn("Mouse", mouseArray);
	Table.setColumn("Synapse", synapseArray);
	Table.setColumn("Vamp", vampArray);
	
	// Merging all of the threshold and channel results tables into one table for saving
	if (initialize) {
		Table.rename("New", "Mouse-Results");
	} else {
		Table.create("Data");
		extras = newArray("Threshold", "Data Type", "Category", "Mouse", "Synapse", "Vamp");
		rowHeaders = Array.concat(colNames,extras);
		for (j = 0; j < rowHeaders.length; j++) {
			transfer = Table.getColumn(rowHeaders[j], "New");
			current = Table.getColumn(rowHeaders[j], "Mouse-Results");
			current = Array.concat(current,transfer);
			Table.setColumn(rowHeaders[j], current, "Data");
		}
		close("New");
		close("Mouse-Results");
		Table.rename("Data", "Mouse-Results");
	}
	

	if (imgBool) {
		selectWindow("Masked image for boutons redirect to "+redirectTitle);
		setThreshold(5, 65535);
		run("Convert to Mask", "method=Default background=Dark black");
		run("Analyze Particles...", "clear add stack");
		selectWindow(redirectTitle);
		run("Duplicate...", "title=pic duplicate");
		roiManager("Set Color", "blue");
		roiManager("Show All without labels");
		run("From ROI Manager");
		close("ROI Manager");
		run("Flatten", "stack");
		saveAs("TIFF", pathSave+"-"+dataType+"-"+thresh+".tif");
		close("Masked image for boutons redirect to "+redirectTitle);
		close("pic");
	}
}



// Main script

// Getting the folder where images are located
imageDir=getDirectory("Select Spinal Cord Folder for Image Processing");
roiDir=imageDir+"Analysis_Results\\";
imagelist = getFileList(imageDir);
imageArray = newArray;
File.makeDirectory(roiDir);
saveSettings();
setOption("BlackBackground", true);
run("Set Measurements...", "area mean standard modal min centroid perimeter shape integrated median skewness stack redirect=None decimal=3");
setBatchMode(true);

start = getTime();
// Grabbing each image individually for analysis
for (image=0;image<imagelist.length;image++) {

// Opening each image that is a 60x image and getting the image title
// This can be modified to whatever file ending is needed (could just be .oib)
	if (endsWith(imagelist[image], "60x2.oib")) {
		imgName=imagelist[image];
		baseNameEnd=indexOf(imgName, ".oib"); 
        baseName=substring(imgName, 0, baseNameEnd);
        maskName = baseName + "-mask.tif";
        savePath = roiDir+baseName;
		fpath = imageDir+imgName;
		mpath = imageDir+maskName;
		run("Bio-Formats Importer", "open=["+fpath+"] color_mode=Grayscale rois_import=[ROI manager] split_channels view=Hyperstack stack_order=XYCZT");
		imageGroups = split(baseName, "_");
		channelTypes = split(imageGroups[2], "-");
		neunResult = 0;
		nucleiResult = 0;

		programStart = getTime();
		// Cycling through the images to rename them
		for (channel = 0; channel < nImages; channel++) {	
			selectImage(channel+1);
			rename("ch"+toString(channel+1));
		}

		//print(mpath);
		open(mpath);
		rename("mmask");
		
	// Checking to see if there are any neurons at all
		selectWindow("mmask");
		run("Analyze Particles...", "size=50-5000 circularity=0.00-1.00 show=Masks clear stack");
		rename("neurons");
		neuronResult = nResults;
		print(neuronResult);

		if (neuronResult > 0) {
		
		// Creating the outline around the cell body for bouton analysis
			// Getting the inside of the doughnut (the core of the neurons)
			selectWindow("ch2");
			run("Duplicate...", "title=boutons duplicate");
			selectWindow("neurons");
			run("Duplicate...", "title=inner2 duplicate");
			run("Duplicate...", "title=inner duplicate");
			run("EDM Binary Operations", "iterations=45 operation=open stack");
			run("EDM Binary Operations", "iterations=8 operation=erode stack");
			run("Analyze Particles...", "size=0-infinity circularity=0.00-1.00 show=Nothing clear stack");
			innerResults = nResults;
			selectWindow("inner2");
			run("Analyze Particles...", "size=0-infinity circularity=0.00-1.00 show=Nothing clear stack");
			inner2Results = nResults;
			if (innerResults == inner2Results) {
				close("inner2");
			} else {
				close("inner");
				selectWindow("inner2");
				rename("inner");
				run("EDM Binary Operations", "iterations=30 operation=open stack");
				run("EDM Binary Operations", "iterations=8 operation=erode stack");
			}
			
			
			//setBatchMode("exit and display");
	
			// Back to finding the boutons
			// Getting the outer edge of where the boutons should be
			print("getting cell edges for boutons");
			selectWindow("neurons");
			run("Duplicate...", "title=outter duplicate");
			run("EDM Binary Operations", "iterations=30 operation=open stack");
			run("BinaryDilateNoMerge8 ", "iterations=25 white");
			run("Analyze Particles...", "size=0-infinity circularity=0.00-1.00 show=Nothing display clear add stack");
			cResult = nResults;
			nerveLocation = Table.getColumn("Slice");
			totalnerves = nerveLocation.length;
			close("Results");
			roiManager("Show None");
			for (k = 0; k < cResult; k++) {
				roiManager("Select", k);
				run("Convex Hull");
				run("Fit Spline");
				run("Fill", "slice");
				run("Select None");
			}
			close("ROI Manager");
			close("neurons");
			
			// For debugging purposes
			//run("Analyze Particles...", "size=0-infinity circularity=0.00-1.00 show=Nothing display clear stack");
			//cResult = nResults;
			//nerveLocation = Table.getColumn("Slice");
			//totalnerves = nerveLocation.length;
			//close("Results");
			
			lastslice = 0;
			intSlice = 0;
			for (intNerve = 0; intNerve < totalnerves; intNerve++) {
			    if (lastslice != nerveLocation[intNerve]) {
			    	selectWindow("outter");
			    	setSlice(nerveLocation[intNerve]);
					run("Create Selection");
					roiManager("Add");
					selectWindow("boutons");
					roiManager("Select", intSlice);
					run("Clear Outside", "slice");
					intSlice++;
					selectWindow("inner");
			    	setSlice(nerveLocation[intNerve]);
					run("Create Selection");
					roiManager("Add");
					selectWindow("boutons");
					roiManager("Select", intSlice);
					run("Clear", "slice");
					intSlice++;
					for (i = lastslice+1; i < nerveLocation[intNerve]; i++) {
						setSlice(i);
						run("Select All");
						run("Clear", "slice");
						run("Select None");
					}
					lastslice = nerveLocation[intNerve];
			    }
			}
			selectWindow("boutons");
			totalSlices = nSlices;
			for (i = lastslice+1; i <= totalSlices; i++) {
				setSlice(i);
				run("Select All");
				run("Clear", "slice");
				run("Select None");
			}
			run("Select None");
			roiManager("Deselect");
			roiManager("Delete");
			close("ROI Manager");
			close("outter");
			close("inner");
			savePics = true;
		} else {
			print("no cells");
			close("boutons");
			selectWindow("neurons");
			rename("boutons");
			savePics = false;
		}

		//setBatchMode("exit and display");

		print("running bouton check");
		// Checking to make sure that particles can be found
		boutonCheck("boutons", 1500, channelTypes, imageGroups, savePath, baseName, savePics);
		boutonCheck("boutons", 2000, channelTypes, imageGroups, savePath, baseName, savePics);
		boutonCheck("boutons", 2500, channelTypes, imageGroups, savePath, baseName, savePics);

		Table.rename("Mouse-Results", "Results-"+image);
		imageArray[image] = "Results-"+toString(image);
			
		run("Close All");
		print("Finished processing the image. Took " + ((getTime()-programStart)/1000)/60 + "minutes");
		//exit;
	} 
}

//imageArray = newArray(16);
//for (i=0;i<imageArray.length;i++) {
//	imageArray[i] = "Results-"+i;
//}

imageArray = Array.deleteValue(imageArray, NaN);
headings = split(Table.headings(imageArray[0]), "\t");
for (k=1; k<imageArray.length; k++) {
	if (imageArray[k] != "") {
		Table.create("Data");
		for (j = 0; j < headings.length; j++) {
			transfer = Table.getColumn(headings[j], imageArray[k]);
			current = Table.getColumn(headings[j], imageArray[0]);
			current = Array.concat(current,transfer);
			Table.setColumn(headings[j], current, "Data");
		}
		close(imageArray[0]);
		close(imageArray[k]);
		Table.rename("Data", imageArray[0]);
	}
}
Table.rename(imageArray[0], "Results");

//imageDir=getDirectory("Select Spinal Cord Folder for Image Processing");
//roiDir=imageDir+"Analysis_Results\\";
//imagelist = getFileList(imageDir);
//for (i=0;i<imagelist.length;i++) {
//	print(imagelist[i]);
//}

saveAs("Results", roiDir+"All-stats.csv");
close("Results");

print("Finished processing all images in "+imageDir+".");
print("Took " + ((getTime()-programStart)/1000)/60 + "minutes");
setBatchMode(false);