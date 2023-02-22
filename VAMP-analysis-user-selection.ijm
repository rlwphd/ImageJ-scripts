// Getting the file or folder where images are located
items = newArray("file", "folder");
Dialog.create("File or Folder Analysis");
Dialog.addChoice("Would you like to analyze a whole folder or an individual file?", items, "file");
Dialog.addMessage("When processing a whole folder, images need to have the same markers in the channels.\nOtherwise the analysis may not work!!");
Dialog.addMessage("This program assumes that you have already ran the VAMP-analysis program.\nIf you have not, please cancel this one and run the other program first.");
Dialog.show();
process = Dialog.getChoice();
if (process == "file") {
	path = File.openDialog("Select File for Image Processing");
	imageDir = File.getParent(path);
	roiDir=imageDir+"\\Analysis_Results\\";
	imagelist = newArray(File.getName(path));
	File.makeDirectory(roiDir);
} else {
	imageDir=getDirectory("Select Spinal Cord Folder for Image Processing");
	roiDir=imageDir+"Analysis_Results\\";
	imagelist = getFileList(imageDir);
	imageArray = newArray;
	File.makeDirectory(roiDir);
	saveSettings();
}

channel=3;
Dialog.create("Fluorescence Analysis");
Dialog.addNumber("Number of Channels in Image:", channel);
Dialog.show();
channel=Dialog.getNumber();
types = newArray("Fast Blue", "VGLUT1", "VGLUT2", "VGAT/VIAAT", "VAMP1", "VAMP2", "NeuN", "ChAT", "Other");
Dialog.create("Channel Fluorescence");
for (i=1; i<=channel; i++) {
	Dialog.addChoice("Channel "+i+" Fluorescent Marker:", types);
}
Dialog.show();
chanArray = newArray(channel);
for (i=0; i<channel; i++) {
	chanArray[i] = Dialog.getChoice();
}
setOption("BlackBackground", true);
setForegroundColor(255, 255, 255);
setBackgroundColor(0, 0, 0);

// Grabbing each image individually for analysis
for (image=0; image<imagelist.length; image++) {

// Opening the desired images (that should already have been processed)
	imgName=imagelist[image];
	baseNameEnd=indexOf(imgName, ".oib"); 
    baseName=substring(imgName, 0, baseNameEnd);
	fpath = imageDir+imgName;

	// Opening all of the necessary images for analysis
	flist = getFileList(roiDir);
	for (f=0; f<flist.length; f++) {
		if (flist[f].contains(baseName) && endsWith(flist[f], "-ROIs.tif")) {
			fileEndIndex = indexOf(flist[f], "-ROIs.tif");
			secondfile = substring(flist[f], 0, fileEndIndex);
			secondfile = secondfile + ".tif";
			open(roiDir+flist[f]);
			open(roiDir+secondfile);
			selectWindow(flist[f]);
			run("Duplicate...", "title="+secondfile+"-UserROIs duplicate");
			setThreshold(1,255);
			run("Convert to Mask", "method=Default background=Dark black");
		}
	}
	
	// Would you like to open the original image as well?
	Dialog.create("Open Original Image?");
	Dialog.addMessage("Would you like to open the original image?\nSelect yes or no and click okay.\nIf you click cancel that will stop this macro.");
	Dialog.addChoice("Open Original Image: ", newArray("yes", "no"), "no");
	Dialog.show();
	original = Dialog.getChoice();
	if (original == "yes") {
		run("Bio-Formats Importer", "open=["+fpath+"] color_mode=Grayscale rois_import=[ROI manager] split_channels view=Hyperstack stack_order=XYCZT");
		open(imageDir+baseName+"-cells.tif");
		open(imageDir+baseName+"-celledges.tif");
	}

	msg = "All of the necessary files should be open.\nIt is now up to you, to use the \"Wand\" tool to select the objects that need to be DELETED\nfrom the images that are labelled \"UserROIs\".\nOnce you are done deleting all of the objects that you wish to remove, click \"OK\".";
	waitForUser("User Input", msg);
	
	numCells=0;
	imgArray = getList("image.titles");
	for (i=0; i<imgArray.length; i++) {
		if (imgArray[i].contains("UserROIs")) {
			lindex = imgArray[i].lastIndexOf("-Cell-");
			numCells = maxOf(numCells, parseInt(imgArray[i].charAt(lindex+1)));
		}
	}
	
	// Create ROI analysis overlap between the Vesicular and Vamp channels
	// Analyze the overlap and store the tables
	run("3D OC Options", "volume surface centroid mean_distance_to_surface std_dev_distance_to_surface bounding_box dots_size=20 font_size=50 store_results_within_a_table_named_after_the_image_(macro_friendly) redirect_to=none");
	run("Set Measurements...", "area mean standard modal min area_fraction stack redirect=None decimal=3");
	Table.create("Results");
	for (j=1; j<=numCells; j++) {
		for (i=0; i<imgArray.length; i++) {
			if (imgArray[i].contains("UserROI") && imgArray[i].contains("-Cell-"+j+"-VG")) {
				vesTitle = imgArray[i];
			} else if (imgArray[i].contains("UserROI") && imgArray[i].contains("-Cell-"+j+"-VA")) {
				vampTitle = imgArray[i];
			}
		}
		imageCalculator("AND create stack", vesTitle, vampTitle);
		rename("cell"+j+"-overlap");

		selectWindow(vesTitle);
		vesTablename = "Statistics for "+vesTitle;
		vesRoiImage = "Objects map of "+vesTitle;
		run("3D Objects Counter", "threshold=128 slice=1 min.=10 max.=4129270 objects statistics");
		roiCount = Table.size(vesTablename);
		for (i=1; i<=roiCount; i++) {
			selectWindow(vesRoiImage);
			setThreshold(i,i);
			run("Analyze Particles...", "size=0-infinity add stack");
			roiManager("show none");
			Table.reset("Results");
			selectWindow("cell"+j+"-overlap");
			roiManager("Select", Array.getSequence(roiManager("count")));
			roiManager("Measure");
			s = Table.getColumn("Slice", "Results");
			a = Table.getColumn("%Area", "Results");
			Table.reset("Results");
			for (k=1; k<s[0]; k++) {
				Table.set("Slice "+k, i, 0, vesTablename);
			}
			for (k=s[0]; k<=s[s.length-1]; k++) {
				val = 0;
				denom = 0;
				for (m=0; m<s.length; m++) {
					if (s[m] == k) {
						val = val + a[m];
						denom++;
					}
				}
				if (denom > 0) {
					Table.set("Slice "+k, i, val/denom, vesTablename);
				} else {
					Table.set("Slice "+k, i, 0, vesTablename);
				}
			}
			for (k=s[s.length-1]+1; k<=slices; k++) {
				Table.set("Slice "+k, i, 0, vesTablename);
			}
		}
		selectWindow(vesRoiImage);
		resetThreshold;
		close("cell"+j+"-overlap");
		close(vesTitle);
		// Saving the analysis table for the vesicular channel
		selectWindow(vesTablename);
		saveName = vesTitle.substring(0,vesTitle.indexOf("ROI"));
		saveAs("results", roiDir+saveName+".csv");
		close(Table.title);
		// Saving the objects image for the vesicular channel
		selectWindow(vesRoiImage);
		saveAs("tiff", roiDir+vesTitle+".tif");
		rename(vesRoiImage);

		// Figuring out which vesicular ROIs the vamp ROIs line up with
		selectWindow(vampTitle);
		vampTablename = "Statistics for "+vampTitle;
		vampRoiImage = "Objects map of "+vampTitle;
		run("3D Objects Counter", "threshold=128 slice=1 min.=10 max.=4129270 objects statistics");
		roiCount = Table.size(vampTablename);
		for (i=1; i<=roiCount; i++) {
			selectWindow(vampRoiImage);
			setThreshold(i,i);
			run("Analyze Particles...", "size=0-infinity add stack");
			roiManager("show none");
			Table.reset("Results");
			selectWindow(vesRoiImage);
			roiManager("Select", Array.getSequence(roiManager("count")));
			roiManager("Measure");
			s = Table.getColumn("Slice", "Results");
			a = Table.getColumn("Mode", "Results");
			b = Table.getColumn("Max", "Results");
			Table.reset("Results");
			for (k=1; k<s[0]; k++) {
				Table.set("Slice "+k, i, 0, vampTablename);
			}
			for (k=s[0]; k<=s[s.length-1]; k++) {
				val = 0;
				denom = 0;
				for (m=0; m<s.length; m++) {
					if (s[m] == k) {
						if (a[m] > 0) {
							val = val + a[m];
							denom++;
						} else {
							val = val + b[m];
							denom++;
						}
					}
				}
				if (denom > 0) {
					Table.set("Slice "+k, i, val/denom, vampTablename);
				} else {
					Table.set("Slice "+k, i, 0, vampTablename);
				}
			}
			for (k=s[s.length-1]+1; k<=slices; k++) {
				Table.set("Slice "+k, i, 0, vampTablename);
			}
			roiManager("reset");
		}
		close("ROI Manager");
		selectWindow(vampRoiImage);
		resetThreshold;
		close(vampTitle);
		// Saving the analysis table for the vamp channel
		selectWindow(vampTablename);
		saveName = vampTitle.substring(0,vampTitle.indexOf("ROI"));
		saveAs("results", roiDir+saveName+".csv");
		close(Table.title);
		// Saving the objects image for the vamp channel
		selectWindow(vampRoiImage);
		saveAs("tiff", roiDir+vampTitle+".tif"););
		rename(vampRoiImage);
		close(vampRoiImage);
		close(vesRoiImage);
	}
	close("Results");
	close("*");
	close
}


