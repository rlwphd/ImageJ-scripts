// Defining the merge box outline function
function mergebox(b,i,j,mm){
	val = i;
	if (i>j) {
		val = j;
	}
	if (mm == "min") {
		b[val] = minOf(b[i], b[j]);
	}
	if (mm == "max") {
		b[val] = maxOf(b[i], b[j]);
	}
	
	return b;
}

// Defining the Endplate Analysis function
function endplate(img, saveDir){
	selectWindow(img);
	run("Z Project...", "projection=[Max Intensity]");
	rename("Z");
	close(img);
	
	// Establishing the correct criteria
	selectWindow("Z");
	//run("Set Scale...", "distance=0 known=0 unit=pixel"); // Option to force image to pixels
	part_size = 20;
	area_limit = 330;
	minor_limit = 15;
	getPixelSize(unit, pw, ph); // Grabbing the units and pixel size of the image
	if (unit != "pixel") {
		part_size = 20*pw*ph;
		area_limit = 330*pw*ph;
		minor_limit = 15*pw;
	}
	
	selectWindow("Z");
	run("Split Channels");
	
	// Cycling through the images to rename them
	// This assumes that the VAChT channel is the 3rd channel 
	// and the Bungarotoxin is the 4th channel
	for (channel = 0; channel < nImages; channel++) {	
		selectImage(channel+1);
		rename("ch"+toString(channel+1));
	}
	
	// Grabbing where the Bungarotoxin staining is plus a bit 
	selectWindow("ch4");
	run("Duplicate...", "title=background");
	run("Enhance Contrast...", "saturated=1.8");
	setAutoThreshold("Triangle");
	run("Convert to Mask");
	run("Despeckle");
	//run("Invert LUT");
	run("Invert");
	

	// Finding the areas of interest
	run("Analyze Particles...", "size="+part_size+"-Infinity show=Nothing clear");
	area = Table.getColumn("Area");
	bx = Table.getColumn("BX");
	by = Table.getColumn("BY");
	bw = Table.getColumn("Width");
	bh = Table.getColumn("Height");
		
	// Adding the width/height of the box to the x/y coordinate to get the coordinates of the box
	for (i=0; i<bx.length; i++) {
		bx[i] = bx[i]/pw; //converting to pixels
		by[i] = by[i]/ph; //converting to pixels
		bw[i] = bw[i]/pw + bx[i]; //converting to pixels
		bh[i] = bh[i]/ph + by[i]; //converting to pixels
	}
	
	// Determining if the any boxes overlap and if they do, merging them
	d = newArray(bx.length);
	Array.fill(d, -1);
	for (i=0; i<bx.length; i++) {
		for (j=0; j<bx.length; j++) {
			if (bw[i]>bx[j] && bx[j]>bx[i] && bh[i]>by[j] && by[j]>by[i]) {
				if (area[i] < area_limit || area[j] < area_limit) {
					bx = mergebox(bx,i,j,"min");
					bw = mergebox(bw,i,j,"max");
					by = mergebox(by,i,j,"min");
					bh = mergebox(bh,i,j,"max");
					area[i] = area[i] + area[j];
					area[j] = area[i];
					if (i>j) {
						d[i] = j;
					} else {
						d[j] = i;
					}
				}
			} else if (bw[i]>bx[j] && bx[j]>bx[i] && bh[i]>bh[j] && bh[j]>by[i]) {
				if (area[i] < area_limit || area[j] < area_limit) {
					bx = mergebox(bx,i,j,"min");
					bw = mergebox(bw,i,j,"max");
					by = mergebox(by,i,j,"min");
					bh = mergebox(bh,i,j,"max");
					area[i] = area[i] + area[j];
					area[j] = area[i];
					if (i>j) {
						d[i] = j;
					} else {
						d[j] = i;
					}
				}
			} else if (bw[i]>bw[j] && bw[j]>bx[i] && bh[i]>by[j] && by[j]>by[i]) {
				if (area[i] < area_limit || area[j] < area_limit) {
					bx = mergebox(bx,i,j,"min");
					bw = mergebox(bw,i,j,"max");
					by = mergebox(by,i,j,"min");
					bh = mergebox(bh,i,j,"max");
					area[i] = area[i] + area[j];
					area[j] = area[i];
					if (i>j) {
						d[i] = j;
					} else {
						d[j] = i;
					}
				}
			} else if (bw[i]>bw[j] && bw[j]>bx[i] && bh[i]>bh[j] && bh[j]>by[i]) {
				if (area[i] < area_limit || area[j] < area_limit) {
					bx = mergebox(bx,i,j,"min");
					bw = mergebox(bw,i,j,"max");
					by = mergebox(by,i,j,"min");
					bh = mergebox(bh,i,j,"max");
					area[i] = area[i] + area[j];
					area[j] = area[i];
					if (i>j) {
						d[i] = j;
					} else {
						d[j] = i;
					}
				}
			}
			for (k=0; k<d.length; k++) {
				if (d[k]>0) {
					bx = Array.deleteIndex(bx, k);
					bw = Array.deleteIndex(bw, k);
					by = Array.deleteIndex(by, k);
					bh = Array.deleteIndex(bh, k);
					area = Array.deleteIndex(area, k);
					
					i = d[k]-1;
					d = Array.deleteIndex(d,k);
				}
			}
		}
	}
	//Returning the box dimensions to x/y and width/height
	for (i=0; i<bx.length; i++) {
		bw[i] = bw[i]-bx[i]; //value is in pixels
		bh[i] = bh[i]-by[i]; //value is in pixels
	}
	
	// Cutting out where the Bungarotoxin staining is in the VAChT & Btx channel
	imageCalculator("AND create", "background","ch3");
	selectWindow("Result of background");
	rename("VAChT");
	imageCalculator("AND create", "background","ch4");
	selectWindow("Result of background");
	rename("Bung");
	selectWindow("background");
	run("Invert");
	imageCalculator("AND create", "background","ch3");
	selectWindow("Result of background");
	rename("VAChT background");
	imageCalculator("AND", "background","ch4");
	rename("Bung background");
	
	// Using the Bungarotoxin image to grab average values
	// for the holes created in the background images
	selectWindow("Bung");
	run("Duplicate...", "title=Outline");
	setThreshold(1, 255);
	run("Convert to Mask");
	run("Create Selection");
	close("Outline");
	
	selectWindow("Bung");
	run("Restore Selection");
	add_val = getValue("Median");
	run("Select None");
	
	selectWindow("Bung background");
	run("Restore Selection");
	run("Add...", "value="+add_val);
	run("Select None");
	
	selectWindow("VAChT");
	run("Restore Selection");
	add_val = 0.5*getValue("Median");
	run("Select None");
	
	selectWindow("VAChT background");
	run("Restore Selection");
	run("Add...", "value="+add_val);
	run("Select None");
	
	// Cutting out the Endplates + 3x the box for background analysis
	results_array = Array.copy(bx);
	Array.fill(results_array, 0);
	for (i=0; i<bx.length; i++) {
		topx = bx[i]-1.5*bw[i];
		topy = by[i]-1.5*bh[i];
		width = 3*bw[i];
		height = 3*bh[i];
		
		// Grabbing the desired Endplate from the Bungarotoxin windows
		selectWindow("Bung");
		makeRectangle(bx[i], by[i], bw[i], bh[i]); //Need everything to be in pixels
		run("Duplicate...", "title=[Bung Endplate-"+i+"]");
		selectWindow("Bung");
		run("Select None");
		
		selectWindow("Bung background");
		makeRectangle(topx, topy, width, height);
		run("Duplicate...", "title=[Bung background-"+i+"]");
		selectWindow("Bung background");
		run("Select None");
		
		// Getting the bungarotoxin's mean value and standard deviation
		selectWindow("Bung background-"+i);
		bung_mean = getValue("Mean");
		bung_std = getValue("StdDev");
		bung_cutoff = 2*bung_std + bung_mean;
		close("Bung background-"+i);
	
		// Grabbing only the Bung labels that are 2x std above the mean
		selectWindow("Bung Endplate-"+i);
		run("Duplicate...", "title=[Bung Mask-"+i+"]");
		setThreshold(bung_cutoff, 255);
		run("Convert to Mask");
		//run("Invert LUT");
		run("Create Selection");
		minor = getValue("Minor");
		run("Select None");
		

		if (minor>minor_limit) {
			selectWindow("VAChT");
			makeRectangle(bx[i], by[i], bw[i], bh[i]);
			run("Duplicate...", "title=[VAChT Endplate-"+i+"]");
			selectWindow("VAChT");
			run("Select None");
			
			selectWindow("VAChT background");
			makeRectangle(topx, topy, width, height);
			run("Duplicate...", "title=[VAChT background-"+i+"]");
			selectWindow("VAChT background");
			run("Select None");

			// Getting the background's mean value and standard deviation
			selectWindow("VAChT background-"+i);
			back_mean = getValue("Mean");
			back_std = getValue("StdDev");
			vacht_cutoff = 2*back_std + back_mean;
			close("VAChT background-"+i);
		
			// Grabbing only the VAChT labels that are 2x std above the mean
			selectWindow("VAChT Endplate-"+i);
			run("Duplicate...", "title=[VAChT Mask-"+i+"]");
			setThreshold(vacht_cutoff, 255);
			run("Convert to Mask");
			//run("Invert LUT");
			
			// Determining Endplate overlap
			selectWindow("Bung Mask-"+i);
			run("Create Selection");
			run("Select None");
			selectWindow("VAChT Mask-"+i);
			run("Restore Selection");
			overlap = getValue("%Area");
			run("Select None");
		
			// Save everything and close windows
			results_array[i] = overlap;
			
			selectWindow("Bung Mask-"+i);
			save_img = getTitle();
			org_img = img.substring(0,img.length-4);
			spath = saveDir+org_img+"-"+save_img+".tif";
			saveAs("Tiff", spath);
			close(org_img+"-"+save_img+".tif");
			selectWindow("VAChT Mask-"+i);
			save_img = getTitle();
			spath = saveDir+org_img+"-"+save_img+".tif";
			saveAs("Tiff", spath);
			close(org_img+"-"+save_img+".tif");
			selectWindow("Bung Endplate-"+i);
			save_img = getTitle();
			spath = saveDir+org_img+"-"+save_img+".tif";
			saveAs("Tiff", spath);
			close(org_img+"-"+save_img+".tif");
			selectWindow("VAChT Endplate-"+i);
			save_img = getTitle();
			spath = saveDir+org_img+"-"+save_img+".tif";
			saveAs("Tiff", spath);
			close(org_img+"-"+save_img+".tif");

		} else {
			close("Bung Endplate-"+i);
			close("Bung Mask-"+i);
		}
	}
	
	// Saving the results array
	Array.show("Results",results_array,bx,by,bw,bh);
	spath = saveDir+org_img+"-Results.csv";
	saveAs("Results", spath);
	run("Close");
	close("*");
}


// ***** MAIN *****

// Allowing user to choose a file or a folder of files to analyze
ff = getBoolean("Would you like to process a single file or a whole folder of files?", "File", "Folder");
// General setup
saveSettings();
setOption("BlackBackground", true);
run("Set Measurements...", "area mean standard perimeter bounding fit median area_fraction redirect=None decimal=3");
//setBatchMode(true);

start = getTime();
if (ff == 1) {
	fpath = File.openDialog("Select a File");
	imageDir = File.getParent(fpath);
	name = File.getName(fpath);
} else {
	// Getting the folder where images are located
	imageDir=getDirectory("Select Endplate Counting Folder for Image Processing");
	imagelist = getFileList(imageDir);
}

roiDir=imageDir+"\\Analysis_Results\\";
File.makeDirectory(roiDir);

if (ff == 1) {
	run("Bio-Formats Importer", "open=["+fpath+"] color_mode=Grayscale rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	endplate(name, roiDir);
} else {
	// Grabbing each image individually for analysis
	for (image=0;image<imagelist.length;image++) {

	// Opening each image that is a oib image and getting the image title
	// This can be modified to whatever file ending is needed (could just be .oib)
		if (endsWith(imagelist[image], ".oib")) {
			imgName=imagelist[image];
			fpath = imageDir+imgName;
			run("Bio-Formats Importer", "open=["+fpath+"] color_mode=Grayscale rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
			endplate(name, roiDir);
		}
	}
}
endTime = (getTime()-start)/1000/60;
print("Analysis is finished! Took "+endTime+" minutes to run.");

