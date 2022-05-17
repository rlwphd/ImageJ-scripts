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
//function endplate(img){
	img="this";
	selectWindow(img);
	run("Z Project...", "projection=[Max Intensity]");
	rename("Z");
	close(img);
	
	selectWindow("Z");
	run("Split Channels");
	
	// Cycling through the images to rename them
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
	run("Invert LUT");
	run("Invert");
	
	// Finding the areas of interest
	run("Analyze Particles...", "size=25-Infinity show=Nothing clear");
	area = Table.getColumn("Area");
	bx = Table.getColumn("BX");
	by = Table.getColumn("BY");
	bw = Table.getColumn("Width");
	bh = Table.getColumn("Height");
		
	// Adding the width/height of the box to the x/y coordinate to get the coordinates of the box
	for (i=0; i<bx.length; i++) {
		bw[i] = bx[i]+bw[i];
		bh[i] = by[i]+bh[i];
	}
	
	// Determining if the any boxes overlap and if they do, merging them
	d = newArray(bx.length);
	Array.fill(d, -1);
	for (i=0; i<bx.length; i++) {
		for (j=0; j<bx.length; j++) {
			if (bw[i]>bx[j] && bx[j]>bx[i] && bh[i]>by[j] && by[j]>by[i]) {
				if (area[i] < 330 || area[j] < 330) {
					bx = mergebox(bx,i,j,"min");
					bw = mergebox(bw,i,j,"max");
					by = mergebox(by,i,j,"min");
					bh = mergebox(bh,i,j,"max");
					if (i>j) {
						d[i] = j;
					} else {
						d[j] = i;
					}
				}
			} else if (bw[i]>bx[j] && bx[j]>bx[i] && bh[i]>bh[j] && bh[j]>by[i]) {
				if (area[i] < 330 || area[j] < 330) {
					bx = mergebox(bx,i,j,"min");
					bw = mergebox(bw,i,j,"max");
					by = mergebox(by,i,j,"min");
					bh = mergebox(bh,i,j,"max");
					if (i>j) {
						d[i] = j;
					} else {
						d[j] = i;
					}
				}
			} else if (bw[i]>bw[j] && bw[j]>bx[i] && bh[i]>by[j] && by[j]>by[i]) {
				if (area[i] < 330 || area[j] < 330) {
					bx = mergebox(bx,i,j,"min");
					bw = mergebox(bw,i,j,"max");
					by = mergebox(by,i,j,"min");
					bh = mergebox(bh,i,j,"max");
					if (i>j) {
						d[i] = j;
					} else {
						d[j] = i;
					}
				}
			} else if (bw[i]>bw[j] && bw[j]>bx[i] && bh[i]>bh[j] && bh[j]>by[i]) {
				if (area[i] < 330 || area[j] < 330) {
					bx = mergebox(bx,i,j,"min");
					bw = mergebox(bw,i,j,"max");
					by = mergebox(by,i,j,"min");
					bh = mergebox(bh,i,j,"max");
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
	
	//Returning the box dimensions to x/y and width/height
	for (i=0; i<bx.length; i++) {
		bw[i] = bw[i]-bx[i];
		bh[i] = bh[i]-by[i];
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
	
	// Cutting out the Endplates + 2x the box for background analysis
	//for (i=0, i<bx.length, i++) {
		//makeRectangle(bx[i], by[i], bw[i], bh[i]);
		//run("Duplicate...", "title=[Endplate-"+i+"]");
	//}
	
	// Getting the background's mean value and standard deviation
	selectWindow("VAChT background");
	List.setMeasurements;
	//print(List.getList);
	back_mean = List.get("Mean");
	back_std = List.get("StdDev");
	Vacht_cutoff = 2*back_std + back_mean;

	// Grabbing only the VAChT labels that are 2x std above the mean
	selectWindow("VAChT");
	setThreshold(Vacht_cutoff, 255);
	run("Convert to Mask");
	//run("Invert LUT");
	
	// Checking the 2x std on the background
	selectWindow("VAChT background");
	setThreshold(Vacht_cutoff, 255);
	run("Convert to Mask");
	//run("Invert LUT");	

	// Getting the bungarotoxin's mean value and standard deviation
	selectWindow("Bung background");
	List.setMeasurements;
	bung_mean = List.get("Mean");
	bung_std = List.get("StdDev");
	bung_cutoff = 2*bung_std + bung_mean;

	// Grabbing only the VAChT labels that are 2x std above the mean
	selectWindow("Bung");
	setThreshold(bung_cutoff, 255);
	run("Convert to Mask");
	//run("Invert LUT");
	
	selectWindow("VAChT");
	selectWindow("Bung");
	run("Colocalization Threshold", "channel_1=VAChT channel_2=Bung use=None channel=[Red : Green] show include");
	
}


// ***** MAIN *****

// Allowing user to choose a file or a folder of files to analyze
ff = getBoolean("Would you like to process a single file or a whole folder of files?", "File", "Folder");
// General setup
saveSettings();
setOption("BlackBackground", true);
run("Set Scale...", "distance=0 known=0 unit=pixel");
run("Set Measurements...", "area mean standard perimeter bounding median area_fraction redirect=None decimal=3");
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

roiDir=imageDir+"Analysis_Results\\";
File.makeDirectory(roiDir);

if (ff == 1) {
	run("Bio-Formats Importer", "open=["+fpath+"] color_mode=Grayscale rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	endplate(name);
} else {
	// Grabbing each image individually for analysis
	for (image=0;image<imagelist.length;image++) {

	// Opening each image that is a oib image and getting the image title
	// This can be modified to whatever file ending is needed (could just be .oib)
		if (endsWith(imagelist[image], ".oib")) {
			imgName=imagelist[image];
			fpath = imageDir+imgName;
			run("Bio-Formats Importer", "open=["+fpath+"] color_mode=Grayscale rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
			endplate(name);
		}
	}
}
endTime = (getTime()-start)/1000/60;
print("Analysis is finished! Took "+endTime+" minutes to run.");

		
