// Defining the merge indicies function
function mergeindex(b,i,j,mm){
	if (mm == "min") {
		b[i] = minOf(b[i], b[j]);
		b[j] = minOf(b[i], b[j]);		
	}
	if (mm == "max") {
		b[i] = maxOf(b[i], b[j]);
		b[j] = maxOf(b[i], b[j]);
	}
	return b;
}


// Defining the Endplate Analysis function
function endplate(image_num, img, roiDir, saveDir) {
	// Creating the name of the image for saving
	name = split(img, "_-");
	name_split = Array.slice(name,1,3);
	org_img = String.join(name_split, "_");
	org_img += "_"+image_num;
	
	selectWindow(img);
	run("Z Project...", "projection=[Max Intensity]");
	rename("Z");
	close(img);
	
	// Establishing the correct criteria
	selectWindow("Z");
	//run("Set Scale...", "distance=0 known=0 unit=pixel"); // Option to force image to pixels
	part_size = 20;
	area_limit = 350;
	minor_limit = 15;
	getPixelSize(unit, pw, ph); // Grabbing the units and pixel size of the image
	if (unit != "pixel") {
		part_size = 20*pw*ph;
		area_limit = 350*pw*ph;
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
	getRawStatistics(nPixels, mean, min, max);
	getHistogram(values, counts, 255);
	Array.reverse(values);
	Array.reverse(counts);
	sum = Array.copy(counts);
	Array.fill(sum, 0);
	for (i=1; i<counts.length; i++) {
		counts[i] += counts[i-1];
		sum[i] = counts[i]*100/nPixels;
		if (sum[i] > 2) {
			stop = i-1;
			i = counts.length;
		}
	}
	setThreshold(values[stop], max);
	run("Convert to Mask");
	run("Analyze Particles...", "size=20-5000 pixel show=Masks clear");
	run("Invert LUT");
	close("background");
	selectWindow("Mask of background");
	rename("160");
	run("Smooth");
	run("Duplicate...", "title=128");
	selectWindow("160");
	setThreshold(160, 255);
	run("Convert to Mask");
	selectWindow("128");
	setThreshold(128, 255);
	run("Convert to Mask");
	
	// Use 128 to get bounding box but use 160 for measurements
	run("Analyze Particles...", "size=0-10000 pixel show=Nothing clear display");
	bx = Table.getColumn("BX");
	by = Table.getColumn("BY");
	bw = Table.getColumn("Width");
	bh = Table.getColumn("Height");
	x = Table.getColumn("X");
	y = Table.getColumn("Y");
	run("Clear Results");
	
	// Adding the width/height of the box to the x/y coordinate to get the coordinates of the box
	getPixelSize(unit, pw, ph); // Grabbing the units and pixel size of the image
	if (unit != "pixel") {
		for (i=0; i<bx.length; i++) {
			bx[i] = bx[i]/pw; //converting to pixels
			by[i] = by[i]/ph; //converting to pixels
			bw[i] = bw[i]/pw + bx[i]; //converting to pixels
			bh[i] = bh[i]/ph + by[i]; //converting to pixels
			x[i] = x[i]/pw; //converting to pixels
			y[i] = y[i]/ph; //converting to pixels
			
		}
	} else {
		for (i=0; i<bx.length; i++) {
			bw[i] = bw[i] + bx[i];
			bh[i] = bh[i] + by[i];
		}	
	}
	
	// Determining if the any boxes overlap and if they do, merging them
	// Defining the merge indicies function
	function mergeindex(b,i,j,mm){
		if (mm == "min") {
			b[i] = minOf(b[i], b[j]);
			b[j] = minOf(b[i], b[j]);		
		}
		if (mm == "max") {
			b[i] = maxOf(b[i], b[j]);
			b[j] = maxOf(b[i], b[j]);
		}
		return b;
	}
	// Merging the boxes if there is overlap between the center and the box with some wiggle room
	for (i=0; i<bx.length-1; i++) {
		for (j=i+1; j<bx.length; j++) {
			if (x[j]>bx[i]-5 && bw[i]+5>x[j] && y[j]>by[i]-5 && bh[i]+5>y[j]) {
				bx = mergeindex(bx,i,j,"min");
				bw = mergeindex(bw,i,j,"max");
				by = mergeindex(by,i,j,"min");
				bh = mergeindex(bh,i,j,"max");
			}
		}
	}
	// Removing the duplicates from the arrays
	end = bx.length;
	current_end = bx.length-1;
	for (k=0; k<end; k++) {
		for (j=k+1; j<current_end+1; j++) {
			if (bx[k]==bx[j] && bw[k]==bw[j]) {
				bx = Array.deleteIndex(bx, j);
				bw = Array.deleteIndex(bw, j);
				by = Array.deleteIndex(by, j);
				bh = Array.deleteIndex(bh, j);
				current_end = current_end - 1;
				j = j - 1;
			}
		}
	}
	
	for (i=0; i<bx.length; i++) {
		bw[i] = bw[i] - bx[i];
		bh[i] = bh[i] - by[i];
	}	
	
	selectWindow("128");
	for (i=0; i<bx.length; i++) {
		makeRectangle(bx[i], by[i], bw[i], bh[i]);
		roiManager("Add");	
	}
	roiManager("deselect");
	roiManager("Measure");
	roiManager("deselect");
	roiManager("delete");
	
	
	// Running through a second time just to make sure
	bx = Table.getColumn("BX");
	by = Table.getColumn("BY");
	bw = Table.getColumn("Width");
	bh = Table.getColumn("Height");
	x = Table.getColumn("X");
	y = Table.getColumn("Y");
	run("Clear Results");
	
	// Adding the width/height of the box to the x/y coordinate to get the coordinates of the box
	getPixelSize(unit, pw, ph); // Grabbing the units and pixel size of the image
	if (unit != "pixel") {
		for (i=0; i<bx.length; i++) {
			bx[i] = bx[i]/pw; //converting to pixels
			by[i] = by[i]/ph; //converting to pixels
			bw[i] = bw[i]/pw + bx[i]; //converting to pixels
			bh[i] = bh[i]/ph + by[i]; //converting to pixels
			x[i] = x[i]/pw; //converting to pixels
			y[i] = y[i]/ph; //converting to pixels
			
		}
	} else {
		for (i=0; i<bx.length; i++) {
			bw[i] = bw[i] + bx[i];
			bh[i] = bh[i] + by[i];
		}	
	}
	
	// Merging the boxes if there is overlap between the center and the box with some wiggle room
	for (i=0; i<bx.length-1; i++) {
		for (j=i+1; j<bx.length; j++) {
			if (x[i]>bx[j]-5 && bw[j]+5>x[i] && y[i]>by[j]-5 && bh[j]+5>y[i]) {
				bx = mergeindex(bx,i,j,"min");
				bw = mergeindex(bw,i,j,"max");
				by = mergeindex(by,i,j,"min");
				bh = mergeindex(bh,i,j,"max");
			}
		}
	}
	// Removing the duplicates from the arrays
	end = bx.length;
	current_end = bx.length-1;
	for (k=0; k<end; k++) {
		for (j=k+1; j<current_end+1; j++) {
			if (bx[k]==bx[j] && bw[k]==bw[j]) {
				bx = Array.deleteIndex(bx, j);
				bw = Array.deleteIndex(bw, j);
				by = Array.deleteIndex(by, j);
				bh = Array.deleteIndex(bh, j);
				current_end -= 1;
				j -= 1;
			}
		}
	}
	
	
	for (i=0; i<bx.length; i++) {
		bw[i] = bw[i] - bx[i];
		bh[i] = bh[i] - by[i];
	}	
	
	selectWindow("128");
	for (i=0; i<bx.length; i++) {
		makeRectangle(bx[i], by[i], bw[i], bh[i]);
		roiManager("Add");	
	}
	roiManager("deselect");
	roiManager("Measure");
	
	// Dropping the selections that don't fit the area requirement for the Endplate
	eparea = Table.getColumn("Area");
	end = eparea.length;
	for (i=0; i<end; i++) {
		if (eparea[i] < 350) {
			roiManager("select", i);
			roiManager("delete");
			eparea = Array.deleteIndex(eparea, i);
			i -= 1;
			end -= 1;
		}
	}
	run("Select None");
	run("Clear Results");
	
	// Removing the extraneous spots
	selectWindow("128");
	roiManager("Combine");
	run("Clear Outside");
	run("Select None");
	roiManager("deselect");
	selectWindow("160");
	roiManager("Combine");
	run("Clear Outside");
	run("Select None");
	roiManager("deselect");
	
	// Setting up the bounding boxes for grabbing the background and storing the original bouding box for later
	selectWindow("ch4");
	run("Select None");
	roiManager("List");
	bx = Table.getColumn("X");
	by = Table.getColumn("Y");
	bw = Table.getColumn("Width");
	bh = Table.getColumn("Height");
	bbx = Table.getColumn("X");
	bby = Table.getColumn("Y");
	bbw = Table.getColumn("Width");
	bbh = Table.getColumn("Height");
	run("Clear Results");
	close("Overlay Elements of ch4");
	
	for (i=0; i<bx.length; i++) {
		bw[i] = bw[i]*3; //converting to pixels and expanding box to be 3 times as big
		bh[i] = bh[i]*3; //converting to pixels	
		bx[i] = bx[i] - bw[i]/2; //converting to pixels
		by[i] = by[i] - bh[i]/2; //converting to pixels
		if (bx[i]<0) {
			bx[i] = 0;
		}
		if(by[i]<0) {
			by[i] = 0;
		}
	}
	
	// Grabbing the background of the VAChT channel
	vmean = Array.copy(bx);
	Array.fill(vmean, -1);
	vstd = Array.copy(vmean);
	selectWindow("ch3");
	for (i = 0; i < vmean.length; i++) {
		makeRectangle(bx[i], by[i], bw[i], bh[i]);
		getStatistics(area, mean, min, max, std);
		vmean[i] = mean;
		vstd[i] = std;
	}
	run("Select None");
	
	// Grabbing the background of the Bungaro channel
	bmean = Array.copy(bx);
	Array.fill(bmean, -1);
	bstd = Array.copy(bmean);
	selectWindow("ch4");
	getStatistics(area, mean, min, max, std);
	gstd = std;
	gmax = max;
	for (i = 0; i < bmean.length; i++) {
		makeRectangle(bx[i], by[i], bw[i], bh[i]);
		getStatistics(area, mean, min, max, std);
		bmean[i] = mean;
		if (std > 2*gstd) {
			bstd[i] = 1.5*gstd;
		} else {
			bstd[i] = std;
		}
	}
	run("Select None");
	
	// Duplicating the desired channels
	selectWindow("ch2");
	run("Duplicate...", "title=NFH");
	selectWindow("ch3");
	run("Duplicate...", "title=VAChT");
	selectWindow("ch4");
	run("Duplicate...", "title=Bung");
	
	// Cutting out the desired areas from the desired images
	roiManager("deselect");
	roiManager("delete");
	selectWindow("160");
	run("Analyze Particles...", "size=0-10000 pixel clear add");
	
	selectWindow("NFH");
	roiManager("combine");
	run("Clear Outside");
	run("Select None");
	roiManager("deselect");
	
	selectWindow("VAChT");
	roiManager("combine");
	run("Clear Outside");
	run("Select None");
	roiManager("deselect");
	
	selectWindow("Bung");
	roiManager("combine");
	run("Clear Outside");
	run("Select None");
	roiManager("deselect");
	roiManager("delete");
	
	// Merging the NFH and VAChT channel since either is sufficient for analysis
	imageCalculator("Max", "VAChT", "NFH");
	close("NFH");
	close("128");
	close("160");
	
	// Re-estabilishing the bounding boxes
	selectWindow("ch4");
	for (i = 0; i<bbx.length; i++) {
		makeRectangle(bbx[i], bby[i], bbw[i], bbh[i]);
		roiManager("add");
	}
	run("Select None");
	
	// Thresholding the localized bounding box for VAChT
	// Grabbing both the first and second std
	varea1 = newArray(roiManager("count"));
	varea2 = Array.copy(varea1);
	selectWindow("VAChT");
	for (i=0; i<roiManager("count"); i++) {
		roiManager("select", i);
		setThreshold(vmean[i]+vstd[i], gmax);
		roiManager("measure");
		val = Table.getColumn("%Area");
		varea1[i] = val[0];
		run("Clear Results");
		run("Select None");
		roiManager("deselect");
		resetThreshold();
		roiManager("select", i);
		setThreshold(vmean[i]+2*vstd[i], gmax);
		roiManager("measure");
		val = Table.getColumn("%Area");
		varea2[i] = val[0];
		run("Clear Results");
		run("Select None");
		roiManager("deselect");
		resetThreshold();
	}
	
	// Thresholding the localized bounding box for Bungaro
	selectWindow("Bung");
	barea1 = Array.copy(varea1);
	barea2 = Array.copy(varea1);
	for (i=0; i<roiManager("count"); i++) {
		roiManager("select", i);
		setThreshold(bmean[i]+bstd[i], gmax);
		roiManager("measure");
		val = Table.getColumn("%Area");
		barea1[i] = val[0];
		run("Clear Results");
		run("Select None");
		roiManager("deselect");
		resetThreshold();
		roiManager("select", i);
		setThreshold(bmean[i]+2*bstd[i], gmax);
		roiManager("measure");
		val = Table.getColumn("%Area");
		barea2[i] = val[0];
		run("Clear Results");
		run("Select None");
		roiManager("deselect");
		resetThreshold();
	}
	
	// Checking the VAChT against the Bungaro
	imageCalculator("Min create", "VAChT", "Bung");
	rename("Merge");
	marea1 = Array.copy(varea1);
	marea2 = Array.copy(varea1);
	for (i=0; i<roiManager("count"); i++) {
		roiManager("select", i);
		setThreshold(vmean[i]+vstd[i], gmax);
		roiManager("measure");
		val = Table.getColumn("%Area");
		marea1[i] = val[0];
		run("Clear Results");
		run("Select None");
		roiManager("deselect");
		resetThreshold();
		roiManager("select", i);
		setThreshold(vmean[i]+2*vstd[i], gmax);
		roiManager("measure");
		val = Table.getColumn("%Area");
		marea2[i] = val[0];
		run("Clear Results");
		run("Select None");
		roiManager("deselect");
		resetThreshold();
	}
	
	// Determining the coverage of the VAChT
	vcov11 = Array.copy(vmean);
	Array.fill(vcov11, -1);
	vcov12 = Array.copy(vcov11);
	vcov22 = Array.copy(vcov11);
	mcov11 = Array.copy(vcov11);
	mcov22 = Array.copy(vcov11);
	
	for (i=0; i<roiManager("count"); i++) {
		vcov11[i] = varea1[i]/barea1[i];
		vcov12[i] = varea1[i]/barea2[i];
		vcov22[i] = varea2[i]/barea2[i];
		mcov11[i] = marea1[i]/barea1[i];
		mcov22[i] = marea2[i]/barea2[i];	
	}
	
	// Displaying the results and saving
	Table.showArrays("Results", varea1, varea2, barea1, barea2, marea1, marea2, vcov11, vcov12, vcov22, mcov11, mcov22);
	spath = saveDir+org_img+"_Results.csv";
	saveAs("Results", spath);
	
	selectWindow("ch3");
	roiManager("show all");
	run("Flatten");
	rename(org_img);
	newImage("Stitch", "8-bit black", 2048, 2048, 1);
	run("Add Image...", "image="+org_img+" x=0 y=0 opacity=100");
	run("Add Image...", "image=ch4 x=1024 y=0 opacity=100");
	run("Add Image...", "image=VAChT x=0 y=1024 opacity=100");
	run("Add Image...", "image=Bung x=1024 y=1024 opacity=100");
	run("Flatten");
	rename("Combined");
	spath = roiDir+org_img+".tiff";
	saveAs("tiff", spath);
	
	close("*");
	roiManager("deselect");
	roiManager("delete");
	run("Clear Results");

}


// ***** MAIN *****

// Allowing user to choose a file or a folder of files to analyze
ff = getBoolean("Would you like to process a single file or a whole folder of files?", "File", "Folder");
// General setup
saveSettings();
setOption("BlackBackground", true);
run("Set Measurements...", "area mean standard centroid perimeter bounding fit shape median area_fraction redirect=None decimal=3");
setBatchMode(true);

start = getTime();
if (ff == 1) {
	fpath = File.openDialog("Select a File");
	imageDir = File.getParent(fpath);
	imgName = File.getName(fpath);
} else {
	// Getting the folder where images are located
	imageDir=getDirectory("Select Endplate Counting Folder for Image Processing");
	imagelist = getFileList(imageDir);
}

resultsDir=imageDir+"\\Analysis_Results\\";
roiDir=imageDir+"\\Individual_Images\\";
File.makeDirectory(resultsDir);
File.makeDirectory(roiDir);

if (ff == 1) {
	run("Bio-Formats Importer", "open=["+fpath+"] color_mode=Grayscale rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
	endplate(1, imgName, roiDir, resultsDir);
} else {
	// Grabbing each image individually for analysis
	for (image=0;image<imagelist.length;image++) {

	// Opening each image that is a oib image and getting the image title
	// This can be modified to whatever file ending is needed (could just be .oib)
		if (endsWith(imagelist[image], ".oib")) {
			imgName=imagelist[image];
			fpath = imageDir+imgName;
			run("Bio-Formats Importer", "open=["+fpath+"] color_mode=Grayscale rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
			endplate(image, imgName, roiDir, resultsDir);
		}
	}
}
endTime = (getTime()-start)/1000/60;
print("Analysis is finished! Took "+endTime+" minutes to run.");
setBatchMode(false);
close("Results");
close("ROI Manager");
