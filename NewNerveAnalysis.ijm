// Getting the folder where images are located
imageDir=getDirectory("Select Nerve Image Processing Folder");
roiDir=imageDir+"AxonList\\";
imagelist = getFileList(imageDir);
File.makeDirectory(roiDir);

for (z=0;z<imagelist.length;z++) {

// Opening each image that is a jpeg and getting the image title
// This can be modified to whatever file format you have	
	if (endsWith(imagelist[z], ".JPG")) {
		open(imageDir+imagelist[z]);
		imgName=getTitle();
		baseNameEnd=indexOf(imgName, ".JPG"); 
        baseName=substring(imgName, 0, baseNameEnd);
	}

// Dialog box for getting the values that you need to adjust
	Dialog.create("Image Scale");
	Dialog.addNumber("Known Distance in Pixels:", 641.0382);
	Dialog.addNumber("Known Distance in Micrometers:", 50);
	Dialog.show();
	
	scalebarpixellength = Dialog.getNumber();
	scalebarmicrometer = Dialog.getNumber();
	
// Defining the min and max axon sizes can adjust this as necessary
	minaxon = 1; // 1 um
	tinyaxon = 3; // axons less than 3 are very small
	smallaxon = 5; // this is used as a criterion for one of the sorting ratios
	largeaxon = 15; // axons greater than 15 are very big
	maxaxon = 25; // 25 um
	paratio = 0.4; // enforced on all axons greater than small
	
// Running the image analysis portion
// ---------------------------------------------
	if (z==0) {
	setOption("BlackBackground", false);
	run("Set Scale...", "distance="+scalebarpixellength+" known="+scalebarmicrometer+" unit=um");
	//run("Set Scale...", "distance=641.0382 known=50 unit=um");
	run("Set Measurements...", "area mean standard centroid perimeter fit shape stack display redirect=None decimal=3");
	}
	
// This section makes a copy of the image to work on
// and splits it into it's red, green, and blue channels
// and finds the edges -------------------------------
	selectWindow(imgName);
	run("Duplicate...", "title=Test");
	selectWindow("Test");
	run("Despeckle");
	run("Despeckle");
	run("Smooth");
	run("Enhance Local Contrast (CLAHE)", "blocksize=127 histogram=256 maximum=3 mask=*None*");
	run("Duplicate...", "title=ROIaxons");
	close(imgName);
	selectWindow("Test");
	run("Split Channels");
	selectWindow("Test (green)");
	rename("0");
	run("FeatureJ Edges", "compute smoothing=2 lower=[] higher=[]");
	run("8-bit");
	close("0");
	selectWindow("Test (red)");
	rename("1");
	run("FeatureJ Edges", "compute smoothing=2 lower=[] higher=[]");
	run("8-bit");
	close("1");
	selectWindow("Test (blue)");
	rename("2");
	run("FeatureJ Edges", "compute smoothing=2 lower=[] higher=[]");
	run("8-bit");
	close("2");

// Creating multiple axon images to ensure we have captured all of the axons
	for (k=0;k<3;k++){
		for (j=0;j<17;j++){
	  		val=180-j*10;
	  		selectWindow(k+" edges");
	  		run("Duplicate...", "title="+k+"A");
	  		selectWindow(k+"A");
	  		setThreshold(val, 255);
	  		run("Convert to Mask");
	  		rename("Axon"+k+"_"+j);
	  		selectWindow("Axon"+k+"_"+j);
	  		run("Invert LUT");
	  		run("Invert");
	  		run("Fill Holes");
	  		run("Analyze Particles...", "size="+minaxon+"-"+smallaxon+" circularity=0.60-1.00 show=Masks exclude clear record");
	  		rename("Axon"+k+" Result"+j);
	  		run("Invert LUT");
	  		run("Invert");
	  		selectWindow("Axon"+k+"_"+j);
  			Table.reset("Results");
	  		run("Analyze Particles...", "size="+smallaxon+"-"+maxaxon+" circularity=0.30-1.00 show=Masks exclude record");
	  		rename("Axon"+k+" Result");
	  		run("Invert LUT");
	  		run("Invert");
	  		for (i=0; i<Table.size; i++) {
				p = Table.get("Perim", i);
				a = Table.get("Area", i);
				Table.set("PA", i, p/a);
			}
			Table.update;
			e = Table.size;
			if(e>0) {
				run("Classify Particles", "class[1]=PA operator[1]=< value[1]="+paratio+" class[2]=-empty- operator[2]=-empty- value[2]=0.0000 class[3]=-empty- operator[3]=-empty- value[3]=0.0000 class[4]=-empty- operator[4]=-empty- value[4]=0.0000 combine=[AND (match all)] output=[Keep members] white");
			}
	  		
	  		imageCalculator("AND", "Axon"+k+" Result"+j,"Axon"+k+" Result");
	  		close("Axon"+k+" Result");
	  		close(k+"A");
	  		close("Axon"+k+"_"+j);
	  		imageCalculator("AND", ""+k+" edges", "Axon"+k+" Result"+j);
		}
	
		close(k+" edges");
		
	// This merges the resulting axon images together creating a no overlap
	// and a overlap set of axons ------------------------------
		
		for (j=1;j<17;j++){
			imageCalculator("OR create", "Axon"+k+" Result0","Axon"+k+" Result"+j);
			rename("Axon"+k+" ResultB"+j);
			imageCalculator("XOR", "Axon"+k+" Result0","Axon"+k+" Result"+j);
			close("Axon"+k+" Result"+j);
			selectWindow("Axon"+k+" Result0");
			run("Invert");
			run("Despeckle");
			run("Despeckle");
			run("Fill Holes");
	  		run("Analyze Particles...", "size="+minaxon+"-"+smallaxon+" circularity=0.70-1.00 show=Masks exclude clear");
	  		rename("Axon"+k+" Result"+j);
	  		run("Invert LUT");
	  		run("Invert");
	  		selectWindow("Axon"+k+" Result0");
	  		run("Analyze Particles...", "size="+smallaxon+"-"+maxaxon+" circularity=0.40-1.00 show=Masks exclude clear");
			rename("Axon 0");
			run("Invert LUT");
	  		run("Invert");
	  		imageCalculator("AND", "Axon 0","Axon"+k+" Result"+j);
	  		close("Axon"+k+" Result"+j);
	  		close("Axon"+k+" Result0");
	  		selectWindow("Axon 0");
	  		rename("Axon"+k+" Result0");
	  		imageCalculator("OR create", "Axon"+k+" Result0","Axon"+k+" ResultB"+j);
			rename("Axon"+k+"B"+j);
			run("Despeckle");
			run("Despeckle");
	  		imageCalculator("XOR", "Axon"+k+" Result0","Axon"+k+" ResultB"+j);
	  		close("Axon"+k+" ResultB"+j);
			selectWindow("Axon"+k+" Result0");
			run("Invert");
			run("Despeckle");
			run("Despeckle");
			if (j>1) {
				imageCalculator("AND", "Axon"+k+"B1","Axon"+k+"B"+j);
				close("Axon"+k+"B"+j);
			}
		}		
	}


// This merges the resulting axon images for each of the color channels
// together for the no overlap and overlap set of axons -------------
	selectWindow("Axon0B1");
	rename("AxonB");
	selectWindow("Axon0 Result0");
	rename("Axon");
	run("Fill Holes");
	run("Erode");
	run("Erode");
	run("Dilate");
	run("Dilate");
	run("Erode");
	run("Erode");
	run("Dilate");
	run("Dilate");	
	
	for (k=1;k<3;k++) {
		imageCalculator("AND", "AxonB","Axon"+k+"B1");
		close("Axon"+k+"B1");
		selectWindow("Axon"+k+" Result0");
		run("Fill Holes");
		run("Erode");
		run("Erode");
		run("Dilate");
		run("Dilate");
		run("Erode");
		run("Erode");
		run("Dilate");
		run("Dilate");
		imageCalculator("OR create", "Axon","Axon"+k+" Result0");
		rename("Axon"+k+"B");
		run("Despeckle");
		run("Despeckle");
		run("Analyze Particles...", "size="+minaxon+"-"+smallaxon+" circularity=0.70-1.00 show=Masks exclude clear");
	  	rename("AxonB 0");
	  	run("Invert LUT");
	  	run("Invert");
	  	selectWindow("Axon"+k+"B");
	  	run("Analyze Particles...", "size="+smallaxon+"-"+maxaxon+" circularity=0.40-1.00 show=Masks exclude clear");
		rename("AxonB 1");
		run("Invert LUT");
	  	run("Invert");
	  	imageCalculator("AND", "AxonB 0","AxonB 1");
	  	close("AxonB 1");
	  	close("Axon"+k+"B");
	  	selectWindow("AxonB 0");
	  	rename("Axon"+k+"B");
			
		imageCalculator("XOR", "Axon","Axon"+k+" Result0");
		close("Axon"+k+" Result0");
		selectWindow("Axon");
		run("Invert");
		run("Despeckle");
		run("Despeckle");
		run("Fill Holes");
		run("Analyze Particles...", "size="+minaxon+"-"+smallaxon+" circularity=0.70-1.00 show=Masks exclude clear");
	  	rename("Axon0");
	  	run("Invert LUT");
	  	run("Invert");
	  	selectWindow("Axon");
	  	run("Analyze Particles...", "size="+smallaxon+"-"+maxaxon+" circularity=0.40-1.00 show=Masks exclude clear");
		rename("Axon1");
		run("Invert LUT");
	  	run("Invert");
	  	imageCalculator("AND", "Axon0","Axon1");
	  	close("Axon1");
	  	close("Axon");
	  	selectWindow("Axon0");
	  	rename("Axon");
	}

// This merges the axon images down to a single no overlap
// and a single overlap axon image -----------------------	
	imageCalculator("OR create", "Axon1B","Axon2B");
	rename("AxonL");
	run("Despeckle");
	run("Despeckle");
	run("Analyze Particles...", "size="+minaxon+"-"+smallaxon+" circularity=0.70-1.00 show=Masks exclude clear");
	rename("AxonL1");
	run("Invert LUT");
	run("Invert");
	selectWindow("AxonL");
	run("Analyze Particles...", "size="+smallaxon+"-"+maxaxon+" circularity=0.40-1.00 show=Masks exclude clear");
	rename("AxonL2");
	run("Invert LUT");
	run("Invert");
	imageCalculator("AND", "AxonL1","AxonL2");
	close("AxonL2");
	close("AxonL");
	selectWindow("AxonL1");
	rename("AxonL");
	imageCalculator("AND", "AxonB","AxonL");
	close("AxonL");

	imageCalculator("XOR", "Axon1B","Axon2B");
	close("Axon2B");
	selectWindow("Axon1B");
	run("Invert");
	run("Despeckle");
	run("Despeckle");
	run("Fill Holes");
	run("Analyze Particles...", "size="+minaxon+"-"+smallaxon+" circularity=0.70-1.00 show=Masks exclude clear");
	rename("Axon0");
	run("Invert LUT");
	run("Invert");
	selectWindow("Axon1B");
	run("Analyze Particles...", "size="+smallaxon+"-"+maxaxon+" circularity=0.40-1.00 show=Masks exclude clear");
	rename("Axon1");
	run("Invert LUT");
	run("Invert");
	imageCalculator("AND", "Axon0","Axon1");
	close("Axon1");
	close("Axon1B");
	imageCalculator("AND", "AxonB","Axon0");
	close("Axon0");

	imageCalculator("Subtract create", "Axon","AxonB");
	selectWindow("Result of Axon");
	run("Invert");
	run("Despeckle");
	run("Despeckle");
	run("Despeckle");
	imageCalculator("AND", "Axon","Result of Axon");
	selectWindow("Result of Axon");
	imageCalculator("XOR", "AxonB","Result of Axon");
	close("Result of Axon");
	selectWindow("AxonB");
	run("Invert");

	
// This saves the resulting selected axons
// and the resulting image with the axons overlayed	
	selectWindow("AxonB");
	run("Analyze Particles...", "size="+minaxon+"-"+maxaxon+" circularity=0.30-1.00 show=Masks exclude clear add");
	rename("ABSave");
	roiManager("Show None");
	roiManager("Set Color", "black");
	roiManager("Set Line Width", 0);
	selectWindow("ROIaxons");
	roiManager("Show All with labels");
	run("Flatten");
	rename("AxonB Outlines");
	selectWindow("ROIaxons");
	roiManager("Show None");
	roiManager("Save", roiDir+baseName+"_RoiSet_Axon_B.zip");
	roiManager("Delete");
	selectWindow("AxonB Outlines");
	saveAs("jpg", roiDir+baseName+"_AxonB_Outlines.jpg");
	selectWindow("ABSave");
	saveAs("jpg", roiDir+baseName+"_AxonB.jpg");
	close("ABSave");
	close("AxonB");
	
	selectWindow("Axon");
	run("Analyze Particles...", "size="+minaxon+"-"+tinyaxon+" circularity=0.30-1.00 show=Masks exclude clear add");
	rename("Tiny Axon");
	roiManager("Show None");
	roiManager("Set Color", "black");
	roiManager("Set Line Width", 0);
	selectWindow("ROIaxons");
	roiManager("Show All with labels");
	run("Flatten");
	rename("Tiny Axon Outlines");
	selectWindow("ROIaxons");
	roiManager("Show None");
	roiManager("Save", roiDir+baseName+"_RoiSet_Tiny_Axon.zip");
	roiManager("Delete");
	selectWindow("Tiny Axon Outlines");
	saveAs("jpg", roiDir+baseName+"_Tiny_Axon_Outlines.jpg");
	selectWindow("Tiny Axon");
	saveAs("jpg", roiDir+baseName+"_Tiny_Axon.jpg");
	
	selectWindow("Axon");
	run("Analyze Particles...", "size="+tinyaxon+"-"+smallaxon+" circularity=0.30-1.00 show=Masks exclude clear add");
	rename("Small Axon");
	roiManager("Show None");
	roiManager("Set Color", "black");
	roiManager("Set Line Width", 0);
	selectWindow("ROIaxons");
	roiManager("Show All with labels");
	run("Flatten");
	rename("Small Axon Outlines");
	selectWindow("ROIaxons");
	roiManager("Show None");
	roiManager("Save", roiDir+baseName+"_RoiSet_Small_Axon.zip");
	roiManager("Delete");
	selectWindow("Small Axon Outlines");
	saveAs("jpg", roiDir+baseName+"_Small_Axon_Outlines.jpg");
	selectWindow("Small Axon");
	saveAs("jpg", roiDir+baseName+"_Small_Axon.jpg");
	
	selectWindow("Axon");
	run("Analyze Particles...", "size="+smallaxon+"-"+largeaxon+" circularity=0.30-1.00 show=Masks exclude clear add");
	rename("N Axon");
	roiManager("Show None");
	roiManager("Set Color", "black");
	roiManager("Set Line Width", 0);
	selectWindow("ROIaxons");
	roiManager("Show All with labels");
	run("Flatten");
	rename("Axon Outlines");
	selectWindow("ROIaxons");
	roiManager("Show None");
	roiManager("Save", roiDir+baseName+"_RoiSet_Axon.zip");
	roiManager("Delete");
	selectWindow("Axon Outlines");
	saveAs("jpg", roiDir+baseName+"_Axon_Outlines.jpg");
	selectWindow("N Axon");
	saveAs("jpg", roiDir+baseName+"_Axon.jpg");

	selectWindow("Axon");
	run("Analyze Particles...", "size="+largeaxon+"-"+maxaxon+" circularity=0.30-1.00 show=Masks exclude clear add");
	rename("Large Axon");
	roiManager("Show None");
	roiManager("Set Color", "black");
	roiManager("Set Line Width", 0);
	selectWindow("ROIaxons");
	roiManager("Show All with labels");
	run("Flatten");
	rename("Large Axon Outlines");
	selectWindow("ROIaxons");
	roiManager("Show None");
	roiManager("Save", roiDir+baseName+"_RoiSet_Large_Axon.zip");
	roiManager("Delete");
	selectWindow("Large Axon Outlines");
	saveAs("jpg", roiDir+baseName+"_Large_Axon_Outlines.jpg");
	selectWindow("Large Axon");
	saveAs("jpg", roiDir+baseName+"_Large_Axon.jpg");

	selectWindow("Axon");
	run("Analyze Particles...", "size="+minaxon+"-"+maxaxon+" circularity=0.30-1.00 show=Masks exclude clear add");
	rename("MM Axon");
	roiManager("Show None");
	roiManager("Set Color", "black");
	roiManager("Set Line Width", 0);
	selectWindow("ROIaxons");
	roiManager("Show All with labels");
	run("Flatten");
	rename("MinMax Axon Outlines");
	selectWindow("ROIaxons");
	roiManager("Show None");
	roiManager("Save", roiDir+baseName+"_RoiSet_MinMax_Axon.zip");
	roiManager("Delete");
	selectWindow("MinMax Axon Outlines");
	saveAs("jpg", roiDir+baseName+"_MinMax_Axon_Outlines.jpg");
	selectWindow("MM Axon");
	saveAs("jpg", roiDir+baseName+"_MinMax_Axon.jpg");
	
	run("Close All");
	print("Finished analyzing image: "+baseName);
}


	