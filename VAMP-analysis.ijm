// Getting the file or folder where images are located
items = newArray("file", "folder");
Dialog.create("File or Folder Analysis");
Dialog.addChoice("Would you like to analyze a whole folder or an individual file?", items, "folder");
Dialog.addMessage("When processing a whole folder, images need to have the same markers in the channels.\nOtherwise the analysis may not work!!");
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
	//setBatchMode(true);
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
run("Set Measurements...", "area mean standard centroid fit median stack redirect=None decimal=3");	

start = getTime();
// Grabbing each image individually for analysis
//for (image=0;image<imagelist.length;image++) {
for (image=0; image<1; image++) {

// This assumes that the image file ending is an .oib
// And that all images in the selected folder should be analyzed
	imgName=imagelist[image];
	baseNameEnd=indexOf(imgName, ".oib"); 
    baseName=substring(imgName, 0, baseNameEnd);
    savePath = roiDir+baseName;
	fpath = imageDir+imgName;
	run("Bio-Formats Importer", "open=["+fpath+"] color_mode=Grayscale rois_import=[ROI manager] split_channels view=Hyperstack stack_order=XYCZT");

	chanArray = newArray("Fast Blue","VGAT/VIAAT","VAMP1");
	for (i = 1; i <= nImages; i++) {	
		selectImage(i);
		if (chanArray[i-1] == "Fast Blue") {
			rename("FB");
		} else if (chanArray[i-1] == "NeuN" || chanArray[i-1] == "ChAT") {
			rename("NC");
		} else if (chanArray[i-1] == "VGLUT1" || chanArray[i-1] == "VGLUT2" || chanArray[i-1] == "VGAT/VIAAT") {
			rename("Vesicular");
		} else if (chanArray[i-1] == "Other") {
			rename("Other");
		} else {
			rename("Vamp");
		}
		slices = nSlices;
		for (j=1; j<=slices; j++) {
	   		setSlice(j);
	   		getRawStatistics(nPixels, mean, min, max, std, histogram);
	   		percent = Array.copy(histogram);
	    	percent[0] = histogram[0]/nPixels;
	    	bottom = 0;
	    	multiply = 1;
	    	for (k=1; k<histogram.length; k++) {
	    		histogram[k] = histogram[k-1]+histogram[k];
	    		percent[k] = histogram[k]*100/nPixels;
	    		if (percent[k]<=10.1) {
	    			bottom = k;
	    		}
	    		if (percent[k]>94 && percent[k]<=97 && k>300 && chanArray[i-1]=="Fast Blue" || chanArray[i-1]=="NeuN" || chanArray[i-1]=="ChAT") {
	    			multiply = 4095/(k-bottom);
	    		} else if (percent[k]>94 && percent[k]<=97 && k<300 && chanArray[i-1]=="Fast Blue" || chanArray[i-1]=="NeuN" || chanArray[i-1]=="ChAT") {
	    			multiply = 0.5;
	    		} else if (percent[k]>90 && percent[k]<=93 && chanArray[i-1]!="Fast Blue" && chanArray[i-1]!="NeuN" && chanArray[i-1]!="ChAT") {
	    			multiply = 4095/(k-bottom);
	    		}
	    	}
	    	run("Subtract...", "value="+bottom+" slice");
	    	run("Multiply...", "value="+multiply+" slice");
		}
		run("Max...", "value=4095 stack");			
		run("8-bit");
	}
	
	xScaled = 1;
	yScaled = 1;
	zScaled = 1;
	toScaled(xScaled,yScaled,zScaled);
	if (xScaled == yScaled) {
		run("Set Scale...", "distance="+(1/xScaled)+" known=1 unit=um global");
	} else {
		run("Set Scale...", "distance="+(1/xScaled)+" known=1 pixel="+(yScaled/xScaled)+" unit=um global");
	}

	imageCalculator("Add create stack", "Vesicular","Vamp");
	rename("nucli");
	setThreshold(0,0);
	run("Convert to Mask", "method=MinError background=Dark black");
	run("Duplicate...", "title=copy duplicate");
	run("Options...", "iterations=2 count=4 black do=Erode stack");
	run("Options...", "iterations=3 count=3 black do=Open stack");
	run("Options...", "iterations=5 count=1 black do=Dilate stack");
	run("Options...", "iterations=5 count=2 black do=Dilate stack");
	run("Options...", "iterations=10 count=3 black do=Dilate stack");
	run("Options...", "iterations=20 count=4 black do=Dilate stack");
	imageCalculator("AND stack", "nucli","copy");
	close("copy");
	run("Options...", "iterations=3 count=3 black do=Close stack");
	run("Analyze Particles...", "size=50-250 circularity=0.2-1.00 show=[Bare Outlines] clear stack");
	r = Table.getColumn("Round");
	s = Table.getColumn("Slice");
	x = Table.getColumn('X');
	y = Table.getColumn('Y');
	toUnscaled(x, y);
	for (i = 1; i < s[0]; i++) {
		selectWindow("nucli");
		setSlice(i);
		run("Select All");
		run("Clear", "slice");
		run("Select None");
	}
	//slices=15;
	for (i=s[s.length-1]+1; i<=slices; i++) {
		selectWindow("nucli");
		setSlice(i);
		run("Select All");
		run("Clear", "slice");
		run("Select None");						
	}
	for (i = 0; i < r.length; i++) {
    	selectWindow("Drawing of nucli");
    	run("Select None");
    	setSlice(s[i]);
		k = i;
	    for (j=i; j<r.length; j++) {
	    	if (s[j]==s[i] && r[j] >= 0.5) {
    			setKeyDown("shift");
    			doWand(x[j],y[j]);
    			setKeyDown("none");
    			k=j;
	    	} else if (s[j]==s[i]) {
	    		k=j;
	    	}
	    }
    	selectWindow("nucli");
    	setSlice(s[i]);
    	run("Restore Selection");
    	wait(10);
    	run("Clear Outside", "slice");
    	wait(10);
    	run("Select None");
    	i = k;
	}
	close("Drawing of nucli");
	selectWindow("nucli");
	run("Fill Holes", "stack");
	run("Mean...", "radius=15 stack");
	run("Convert to Mask", "method=Minimum background=Dark calculate black");
	close("Log");			
	run("Options...", "iterations=10 count=1 black do=Dilate stack");
	run("Options...", "iterations=20 count=2 black do=Dilate stack");
	run("Options...", "iterations=10 count=3 black do=Dilate stack");
	run("Analyze Particles...", "size=100-450 circularity=0.6-1.00 show=Ellipses display exclude clear stack");
	s = Table.getColumn("Slice");
	x = Table.getColumn('X');
	y = Table.getColumn('Y');
	toUnscaled(x, y);
	c = newArray(nResults);
	run("Invert", "stack");
	run("Fill Holes", "stack");
	run("Options...", "iterations=10 count=2 black do=Erode stack");
	rename("nucleus");
	close("nucli");
	close("Results");
	run("Duplicate...", "title=blank");
	run("Select All");
	run("Clear");
	run("Select None");
	c[0] = 1;
	setForegroundColor(1,1,1);
	makeRectangle(x[0]-50, y[0]-50, 100, 100);
	run("Fill", "slice");
	numCells = 1;
	updateDisplay();
	for (i = 1; i < c.length; i++) {
		if (getPixel(x[i], y[i]) != 0) {
			numCells = getPixel(x[i],y[i]);
		} else {
			Array.getStatistics(c, min, max);
			numCells = max+1;
		}
	    setForegroundColor(numCells,numCells,numCells);
		makeRectangle(x[i]-50, y[i]-50, 100, 100);
		run("Fill", "slice");
		c[i] = numCells;
		updateDisplay();
	}
	close("blank");
	Array.sort(c,s,x,y);
	setForegroundColor(255,255,255);
	for (i = 1; i < c.length; i++) {
		if (c[i] == c[i-1] && s[i] != s[i-1]+1) {
			selectWindow("nucleus");
			setSlice(s[i-1]);
			doWand(x[i-1], y[i-1]);
			run("Fit Ellipse");
			for (j=s[i-1]+1; j<s[i]; j++) {
				setSlice(j);
				run("Restore Selection");
				run("Fill", "slice");
				run("Select None");
			}
		}
	}
	run("Select None");
	run("Duplicate...", "title=nucli duplicate");
	run("Options...", "iterations=5 count=1 black do=Dilate stack");
	run("Options...", "iterations=5 count=2 black do=Dilate stack");
	run("Options...", "iterations=5 count=3 black do=Dilate stack");
	run("Options...", "iterations=5 count=4 black do=Dilate stack");	
	imageCalculator("Subtract stack", "Vesicular","nucli");
	imageCalculator("Subtract stack", "Vamp","nucli");
	
	fastblue=false;
	if (isOpen("FB")) {
		setForegroundColor(255, 255, 255);
		setBackgroundColor(0, 0, 0);
		selectWindow("FB");
		run("Duplicate...", "title=Finding duplicate");
		run("Convert to Mask", "method=Otsu background=Dark calculate black");
		run("Options...", "iterations=3 count=3 black do=Open stack");
		run("Fill Holes", "stack");
		run("Analyze Particles...", "size=500-2000 show=Masks clear stack");
		run("Invert LUT");
		rename("FBcells");
		close("Finding");
		if (nResults >= 5) {
			s = Table.getColumn('Slice');
			x = Table.getColumn('X');
			y = Table.getColumn('Y');
			toUnscaled(x, y);
			//slices = nSlices;
			for (j = 1; j < s[0]; j++) {
				selectWindow("FBcells");
				setSlice(j);
				run("Duplicate...", "title=Ave"+j);
			}
			for (j = s[0]; j <= s[s.length-1]; j++) {
				selectWindow("FBcells");
				if (j < 2) {
			    	run("Z Project...", "start="+j-1+" stop="+j+7+" projection=[Average Intensity]");
				} else if (j < 5) {
			    	run("Z Project...", "start="+j-2+" stop="+j+6+" projection=[Average Intensity]");
				} else if (j < s.length-5) {
					run("Z Project...", "start="+j-4+" stop="+j+4+" projection=[Average Intensity]");
				} else if (j < s.length-3) {
					run("Z Project...", "start="+j-6+" stop="+j+2+" projection=[Average Intensity]");
				} else {
					run("Z Project...", "start="+j-7+" stop="+j+1+" projection=[Average Intensity]");
				}
				rename("Ave"+j);
			}
			for (j = s[s.length-1]+1; j <= slices; j++) {
				selectWindow("FBcells");
				setSlice(j);
				run("Duplicate...", "title=Ave"+j);
			}
			close("FBcells");
			run("Images to Stack", "name=FBcells title=Ave");
			setThreshold(140,255);
			run("Convert to Mask", "method=Default background=Dark black");
			for (j = 0; j < s.length; j++) {
				selectWindow("FBcells");
				setSlice(s[j]);
				doWand(x[j],y[j]);
				k = j;
				for (i = j+1; i < s.length; i++) {
					if (s[i] == s[j]) {
						setKeyDown("shift");
						doWand(x[i], y[i]);
						k = i;
					}
				}
				run("Enlarge...", "enlarge=-4");
				wait(50);
				run("Enlarge...", "enlarge=4.5");
				wait(50);
				run("Clear Outside", "slice");
				run("Fill Holes", "slice");
				run("Select None");
				doWand(x[j],y[j]);
				for (i = j+1; i < s.length; i++) {
					if (s[i] == s[j]) {
						setKeyDown("shift");
						doWand(x[i], y[i]);
					}
				}
				run("Enlarge...", "enlarge=-2");
				wait(50);
				run("Enlarge...", "enlarge=1");
				wait(50);
				run("Clear Outside", "slice");
				run("Fill Holes", "slice");
				run("Select None");
				j = k;
			}
			imageCalculator("AND create stack", "FBcells", "nucli");
			rename("selected");
			run("Analyze Particles...", "size=100-1000 show=Nothing clear stack");
			s = Table.getColumn('Slice');
			x = Table.getColumn('X');
			y = Table.getColumn('Y');
			toUnscaled(x, y);
			for (j = 1; j < s[0]; j++) {
				selectWindow("FBcells");
				setSlice(j);
				run("Select All");
				run("Clear", "slice");
				run("Select None");
			}
			for (j = 0; j < s.length; j++) {
				selectWindow("FBcells");
				setSlice(s[j]);
				doWand(x[j],y[j]);
				k = j;
				for (i = j+1; i < s.length; i++) {
					if (s[i] == s[j]) {
						setKeyDown("shift");
						doWand(x[i], y[i]);
						setKeyDown("none");
						k = i;
					}
				}
				run("Clear Outside", "slice");
				run("Select None");
				j = k;
			}
			for (j = s[s.length-1]+1; j <= slices; j++) {
				selectWindow("FBcells");
				setSlice(j);
				run("Select All");
				run("Clear", "slice");
				run("Select None");
			}
			close("selected");
			selectWindow("FBcells");
			run("Duplicate...", "title=celledges duplicate");
			run("Analyze Particles...", "size=500-2500 show=Nothing display clear stack");
			s = Table.getColumn('Slice'); 
			x = Table.getColumn('X');
			y = Table.getColumn('Y');
			toUnscaled(x, y);
			close("Results");
			for (j = 0; j < s.length; j++) {
				selectWindow("celledges");
				setSlice(s[j]);
				doWand(x[j],y[j]);
				k = j;
				for (i = j+1; i < s.length; i++) {
					if (s[i] == s[j]) {
						setKeyDown("shift");
						doWand(x[i], y[i]);
						setKeyDown("none");
						k = i;
					}
				}
				run("Enlarge...", "enlarge=3");
				wait(50);
				setForegroundColor(255, 255, 255);
				run("Fill", "slice");
				setBackgroundColor(0, 0, 0);
				run("Clear Outside", "slice");
				run("Select None");
				j = k;
			}
			for (j = 0; j < s.length; j++) {
				selectWindow("FBcells");
				setSlice(s[j]);
				doWand(x[j],y[j]);
				k = j;
				for (i = j+1; i < s.length; i++) {
					if (s[i] == s[j]) {
						setKeyDown("shift");
						doWand(x[i], y[i]);
						setKeyDown("none");
						k = i;
					}
				}
				run("Enlarge...", "enlarge=-1");
				wait(50);
				run("Clear Outside", "slice");
				run("Select None");
				j = k;
			}
			rename("cells");
			run("Fill Holes", "stack");
			imageCalculator("Subtract stack", "celledges", "cells");
			close("nucli");
			close("nucleus");
			fastblue = true;
		} else {
			close("FBcells");
		}
	}

	if (!fastblue) {
		setForegroundColor(255, 255, 255);
		setBackgroundColor(0, 0, 0);			
		selectWindow("Vamp");
		//slices = nSlices;
		run("Duplicate...", "duplicate");
		run("Convert to Mask", "method=Default background=Dark calculate black");
		run("Analyze Particles...", "size=0-0.6 circularity=0.0-1.00 show=Masks clear stack");
		run("Invert LUT");
		run("Options...", "iterations=1 count=2 black do=Dilate stack");
		imageCalculator("Subtract stack", "Vamp","Mask of Vamp-1");
		close("Vamp-1");
		close("Mask of Vamp-1");
		selectWindow("Vesicular");
		run("Duplicate...", "duplicate");
		run("Convert to Mask", "method=Default background=Dark calculate black");
		run("Analyze Particles...", "size=0.00-0.60 show=Masks clear stack");
		run("Invert LUT");
		run("Options...", "iterations=1 count=2 black do=Dilate stack");
		imageCalculator("Subtract stack", "Vesicular","Mask of Vesicular-1");
		close("Vesicular-1");
		close("Mask of Vesicular-1");
		
		selectWindow("Vamp");
		run("Duplicate...", "title=raw duplicate");
		run("Gaussian Blur...", "sigma=15 stack");
		run("Find Edges", "stack");
		run("Convert to Mask", "method=Huang background=Dark calculate black");
		run("Options...", "iterations=10 count=3 black do=Close stack");
		run("Skeletonize", "stack");
		run("Options...", "iterations=5 count=3 black do=Dilate stack");
		run("Options...", "iterations=1 count=3 black do=Open stack");
		run("Analyze Particles...", "size=2-Infinity show=Masks clear stack");
		run("Invert LUT");
		rename("rawline");
		close("raw");
		imageCalculator("Subtract stack", "rawline","nucli");
		selectWindow("rawline");
		run("Duplicate...", "title=minline duplicate");
		run("Options...", "iterations=20 count=4 black do=Close stack");
		run("Skeletonize", "stack");
		run("Options...", "iterations=3 count=3 black do=Dilate stack");
		imageCalculator("Add stack", "Vesicular","minline");
		close("minline");
		selectWindow("rawline");
		run("Duplicate...", "title=line2 duplicate");
		run("Options...", "iterations=3 count=3 black do=Dilate stack");
		run("Options...", "iterations=2 count=2 black do=Dilate stack");
		run("Options...", "iterations=2 count=1 black do=Dilate stack");
		run("Invert", "stack");
		run("Fill Holes", "stack");
		run("Analyze Particles...", "size=200-infinity show=Masks clear stack");
		run("Invert LUT");
		rename("celled1");
		close("line2");
		selectWindow("rawline");
		run("Options...", "iterations=10 count=3 black do=Dilate stack");
		imageCalculator("Subtract stack", "rawline","nucli");		
		
		selectWindow("Vamp");
		run("Duplicate...", "title=1stpass duplicate");
		run("Convert to Mask", "method=RenyiEntropy background=Dark calculate black");
		run("Options...", "iterations=3 count=4 black do=Erode stack");
		run("Analyze Particles...", "size=0.00-1.00 show=Masks clear stack");
		run("Invert LUT");
		imageCalculator("Subtract stack", "1stpass","Mask of 1stpass");
		close("Mask of 1stpass");
		selectWindow("1stpass");
		run("Options...", "iterations=12 count=3 black do=Dilate stack");
		run("Skeletonize", "stack");
		run("Options...", "iterations=5 count=3 black do=Dilate stack");
		run("Options...", "iterations=1 count=3 black do=Open stack");
		run("Analyze Particles...", "size=2-Infinity show=Masks clear stack");
		run("Invert LUT");
		rename("3line");
		close("1stpass");
		imageCalculator("Subtract stack", "3line","nucli");
		selectWindow("3line");
		run("Duplicate...", "title=minline duplicate");
		run("Options...", "iterations=20 count=4 black do=Close stack");
		run("Skeletonize", "stack");
		run("Options...", "iterations=3 count=3 black do=Dilate stack");
		imageCalculator("Add stack", "Vamp","minline");
		close("minline");
		selectWindow("3line");
		run("Duplicate...", "title=line3 duplicate");
		run("Options...", "iterations=3 count=3 black do=Dilate stack");
		run("Options...", "iterations=2 count=2 black do=Dilate stack");
		run("Options...", "iterations=2 count=1 black do=Dilate stack");
		run("Invert", "stack");
		run("Fill Holes", "stack");
		run("Analyze Particles...", "size=200-infinity show=Masks clear stack");
		run("Invert LUT");
		rename("celled2");
		imageCalculator("Max create stack", "Vesicular","Vamp");
		rename("max");
		selectWindow("3line");
		run("Options...", "iterations=10 count=3 black do=Dilate stack");
		imageCalculator("Subtract stack", "3line","nucli");					
		imageCalculator("AND stack", "3line", "rawline");
		imageCalculator("Add stack", "max", "3line");
		close("rawline");
		close("3line");
		close("line3");
		
		imageCalculator("Min create stack", "Vesicular","Vamp");
		rename("min");
		run("Morphological Filters (3D)", "operation=Closing element=Octagon x-radius=5 y-radius=5 z-radius=2");
		rename("filter");
		run("Duplicate...", "title=dup duplicate");
		run("Duplicate...", "title=dup2 duplicate");
		run("Multiply...", "value=2.000 stack");
		run("Gaussian Blur...", "sigma=10 stack");
		selectWindow("dup");
		run("Subtract...", "value=25 stack");
		run("Multiply...", "value=2.000 stack");
		run("Mean...", "radius=10 stack");
		imageCalculator("Average stack", "dup","dup2");
		selectWindow("filter");
		run("Mean...", "radius=10 stack");
		run("Subtract...", "value=25 stack");
		run("Multiply...", "value=1.25 stack");
		imageCalculator("Average stack","filter", "dup");
		close("dup");
		close("dup2");
		run("Duplicate...", "title=filteredge duplicate");
		run("Find Edges", "stack");
		
		imageCalculator("Average create stack", "Vesicular","Vamp");
		rename("ave");
		run("Morphological Filters (3D)", "operation=Gradient element=Ball x-radius=3 y-radius=3 z-radius=2");
		run("Mean...", "radius=10 stack");
		run("Find Edges", "stack");
		rename("gradedge");
		imageCalculator("Add stack", "gradedge","filteredge");
		close("filteredge");
		run("Convert to Mask", "method=Mean background=Dark calculate black");
		run("Options...", "iterations=5 count=3 black do=Open stack");
		run("Options...", "iterations=10 count=3 black do=Dilate stack");
		run("Skeletonize", "stack");
		run("Options...", "iterations=5 count=3 black do=Dilate stack");
		run("Options...", "iterations=1 count=3 black do=Open stack");
		imageCalculator("Subtract stack", "gradedge","nucli");
		run("Options...", "iterations=3 count=3 black do=Dilate stack");
		run("Options...", "iterations=2 count=2 black do=Dilate stack");
		run("Options...", "iterations=2 count=1 black do=Dilate stack");
		run("Invert", "stack");
		run("Fill Holes", "stack");
		run("Analyze Particles...", "size=200-infinity show=Masks clear stack");
		run("Invert LUT");
		rename("celled3");
		close("gradedge");
						
		selectWindow("ave");
		run("Morphological Filters (3D)", "operation=Laplacian element=Octagon x-radius=10 y-radius=10 z-radius=3");
		rename("laplacian");
		run("Duplicate...", "title=lapedge duplicate");
		run("Mean...", "radius=10 stack");
		run("Find Edges", "stack");
		run("Convert to Mask", "method=Default background=Dark calculate black");
		run("Options...", "iterations=2 count=3 black do=Open stack");
		run("Options...", "iterations=10 count=3 black do=Dilate stack");
		run("Skeletonize", "stack");
		run("Options...", "iterations=5 count=3 black do=Dilate stack");
		run("Options...", "iterations=1 count=3 black do=Open stack");
		imageCalculator("Subtract stack", "lapedge","nucli");		
		run("Options...", "iterations=3 count=3 black do=Dilate stack");
		run("Options...", "iterations=2 count=2 black do=Dilate stack");
		run("Options...", "iterations=2 count=1 black do=Dilate stack");
		run("Invert", "stack");
		run("Fill Holes", "stack");
		run("Analyze Particles...", "size=200-infinity show=Masks clear stack");
		run("Invert LUT");
		rename("celled4");
		close("lapedge");
		
		selectWindow("laplacian");
		run("Invert", "stack");
		imageCalculator("Subtract stack", "laplacian","nucleus");
		run("Duplicate...", "duplicate");
		run("Subtract...", "value=50 stack");
		imageCalculator("Average stack", "laplacian-1","filter");
		imageCalculator("Average stack", "laplacian-1","filter");
		run("Multiply...", "value=1.5 stack");
		selectWindow("laplacian");
		run("Convert to Mask", "method=Shanbhag background=Dark calculate black");
		run("Despeckle", "stack");
		run("Despeckle", "stack");
		run("Options...", "iterations=5 count=4 black do=Dilate stack");
		run("Analyze Particles...", "size=0-80 pixel show=Masks clear stack");
		run("Invert LUT");
		imageCalculator("Subtract stack", "laplacian","Mask of laplacian");
		close("Mask of laplacian");
		selectWindow("laplacian");
		run("Duplicate...", "title=celled5 duplicate");
		run("Invert", "stack");
		selectWindow("laplacian");
		run("Subtract...", "value=150 stack");
		imageCalculator("Add stack", "laplacian-1","laplacian");
		selectWindow("laplacian");
		run("Multiply...", "value=2 stack");
		imageCalculator("Add stack", "max","laplacian");
		imageCalculator("Add stack", "max","min");
		close("laplacian");
		close("min");
		
		selectWindow("max");
		run("Duplicate...", "duplicate");
		run("Convert to Mask", "method=Huang background=Dark calculate black");
		run("Analyze Particles...", "size=0.00-0.80 show=Masks clear stack");
		run("Invert LUT");
		run("Subtract...", "value=127 stack");
		imageCalculator("Subtract stack", "max","Mask of max-1");
		selectWindow("Mask of max-1");
		run("Subtract...", "value=47 stack");
		imageCalculator("Subtract stack", "laplacian-1","Mask of max-1");
		close("max-1");
		close("Mask of max-1");
		
		selectWindow("laplacian-1");
		run("Duplicate...", "duplicate");
		run("Convert to Mask", "method=Percentile background=Dark calculate black");
		run("Despeckle", "stack");		
		run("Invert", "stack");
		run("Options...", "iterations=5 count=4 black do=Close stack");
		run("Options...", "iterations=5 count=3 black do=Erode stack");
		run("Options...", "iterations=3 count=3 black do=Open stack");
		run("Options...", "iterations=5 count=3 black do=Dilate stack");
		imageCalculator("Add stack", "laplacian-2","nucli");
		run("Duplicate...", "title=celled6 duplicate");
		selectWindow("laplacian-2");		
		run("Outline", "stack");
		run("Options...", "iterations=2 count=3 black do=Dilate stack");
		run("Analyze Particles...", "size=5-infinity show=[Bare Outlines] clear stack");
		run("Invert", "stack");
		run("Options...", "iterations=2 count=3 black do=Dilate stack");
		rename("laplacian");
		close("laplacian-2");
				
		selectWindow("max");
		run("Morphological Filters (3D)", "operation=Dilation element=Ball x-radius=4 y-radius=4 z-radius=2");
		run("Mean...", "radius=10 stack");
		rename("edgemax");
		imageCalculator("Average stack", "filter","laplacian-1");
		imageCalculator("Average stack", "max","ave");
		imageCalculator("Average stack", "max","filter");
		close("laplacian-1");
		close("ave");
		close("filter");	
		selectWindow("max");
		run("Morphological Filters (3D)", "operation=Dilation element=Ball x-radius=4 y-radius=4 z-radius=2");
		run("Mean...", "radius=10 stack");
		rename("edgemax2");
		imageCalculator("Average stack", "edgemax","edgemax2");
		close("edgemax2");
		selectWindow("edgemax");
		run("Convert to Mask", "method=Triangle background=Dark calculate black");
		run("Duplicate...", "title=celled7 duplicate");
		run("Invert", "stack");
		run("Options...", "iterations=10 count=4 black do=Close stack");
		close("edgemax");
		
		selectWindow("max");
		run("Convert to Mask", "method=Minimum background=Dark calculate black");
		imageCalculator("Add stack", "max","laplacian");
		run("Options...", "iterations=15 count=3 black do=Dilate stack");
		run("Duplicate...", "title=celled8 duplicate");
		run("Invert", "stack");
		close("max");
		close("laplacian");
		
		selectWindow("Vamp");
		run("Duplicate...", "title=entropy duplicate");
		run("Convert to Mask", "method=MaxEntropy background=Dark black");
		run("Despeckle", "stack");
		run("Duplicate...", "title=celled9 duplicate");
		run("Invert", "stack");
		close("entropy");	
				
		imageCalculator("Average stack", "celled1","celled4");
		imageCalculator("Average stack", "celled2","celled3");
		imageCalculator("Average stack", "celled7","celled9");
		imageCalculator("Average stack", "celled6","celled8");
		imageCalculator("Average stack", "celled1","celled2");
		imageCalculator("Average stack", "celled7","celled5");
		imageCalculator("Average stack", "celled6","celled7");
		imageCalculator("Average stack", "celled1","celled6");
		imageCalculator("Add stack", "celled1", "nucleus");
		close("celled2");
		close("celled3");
		close("celled4");
		close("celled5");
		close("celled6");
		close("celled7");
		close("celled8");
		close("celled9");
		selectWindow("celled1");
		rename("default");
		//run("Duplicate...", "title=default duplicate");
		run("Duplicate...", "title=moments duplicate");
		run("Duplicate...", "title=tick1 duplicate");
		run("Duplicate...", "title=tick0 duplicate");
		setThreshold(255, 255);
		run("Convert to Mask", "method=Default background=Dark black");
		run("Options...", "iterations=10 count=4 black do=Close stack");
		run("Options...", "iterations=5 count=3 black do=Dilate stack");
		run("Options...", "iterations=5 count=4 black do=Dilate stack");
		run("Analyze Particles...", "size=300-infinity show=Masks clear stack");
		run("Invert LUT");
		rename("selected0");
		close("tick0");
		selectWindow("tick1");
		setThreshold(200, 255);
		run("Convert to Mask", "method=Default background=Dark black");
		run("Options...", "iterations=10 count=4 black do=Close stack");
		run("Options...", "iterations=5 count=3 black do=Dilate stack");
		run("Options...", "iterations=5 count=4 black do=Dilate stack");
		run("Analyze Particles...", "size=300-1200 show=Masks clear stack");
		run("Invert LUT");
		rename("selected1");
		close("tick1");
		selectWindow("moments");
		run("Convert to Mask", "method=Moments background=Dark calculate black");
		run("Outline", "stack");
		run("Duplicate...", "duplicate");
		run("Options...", "iterations=10 count=1 black do=Dilate stack");
		run("Invert", "stack");
		run("Options...", "iterations=5 count=4 black do=Open stack");
		run("Fill Holes", "stack");
		run("Analyze Particles...", "size=300-1100 show=Masks clear stack");
		run("Invert LUT");
		rename("selected2");
		close("moments-1");
		selectWindow("moments");
		run("Options...", "iterations=5 count=4 black do=Dilate stack");
		run("Invert", "stack");
		run("Analyze Particles...", "size=300-1200 show=Masks clear stack");
		run("Invert LUT");
		run("Fill Holes", "stack");
		rename("selected3");
		close("moments");	
		selectWindow("default");
		run("Convert to Mask", "method=Default background=Dark calculate black");
		run("Outline", "stack");
		run("Options...", "iterations=10 count=1 black do=Dilate stack");
		run("Invert", "stack");
		run("Options...", "iterations=5 count=4 black do=Open stack");
		run("Fill Holes", "stack");
		run("Analyze Particles...", "size=300-1200 show=Masks clear stack");
		run("Invert LUT");
		rename("selected4");
		close("default");
		
		selectWindow("nucli");
		run("Z Project...", "projection=[Sum Slices]");
		run("Convert to Mask");
		rename("center");
		imageCalculator("AND create stack", "selected0","center");
		rename("se0");
		run("Analyze Particles...", "size=100-1000 show=Masks clear stack");
		run("Invert LUT");
		rename("seed0");
		close("se0");
		close("center");
		close("nucleus");
		selectWindow("selected0");
		run("Stack to Images");
		selectWindow("seed0");
		run("Stack to Images");
		//slices=25;
		for (i = 1; i < slices+1; i++) {
			if (i < 10) {
				run("BinaryReconstruct ", "mask=selected0-000"+i+" seed=seed0-000"+i+" white");
				close("selected0-000"+i);
			} else {
				run("BinaryReconstruct ", "mask=selected0-00"+i+" seed=seed0-00"+i+" white");
				close("selected0-00"+i);
			}
		}
		run("Images to Stack", "name=selected0 title=seed0");
		
//			if (xScaled == yScaled) {
//				run("Set Scale...", "distance="+(1/xScaled)+" known=1 unit=um global");
//			} else {
//				run("Set Scale...", "distance="+(1/xScaled)+" known=1 pixel="+(yScaled/xScaled)+" unit=um global");
//			}
				
		selectWindow("Vamp");
		run("Invert", "stack");
		imageCalculator("AND create stack", "selected0","Vamp");
		rename("fake1");
		imageCalculator("Subtract stack", "fake1", "nucli");
		close("nucli");
		selectWindow("fake1");
		setThreshold(255, 255);
		run("Convert to Mask", "method=Default background=Dark black");
		run("Despeckle", "stack");
		run("Fill Holes", "stack");
		run("Analyze Particles...", "size=20-70 circularity=0.10-1.00 show=Masks clear stack");
		run("Invert LUT");
		rename("fake2");
		close("fake1");
		run("Options...", "iterations=10 count=4 black do=Dilate stack");
		run("Analyze Particles...", "size=0-infinity circularity=0.45-1.00 show=Masks clear stack");
		run("Invert LUT");
		rename("fake");
		close("fake2");
				
		imageCalculator("AND create stack", "selected0","Vamp");
		rename("nucli");
		setThreshold(255, 255);
		run("Convert to Mask", "method=Default background=Dark black");
		run("Despeckle", "stack");
		run("Fill Holes", "stack");
		run("Analyze Particles...", "size=45-infinity circularity=0.10-1.00 show=Masks display clear stack");
		s = Table.getColumn('Slice');
		x = Table.getColumn('X');
		y = Table.getColumn('Y');
		toUnscaled(x, y);
		c = newArray(nResults);
		run("Invert LUT");
		rename("nucleus");
		close("nucli");
		close("Results");
		run("Duplicate...", "title=blank");
		run("Select All");
		run("Clear", "slice");
		run("Select None");
		c[0] = 1;
		setForegroundColor(1,1,1);
		makeRectangle(x[0]-75, y[0]-75, 150, 150);
		run("Fill", "slice");
		numCells = 1;
		updateDisplay();
		for (i = 1; i < c.length; i++) {
			if (getPixel(x[i], y[i]) != 0) {
				numCells = getPixel(x[i],y[i]);
			} else {
				Array.getStatistics(c, min, max);
				numCells = max+1;
			}
		    setForegroundColor(numCells,numCells,numCells);
			makeRectangle(x[i]-75, y[i]-75, 150, 150);
			run("Fill", "slice");
			c[i] = numCells;
			updateDisplay();
		}
		close("blank");
		Array.sort(c,s,x,y);
		skip = Array.copy(c);
		Array.fill(skip, 0);
		nucStart = 0;
		for (i = 1; i < c.length; i++) {
			if (c[i] != c[i-1]) {
				nucStop = i;
				if (nucStop - nucStart < 5) {
					for (j=nucStart; j < nucStop; j++) {
						skip[j] = 1; 
					}
				}
				nucStart = nucStop;
			}
		}
		nucStop = c.length;
		if (nucStop - nucStart < 5) {
			for (j=nucStart; j < nucStop; j++) {
				skip[j] = 1; 
			}
		}
		setForegroundColor(255,255,255);
		selectWindow("nucleus");
		run("Options...", "iterations=2 count=4 black do=Erode stack");
		for (i = 0; i < c.length; i++) {
			selectWindow("nucleus");
			setSlice(s[i]);
			doWand(x[i], y[i]);
			if (skip[i] == 0) {
				run("Fit Ellipse");
				run("Fill", "slice");
				if (getValue("Area") < 60) {
					run("Make Band...", "band=1");
					run("Fill", "slice");
				}
				run("Select None");
			} else {
				run("Select Bounding Box");
				run("Clear", "slice");
				run("Select None");
			}
		}
		for (i = 0; i < c.length; i++) {
			if (i > 0 && skip[i] == 0) {				
				if (c[i] == c[i-1] && s[i] != s[i-1]+1) {
					selectWindow("nucleus");
					setSlice(s[i-1]);
					doWand(x[i-1], y[i-1]);
					run("Fit Ellipse");
					for (j=s[i-1]+1; j<s[i]; j++) {
						setSlice(j);
						run("Restore Selection");
						run("Fill", "slice");
						run("Select None");
					}
				}
			}
		}
		
		
		imageCalculator("Subtract stack", "fake","nucleus");
		selectWindow("fake");
		run("Options...", "iterations=3 count=1 black do=Dilate stack");
		imageCalculator("Subtract stack", "selected0", "fake");
		imageCalculator("Subtract stack", "selected1", "fake");
		imageCalculator("Subtract stack", "selected2", "fake");
		imageCalculator("Subtract stack", "selected3", "fake");
		imageCalculator("Subtract stack", "selected4", "fake");
		close("fake");
		
		imageCalculator("AND create stack", "selected0","nucleus");
		rename("se0");
		run("Analyze Particles...", "size=60-1000 show=Masks clear stack");
		run("Invert LUT");
		rename("seed0");
		close("se0");
		imageCalculator("AND create stack", "selected1","nucleus");
		rename("se1");
		run("Analyze Particles...", "size=60-1000 show=Masks clear stack");
		run("Invert LUT");
		rename("seed1");
		close("se1");
		imageCalculator("AND create stack", "selected2","nucleus");
		rename("se2");
		run("Analyze Particles...", "size=60-1000 show=Masks clear stack");
		run("Invert LUT");
		rename("seed2");
		close("se2");
		imageCalculator("AND create stack", "selected3","nucleus");
		rename("se3");
		run("Analyze Particles...", "size=60-1000 show=Masks clear stack");
		run("Invert LUT");
		rename("seed3");
		close("se3");
		imageCalculator("AND create stack", "selected4","nucleus");
		rename("se4");
		run("Analyze Particles...", "size=60-1000 show=Masks clear stack");
		run("Invert LUT");
		rename("seed4");
		close("se4");
		close("nucleus");
		
		for (j=0; j<5; j++) {
			selectWindow("selected"+j);
			run("Stack to Images");
			selectWindow("seed"+j);
			run("Stack to Images");
			//slices=25;
			for (i = 1; i < slices+1; i++) {
				if (i < 10) {
					run("BinaryReconstruct ", "mask=selected"+j+"-000"+i+" seed=seed"+j+"-000"+i+" white");
					close("selected"+j+"-000"+i);
				} else {
					run("BinaryReconstruct ", "mask=selected"+j+"-00"+i+" seed=seed"+j+"-00"+i+" white");
					close("selected"+j+"-00"+i);
				}
			}
			run("Images to Stack", "name=selected"+j+" title=seed"+j);
		}

		
		imageCalculator("Average create stack", "selected2","selected4");
		rename("merge");
		imageCalculator("Average stack", "selected1","selected3");
		imageCalculator("Add stack", "selected1","selected0");
		imageCalculator("Average stack", "merge","selected1");
		close("selected0");
		close("selected1");
		close("selected2");
		close("selected3");
		close("selected4");
		selectWindow("merge");
		setThreshold(100, 255);
		run("Convert to Mask", "method=Default background=Dark black");
		run("Options...", "iterations=5 count=4 black do=Close stack");
		run("Analyze Particles...", "size=100-infinity show=Nothing display clear stack");
		s = Table.getColumn('Slice');
		close("Results");
		//slices=25;
		for (j = 1; j < s[0]; j++) {
			selectWindow("merge");
			setSlice(j);
			run("Duplicate...", "title=Ave"+j);
		}
		for (j = s[0]; j <= slices; j++) {
			selectWindow("merge");
			if (j < 2) {
		    	run("Z Project...", "start="+j-1+" stop="+j+9+" projection=[Average Intensity]");
			} else if (j < 5) {
		    	run("Z Project...", "start="+j-2+" stop="+j+8+" projection=[Average Intensity]");
			} else if (j < s.length-5) {
				run("Z Project...", "start="+j-4+" stop="+j+4+" projection=[Average Intensity]");
			} else if (j < s.length-3) {
				run("Z Project...", "start="+j-8+" stop="+j+2+" projection=[Average Intensity]");
			} else {
				run("Z Project...", "start="+j-9+" stop="+j+1+" projection=[Average Intensity]");
			}
			rename("Ave"+j);
		}
		close("merge");
		run("Images to Stack", "name=merge title=Ave");
		setThreshold(100,255);
		run("Convert to Mask", "method=Default background=Dark black");
		run("Analyze Particles...", "size=500-1500 show=Nothing display clear stack");
		s = Table.getColumn('Slice'); 
		x = Table.getColumn('X');
		y = Table.getColumn('Y');
		toUnscaled(x, y);
		close("Results");
		for (j = 0; j < s.length; j++) {
			selectWindow("merge");
			setSlice(s[j]);
			doWand(x[j],y[j]);
			k = j;
			for (i = j+1; i < s.length; i++) {
				if (s[i] == s[j]) {
					setKeyDown("shift");
					doWand(x[i], y[i]);
					k = i;
				}
			}
			run("Enlarge...", "enlarge=-4");
			wait(50);
			run("Enlarge...", "enlarge=4.5");
			wait(50);
			run("Clear Outside", "slice");
			run("Fill Holes", "slice");
			run("Select None");
			doWand(x[j],y[j]);
			for (i = j+1; i < s.length; i++) {
				if (s[i] == s[j]) {
					setKeyDown("shift");
					doWand(x[i], y[i]);
				}
			}
			run("Enlarge...", "enlarge=-2");
			wait(50);
			run("Enlarge...", "enlarge=2");
			wait(50);
			run("Clear Outside", "slice");
			run("Fill Holes", "slice");
			run("Select None");
			j = k;
		}
		run("Analyze Particles...", "size=100-infinity show=Nothing display clear stack");
		s = Table.getColumn('Slice');
		close("Results");
		//slices=25;
		for (j = 1; j < s[0]; j++) {
			selectWindow("merge");
			setSlice(j);
			run("Duplicate...", "title=Ave"+j);
		}
		for (j = s[0]; j <= slices; j++) {
			selectWindow("merge");
			if (j < 2) {
		    	run("Z Project...", "start="+j-1+" stop="+j+9+" projection=[Average Intensity]");
			} else if (j < 5) {
		    	run("Z Project...", "start="+j-2+" stop="+j+8+" projection=[Average Intensity]");
			} else if (j < s.length-5) {
				run("Z Project...", "start="+j-4+" stop="+j+4+" projection=[Average Intensity]");
			} else if (j < s.length-3) {
				run("Z Project...", "start="+j-8+" stop="+j+2+" projection=[Average Intensity]");
			} else {
				run("Z Project...", "start="+j-9+" stop="+j+1+" projection=[Average Intensity]");
			}
			rename("Ave"+j);
		}
		close("merge");
		run("Images to Stack", "name=merge title=Ave");
		setThreshold(100,255);
		run("Convert to Mask", "method=Default background=Dark black");
		run("Duplicate...", "title=celledges duplicate");
		run("Analyze Particles...", "size=500-1500 show=Nothing display clear stack");
		s = Table.getColumn('Slice'); 
		x = Table.getColumn('X');
		y = Table.getColumn('Y');
		toUnscaled(x, y);
		close("Results");
		for (j = 0; j < s.length; j++) {
			selectWindow("celledges");
			setSlice(s[j]);
			doWand(x[j],y[j]);
			k = j;
			for (i = j+1; i < s.length; i++) {
				if (s[i] == s[j]) {
					setKeyDown("shift");
					doWand(x[i], y[i]);
					setKeyDown("none");
					k = i;
				}
			}
			run("Enlarge...", "enlarge=3");
			wait(50);
			setForegroundColor(255, 255, 255);
			run("Fill", "slice");
			setBackgroundColor(0, 0, 0);
			run("Clear Outside", "slice");
			run("Select None");
			j = k;
		}
		for (j = 0; j < s.length; j++) {
			selectWindow("merge");
			setSlice(s[j]);
			doWand(x[j],y[j]);
			k = j;
			for (i = j+1; i < s.length; i++) {
				if (s[i] == s[j]) {
					setKeyDown("shift");
					doWand(x[i], y[i]);
					setKeyDown("none");
					k = i;
				}
			}
			run("Enlarge...", "enlarge=-1");
			wait(50);
			run("Clear Outside", "slice");
			run("Select None");
			j = k;
		}
		rename("cells");
		run("Fill Holes", "stack");
		imageCalculator("Subtract stack", "celledges", "cells");
	}
	close("FB");
	close("Vesicular");
	close("Vamp");
	selectWindow("cells");
	//saveAs("tiff", savePath+"-cells.tif");
	rename("cells");
	selectWindow("celledges");
	//saveAs("tiff", savePath+"-celledges.tif");
	rename("celledges");

	//open the original file again
	run("Bio-Formats Importer", "open=["+fpath+"] color_mode=Grayscale rois_import=[ROI manager] split_channels view=Hyperstack stack_order=XYCZT");
	//chanArray = newArray("Fast Blue","VGAT/VIAAT","VAMP1");
	chan = 0;
	for (i = 1; i < nImages+1; i++) {	
		selectImage(i);
		img = getTitle();			
		if (!img.contains("cell")) {
			if (chanArray[chan] == "Fast Blue") {
				rename("FB");
			} else if (chanArray[chan] == "NeuN" || chanArray[chan] == "ChAT") {
				rename("NC");
			} else if (chanArray[chan] == "VGLUT1" || chanArray[chan] == "VGLUT2" || chanArray[chan] == "VGAT/VIAAT") {
				rename("Vesicular");
				vesChannel = chan;
			} else if (chanArray[chan] == "Other") {
				rename("Other");
			} else {
				rename("Vamp");
				vampChannel = chan;
			}
			chan++;
			//slices = nSlices;	
			for (j=1; j<=slices; j++) {
		   		setSlice(j);
		   		getRawStatistics(nPixels, mean, min, max, std, histogram);
		   		percent = Array.copy(histogram);
		    	percent[0] = histogram[0]/nPixels;
		    	bottom = 0;
		    	multiply = 1;
		    	for (k=1; k<histogram.length; k++) {
		    		histogram[k] = histogram[k-1]+histogram[k];
		    		percent[k] = histogram[k]*100/nPixels;
		    		if (percent[k]<=1.1) {
		    			bottom = k;
		    		}
		    		if (percent[k]>95 && percent[k]<=99) {
		    			multiply = 4095/(k-bottom);
		    		}
		    	}
		    	run("Subtract...", "value="+bottom+" slice");
		    	run("Multiply...", "value="+multiply+" slice");
			}
			run("Max...", "value=4095 stack");		
		}
	}

	selectWindow("celledges");
	run("Duplicate...", "title=box duplicate");
	run("Set Measurements...", "centroid bounding stack redirect=None decimal=3");
	run("Analyze Particles...", "size=0.00-infinity show=Nothing display clear add stack");
	s = Table.getColumn('Slice');
	x = Table.getColumn('X');
	y = Table.getColumn('Y');
	bx = Table.getColumn('BX');
	by = Table.getColumn('BY');
	bw = Table.getColumn('Width');
	bh = Table.getColumn('Height');
	toUnscaled(x, y);
	c = newArray(nResults);
	close("Results");
	numCells = 1;
	c[0] = numCells;
	roiManager("Show None");
	n = roiManager("count");
	setForegroundColor(255, 255, 255);
	for (i = 0; i < n; i++) {
	    roiManager("select", i);
		run("Enlarge...", "enlarge=1");
		wait(100);
		run("Fill", "slice");
		run("To Bounding Box");
		bx[i] = getValue("BX");
		by[i] = getValue("BY");
		bw[i] = getValue("Width");
		bh[i] = getValue("Height");
		run("Select None");
		roiManager("select", i);
		run("Clear", "slice");
		run("Select None");
	}
	close("ROI Manager");
	toUnscaled(bx, by);
	toUnscaled(bw, bh);
	run("Duplicate...", "title=blank");
	run("Select All");
	run("Clear", "slice");
	run("Select None");
	setForegroundColor(numCells, numCells, numCells);
	makeRectangle(x[0]-75, y[0]-75, 150, 150);
	run("Fill", "slice");
	updateDisplay();
	for (i = 1; i < nResults; i++) {
		if (getPixel(x[i], y[i]) != 0) {
			numCells = getPixel(x[i],y[i]);
		} else {
			Array.getStatistics(c, min, max);
			numCells = max+1;
		}
	    setForegroundColor(numCells, numCells, numCells);
		makeRectangle(x[i]-75, y[i]-75, 150, 150);
		run("Fill", "slice");
		c[i] = numCells;
		updateDisplay();
	}
	close("blank");
	close("box");
	Array.sort(c,s,x,y,bx,by,bw,bh);
	cropParam = newArray(4*numCells);
	bxmin = bx[0];
	bymin = by[0];
	bwmax = bw[0];
	bhmax = bh[0];
	cPstart = 0;
	for (i = 1; i < c.length; i++) {
		if (c[i] == c[i-1]) {
			bxmin = minOf(bxmin, bx[i]);
			bymin = minOf(bymin, by[i]);
			bwmax = maxOf(bwmax, bw[i]);
			bhmax = maxOf(bhmax, bh[i]);
		} else {
			cropParam[cPstart] = bxmin;
			cropParam[cPstart+1] = bymin;
			cropParam[cPstart+2] = bwmax;
			cropParam[cPstart+3] = bhmax;
			cPstart = cPstart+4;
			bxmin = bx[i];
			bymin = by[i];
			bwmax = bw[i];
			bhmax = bh[i];
		}
	}
	cropParam[cPstart] = bxmin;
	cropParam[cPstart+1] = bymin;
	cropParam[cPstart+2] = bwmax;
	cropParam[cPstart+3] = bhmax;
	
	for (i = 2; i <= (cropParam.length/4); i++) {
		selectWindow("cells");
		run("Duplicate...", "title=cell"+i+" duplicate");
		makeRectangle(cropParam[4*(i-1)], cropParam[4*(i-1)+1], cropParam[4*(i-1)+2], cropParam[4*(i-1)+3]);
		run("Crop");
		run("Select None");
		selectWindow("celledges");
		run("Duplicate...", "title=celledge"+i+" duplicate");
		makeRectangle(cropParam[4*(i-1)], cropParam[4*(i-1)+1], cropParam[4*(i-1)+2], cropParam[4*(i-1)+3]);
		run("Crop");			
		run("Select None");
		selectWindow("Vesicular");
		run("Duplicate...", "title=Vesicular-cell"+i+" duplicate");
		makeRectangle(cropParam[4*(i-1)], cropParam[4*(i-1)+1], cropParam[4*(i-1)+2], cropParam[4*(i-1)+3]);
		run("Crop");
		run("Select None");
		selectWindow("Vamp");
		run("Duplicate...", "title=Vamp-cell"+i+" duplicate");
		makeRectangle(cropParam[4*(i-1)], cropParam[4*(i-1)+1], cropParam[4*(i-1)+2], cropParam[4*(i-1)+3]);
		run("Crop");			
		run("Select None");			
	}
	selectWindow("cells");
	rename("cell1");
	makeRectangle(cropParam[0], cropParam[1], cropParam[2], cropParam[3]);
	run("Crop");
	run("Select None");
	selectWindow("celledges");
	rename("celledge1");
	makeRectangle(cropParam[0], cropParam[1], cropParam[2], cropParam[3]);
	run("Crop");
	run("Select None");	
	selectWindow("Vesicular");
	rename("Vesicular-cell1");
	makeRectangle(cropParam[0], cropParam[1], cropParam[2], cropParam[3]);
	run("Crop");
	run("Select None");
	selectWindow("Vamp");
	rename("Vamp-cell1");
	makeRectangle(cropParam[0], cropParam[1], cropParam[2], cropParam[3]);
	run("Crop");			
	run("Select None");			
	
	run("Set Measurements...", "area mean standard centroid bounding shape feret's stack redirect=None decimal=3");
	numCells=i-1;
	backgroundTable = "Background";
	Table.create(backgroundTable);
	for (j = 1; j <= numCells; j++) {
		selectWindow("celledge"+j);
		run("Analyze Particles...", "size=100-infinity show=Nothing display clear add stack");
		s = Table.getColumn('Slice',"Results");
		roiManager("Show None");
		n = roiManager("count");
		
		selectWindow("Vesicular-cell"+j);
		run("Duplicate...", "title=Vesicular-celledge"+j+" duplicate");
		setBackgroundColor(0, 0, 0);
		for (k = 1; k < s[0]; k++) {
			setSlice(k);
			run("Select All");
			run("Clear", "slice");
			run("Select None");
		}
		for (k = 0; k < n; k++) {
	    	roiManager("select", k);
	    	run("Clear Outside", "slice");
	    	run("Select None");
		}
		for (k = s[s.length-1]+1; k <= slices; k++) {
			setSlice(k);
			run("Select All");
			run("Clear", "slice");
			run("Select None");
		}
		selectWindow("Vamp-cell"+j);
		run("Duplicate...", "title=Vamp-celledge"+j+" duplicate");
		setBackgroundColor(0, 0, 0);
		for (k = 1; k < s[0]; k++) {
			setSlice(k);
			run("Select All");
			run("Clear", "slice");
			run("Select None");
		}
		for (k = 0; k < n; k++) {
	    	roiManager("select", k);
	    	run("Clear Outside", "slice");
	    	run("Select None");
		}
		for (k = s[s.length-1]+1; k <= slices; k++) {
			setSlice(k);
			run("Select All");
			run("Clear", "slice");
			run("Select None");
		}
		close("ROI Manager");
		close("Results");		
		
		selectWindow("cell"+j);
		run("Analyze Particles...", "size=100-infinity show=Nothing display clear add stack");
		s = Table.getColumn('Slice',"Results");
		roiManager("Show None");
		n = roiManager("count");
		
		selectWindow("Vesicular-cell"+j);
		setBackgroundColor(0, 0, 0);
		for (k = 1; k < s[0]; k++) {
			setSlice(k);
			run("Select All");
			run("Clear", "slice");
			run("Select None");
		}
		Table.reset("Results");
		for (k = 0; k < n; k++) {
	    	roiManager("select", k);
	    	run("Measure");
	    	run("Clear Outside", "slice");
	    	run("Select None");
		}
		bmean = Table.getColumn("Mean","Results");
		bstd = Table.getColumn("StdDev", "Results");
		for (k = 0; k < n; k++) {
		    Table.set("Vesicular-Cell-"+j+" Mean", k, bmean[k], backgroundTable);
			Table.set("Vesicular-Cell-"+j+" StdDev", k, bstd[k], backgroundTable);
		}
		for (k = s[s.length-1]+1; k <= slices; k++) {
			setSlice(k);
			run("Select All");
			run("Clear", "slice");
			run("Select None");
		}
			
		selectWindow("Vesicular-celledge"+j);
		setBackgroundColor(0, 0, 0);
		for (k = 0; k < n; k++) {
	    	roiManager("select", k);
	    	run("Clear", "slice");
	    	run("Select None");
	    	run("Subtract...", "value="+(bmean[k]+2*bstd[k])+" slice");
		}
					
		selectWindow("Vamp-cell"+j);
		setBackgroundColor(0, 0, 0);
		for (k = 1; k < s[0]; k++) {
			setSlice(k);
			run("Select All");
			run("Clear", "slice");
			run("Select None");
		}
		Table.reset("Results");
		for (k = 0; k < n; k++) {
	    	roiManager("select", k);
	    	run("Measure");
	    	run("Clear Outside", "slice");
	    	run("Select None");
		}
		bmean = Table.getColumn("Mean","Results");
		bstd = Table.getColumn("StdDev", "Results");
		for (k = 0; k < n; k++) {
		    Table.set("Vesicular-Cell-"+j+" Mean", k, bmean[k], backgroundTable);
			Table.set("Vesicular-Cell-"+j+" StdDev", k, bstd[k], backgroundTable);
		}
		for (k = s[s.length-1]+1; k <= slices; k++) {
			setSlice(k);
			run("Select All");
			run("Clear", "slice");
			run("Select None");
		}

		selectWindow("Vamp-celledge"+j);
		setBackgroundColor(0, 0, 0);
		for (k = 0; k < n; k++) {
	    	roiManager("select", k);
	    	run("Clear", "slice");
	    	run("Select None");
	    	run("Subtract...", "value="+(bmean[k]+2*bstd[k])+" slice");
		}
		close("ROI Manager");
		close("Results");
	}
	//saveAs(backgroundTable, savePath+"-Background.csv");
	close(backgroundTable);

	run("3D OC Options", "volume surface centroid mean_distance_to_surface std_dev_distance_to_surface bounding_box dots_size=20 font_size=50 redirect_to=none");
	setForegroundColor(255,255,255);
	//chanArray = newArray("Fast Blue","VGAT/VIAAT","VAMP1");
	for (i=1; i<=2; i++) {			
		sizeMin = 0.1;
		sizeMax = 4;
		if (i>1) {
			begin = "Vamp";
		} else {
			begin = "Vesicular";
			check = Array.filter(chanArray, "VGLUT1");
			if (check.length > 0) {
				sizeMin = 0.5;
				sizeMax = 6;
			}
		}
		for (j=1; j<=numCells; j++) {
			selectWindow(begin+"-celledge"+j);
			setThreshold(1,4095);
			run("Analyze Particles...", "size=0-infinity show=Nothing clear stack");
			s = Table.getColumn('Slice');
			//slices=nSlices;
			Array.getStatistics(s, min, max);
			resetThreshold();
			for (k=1; k<min; k++) {
				selectWindow(begin+"-celledge"+j);
				setSlice(k);
				run("Duplicate...", "title=selected"+k);
				run("8-bit");
			}
			for (k=min; k<=max; k++) {
				selectWindow(begin+"-celledge"+j);
				setSlice(k);
				pixelMax = 4096;
				pixelMin = maxOf(getValue("Min"),1);
				for (m=0; m<15; m++) {
					selectWindow(begin+"-celledge"+j);
					run("Select None");
					skip = false;
					if (m==0) {
						lowValue = getValue("Max")-100;
						setThreshold(lowValue, pixelMax);
						run("Analyze Particles...", "size="+sizeMin+"-"+sizeMax+" show=Masks clear slice");
						rename("possible");
						run("Invert LUT");
						run("Options...", "iterations=1 count=4 black do=Open");
						run("Options...", "iterations=1 count=4 black do=Dilate");
						run("Analyze Particles...", "size="+sizeMin+"-"+sizeMax+" show=Masks clear slice");
						rename(begin+"-cell"+j+"-selection"+k);
						run("Invert LUT");
						close("possible");
						skip = true;
					} else if (m<3) {
						lowValue = (getValue("Max")-100)*(4-m)/4;
					} else {
						lowValue = (getValue("Max")-100)/2 - m*((getValue("Max")-100)/2-pixelMin)/14;
					}
					if (!skip) {
						setThreshold(lowValue, pixelMax);
						run("Analyze Particles...", "size="+sizeMin+"-"+sizeMax+" show=Masks clear slice");
						rename("possible");
						run("Invert LUT");
						run("Options...", "iterations=1 count=4 black do=Open");
						run("Options...", "iterations=1 count=4 black do=Dilate");
						run("Analyze Particles...", "size="+sizeMin+"-"+sizeMax+" show=Nothing clear slice");
						if (nResults>=1) {
							f = Table.getColumn("Feret");
							r = Table.getColumn("Round");
							x = Table.getColumn('X');
							y = Table.getColumn('Y');
							toUnscaled(x, y);
							for (n=0; n<nResults; n++) {
								if (f[n] >= 1 && f[n] <= 3 && r[n] > 0.45) {
									selectWindow("possible");
									doWand(x[n],y[n]);
									wait(10);
									selectWindow(begin+"-cell"+j+"-selection"+k);
									run("Restore Selection");
									wait(10);
									run("Fill", "slice");
									run("Select None");
								}
							}
						}
						close("Results");
						close("possible");
					}
					selectWindow(begin+"-celledge"+j);
					run("Select None");
					resetThreshold();
				}
				selectWindow(begin+"-cell"+j+"-selection"+k);
				rename("selected"+k);
			}
			for (k=max+1; k<=slices; k++) {
				selectWindow(begin+"-celledge"+j);
				setSlice(k);
				run("Duplicate...", "title=selected"+k);
				run("8-bit");
			}
			run("Images to Stack", "name="+begin+"-cell"+j+"-analysis title=selected");
			close("Results");
			if (i==1) {
				run("3D Objects Counter", "threshold=128 slice=1 min.=10 max.=4129270 statistics");
				bdepth = Table.getColumn("B-depth");
				bz = Table.getColumn("BZ");
				x = Table.getColumn('X');
				y = Table.getColumn('Y');
				for (n=0; n<bdepth.length; n++) {
					if (bdepth[n] == 1) {
						selectWindow(begin+"-cell"+j+"-analysis");
						setSlice(bz[n]);
						doWand(x[n], y[n]);
						wait(10);
						run("Clear", "slice");
						wait(10);
						run("Select None");
					}
				}
				close("Results");
			}
			
			selectWindow(begin+"-cell"+j);
			if (i>1) {
				//saveAs("tiff", savePath+"-Cell-"+j+"-"+chanArray[vampChannel]+".tif");
			} else {
				//saveAs("tiff", savePath+"-Cell-"+j+"-"+chanArray[vesChannel]+".tif");
			}
			rename(begin+"-cell"+j);
			close(begin+"-cell"+j);
			selectWindow(begin+"-celledge"+j);
			if (i>1) {
				//saveAs("tiff", savePath+"-Celledge-"+j+"-"+chanArray[vampChannel]+".tif");
			} else {
				//saveAs("tiff", savePath+"-Celledge-"+j+"-"+chanArray[vesChannel]+".tif");
			}
			rename(begin+"-cell"+j);
			close(begin+"-cell"+j);
		}
	}
	
	// Create ROI analysis overlap between the Vesicular and Vamp channels
	for (j=1; j<=numCells; j++) {
		imageCalculator("AND create stack", "Vesicular-cell"+j+"-analysis", "Vamp-cell"+j+"-analysis");
		rename("cell"+j+"-overlap");
	}
	
	// Analyze the overlap and store the tables
	run("3D OC Options", "volume surface centroid mean_distance_to_surface std_dev_distance_to_surface bounding_box dots_size=20 font_size=50 store_results_within_a_table_named_after_the_image_(macro_friendly) redirect_to=none");
	run("Set Measurements...", "area mean standard modal min area_fraction stack redirect=None decimal=3");
	Table.create("Results");
	for (j=1; j<=numCells; j++) {
		selectWindow("Vesicular-cell"+j+"-analysis");
		vesTitle = getTitle();
		vesTablename = "Statistics for "+vesTitle;
		vesRoiImage = "Objects map of "+vesTitle;
		run("3D Objects Counter", "threshold=128 slice=1 min.=10 max.=4129270 objects statistics");
		tableCount = Table.size(vesTablename);
		for (i=1; i<=tableCount; i++) {
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
				Table.set("Slice "+k, i-1, 0, vesTablename);
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
					Table.set("Slice "+k, i-1, val/denom, vesTablename);
				} else {
					Table.set("Slice "+k, i-1, 0, vesTablename);
				}
			}
			for (k=s[s.length-1]+1; k<=slices; k++) {
				Table.set("Slice "+k, i-1, 0, vesTablename);
			}
			roiManager("reset");
		}
		close("ROI Manager");
		selectWindow(vesRoiImage);
		resetThreshold;
		//close("cell"+j+"-overlap");
		//close(vesTitle);
		// Saving the analysis table for the vesicular channel
		Table.update(vesTablename);
		selectWindow(vesTablename);
		//saveAs("results", savePath+"-Cell-"+j+"-"+chanArray[vesChannel]+".csv");
		//close(Table.title);
		// Saving the objects image for the vesicular channel
		selectWindow(vesRoiImage);
		//saveAs("tiff", savePath+"-Cell-"+j+"-"+chanArray[vesChannel]+"-ROIs.tif");
		//rename(vesRoiImage);

		// Figuring out which vesicular ROIs the vamp ROIs line up with
		selectWindow("Vamp-cell"+j+"-analysis");
		vampTitle = getTitle();
		vampTablename = "Statistics for "+vampTitle;
		vampRoiImage = "Objects map of "+vampTitle;
		run("3D Objects Counter", "threshold=128 slice=1 min.=10 max.=4129270 objects statistics");
		tableCount = Table.size(vampTablename);
		for (i=1; i<=tableCount; i++) {
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
				Table.set("Slice "+k, i-1, 0, vampTablename);
			}
			for (k=s[0]; k<=s[s.length-1]; k++) {
				val = 0;
				denom = 0;
				for (m=0; m<s.length; m++) {
					if (s[m] == k) {
						if (a[m] > 0 && a[m]%1 == 0) {
							val = val + a[m];
							denom++;
						} else {
							val = val + b[m];
							denom++;
						}
					}
				}
				if (denom > 0) {
					Table.set("Slice "+k, i-1, val/denom, vampTablename);
				} else {
					Table.set("Slice "+k, i-1, 0, vampTablename);
				}
			}
			for (k=s[s.length-1]+1; k<=slices; k++) {
				Table.set("Slice "+k, i-1, 0, vampTablename);
			}
			roiManager("reset");
		}
		close("ROI Manager");
		selectWindow(vampRoiImage);
		resetThreshold;
		//close(vampTitle);
		// Saving the analysis table for the vamp channel
		Table.update(vampTablename);
		selectWindow(vampTablename);
		//saveAs("results", savePath+"-Cell-"+j+"-"+chanArray[vampChannel]+".csv");
		//close(Table.title);
		// Saving the objects image for the vamp channel
		selectWindow(vampRoiImage);
		//saveAs("tiff", savePath+"-Cell-"+j+"-"+chanArray[vampChannel]+"-ROIs.tif");
		//rename(vampRoiImage);
		close(vampRoiImage);
		close(vesRoiImage);
	}
	close("Results");
	close("*");
	close
}

print("Finished processing all images in "+imageDir+".");
print("Took " + ((getTime()-programStart)/1000)/60 + "minutes");
//setBatchMode(false);



