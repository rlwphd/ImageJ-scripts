// Getting the folder where images are located
imageDir=getDirectory("Select Spinal Cord Folder for Image Processing");
roiDir=imageDir+"Analysis_Results\\";
imagelist = getFileList(imageDir);
imageArray = newArray;
File.makeDirectory(roiDir);
saveSettings();
setOption("BlackBackground", true);
setBackgroundColor(0, 0, 0);
run("Set Measurements...", "area mean standard fit median stack redirect=None decimal=3");
//setBatchMode(true);

start = getTime();
// Grabbing each image individually for analysis
//for (image=0;image<imagelist.length;image++) {
for (image=0;image<1;image++) {


// Opening each image that is a 60x image and getting the image title
// This can be modified to whatever file ending is needed (could just be .oib)
	if (endsWith(imagelist[image], ".oib")) {
		imgName=imagelist[image];
		baseNameEnd=indexOf(imgName, ".oib"); 
        baseName=substring(imgName, 0, baseNameEnd);
        maskName = baseName + "-mask-";
        savePath = roiDir+baseName;
		fpath = imageDir+imgName;
		mpath = roiDir+maskName;
		run("Bio-Formats Importer", "open=["+fpath+"] color_mode=Grayscale rois_import=[ROI manager] split_channels view=Hyperstack stack_order=XYCZT");

		for (i = 0; i < nImages; i++) {	
			selectImage(i+1);
			rename("ch"+toString(i+1));
			slices = nSlices;
			for (j=1; j<=nSlices; j++) {
		   		setSlice(j);
		   		getRawStatistics(nPixels, mean, min, max, std, histogram);
		   		percent = Array.copy(histogram);
		    	percent[0] = histogram[0]/nPixels;
		    	bottom = 0;
		    	for (k=1; k<histogram.length; k++) {
		    		histogram[k] = histogram[k-1]+histogram[k];
		    		percent[k] = histogram[k]*100/nPixels;
		    		if (percent[k]<=10.1) {
		    			bottom = k;
		    		}
		    		if (percent[k]>94 && percent[k]<=97 && k>300 && i==0) {
		    			multiply = 4095/(k-bottom);
		    		} else if (percent[k]>94 && percent[k]<=97 && k<300 && i==0) {
		    			multiply = 0.5;
		    		} else if (percent[k]>90 && percent[k]<=93 && i!=0) {
		    			multiply = 4095/(k-bottom);
		    		}
		    	}
		    	run("Subtract...", "value="+bottom+" slice");
		    	run("Multiply...", "value="+multiply+" slice");
			}			
			run("8-bit");
		}
		
		selectWindow("ch1");
		
		
		
		
		
		selectWindow("ch3");
		run("Duplicate...", "duplicate");
		run("Convert to Mask", "method=Default background=Dark calculate black");
		run("Analyze Particles...", "size=0-0.6 circularity=0.0-1.00 show=Masks clear stack");
		run("Invert LUT");
		run("Options...", "iterations=1 count=2 black do=Dilate stack");
		imageCalculator("Subtract stack", "ch3","Mask of ch3-1");
		close("ch3-1");
		close("Mask of ch3-1");
		selectWindow("ch2");
		run("Duplicate...", "duplicate");
		run("Convert to Mask", "method=Default background=Dark calculate black");
		run("Analyze Particles...", "size=0.00-0.60 show=Masks clear stack");
		run("Invert LUT");
		run("Options...", "iterations=1 count=2 black do=Dilate stack");
		imageCalculator("Subtract stack", "ch2","Mask of ch2-1");
		close("ch2-1");
		close("Mask of ch2-1");
		
		imageCalculator("Max create stack", "ch2","ch3");
		rename("max");
		run("Statistical Region Merging", "q=20 showaverages 3d");
		rename("filternucleus");
		run("8-bit");
		run("Subtract...", "value=18 stack");
		setThreshold(0, 0);
		run("Convert to Mask", "method=MinError background=Dark black");
		run("Despeckle", "stack");
		run("Fill Holes", "stack");
		run("Outline", "stack");
		run("Options...", "iterations=5 count=2 black do=Dilate stack");
		run("Fill Holes", "stack");
		run("Options...", "iterations=5 count=1 black do=Erode stack");
		run("Analyze Particles...", "size=4500-30000 circularity=0.10-1.00 show=Masks clear stack");
		run("Invert LUT");
		run("Fill Holes", "stack");
		rename("nucli");
		close("filternucleus");
		selectWindow("nucli");
		run("Mean...", "radius=25 stack");
		run("Convert to Mask", "method=Minimum background=Dark calculate black");
		run("Options...", "iterations=30 count=3 black do=Dilate stack");
		run("Analyze Particles...", "size=10000-40000 circularity=0.65-1.00 show=Ellipses exclude clear stack");
		run("Invert", "stack");
		run("Fill Holes", "stack");
		rename("nucleus");
		close("nucli");
		imageCalculator("Subtract stack", "ch2","nucleus");
		imageCalculator("Subtract stack", "ch3","nucleus");
		
		selectWindow("ch3");
		run("Duplicate...", "title=raw duplicate");
		run("Gaussian Blur...", "sigma=15 stack");
		run("Find Edges", "stack");
		selectWindow("ch3");
		run("Duplicate...", "title=1stpass duplicate");
		setAutoThreshold("RenyiEntropy dark");
		setOption("BlackBackground", true);
		run("Convert to Mask", "method=RenyiEntropy background=Dark calculate black");
		run("Options...", "iterations=3 count=4 black do=Erode stack");
		run("Analyze Particles...", "size=0.00-1.00 show=Masks clear stack");
		run("Invert LUT");
		imageCalculator("Subtract stack", "1stpass","Mask of 1stpass");
		close("Mask of 1stpass");
		selectWindow("1stpass");
		run("Options...", "iterations=12 count=3 black do=Dilate stack");
		run("Invert", "stack");
		run("Gaussian Blur...", "sigma=15 stack");
		run("Find Edges", "stack");
		imageCalculator("Add stack", "1stpass", "raw");
		close("raw");		
		
		imageCalculator("Min create stack", "ch2","ch3");
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
		
		imageCalculator("Average create stack", "ch2","ch3");
		rename("ave");
		run("Morphological Filters (3D)", "operation=Gradient element=Ball x-radius=3 y-radius=3 z-radius=2");
		run("Mean...", "radius=10 stack");
		run("Find Edges", "stack");
		rename("gradedge");
		selectWindow("ave");
		run("Morphological Filters (3D)", "operation=Laplacian element=Octagon x-radius=10 y-radius=10 z-radius=3");
		rename("laplacian");
		run("Duplicate...", "title=lapedge duplicate");
		run("Mean...", "radius=10 stack");
		run("Find Edges", "stack");
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
		run("Convert to Mask", "method=Shanbhag background=Dark calculate black");
		run("Invert", "stack");
		run("Fill Holes", "stack");
		run("Outline", "stack");
		selectWindow("nucleus");
		run("Duplicate...", "duplicate");
		run("Options...", "iterations=10 count=1 black do=Dilate stack");
		run("Options...", "iterations=10 count=2 black do=Dilate stack");
		imageCalculator("Subtract stack", "laplacian-2","nucleus-1");
		close("nucleus-1");
		selectWindow("laplacian-2");
		run("Fill Holes", "stack");
		run("Options...", "iterations=3 count=4 black do=Close stack");
		run("Fill Holes", "stack");
		run("Options...", "iterations=10 count=4 black do=Erode stack");
		run("Options...", "iterations=10 count=4 black do=Dilate stack");
		run("Analyze Particles...", "size=15-Infinity show=Masks clear stack");
		run("Invert LUT");
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
		selectWindow("max");
		run("Duplicate...", "duplicate");
		run("Morphological Filters (3D)", "operation=Dilation element=Ball x-radius=4 y-radius=4 z-radius=2");
		run("Mean...", "radius=10 stack");
		rename("edgemax2");
		imageCalculator("Average stack", "edgemax","edgemax2");
		close("edgemax2");
		selectWindow("edgemax");
		run("Find Edges", "stack");
		selectWindow("max-1");
		run("Convert to Mask", "method=Minimum background=Dark calculate black");
		imageCalculator("Add stack", "max-1","laplacian");
		close("laplacian");
		run("Duplicate...", "title=edges duplicate");
		run("Gaussian Blur...", "sigma=10 stack");
		run("Find Edges", "stack");
		
		selectWindow("max-1");
		run("Options...", "iterations=15 count=3 black do=Dilate stack");
		run("Gaussian Blur...", "sigma=20 stack");
		run("Morphological Filters (3D)", "operation=Erosion element=Ball x-radius=5 y-radius=5 z-radius=2");
		run("Find Edges", "stack");
		rename("erosion");
		close("max-1");
		
		imageCalculator("Add create stack", "gradedge","filteredge");
		rename("addedge");
		imageCalculator("Add stack", "edges","edgemax");
		imageCalculator("Add stack", "lapedge","erosion");
		imageCalculator("Add stack", "addedge","edges");
		imageCalculator("Add stack", "addedge","lapedge");
		imageCalculator("Add stack", "addedge","1stpass");
		run("Convert to Mask", "method=Percentile background=Dark calculate black");
		imageCalculator("Average create stack", "gradedge","filteredge");
		rename("aveedge");
		imageCalculator("Average stack", "edges","edgemax");
		imageCalculator("Average stack", "lapedge","erosion");
		imageCalculator("Average stack", "aveedge","edges");
		imageCalculator("Average stack", "aveedge","lapedge");
		imageCalculator("Average stack", "aveedge","1stpass");
		run("Convert to Mask", "method=Li background=Dark calculate black");
		imageCalculator("Max create stack", "gradedge","filteredge");
		rename("maxedge");
		imageCalculator("Max stack", "edges","edgemax");
		imageCalculator("Max stack", "lapedge","erosion");
		imageCalculator("Max stack", "maxedge","edges");
		imageCalculator("Max stack", "maxedge","lapedge");
		imageCalculator("Max stack", "maxedge","1stpass");
		run("Convert to Mask", "method=Percentile background=Dark calculate black");
		imageCalculator("Average stack", "addedge","aveedge");
		imageCalculator("Average stack", "addedge","maxedge");
		selectWindow("addedge");
		setThreshold(225, 255);
		run("Convert to Mask", "method=Default background=Dark black");
		run("Invert", "stack");
		run("Options...", "iterations=25 count=4 black do=Erosion stack");
		run("Gaussian Blur...", "sigma=20 stack");
		run("Convert to Mask", "method=Default background=Dark calculate black");
		run("Watershed Irregular Features", "erosion=1 convexity_threshold=0.98 separator_size=0-Infinity stack");
		run("Stack to Images");
		selectWindow("nucleus");
		run("Duplicate...", "title=nucli duplicate");
		run("Stack to Images");
		close("lapedge");
		close("filteredge");
		close("gradedge");
		close("edgemax");
		close("edges");
		close("erosion");
		close("1stpass");
		close("aveedge");
		close("maxedge");
		close("max");
		close("filter");
		close("ch1");
		close("ch2");
		close("ch3");
		
		//slices=33;
		for (i = 1; i < slices+1; i++) {
			if (i < 10) {
				run("BinaryReconstruct ", "mask=addedge-000"+i+" seed=nucli-000"+i+" white");
				close("addedge-000"+i);
			} else if (i < 100) {
				run("BinaryReconstruct ", "mask=addedge-00"+i+" seed=nucli-00"+i+" white");
				close("addedge-00"+i);
			} else {
				run("BinaryReconstruct ", "mask=addedge-0"+i+" seed=nucli-0"+i+" white");
				close("addedge-0"+i);
			}
		}
		run("Images to Stack", "name=cellbodies title=nucli");
		run("Options...", "iterations=5 count=4 black do=Dilate stack");
		run("Analyze Particles...", "size=0-Infinity show=Nothing clear add stack");
		close("cellbodies");
		
		//open the original file again
		run("Bio-Formats Importer", "open=["+fpath+"] color_mode=Grayscale rois_import=[ROI manager] split_channels view=Hyperstack stack_order=XYCZT");
		for (i = 0; i < nImages; i++) {	
			selectImage(i+1);
			img = getTitle();
			if (img != "nucleus") {
				rename("ch"+toString(i));
			}
		}
		n = roiManager("count");
		for (i=0; i<n; i++) {
			selectWindow("ch3");
		    roiManager("select", i);
			run("Make Band...", "band=4");
			run("Duplicate...", "title=ch3-roi-"+i);
			run("Clear Outside");
		}
		run("Images to Stack", "name=ch3-roi title=ch3-roi");
		
		for (i=0; i<n; i++) {
			selectWindow("ch2");
		    roiManager("select", i);
			run("Make Band...", "band=4");
			run("Duplicate...", "title=ch2-roi-"+i);
			run("Clear Outside");
		}
		run("Images to Stack", "name=ch2-roi title=ch2-roi");
		run("Duplicate...", "title=synapse duplicate");
		run("Convert to Mask", "method=IJ_IsoData background=Dark calculate black");
		run("Despeckle", "stack");
		run("Options...", "iterations=5 count=4 black do=Open stack");
		run("Options...", "iterations=5 count=4 black do=Dilate stack");
		run("Watershed Irregular Features", "erosion=1 convexity_threshold=0.99 separator_size=0-Infinity stack");
		run("Set Measurements...", "area mean standard fit median stack redirect=None decimal=3");
		run("Analyze Particles...", "size=0.25-5 show=Nothing clear add stack");
		
		run("Clear Results");
		n = roiManager("count");
		for (i=0; i<n; i++) {
			selectWindow("ch2-roi");
		    roiManager("select", i);
			roiManager("Measure");
			run("To Bounding Box");
			horizontal = getValue("Width");
			vertical = getValue("Height");
			if (vertical >= horizontal) {
				setKeyDown("alt");
				profile = getProfile();
				setKeyDown("none");
			} else {
				profile = getProfile();
			}
			for (j=0; j<profile.length; j++) {
				setResult("Value-"+j, i, profile[j]);
			}
			updateResults();
		}
		//save Results table and image
		//saveAs("Results", savePath + "-Results-ch2.csv");
		run("Clear Results");
		//saveAs("TIFF", mpath+"ch2.tif");
		n = roiManager("count");
		for (i=0; i<n; i++) {
			selectWindow("ch3-roi");
		    roiManager("select", i);
			roiManager("Measure");
			run("To Bounding Box");
			horizontal = getValue("Width");
			vertical = getValue("Height");
			if (vertical >= horizontal) {
				setKeyDown("alt");
				profile = getProfile();
				setKeyDown("none");
			} else {
				profile = getProfile();
			}
			for (j=0; j<profile.length; j++) {
				setResult("Value-"+j, i, profile[j]);
			}
			updateResults();
		}
		//save Results table and image
		//saveAs("Results", savePath + "-Results-ch3.csv");
		run("Clear Results");
		//saveAs("TIFF", mpath+"ch3.tif");
		roiManager("reset");
		
		selectWindow("nucleus");
		run("Analyze Particles...", "size=0.00-infinity show=Nothing clear add stack");
		
		run("Clear Results");
		n = roiManager("count");
		chosen = Math.round(random * n);
		selectWindow("ch2");
		for (i = 1; i <= nSlices; i++) {
			setSlice(i);
			roiManager("select", chosen);
		    Roi.setPosition(i);
		    roiManager("Measure");
		    roiManager("select", chosen+1);
		    Roi.setPosition(i);
		    roiManager("Measure");

		}
		//save Results table and image
		//saveAs("Results", savePath + "-Background-ch2.csv");
		//run("Clear Results");

		selectWindow("ch3");
		chosen = 10;
		for (i = 1; i <= nSlices; i++) {
			setSlice(i);
			roiManager("select", chosen);
		    Roi.setPosition(i);
		    roiManager("Measure");
		    roiManager("select", chosen+1);
		    Roi.setPosition(i);
		    roiManager("Measure");
		}
		//save Results table and image
		//saveAs("Results", savePath + "-Background-ch3.csv");
		run("Clear Results");
		roiManager("reset");
		run("Close");
		close("Results");
	}
}

close("Results");

print("Finished processing all images in "+imageDir+".");
print("Took " + ((getTime()-programStart)/1000)/60 + "minutes");
//setBatchMode(false);

