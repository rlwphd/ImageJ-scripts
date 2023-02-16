// Getting the folder where images are located
imageDir=getDirectory("Select Spinal Cord Folder for Image Processing");
roiDir=imageDir+"Analysis_Results\\";
imagelist = getFileList(imageDir);
imageArray = newArray;
File.makeDirectory(roiDir);
saveSettings();
setOption("BlackBackground", true);
setBackgroundColor(0, 0, 0);
run("Set Measurements...", "area mean standard centroid fit median stack redirect=None decimal=3");
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
			run("Max...", "value=4095 stack");			
			run("8-bit");
		}
		xScaled = 1;
		yScaled = 1;
		zScaled = 1;
		toScaled(xScaled,yScaled,zScaled);
		
		fastblue=false
		selectWindow("ch1");
			
		
		
		
		if (fastblue != true) {
			setForegroundColor(255, 255, 255);
			setBackgroundColor(0, 0, 0);
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
			run("Statistical Region Merging", "q=30 showaverages 3d");
			rename("filternucleus");
			run("8-bit");
			run("Subtract...", "value=15 stack");
			setThreshold(0, 0);
			run("Convert to Mask", "method=MinError background=Dark black");
			run("Despeckle", "stack");
			run("Fill Holes", "stack");
			run("Outline", "stack");
			run("Options...", "iterations=5 count=2 black do=Dilate stack");
			run("Invert", "stack");
			run("Options...", "iterations=5 count=3 black do=Erode stack");
			if (xScaled == yScaled) {
				run("Set Scale...", "distance="+(1/xScaled)+" known=1 unit=um global");
			} else {
				run("Set Scale...", "distance="+(1/xScaled)+" known=1 pixel="+(yScaled/xScaled)+" unit=um global");
			}
			run("Analyze Particles...", "size=45-200 circularity=0.2-1.00 show=Masks clear stack");
			run("Invert LUT");
			run("Fill Holes", "stack");
			rename("nucli");
			close("filternucleus");
			selectWindow("nucli");
			run("Mean...", "radius=25 stack");
			run("Convert to Mask", "method=Minimum background=Dark calculate black");
			close("Log");			
			run("Options...", "iterations=10 count=1 black do=Dilate stack");
			run("Options...", "iterations=20 count=2 black do=Dilate stack");
			run("Options...", "iterations=10 count=3 black do=Dilate stack");
			run("Analyze Particles...", "size=90-350 circularity=0.6-1.00 show=Ellipses display exclude clear stack");
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
			run("Duplicate...", "title=nucli duplicate");
			run("Options...", "iterations=5 count=1 black do=Dilate stack");
			run("Options...", "iterations=5 count=2 black do=Dilate stack");
			run("Options...", "iterations=5 count=3 black do=Dilate stack");
			run("Options...", "iterations=5 count=4 black do=Dilate stack");	
			imageCalculator("Subtract stack", "ch2","nucli");
			imageCalculator("Subtract stack", "ch3","nucli");
			
			selectWindow("ch3");
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
			imageCalculator("Add stack", "ch2","minline");
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
			
			selectWindow("ch3");
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
			imageCalculator("Add stack", "ch3","minline");
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
			selectWindow("3line");
			run("Options...", "iterations=10 count=3 black do=Dilate stack");
			imageCalculator("Subtract stack", "3line","nucli");					
			imageCalculator("AND stack", "3line", "rawline");
			imageCalculator("Add stack", "max", "3line");
			close("rawline");
			close("3line");
			close("line3");
			
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
			
			selectWindow("ch3");
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
			
			if (xScaled == yScaled) {
				run("Set Scale...", "distance="+(1/xScaled)+" known=1 unit=um global");
			} else {
				run("Set Scale...", "distance="+(1/xScaled)+" known=1 pixel="+(yScaled/xScaled)+" unit=um global");
			}
					
			selectWindow("ch3");
			run("Invert", "stack");
			imageCalculator("AND create stack", "selected0","ch3");
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
					
			imageCalculator("AND create stack", "selected0","ch3");
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
			close("ch1");
			close("ch2");
			close("ch3");
		}
		
		//open the original file again
		run("Bio-Formats Importer", "open=["+fpath+"] color_mode=Grayscale rois_import=[ROI manager] split_channels view=Hyperstack stack_order=XYCZT");
		chan = 1;
		for (i = 1; i < nImages+1; i++) {	
			selectImage(i);
			img = getTitle();
			if (!img.contains("cell")) {
				rename("ch"+toString(chan));
				chan++;
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
			selectWindow("ch2");
			run("Duplicate...", "title=ch2-cell"+i+" duplicate");
			makeRectangle(cropParam[4*(i-1)], cropParam[4*(i-1)+1], cropParam[4*(i-1)+2], cropParam[4*(i-1)+3]);
			run("Crop");
			run("Select None");
			selectWindow("ch3");
			run("Duplicate...", "title=ch3-cell"+i+" duplicate");
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
		selectWindow("ch2");
		rename("ch2-cell1");
		makeRectangle(cropParam[0], cropParam[1], cropParam[2], cropParam[3]);
		run("Crop");
		run("Select None");
		selectWindow("ch3");
		rename("ch3-cell1");
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
			
			selectWindow("ch2-cell"+j);
			run("Duplicate...", "title=ch2-celledge"+j+" duplicate");
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
			selectWindow("ch3-cell"+j);
			run("Duplicate...", "title=ch3-celledge"+j+" duplicate");
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
			close("ROI Manager");
			close("Results");		
			
			selectWindow("cell"+j);
			run("Analyze Particles...", "size=100-infinity show=Nothing display clear add stack");
			s = Table.getColumn('Slice',"Results");
			roiManager("Show None");
			n = roiManager("count");
			
			selectWindow("ch2-cell"+j);
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
			    Table.set("Ch2-Cell-"+j+" Mean", k, bmean[k], backgroundTable);
				Table.set("Ch2-Cell-"+j+" StdDev", k, bstd[k], backgroundTable);
			}
				
			selectWindow("ch2-celledge"+j);
			setBackgroundColor(0, 0, 0);
			for (k = 0; k < n; k++) {
		    	roiManager("select", k);
		    	run("Clear", "slice");
		    	run("Select None");
		    	run("Subtract...", "value="+(bmean[k]+2*bstd[k])+" slice");
			}
						
			selectWindow("ch3-cell"+j);
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
			    Table.set("Ch2-Cell-"+j+" Mean", k, bmean[k], backgroundTable);
				Table.set("Ch2-Cell-"+j+" StdDev", k, bstd[k], backgroundTable);
			}

			selectWindow("ch3-celledge"+j);
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

		run("3D OC Options", "volume surface centroid mean_distance_to_surface std_dev_distance_to_surface bounding_box dots_size=20 font_size=50 redirect_to=none");
		setForegroundColor(255,255,255);
		for (i=2; i<=3; i++) {
			for (j=1; j<=numCells; j++) {
				selectWindow("ch"+i+"-celledge"+j);
				i=2;
				j=2;
				setThreshold(1,4095);
				run("Analyze Particles...", "size=0-infinity show=Nothing display clear stack");
				s = Table.getColumn('Slice',"Results");
				Array.getStatistics(s, min, max);
				resetThreshold();
				for (k=min; k<=max; k++) {
					selectWindow("ch"+i+"-celledge"+j);
					setSlice(k);
					pixelMax = getValue("Max");
					pixelMin = maxOf(getValue("Min"),1);
					run("Duplicate...", "title=ch"+i+"-cell"+j+"-selection"+k);
					run("Select All");
					wait(10);
					run("Clear", "slice");
					wait(10);
					run("Select None");
					run("8-bit");
					for (m=0; m<20; m++) {
						selectWindow("ch"+i+"-celledge"+j);
						run("Select None");
						skip = false;
						if (m==0 && pixelMax > 2047) {
							lowValue = pixelMax-100;
							setThreshold(lowValue, pixelMax+1);
							run("Analyze Particles...", "size=0.1-5 show=Nothing display clear slice");
							x = Table.getColumn('X');
							y = Table.getColumn('Y');
							toUnscaled(x, y);
							for (n=0; n<nResults; n++) {
								selectWindow("ch"+i+"-celledge"+j);
								doWand(x[n],y[n]);
								selectWindow("ch"+i+"-cell"+j+"-selection"+k);
								run("Restore Selection");
								run("Fill", "slice");
								run("Select None");
							}
							skip = true;
						} else if (m==0 && pixelMax <= 2047) {
							lowValue = pixelMax-100;
						} else if (m<3) {
							lowValue = (pixelMax-100)*(4-m)/4;
						} else {
							lowValue = (pixelMax-100)/2 - m*((pixelMax-100)/2-pixelMin)/20;
						}
						if (!skip) {
							setThreshold(lowValue, pixelMax+1);
							run("Analyze Particles...", "size=0.1-5 show=Nothing display clear slice");
							f = Table.getColumn("Feret");
							r = Table.getColumn("Round");
							x = Table.getColumn('X');
							y = Table.getColumn('Y');
							toUnscaled(x, y);
							for (n=0; n<nResults; n++) {
								if (f[n] >= 1 && f[n] <= 3 && r[n] > 0.45) {
									selectWindow("ch"+i+"-celledge"+j);
									doWand(x[n],y[n]);
									selectWindow("ch"+i+"-cell"+j+"-selection"+k);
									run("Restore Selection");
									run("Fill", "slice");
									run("Select None");
								}
							}
						}
						selectWindow("ch"+i+"-celledge"+j);
						run("Select None");
						resetThreshold();
					}
					selectWindow("ch"+i+"-cell"+j+"-selection"+k);
					run("Options...", "iterations=1 count=4 black do=Open");
					run("Options...", "iterations=1 count=4 black do=Dilate");
					run("Analyze Particles...", "size=0.1-5 show=Masks display clear slice");
					f = Table.getColumn("Feret");
					r = Table.getColumn("Round");
					x = Table.getColumn('X');
					y = Table.getColumn('Y');
					toUnscaled(x, y);
					run("Invert LUT");
					rename("selected"+k);
					for (n=0; n<nResults; n++) {
						if (f[n] < 0.75 || r[n] < 0.4) {
							selectWindow("selected"+k);
							doWand(x[n],y[n]);
							run("Clear", "slice");
							run("Select None");
						}
					}
					close("ch"+i+"-cell"+j+"-selection"+k);
				}
				run("Images to Stack", "name=ch"+i+"-cell"+j+"-analysis title=selected");
				close("Results");
				run("3D Objects Counter", "threshold=128 slice=1 min.=10 max.=4129270 statistics");
				bdepth = Table.getColumn("B-depth");
				bz = Table.getColumn("BZ");
				x = Table.getColumn('X');
				y = Table.getColumn('Y');
				i = j = 2;
				for (n=0; n<bdepth.length; n++) {
					if (bdepth[n] == 1) {
						selectWindow("ch"+i+"-cell"+j+"-analysis");
						setSlice(bz[n]);
						doWand(x[n], y[n]);
						wait(10);
						run("Clear", "slice");
						wait(10);
					}
				}
			}
		}
	
	
	Table.deleteRows(0, 23);
Table.sort("Circ.");
		
		selectWindow("ch3-celledge2");
		slices = nSlices;
		for (s = 15; s <= slices; s++) {
			selectWindow("ch3-celledge2");
		    run("Z Project...", "start="+s-2+" stop="+s+2+" projection=[Average Intensity]");
			rename("Ave"+s);
		}
		run("Images to Stack", "name=5 title=Ave");
		
		
		//perform analysis slice by slice
		//get max value
		//threshold the image at the different thresholds - jump 10-15 at a time starting at max and going to max/2
		//analyze images - create masks for size between 0.25-1 um
		//merge masks into a stack and then potentially sum stack
		// or address differences by shape descriptors
		// need to also address those that are already bigger than 1 um at max value - potentially split

		
		
		
		
		
		len = 2;
		if (len>1) {
			setForegroundColor(255, 255, 255);
			selectWindow("cell-1");
			run("Duplicate...", "title=merge duplicate");
			// need to change this to be imageCalculator and run through all of the cells really fast
			for (i=1; i<len+1; i++) {
				selectWindow("cell-"+i);
				run("Z Project...", "projection=[Sum Slices]");
				rename("zstack-"+i);
				run("8-bit");
				
				
			}
		}
		for (j=1; j<len+1; j++) {
			selectWindow("cell-"+j);
			run("Analyze Particles...", "size=0.00-infinity show=Nothing clear add stack");
			n = roiManager("count");
						
			selectWindow("ch3");
			run("Duplicate...", "title=ch3-cell-"+j+" duplicate");
			run("Duplicate...", "title=ch3-background-"+j+" duplicate");
			selectWindow("ch2");
			run("Duplicate...", "title=ch2-cell-"+j+" duplicate");
			run("Duplicate...", "title=ch2-background-"+j+" duplicate");

			selectWindow("cell-"+j);
			lslice = nSlices+1;
			for (k=1; k<lslice; k++) {
				selectWindow("cell-"+j);
				setSlice(k);
				run("Select None");
				mean = getValue("Mean");
				if (mean < 2) {
					selectWindow("ch3-background-"+j);
					setSlice(k);
					run("Select All");
					run("Clear", "slice");
					run("Select None");
					
					selectWindow("ch2-background-"+j);
					setSlice(k);
					run("Select All");
					run("Clear", "slice");
					run("Select None");
					
					selectWindow("ch3-cell-"+j);
					setSlice(k);
					run("Select All");
					run("Clear", "slice");
					run("Select None");
					
					selectWindow("ch2-cell-"+j);
					setSlice(k);
					run("Select All");
					run("Clear", "slice");
					run("Select None");
				}
			}			
			
			for (i=0; i<n; i++) {
				selectWindow("ch3-background-"+j);
			    roiManager("select", i);
				run("Clear Outside", "slice");
				run("Select None");
				selectWindow("ch3-cell-"+j);
			    roiManager("select", i);
				run("Make Band...", "band=4");
				run("Clear Outside", "slice");
				run("Select None");
				selectWindow("ch2-background-"+j);
			    roiManager("select", i);
			    run("Clear Outside", "slice");
				run("Select None");
				selectWindow("ch2-cell-"+j);
			    roiManager("select", i);
				run("Make Band...", "band=4");
				run("Clear Outside", "slice");
				run("Select None");
				if (len>1 && j==1) {
					selectWindow("merge");
					roiManager("select", i);
					run("Make Band...", "band=4");
					run("Fill", "slice");
				} else if (len>1 && j>1) {
					selectWindow("merge");
					run("Duplicate...", "title=merge"+j+" duplicate");
					roiManager("select", i);
					
			}
			

			if (j>1) {
				imageCalculator("Add stack", "ch3-background-1", "ch3-background-"+j);
				imageCalculator("Add stack", "ch2-background-1", "ch2-background-"+j);
				imageCalculator("Add stack", "ch3-cell-1", "ch3-cell-"+j);
				imageCalculator("Add stack", "ch2-cell-1", "ch2-cell-"+j);
				close("ch3-background-"+j);
				close("ch2-background-"+j);
				close("ch3-cell-"+j);
				close("ch2-cell-"+j);
			}
		}	
			
		selectWindow("ch2-cell-1");
		run("Duplicate...", "title=synapse duplicate");
		imageCalculator("Add stack", "ch3-cell-1", "ch3-background-1");
		imageCalculator("Add stack", "ch2-cell-1", "ch2-background-1");
		close("ch3-background-1");
		close("ch2-background-1");
		
		for (i=1; i<n; i++) {
			selectWindow("ch2");
			roiManager("select", i);
			baseline = getValue("Mean");
			basestd = getValue("StdDev");
			selectWindow("synapse");
			roiManager("select", i);
			run("Select None");
		    val = baseline+2*basestd;
		    run("Subtract...", "value="+val+" slice");
		}
		run("Convert to Mask", "method=IsoData background=Dark calculate black");
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
}

close("Results");

print("Finished processing all images in "+imageDir+".");
print("Took " + ((getTime()-programStart)/1000)/60 + "minutes");
//setBatchMode(false);





		//Don't need? Pretty sure I don't
		imageCalculator("Average stack", "max","edgemax");
		imageCalculator("Average stack", "lapedge","laplacian");
		imageCalculator("Average stack", "gradedge","entropy");
		imageCalculator("Average stack", "max","lapegde");
		imageCalculator("Average stack", "max","gradedge");
		close("edgemax");
		close("laplacian");
		close("lapedge");
		close("gradedge");
		selectWindow("max");
		run("Duplicate...", "duplicate");
		run("Convert to Mask", "method=IsoData background=Dark calculate black");
		selectWindow("max");
		run("Convert to Mask", "method=Intermodes background=Dark calculate black");
		
		//Also Don't need
		run("Gaussian Blur...", "sigma=20 stack");
		run("Convert to Mask", "method=Default background=Dark calculate black");
		run("Watershed Irregular Features", "erosion=1 convexity_threshold=0.98 separator_size=0-Infinity stack");
		run("Stack to Images");
		selectWindow("nucleus");
		run("Stack to Images");
		//close("ch1");
		//close("ch2");
		//close("ch3");

		//Also don't need
		selectWindow("merge");
		run("Duplicate...", "title=bottom");
		run("Select All");
		run("Clear", "slice");
		run("Select None");
		run("Duplicate...", "title=top");
		selectWindow("merge");
		run("Analyze Particles...", "size=100-infinity show=Nothing clear add stack");
		roiManager("Show None");
		n = roiManager("count");
		//slices=33;
		for (i = 0; i < n; i++) {
			b = 0;
			t = 0;
			selectWindow("merge");
			roiManager("select", i);
			s = getSliceNumber();
			x = getValue("X");
			y = getValue("Y");
			if (s > 1) {
				selectWindow("merge");
				setSlice(s-1);
				p = getValue("Mean");
				run("Select None");
				if (p > 5 ) {
					doWand(x,y);
					selectWindow("bottom");
					run("Restore Selection");
					run("Fill", "slice");
					run("Select None");
					b = 1;
				}
			}
			
			if (s < slices) {
				selectWindow("merge");
				setSlice(s+1);
				p = getValue("Mean");
				run("Select None");
				if (p > 5 ) {
					doWand(x,y);
					selectWindow("top");
					run("Restore Selection");
					run("Fill", "slice");
					run("Select None");
					t = 1;
				}
			}
			
			if (b == 1 && t == 1) {
				imageCalculator("AND stack", "bottom","top");
				selectWindow("bottom");
				doWand(x,y);
				selectWindow("merge");
				setSlice(s);
				run("Restore Selection");
				run("Fill", "slice");
				run("Select None");
				print(s, " Worked ", "X&Y ", x, " ", y);				
				selectWindow("bottom");
				run("Select All");
				run("Clear", "slice");
				run("Select None");
				selectWindow("top");
				run("Select All");
				run("Clear", "slice");
				run("Select None");
			} else if (b == 1 && t == 0) {

				selectWindow("bottom");
				run("Select All");
				run("Clear", "slice");
				run("Select None");
				selectWindow("top");
				run("Select All");
				run("Clear", "slice");
				run("Select None");			
			} else if (b == 0 && t == 1) {

				selectWindow("bottom");
				run("Select All");
				run("Clear", "slice");
				run("Select None");
				selectWindow("top");
				run("Select All");
				run("Clear", "slice");
				run("Select None");				
			}
		}
		
				
						
								
												
		//slices=25;
		for (i = 1; i < slices+1; i++) {
			if (i < 10) {
				run("BinaryReconstruct ", "mask=addedge-000"+i+" seed=nucleus-000"+i+" white");
				close("addedge-000"+i);
			} else if (i < 100) {
				run("BinaryReconstruct ", "mask=addedge-00"+i+" seed=nucleus-00"+i+" white");
				close("addedge-00"+i);
			} else {
				run("BinaryReconstruct ", "mask=addedge-0"+i+" seed=nucleus-0"+i+" white");
				close("addedge-0"+i);
			}
		}
		run("Images to Stack", "name=cellbodies title=nucleus");
		run("Options...", "iterations=5 count=2 black do=Dilate stack");
		run("Options...", "iterations=3 count=1 black do=Dilate stack");
		run("Fill Holes", "stack");

		
		
		
		
		// Might need this
		run("3D OC Options", "volume surface nb_of_obj._voxels nb_of_surf._voxels integrated_density mean_gray_value std_dev_gray_value median_gray_value minimum_gray_value maximum_gray_value centroid mean_distance_to_surface std_dev_distance_to_surface median_distance_to_surface centre_of_mass bounding_box dots_size=20 font_size=50 store_results_within_a_table_named_after_the_image_(macro_friendly) redirect_to=none");
		run("3D Objects Counter", "threshold=128 slice=12 min.=10 max.=26214400 objects statistics summary");
		selectWindow("Statistics for cellbodies");
		len = Table.size;
		for (i=1; i<len+1; i++) {
			selectWindow("Objects map of cellbodies");
			run("Duplicate...", "duplicate");
			setThreshold(i, i);
			run("Convert to Mask", "method=Default background=Dark black");
			rename("cell-"+toString(i));
		}
		close("Objects map of cellbodies");
		//close("cellbodies");
		close("Statistics for cellbodies");
		close("Log");
		
		
		
		selectWindow("cellbodies");
		imageCalculator("Subtract create stack", "cellbodies","entropy");
		selectWindow("Result of cellbodies");
		rename("cutout");
		run("Outline");

		selectWindow("ch3");
		run("Duplicate...", "title=entropy duplicate");
		run("Convert to Mask", "method=MaxEntropy background=Dark black");
		run("Despeckle", "stack");		
		selectWindow("entropy");
		run("Duplicate...", "duplicate");
		run("Maximum...", "radius=6 stack");
		run("Options...", "iterations=1 count=2 black do=Dilate stack");
		run("Outline", "stack");
		run("Options...", "iterations=2 count=2 black do=Dilate stack");
		run("Invert", "stack");
		run("Analyze Particles...", "size=0-35 show=Masks clear stack");
		run("Invert LUT");
		rename("particles");
		imageCalculator("AND stack", "particles","entropy");
		close("entropy-1");
		selectWindow("particles");
		run("Analyze Particles...", "size=5-infinity show=Masks clear stack");
		run("Invert LUT");
		rename("particles-1");
		imageCalculator("Subtract stack", "particles","particles-1");
		close("particles-1");
		imageCalculator("Subtract create stack", "entropy","particles");
		rename("entropy-1");
		run("Maximum...", "radius=6 stack");
		run("Options...", "iterations=2 count=2 black do=Dilate stack");
		run("Invert", "stack");
		run("Options...", "iterations=10 count=4 black do=Dilate stack");
		run("Analyze Particles...", "size=200-infinity show=Masks clear stack");
		run("Invert LUT");
		rename("cells");
		close("entropy-1");
		close("particles");
		selectWindow("cells");
		run("Duplicate...", "duplicate");
		run("Outline", "stack");
		run("Options...", "iterations=10 count=1 black do=Dilate stack");
		run("Options...", "iterations=10 count=2 black do=Dilate stack");
		run("Invert", "stack");
		run("Analyze Particles...", "size=100-1500 show=Masks clear stack");
		run("Invert LUT");
		rename("selected");
		close("cells-1");		
		imageCalculator("AND stack", "cells","selected");
		close("selected");
		selectWindow("cells");
		run("Options...", "iterations=10 count=1 black do=Dilate stack");
		run("Options...", "iterations=10 count=2 black do=Dilate stack");
		run("Options...", "iterations=15 count=4 black do=Dilate stack");
		run("Fill Holes", "stack");
		
		
		
		
		
		run("Stack to Images");
		selectWindow("cells-0001");
		run("Duplicate...", "title=down-0001");	
		slices=25;
		for (i = 1; i < slices; i++) {
			j = i+1;
			if (j < 10) {
				imageCalculator("XOR create", "cells-000"+i, "cells-000"+j);
				rename("result");
				run("Options...", "iterations=2 count=2 black do=Open");
				run("Analyze Particles...", "size=30-200 show=Masks clear");
				run("Invert LUT");
				rename("delete");
				imageCalculator("Subtract create", "cells-000"+j, "delete");
				rename("down-000"+j);
				imageCalculator("Subtract create", "cells-000"+i, "delete");
				rename("up-000"+i);
				close("result");
				close("delete");
			} else if (j == 10) {
				imageCalculator("XOR create", "cells-000"+i, "cells-00"+j);
				rename("result");
				run("Options...", "iterations=2 count=2 black do=Open");
				run("Analyze Particles...", "size=30-200 show=Masks clear");
				run("Invert LUT");
				rename("delete");
				imageCalculator("Subtract create", "cells-00"+j, "delete");
				rename("down-00"+j);
				imageCalculator("Subtract create", "cells-000"+i, "delete");
				rename("up-000"+i);
				close("result");
				close("delete");
			} else {
				imageCalculator("XOR create", "cells-00"+i, "cells-00"+j);
				rename("result");
				run("Options...", "iterations=2 count=2 black do=Open");
				run("Analyze Particles...", "size=30-200 show=Masks clear");
				run("Invert LUT");
				rename("delete");
				imageCalculator("Subtract create", "cells-00"+j, "delete");
				rename("down-00"+j);
				imageCalculator("Subtract create", "cells-00"+i, "delete");
				rename("up-00"+i);
				close("result");
				close("delete");
			}
		}
		selectWindow("cells-00"+slices);
		run("Duplicate...", "title=up-00"+slices);
		run("Images to Stack", "name=up title=up");
		run("Images to Stack", "name=down title=down");
		run("Images to Stack", "name=cells title=cells");
		close("cells");
		imageCalculator("AND create stack", "up", "down");
		rename("cells");
		close("up");
		close("down");
		
		selectWindow("entropy");
		run("Duplicate...", "title=cells duplicate");
		run("Stack to Images");
		selectWindow("cells-0001");
		run("Duplicate...", "title=down-0001");	
		slices=25;
		for (i = 1; i < slices; i++) {
			j = i+1;
			if (j < 10) {
				imageCalculator("XOR create", "cells-000"+i, "cells-000"+j);
				rename("result");
				run("Options...", "iterations=2 count=2 black do=Open");
				run("Options...", "iterations=1 count=2 black do=Dilate");
				imageCalculator("Subtract create", "cells-000"+j, "result");
				rename("down-000"+j);
				imageCalculator("Subtract create", "cells-000"+i, "result");
				rename("up-000"+i);
				close("result");
			} else if (j == 10) {
				imageCalculator("XOR create", "cells-000"+i, "cells-00"+j);
				rename("result");
				run("Options...", "iterations=2 count=2 black do=Open");
				run("Options...", "iterations=1 count=2 black do=Dilate");
				imageCalculator("Subtract create", "cells-00"+j, "result");
				rename("down-00"+j);
				imageCalculator("Subtract create", "cells-000"+i, "result");
				rename("up-000"+i);
				close("result");
			} else {
				imageCalculator("XOR create", "cells-00"+i, "cells-00"+j);
				rename("result");
				run("Options...", "iterations=2 count=2 black do=Open");
				run("Options...", "iterations=1 count=2 black do=Dilate");
				imageCalculator("Subtract create", "cells-00"+j, "result");
				rename("down-00"+j);
				imageCalculator("Subtract create", "cells-00"+i, "result");
				rename("up-00"+i);
				close("result");
			}
		}
		selectWindow("cells-00"+slices);
		run("Duplicate...", "title=up-00"+slices);
		run("Images to Stack", "name=up title=up");
		run("Images to Stack", "name=down title=down");
		run("Images to Stack", "name=cells title=cells");
		close("cells");
		imageCalculator("AND create stack", "up", "down");
		rename("cells");
		close("up");
		close("down");

		run("3D OC Options", "centroid centre_of_mass bounding_box dots_size=20 font_size=50 store_results_within_a_table_named_after_the_image_(macro_friendly) redirect_to=none");
		run("3D Objects Counter", "threshold=128 slice=16 min.=10 max.=34603008 objects statistics");
		bz = Table.getColumn("B-depth");
		
		imageCalculator("Add create stack", "ch2","ch3");
		selectWindow("Result of ch2");
		rename("add");
		run("Convert to Mask", "method=MaxEntropy background=Dark black");
		
		
		
		
		

