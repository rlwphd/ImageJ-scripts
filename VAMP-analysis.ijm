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
		run("Invert", "stack");
		run("Options...", "iterations=5 count=3 black do=Erode stack");
		run("Analyze Particles...", "size=5000-30000 circularity=0.2-1.00 show=Masks clear stack");
		run("Invert LUT");
		run("Fill Holes", "stack");
		rename("nucli");
		close("filternucleus");
		selectWindow("nucli");
		run("Mean...", "radius=25 stack");
		run("Convert to Mask", "method=Minimum background=Dark calculate black");
		run("Options...", "iterations=10 count=1 black do=Dilate stack");
		run("Options...", "iterations=20 count=2 black do=Dilate stack");
		run("Options...", "iterations=10 count=3 black do=Dilate stack");
		run("Analyze Particles...", "size=10000-55000 circularity=0.6-1.00 show=Ellipses exclude clear stack");
		run("Invert", "stack");
		run("Fill Holes", "stack");
		run("Options...", "iterations=10 count=2 black do=Erode stack");
		rename("nucleus");
		close("nucli");
		run("Duplicate...", "title=blank");
		s = Table.getColumn('Slice');
		x = Table.getColumn('X');
		y = Table.getColumn('Y');
		c = newArray(nResults);
		c[0] = 1;
		setForegroundColor(1,1,1);
		makeRectangle(x[0]-50, y[0]-50, 100, 100);
		run("Fill", "slice");
		numCells = 1;
		updateDisplay();
		for (i = 1; i < nResults(); i++) {
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
		imageCalculator("Subtract stack", "ch2","nucleus");
		imageCalculator("Subtract stack", "ch3","nucleus");
		
		selectWindow("ch3");
		run("Duplicate...", "title=raw duplicate");
		run("Gaussian Blur...", "sigma=15 stack");
		run("Find Edges", "stack");
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
		
		selectWindow("ch3");
		run("Duplicate...", "title=entropy duplicate");
		run("Convert to Mask", "method=MaxEntropy background=Dark black");
		run("Despeckle", "stack");		
		run("Outline", "stack");
		run("Options...", "iterations=10 count=4 black do=Dilate stack");
		run("Duplicate...", "duplicate");
		run("Options...", "iterations=20 count=4 black do=Dilate stack");
		run("Skeletonize", "stack");
		run("Duplicate...", "duplicate");
		selectWindow("entropy-1");
		run("Options...", "iterations=20 count=3 black do=Dilate stack");
		run("Options...", "iterations=1 count=3 black do=Open stack");
		run("Skeletonize", "stack");
		selectWindow("entropy-2");
		run("Options...", "iterations=20 count=2 black do=Dilate stack");
		run("Skeletonize", "stack");
		imageCalculator("Add stack", "entropy-1","entropy-2");
		close("entropy-2");
		run("Duplicate...", "duplicate");
		selectWindow("entropy-1");
		run("Options...", "iterations=30 count=3 black do=Dilate stack");
		run("Options...", "iterations=1 count=3 black do=Open stack");
		run("Skeletonize", "stack");
		selectWindow("entropy-2");
		run("Options...", "iterations=30 count=2 black do=Dilate stack");
		run("Skeletonize", "stack");
		imageCalculator("Add stack", "entropy-1","entropy-2");
		close("entropy-2");
		run("Options...", "iterations=10 count=3 black do=Dilate stack");
		run("Options...", "iterations=1 count=3 black do=Open stack");
		run("Skeletonize", "stack");
		run("Options...", "iterations=5 count=3 black do=Dilate stack");
		imageCalculator("Add stack", "entropy","entropy-1");
		close("entropy-1");
		run("Analyze Particles...", "size=0-7.5 show=Masks clear stack");
		run("Invert LUT");
		imageCalculator("Subtract stack", "entropy","Mask of entropy");
		close("Mask of entropy");
		run("Skeletonize", "stack");
		run("Options...", "iterations=5 count=3 black do=Dilate stack");
		
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
		imageCalculator("Add stack", "addedge","entropy");
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
		run("Duplicate...", "title=smoother duplicate");
		run("Options...", "iterations=25 count=4 black do=Erosion stack");
		run("Gaussian Blur...", "sigma=20 stack");
		run("Convert to Mask", "method=Default background=Dark calculate black");
		imageCalculator("Add stack", "addedge","smoother");
		close("smoother");
		selectWindow("addedge");
		run("Options...", "iterations=25 count=4 black do=Erosion stack");
		run("Gaussian Blur...", "sigma=20 stack");
		run("Convert to Mask", "method=Default background=Dark calculate black");
		run("Watershed Irregular Features", "erosion=1 convexity_threshold=0.98 separator_size=0-Infinity stack");
		run("Stack to Images");
		selectWindow("nucleus");
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
		close("entropy");
		close("Log");
		//close("ch1");
		//close("ch2");
		//close("ch3");
		
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
		run("Analyze Particles...", "size=0.00-infinity show=Nothing clear add stack");
		n = roiManager("count");
		for (i = 0; i < n; i++) {
		    roiManager("select", i);
			run("Make Band...", "band=50");
			run("Fill", "slice");
			run("Select None");
		}
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
		
		//open the original file again
		run("Bio-Formats Importer", "open=["+fpath+"] color_mode=Grayscale rois_import=[ROI manager] split_channels view=Hyperstack stack_order=XYCZT");
		chan = 1;
		for (i = 1; i < nImages+1; i++) {	
			selectImage(i);
			img = getTitle();
			if (img != "nucleus" && !img.contains("cell")) {
				rename("ch"+toString(chan));
				chan++;
			}
		}
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

