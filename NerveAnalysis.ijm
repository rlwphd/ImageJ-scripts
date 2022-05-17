names = newArray("Nerve1Bottom","Nerve1Middle","Nerve1Top","Nerve2Left","Nerve2Middle","Nerve2Right","Nerve3Bottom","Nerve3Middle","Nerve3Top","Nerve4Left","Nerve4Right","Nerve5left","Nerve5Middle","Nerve5TopRight","Nerve5BottomRight","Nerve6Middle","Nerve6Top","Nerve6Right","Nerve7Left","Nerve7Right","Nerve8Left","Nerve8Right","Nerve9Middle","Nerve9Top","Nerve9Bottom","Nerve10Right","Nerve10Bottom","Nerve10Middle","Nerve11Left","Nerve11Right","Nerve12Left","Nerve12Middle","Nerve12Right","Nerve13Left","Nerve13Right");
path1 = "J:\\Research\\Histology-Imaging\\Images\\Working\\Fibers\\";
path2 = "J:\\Research\\Histology-Imaging\\Images\\Working\\Myelin\\";
path3 = "J:\\Research\\Histology-Imaging\\Images\\Working\\";

for(j=0; j<lengthOf(names); j++) {
	open(path3+names[j]+".jpg");
	run("8-bit");
	run("Enhance Local Contrast (CLAHE)", "blocksize=100 histogram=256 maximum=5 mask=*None*");	
	save(path1+names[j]+".jpg");
}

//Creating the files that contain the axons
for(j=0; j<lengthOf(names); j++) {
	open(path1+names[j]+".jpg");
	rename("Nerve");
	run("Set Scale...", "distance=12.88 known=1 unit=um global");
	run("Set Measurements...", "area perimeter bounding fit shape feret's stack display redirect=None decimal=3");	
	thres = newArray(115,130,140,150,160,170,180,190,200,215);	
	
	for(it=0; it<10; it++) {	
		selectWindow("Nerve");
		run("Duplicate...", " ");	
		setThreshold(thres[it],255);
		setOption("BlackBackground", true);
		run("Convert to Mask");
		rename("Thres"+it);
		if(it<3) {
			selectWindow("Thres"+it);
			run("Duplicate...", " ");
			rename("S");
			run("Invert");
			run("BinaryFilterReconstruct ", "erosions=5 white");
			floodFill(0, 0);
			run("Invert");
			run("BinaryFilterReconstruct ", "erosions=3 white");
			run("Analyze Particles...", "size=.3-4 show=Masks clear");
			run("Invert LUT");
			rename("Small");
			close("S");
			selectWindow("Thres"+it);
			run("Invert");
			run("BinaryFilterReconstruct ", "erosions=8 white");
			floodFill(0, 0);
			run("Invert");
			run("BinaryFilterReconstruct ", "erosions=5 white");
			run("Analyze Particles...", "size=4.01-180 show=Masks clear");
			run("Invert LUT");
			rename("Axon"+it);
			close("Thres"+it);
			imageCalculator("Add", "Axon"+it, "Small");
			close("Small");
		} else{
			run("Invert");
			run("BinaryFilterReconstruct ", "erosions=8 white");
			floodFill(0, 0);
			run("Invert");
			run("BinaryFilterReconstruct ", "erosions=5 white");
			run("Analyze Particles...", "size=.3-180 show=Masks clear");
			run("Invert LUT");
			rename("Axon"+it);
			close("Thres"+it);
		}
	}
	
	for(it=0; it<10; it++) {
		if(it!=5) {
			imageCalculator("Difference create", "Axon5", "Axon"+it);
			rename("Dif"+it);
			run("Invert");
			close("Axon"+it);
		}
	}
	
	for(it=0; it<10; it++) {
		if(it==5) {
			selectWindow("Axon5");
			run("Duplicate...", " ");
			rename("A5");
			run("Analyze Particles...", "size=70-250 show=Nothing clear record");
			a = nResults;
			if(a>0) {
				run("Classify Particles", "class[1]=Solidity operator[1]=< value[1]=.7 class[2]=-empty- operator[2]=-empty- value[2]=0.0000 class[3]=-empty- operator[3]=-empty- value[3]=0.0000 class[4]=-empty- operator[4]=-empty- value[4]=0.0000 combine=[AND (match all)] output=[Keep members] white");
				close("A5");
				selectWindow("Subset");
				rename("A5");
				imageCalculator("Subtract", "Axon5", "A5");
			}
			close("A5");
			selectWindow("Axon5");
			run("Duplicate...", " ");
			rename("P"+it);
			run("Ultimate Points");
			setThreshold(1,255);
			setOption("BlackBackground", true);
			run("Convert to Mask");
		} else {
			selectWindow("Dif"+it);
			run("Analyze Particles...", "size=.3-180 show=Nothing clear record");
			a = nResults;
			if(a>0) {
				run("Classify Particles", "class[1]=Solidity operator[1]=>= value[1]=.85 class[2]=-empty- operator[2]=-empty- value[2]=0.0000 class[3]=-empty- operator[3]=-empty- value[3]=0.0000 class[4]=-empty- operator[4]=-empty- value[4]=0.0000 combine=[AND (match all)] output=[Keep members] white");		
				rename("Axon"+it);
				close("Dif"+it);
				if(it<5) {
					selectWindow("Axon"+it);
					run("Duplicate...", " ");
					rename("P"+it);
					run("Ultimate Points");
					setThreshold(1,255);
					setOption("BlackBackground", true);
					run("Convert to Mask");
				}
			} else{
				selectWindow("Dif"+it);
				run("Fill Holes");
				run("Invert");
				rename("Axon"+it);
				if(it<5) {
					run("Duplicate...", " ");
					rename("P"+it);
				}
			}
		}
	}

	for(it=0; it<5; it++) {
		imageCalculator("AND", "P"+it, "P5");
		imageCalculator("Subtract", "P5", "P"+it);
		run("BinaryReconstruct ", "mask=[Axon"+it+"] seed=[P"+it+"] white");
		close("Axon"+it);
		if(it!=0) {
			imageCalculator("Add", "P0", "P"+it);
			close("P"+it);
		}
	}

	for(it=6; it<10; it++) {
		imageCalculator("And create", "Axon"+it, "Axon5");
		imageCalculator("Subtract", "Axon"+it, "Result of Axon"+it);
		close("Result of Axon"+it);
		if(it!=6) {
			imageCalculator("Add", "Axon6", "Axon"+it);
			close("Axon"+it);
		}
	}

	run("BinaryReconstruct ", "mask=[Axon5] seed=[P5] white");
	close("Axon5");
	imageCalculator("Add", "P5", "Axon6");
	close("Axon6");
	imageCalculator("AND create", "P5", "P0");
	rename("A1");
	run("Analyze Particles...", "size=.3-180 show=Nothing clear");
	a = nResults;
	z = 0;
	if(a>0) {
		run("BinaryReconstruct ", "mask=[P0] seed=[A1] white");
		imageCalculator("Subtract", "P0", "A1");
		selectWindow("A1");
		saveAs("Jpeg", path2+names[j]+"A1.jpg");
		z = 1;
	}
	imageCalculator("Add", "P5", "P0");
	selectWindow("P5");
	run("Analyze Particles...", "size=.3-180 show=Nothing clear record");
	run("Classify Particles", "class[1]=Solidity operator[1]=>= value[1]=.7 class[2]=-empty- operator[2]=-empty- value[2]=0.0000 class[3]=-empty- operator[3]=-empty- value[3]=0.0000 class[4]=-empty- operator[4]=-empty- value[4]=0.0000 combine=[AND (match all)] output=[Grey non-members] white");		
	rename("Axon0");
	run("Duplicate...", " ");	
	setThreshold(120,130);
	setOption("BlackBackground", true);
	run("Convert to Mask");
	rename("Axon1");
	saveAs("Jpeg", path1+names[j]+"Axon1.jpg");
	selectWindow("Axon0");
	setThreshold(250,255);
	setOption("BlackBackground", true);
	run("Convert to Mask");
	saveAs("Jpeg", path1+names[j]+"Axon0.jpg");
	run("Close All");
}


//Analyzing the axons for myelin
for(h=0; h<lengthOf(names); h++) {

	if(z==1) {
		k = 3;
	} else {
		k = 2;
	}
	
	for(g=0; g<k; g++) {

		open(path1+names[h]+".jpg");
		rename("Nerve");
		run("Set Scale...", "distance=12.88 known=1 unit=um global");
		run("Set Measurements...", "area perimeter bounding fit shape feret's stack display redirect=None decimal=3");
		
		if(g!=2) {
			open(path1+names[h]+"Axon"+g+".jpg");
			selectWindow(names[h]+"Axon"+g+".jpg");
			rename("Axons");
			run("Make Binary");
		} else {
			open(path2+names[h]+"A1.jpg");
			selectWindow(names[h]+"A1.jpg");
			rename("Axons");
			run("Make Binary");		
		}
		
		run("Analyze Particles...", "size=.3-1.5 show=Masks clear");
		rename("Groups"+0);
		selectWindow("Axons");
		run("Analyze Particles...", "size=1.51-3 show=Masks clear");
		rename("Groups"+1);
		selectWindow("Axons");
		run("Analyze Particles...", "size=3.01-5 show=Masks clear");
		rename("Groups"+2);
		selectWindow("Axons");
		run("Analyze Particles...", "size=5.01-7.50 show=Masks clear");
		rename("Groups"+3);
		selectWindow("Axons");	
		run("Analyze Particles...", "size=7.51-10 show=Masks clear");
		rename("Groups"+4);
	
		for(i=0; i<5; i++) {
			selectWindow("Groups"+i);
			run("Invert LUT");
			run("Duplicate...", " ");
			rename("Dilate"+i);
			if(i<2) {
				run("Options...", "iterations=5 count=2 black do=Dilate");
				imageCalculator("Subtract", "Dilate"+i, "Groups"+i);
				imageCalculator("AND", "Dilate"+i, "Nerve");
				setThreshold(1,145);
				setOption("BlackBackground", true);
				run("Convert to Mask");
				run("Fill Holes");
				run("Options...", "iterations=3 count=3 black do=Open");			
			} else {
				run("Options...", "iterations=15 count=2 black do=Dilate");
				imageCalculator("Subtract", "Dilate"+i, "Groups"+i);
				imageCalculator("AND", "Dilate"+i, "Nerve");
				setThreshold(1,145);
				setOption("BlackBackground", true);
				run("Convert to Mask");
				run("Fill Holes");
				run("Options...", "iterations=5 count=3 black do=Open");		
			}
			close("Groups"+i);	
		}
	
		for(i=0; i<5; i++) {
			selectWindow("Dilate"+i);
			if(i==0) {
				run("Analyze Particles...", "size=.6-8 show=Nothing clear record");
			}
			if(i==1) {
				run("Analyze Particles...", "size=3-12 show=Nothing clear record");
			}
			if(i==2) {
				run("Analyze Particles...", "size=6-20 show=Nothing clear record");
			}
			if(i==3) {
				run("Analyze Particles...", "size=10-25 show=Nothing clear record");
			}
			if(i==4) {
				run("Analyze Particles...", "size=15-36 show=Nothing clear record");
			}
			a = nResults;
			if(a>0) {
				run("Classify Particles", "class[1]=Solidity operator[1]=> value[1]=.85 class[2]=-empty- operator[2]=-empty- value[2]=0.0000 class[3]=-empty- operator[3]=-empty- value[3]=0.0000 class[4]=-empty- operator[4]=-empty- value[4]=0.0000 combine=[AND (match all)] output=[Keep members] white");
				selectWindow("Subset");
				rename("Myelin"+i);
			} else {
				selectWindow("Dilate"+i);
				rename("Myelin"+i);
			}
			if(a>0) {
				imageCalculator("Subtract", "Dilate"+i, "Myelin"+i);
				if(i>1) {
					selectWindow("Dilate"+i);
					run("Watershed");
					if(i==0) {
						run("Analyze Particles...", "size=.6-8 show=Nothing clear record");
					}
					if(i==1) {
						run("Analyze Particles...", "size=3-12 show=Nothing clear record");
					}
					if(i==2) {
						run("Analyze Particles...", "size=6-20 show=Nothing clear record");
					}
					if(i==3) {
						run("Analyze Particles...", "size=10-25 show=Nothing clear record");
					}
					if(i==4) {
						run("Analyze Particles...", "size=15-36 show=Nothing clear record");
					}
					b = nResults;
					if(b>0) {
						run("Classify Particles", "class[1]=Solidity operator[1]=> value[1]=.85 class[2]=-empty- operator[2]=-empty- value[2]=0.0000 class[3]=-empty- operator[3]=-empty- value[3]=0.0000 class[4]=-empty- operator[4]=-empty- value[4]=0.0000 combine=[AND (match all)] output=[Keep members] white");
						selectWindow("Subset");
						rename("Myelin2"+i);
						close("Dilate"+i);
						imageCalculator("Add", "Myelin"+i, "Myelin2"+i);
						close("Myelin2"+i);
					} 
				}
			}
		}
		
		for(i=0; i<5; i++) {
			selectWindow("Myelin"+i);
			saveAs("Jpeg", path2+names[h]+"MyelinA"+g+"-"+i+".jpg");
			close("Myelin"+i);
			close("MyelinA"+g+"-"+i);
			close(names[h]+"MyelinA"+g+"-"+i+".jpg");	
		}
	
//END OF 0-4
	
		selectWindow("Axons");	
		run("Analyze Particles...", "size=10.01-12.50 show=Masks clear");
		rename("Groups"+5);	
		selectWindow("Axons");
		run("Analyze Particles...", "size=12.51-15 show=Masks clear");
		rename("Groups"+6);
		selectWindow("Axons");
		run("Analyze Particles...", "size=15.01-17 show=Masks clear");
		rename("Groups"+7);
		selectWindow("Axons");
		run("Analyze Particles...", "size=17.01-19.5 show=Masks clear");
		rename("Groups"+8);
		selectWindow("Axons");
		run("Analyze Particles...", "size=19.51-22 show=Masks clear");
		rename("Groups"+9);
	
		
		for(i=5; i<10; i++) {
			selectWindow("Groups"+i);
			run("Invert LUT");
			run("Duplicate...", " ");
			rename("Dilate"+i);
			if(i<7) {
				run("Options...", "iterations=18 count=2 black do=Dilate");
				imageCalculator("Subtract", "Dilate"+i, "Groups"+i);
				imageCalculator("AND", "Dilate"+i, "Nerve");
				setThreshold(1,145);
				setOption("BlackBackground", true);
				run("Convert to Mask");
				run("Fill Holes");
				run("Options...", "iterations=3 count=3 black do=Open");			
			} else {
				run("Options...", "iterations=23 count=2 black do=Dilate");
				imageCalculator("Subtract", "Dilate"+i, "Groups"+i);
				imageCalculator("AND", "Dilate"+i, "Nerve");
				setThreshold(1,145);
				setOption("BlackBackground", true);
				run("Convert to Mask");
				run("Fill Holes");
				run("Options...", "iterations=25 count=3 black do=Open");		
			}
			close("Groups"+i);	
		}
		
		for(i=5; i<10; i++) {
			selectWindow("Dilate"+i);
			if(i==5) {
				run("Analyze Particles...", "size=20-42 show=Nothing clear record");
			}
			if(i==6) {
				run("Analyze Particles...", "size=25-52 show=Nothing clear record");
			}
			if(i==7) {
				run("Analyze Particles...", "size=30-62 show=Nothing clear record");
			}
			if(i==8) {
				run("Analyze Particles...", "size=35-70 show=Nothing clear record");
			}
			if(i==9) {
				run("Analyze Particles...", "size=40-75 show=Nothing clear record");
			}
			a = nResults;
			if(a>0) {
				run("Classify Particles", "class[1]=Solidity operator[1]=> value[1]=.85 class[2]=-empty- operator[2]=-empty- value[2]=0.0000 class[3]=-empty- operator[3]=-empty- value[3]=0.0000 class[4]=-empty- operator[4]=-empty- value[4]=0.0000 combine=[AND (match all)] output=[Keep members] white");
				selectWindow("Subset");
				rename("Myelin"+i);
			} else {
				selectWindow("Dilate"+i);
				rename("Myelin"+i);
			}
			if(a>0) {
				imageCalculator("Subtract", "Dilate"+i, "Myelin"+i);
				selectWindow("Dilate"+i);
				run("Watershed");
				if(i==5) {
					run("Analyze Particles...", "size=20-42 show=Nothing clear record");
				}
				if(i==6) {
					run("Analyze Particles...", "size=25-52 show=Nothing clear record");
				}
				if(i==7) {
					run("Analyze Particles...", "size=30-62 show=Nothing clear record");
				}
				if(i==8) {
					run("Analyze Particles...", "size=35-70 show=Nothing clear record");
				}
				if(i==9) {
					run("Analyze Particles...", "size=40-75 show=Nothing clear record");
				}
				b = nResults;
				if(b>0) {
					run("Classify Particles", "class[1]=Solidity operator[1]=> value[1]=.85 class[2]=-empty- operator[2]=-empty- value[2]=0.0000 class[3]=-empty- operator[3]=-empty- value[3]=0.0000 class[4]=-empty- operator[4]=-empty- value[4]=0.0000 combine=[AND (match all)] output=[Keep members] white");
					selectWindow("Subset");
					rename("Myelin2"+i);
					close("Dilate"+i);
					imageCalculator("Add", "Myelin"+i, "Myelin2"+i);
					close("Myelin2"+i);
				} 
			}
		}
		
		for(i=5; i<10; i++) {
			selectWindow("Myelin"+i);
			saveAs("Jpeg", path2+names[h]+"MyelinA"+g+"-"+i+".jpg");
			close("Myelin"+i);
			close("MyelinA"+g+"-"+i);
			close(names[h]+"MyelinA"+g+"-"+i+".jpg");
		}
		
//END OF 5-9	
		
		selectWindow("Axons");
		run("Analyze Particles...", "size=22.01-24.5 show=Masks clear");
		rename("Groups"+10);
		selectWindow("Axons");
		run("Analyze Particles...", "size=24.51-27 show=Masks clear");
		rename("Groups"+11);
		selectWindow("Axons");
		run("Analyze Particles...", "size=27.01-30 show=Masks clear");
		rename("Groups"+12);
		selectWindow("Axons");
		run("Analyze Particles...", "size=30.01-33 show=Masks clear");
		rename("Groups"+13);
		
		for(i=10; i<14; i++) {
			selectWindow("Groups"+i);
			run("Invert LUT");
			run("Duplicate...", " ");
			rename("Dilate"+i);
			if(i<11) {
				run("Options...", "iterations=25 count=2 black do=Dilate");
				imageCalculator("Subtract", "Dilate"+i, "Groups"+i);
				imageCalculator("AND", "Dilate"+i, "Nerve");
				setThreshold(1,145);
				setOption("BlackBackground", true);
				run("Convert to Mask");
				run("Fill Holes");
				run("Options...", "iterations=3 count=3 black do=Open");			
			} else {
				run("Options...", "iterations=30 count=2 black do=Dilate");
				imageCalculator("Subtract", "Dilate"+i, "Groups"+i);
				imageCalculator("AND", "Dilate"+i, "Nerve");
				setThreshold(1,145);
				setOption("BlackBackground", true);
				run("Convert to Mask");
				run("Fill Holes");
				run("Options...", "iterations=25 count=3 black do=Open");		
			}
			close("Groups"+i);	
		}
	
		for(i=10; i<14; i++) {
			selectWindow("Dilate"+i);
			if(i==10) {
				run("Analyze Particles...", "size=45-80 show=Nothing clear record");
			}
			if(i==11) {
				run("Analyze Particles...", "size=50-85 show=Nothing clear record");
			}
			if(i==12) {
				run("Analyze Particles...", "size=55-95 show=Nothing clear record");
			}
			if(i==13) {
				run("Analyze Particles...", "size=60-105 show=Nothing clear record");
			}
			a = nResults;
			if(a>0) {
				run("Classify Particles", "class[1]=Solidity operator[1]=> value[1]=.85 class[2]=-empty- operator[2]=-empty- value[2]=0.0000 class[3]=-empty- operator[3]=-empty- value[3]=0.0000 class[4]=-empty- operator[4]=-empty- value[4]=0.0000 combine=[AND (match all)] output=[Keep members] white");
				selectWindow("Subset");
				rename("Myelin"+i);
			} else {
				selectWindow("Dilate"+i);
				rename("Myelin"+i);
			}
			if(a>0) {
				imageCalculator("Subtract", "Dilate"+i, "Myelin"+i);
				selectWindow("Dilate"+i);
				run("Watershed");
				if(i==10) {
					run("Analyze Particles...", "size=45-80 show=Nothing clear record");
				}
				if(i==11) {
					run("Analyze Particles...", "size=50-85 show=Nothing clear record");
				}
				if(i==12) {
					run("Analyze Particles...", "size=55-95 show=Nothing clear record");
				}
				if(i==13) {
					run("Analyze Particles...", "size=60-105 show=Nothing clear record");
				}
				b = nResults;
				if(b>0) {
					run("Classify Particles", "class[1]=Solidity operator[1]=> value[1]=.85 class[2]=-empty- operator[2]=-empty- value[2]=0.0000 class[3]=-empty- operator[3]=-empty- value[3]=0.0000 class[4]=-empty- operator[4]=-empty- value[4]=0.0000 combine=[AND (match all)] output=[Keep members] white");
					selectWindow("Subset");
					rename("Myelin2"+i);
					close("Dilate"+i);
					imageCalculator("Add", "Myelin"+i, "Myelin2"+i);
					close("Myelin2"+i);
				} 
			}
		}
		
		for(i=10; i<14; i++) {
			selectWindow("Myelin"+i);
			saveAs("Jpeg", path2+names[h]+"MyelinA"+g+"-"+i+".jpg");
			close("Myelin"+i);
			close("MyelinA"+g+"-"+i);
			close(names[h]+"MyelinA"+g+"-"+i+".jpg");
		}
		
//END OF 10-13	
		
		selectWindow("Axons");
		run("Analyze Particles...", "size=33.01-36 show=Masks clear");
		rename("Groups"+14);
		selectWindow("Axons");
		run("Analyze Particles...", "size=36.01-39 show=Masks clear");
		rename("Groups"+15);
		selectWindow("Axons");
		run("Analyze Particles...", "size=39.01-42 show=Masks clear");
		rename("Groups"+16);
		selectWindow("Axons");
		run("Analyze Particles...", "size=42.01-45 show=Masks clear");
		rename("Groups"+17);
		
		for(i=14; i<18; i++) {
			selectWindow("Groups"+i);
			run("Invert LUT");
			run("Duplicate...", " ");
			rename("Dilate"+i);
			if(i<16) {
				run("Options...", "iterations=30 count=2 black do=Dilate");
				imageCalculator("Subtract", "Dilate"+i, "Groups"+i);
				imageCalculator("AND", "Dilate"+i, "Nerve");
				setThreshold(1,145);
				setOption("BlackBackground", true);
				run("Convert to Mask");
				run("Fill Holes");
				run("Options...", "iterations=3 count=3 black do=Open");			
			} else {
				run("Options...", "iterations=35 count=2 black do=Dilate");
				imageCalculator("Subtract", "Dilate"+i, "Groups"+i);
				imageCalculator("AND", "Dilate"+i, "Nerve");
				setThreshold(1,145);
				setOption("BlackBackground", true);
				run("Convert to Mask");
				run("Fill Holes");
				run("Options...", "iterations=25 count=3 black do=Open");		
			}
			close("Groups"+i);	
		}
		
		for(i=14; i<18; i++) {
			selectWindow("Dilate"+i);
			if(i==14) {
				run("Analyze Particles...", "size=65-115 show=Nothing clear record");
			}
			if(i==15) {
				run("Analyze Particles...", "size=70-135 show=Nothing clear record");
			}
			if(i==16) {
				run("Analyze Particles...", "size=80-145 show=Nothing clear record");
			}
			if(i==17) {
				run("Analyze Particles...", "size=85-155 show=Nothing clear record");
			}
			a = nResults;
			if(a>0) {
				run("Classify Particles", "class[1]=Solidity operator[1]=> value[1]=.85 class[2]=-empty- operator[2]=-empty- value[2]=0.0000 class[3]=-empty- operator[3]=-empty- value[3]=0.0000 class[4]=-empty- operator[4]=-empty- value[4]=0.0000 combine=[AND (match all)] output=[Keep members] white");
				selectWindow("Subset");
				rename("Myelin"+i);
			} else {
				selectWindow("Dilate"+i);
				rename("Myelin"+i);
			}
			if(a>0) {
				imageCalculator("Subtract", "Dilate"+i, "Myelin"+i);
				selectWindow("Dilate"+i);
				run("Watershed");
				if(i==14) {
					run("Analyze Particles...", "size=65-115 show=Nothing clear record");
				}
				if(i==15) {
					run("Analyze Particles...", "size=70-135 show=Nothing clear record");
				}
				if(i==16) {
					run("Analyze Particles...", "size=80-145 show=Nothing clear record");
				}
				if(i==17) {
					run("Analyze Particles...", "size=85-155 show=Nothing clear record");
				}
				b = nResults;
				if(b>0) {
					run("Classify Particles", "class[1]=Solidity operator[1]=> value[1]=.85 class[2]=-empty- operator[2]=-empty- value[2]=0.0000 class[3]=-empty- operator[3]=-empty- value[3]=0.0000 class[4]=-empty- operator[4]=-empty- value[4]=0.0000 combine=[AND (match all)] output=[Keep members] white");
					selectWindow("Subset");
					rename("Myelin2"+i);
					close("Dilate"+i);
					imageCalculator("Add", "Myelin"+i, "Myelin2"+i);
					close("Myelin2"+i);
				} 
			}
		}
	
		for(i=14; i<18; i++) {
			selectWindow("Myelin"+i);
			saveAs("Jpeg", path2+names[h]+"MyelinA"+g+"-"+i+".jpg");
			close("Myelin"+i);
			close("MyelinA"+g+"-"+i);
			close(names[h]+"MyelinA"+g+"-"+i+".jpg");
		}
		
//END OF 14-17	
		
		selectWindow("Axons");
		run("Analyze Particles...", "size=45.01-50 show=Masks clear");
		rename("Groups"+18);
		selectWindow("Axons");
		run("Analyze Particles...", "size=50.01-55 show=Masks clear");
		rename("Groups"+19);
		selectWindow("Axons");
		run("Analyze Particles...", "size=55.01-60 show=Masks clear");
		rename("Groups"+20);
		selectWindow("Axons");
		run("Analyze Particles...", "size=60.01-70 show=Masks clear");
		rename("Groups"+21);
		selectWindow("Axons");
		run("Analyze Particles...", "size=70.01-80 show=Masks clear");
		rename("Groups"+22);
		selectWindow("Axons");
		run("Analyze Particles...", "size=80.01-250 show=Masks clear");
		rename("Groups"+23);
	
		for(i=18; i<24; i++) {
			selectWindow("Groups"+i);
			run("Invert LUT");
			run("Duplicate...", " ");
			rename("Dilate"+i);
			if(i<22) {
				run("Options...", "iterations=40 count=2 black do=Dilate");
				imageCalculator("Subtract", "Dilate"+i, "Groups"+i);
				imageCalculator("AND", "Dilate"+i, "Nerve");
				setThreshold(1,145);
				setOption("BlackBackground", true);
				run("Convert to Mask");
				run("Fill Holes");
				run("Options...", "iterations=3 count=3 black do=Open");			
			} else {
				run("Options...", "iterations=45 count=2 black do=Dilate");
				imageCalculator("Subtract", "Dilate"+i, "Groups"+i);
				imageCalculator("AND", "Dilate"+i, "Nerve");
				setThreshold(1,145);
				setOption("BlackBackground", true);
				run("Convert to Mask");
				run("Fill Holes");
				run("Options...", "iterations=25 count=3 black do=Open");		
			}			
			close("Groups"+i);	
		}
		
		for(i=18; i<24; i++) {
			selectWindow("Dilate"+i);
			if(i==18) {
				run("Analyze Particles...", "size=90-165 show=Nothing clear record");
			}
			if(i==19) {
				run("Analyze Particles...", "size=100-180 show=Nothing clear record");
			}
			if(i==20) {
				run("Analyze Particles...", "size=110-200 show=Nothing clear record");
			}
			if(i==21) {
				run("Analyze Particles...", "size=120-240 show=Nothing clear record");
			}
			if(i==22) {
				run("Analyze Particles...", "size=140-290 show=Nothing clear record");
			}
			if(i==23) {
				run("Analyze Particles...", "size=160-350 show=Nothing clear record");
			}
			a = nResults;
			if(a>0) {
				run("Classify Particles", "class[1]=Solidity operator[1]=> value[1]=.85 class[2]=-empty- operator[2]=-empty- value[2]=0.0000 class[3]=-empty- operator[3]=-empty- value[3]=0.0000 class[4]=-empty- operator[4]=-empty- value[4]=0.0000 combine=[AND (match all)] output=[Keep members] white");
				selectWindow("Subset");
				rename("Myelin"+i);
			} else {
				selectWindow("Dilate"+i);
				rename("Myelin"+i);
			}
			if(a>0) {
				imageCalculator("Subtract", "Dilate"+i, "Myelin"+i);
				selectWindow("Dilate"+i);
				run("Watershed");
				if(i==18) {
					run("Analyze Particles...", "size=90-165 show=Nothing clear record");
				}
				if(i==19) {
					run("Analyze Particles...", "size=100-180 show=Nothing clear record");
				}
				if(i==20) {
					run("Analyze Particles...", "size=110-200 show=Nothing clear record");
				}
				if(i==21) {
					run("Analyze Particles...", "size=120-240 show=Nothing clear record");
				}
				if(i==22) {
					run("Analyze Particles...", "size=140-280 show=Nothing clear record");
				}
				if(i==23) {
					run("Analyze Particles...", "size=160-350 show=Nothing clear record");
				}
				b = nResults;
				if(b>0) {
					run("Classify Particles", "class[1]=Solidity operator[1]=> value[1]=.85 class[2]=-empty- operator[2]=-empty- value[2]=0.0000 class[3]=-empty- operator[3]=-empty- value[3]=0.0000 class[4]=-empty- operator[4]=-empty- value[4]=0.0000 combine=[AND (match all)] output=[Keep members] white");
					selectWindow("Subset");
					rename("Myelin2"+i);
					close("Dilate"+i);
					imageCalculator("Add", "Myelin"+i, "Myelin2"+i);
					close("Myelin2"+i);
				} 
			}
		}
	
		for(i=18; i<24; i++) {
			selectWindow("Myelin"+i);
			saveAs("Jpeg", path2+names[h]+"MyelinA"+g+"-"+i+".jpg");
			close("Myelin"+i);
			close("MyelinA"+g+"-"+i);
			close(names[h]+"MyelinA"+g+"-"+i+".jpg");
		}
		
	//END OF 18-23
	run("Close All");
	}	
}

//Merging Files
for(h=0; h<lengthOf(names); h++) {

	if(z==1) {
		k = 3;
	} else {
		k = 2;
	}
	
	for(g=0; g<k; g++) {
		
		for(i=0; i<24; i++) {
			open(path2+names[h]+"MyelinA"+g+"-"+i+".jpg");
			rename("MyelinA"+g+"-"+i);
			run("Make Binary");
		}
		
		for(i=1; i<24; i++) {
			imageCalculator("AND create", "MyelinA"+g+"-0", "MyelinA"+g+"-"+i);
			rename("M"+g+"-"+i);
			run("Analyze Particles...", "size=.1-350 show=Nothing clear");
			a = nResults;
			if(a>0) {
				run("BinaryReconstruct ", "mask=[MyelinA"+g+"-"+i+"] seed=[M"+g+"-"+i+"] white");
				imageCalculator("Subtract", "MyelinA"+g+"-"+i, "M"+g+"-"+i);
			}
			imageCalculator("Add", "MyelinA"+g+"-0", "MyelinA"+g+"-"+i);
			close("MyelinA"+g+"-"+i);		
		}

		for(i=2; i<24; i++) {
			imageCalculator("AND create", "M"+g+"-1", "M"+g+"-"+i);
			rename("MM"+g+"-"+i);
			run("Analyze Particles...", "size=.1-350 show=Nothing clear");
			a = nResults;
			if(a>0) {
				run("BinaryReconstruct ", "mask=[M"+g+"-"+i+"] seed=[MM"+g+"-"+i+"] white");
				imageCalculator("Subtract", "M"+g+"-"+i, "MM"+g+"-"+i);
			}
			imageCalculator("Add", "M"+g+"-1", "M"+g+"-"+i);
			close("M"+g+"-"+i);		
		}

		for(i=3; i<24; i++) {
			imageCalculator("Add", "MM"+g+"-2", "MM"+g+"-"+i);	
			close("MM"+g+"-"+i);
		}

	}

	for(i=0; i<k; i++) {
		selectWindow("MyelinA"+i+"-0");
		rename("Myelin"+i);
	}

	for(i=0; i<k; i++) {
		selectWindow("M"+i+"-1");
		rename("M"+i);
	}
			
	for(i=0; i<k; i++) {
		selectWindow("MM"+i+"-2");
		rename("MM"+i);
	}

	for(i=1; i<k; i++) {
		imageCalculator("AND create", "Myelin0", "Myelin"+i);
		rename("Mye"+i);
		run("Analyze Particles...", "size=.1-350 show=Nothing clear");
		a = nResults;
		b = 0;
		if(a>0) {
			run("BinaryReconstruct ", "mask=[Myelin"+i+"] seed=[Mye"+i+"] white");
			imageCalculator("Subtract", "Myelin"+i, "Mye"+i);
			if(i==2) {
				imageCalculator("Add", "Mye1", "Mye2")
				selectWindow("Mye1");
				run("Analyze Particles...", "size=.1-350 show=Nothing clear");
				a = nResults;
				if(a<1) {
					close("Mye1");
					b = 1;					
				}
			}
		}
		imageCalculator("Add", "Myelin0", "Myelin"+i);
		close("Myelin"+i);		
	}

	for(i=1; i<k; i++) {
		imageCalculator("AND create", "M0", "M"+i);
		rename("MMye"+i);
		run("Analyze Particles...", "size=.1-350 show=Nothing clear");
		a = nResults;
		c = 0;
		if(a>0) {
			run("BinaryReconstruct ", "mask=[M"+i+"] seed=[MMye"+i+"] white");
			imageCalculator("Subtract", "M"+i, "MMye"+i);
			if(i==2) {
				imageCalculator("Add", "MMye1", "MMye2")
				selectWindow("MMye1");
				run("Analyze Particles...", "size=.1-350 show=Nothing clear");
				a = nResults;
				if(a<1) {
					close("MMye1");
					c = 1;					
				}
			}
		}
		imageCalculator("Add", "M0", "M"+i);
		close("M"+i);		
	}

	for(i=1; i<k; i++) {
		imageCalculator("Add", "MM0", "MM"+i);
		close("MM"+i);
	}

	if(b==1) {
		imageCalculator("Add", "MM0", "Mye1");
		close("Mye1");
	}
	if(c==1) {
		imageCalculator("Add", "MM0", "MMye1");
		close("MMye1");
	}
		
	imageCalculator("AND create", "Myelin0", "M0");
	rename("Mye1");
	run("Analyze Particles...", "size=.1-350 show=Nothing clear");
	a = nResults;
	d = 0;
	if(a>0) {
		run("BinaryReconstruct ", "mask=[M0] seed=[Mye1] white");
		imageCalculator("Subtract", "M0", "Mye1");
		d = 1;
	}
	imageCalculator("Add", "Myelin0", "M0");
	close("M0");

	imageCalculator("AND create", "Myelin0", "MM0");
	rename("Mye3");
	run("Analyze Particles...", "size=.1-350 show=Nothing clear");
	a = nResults;
	e = 0;
	if(a>0) {
		run("BinaryReconstruct ", "mask=[MM0] seed=[Mye3] white");
		imageCalculator("Subtract", "MM0", "Mye3");
		e = 1;
	}
	imageCalculator("Add", "Myelin0", "MM0");
	close("MM0");
	selectWindow("Myelin0");
	saveAs("Jpeg", path1+names[h]+"Myelin.jpg");

	if(d==1) {
		if(e==1) {
			imageCalculator("AND create", "Mye1", "Mye3");
			rename("Mye2");
			run("Analyze Particles...", "size=.1-350 show=Nothing clear");
			a = nResults;
			y = 0;
			if(a>0) {
				run("BinaryReconstruct ", "mask=[Mye3] seed=[Mye2] white");
				imageCalculator("Subtract", "Mye3", "Mye2");
				selectWindow("Mye2");
				saveAs("Jpeg", path1+names[h]+"Mye2.jpg");
				y = 1;
			}
			imageCalculator("Add", "Mye1", "Mye3");
			close("Mye3");
		}
		
		selectWindow("Mye1");
		saveAs("Jpeg", path1+names[h]+"Mye1.jpg");
	}

	run("Close All");
}


names = newArray("Nerve4Left","Nerve4Right","Nerve11Left","Nerve11Right","Nerve12Right");
path1 = "J:\\Research\\Histology-Imaging\\Images\\Working\\Fibers\\";
path3 = "J:\\Research\\Histology-Imaging\\Images\\Working\\";
//Outlining Myelin on original pictures for results
for(j=0; j<lengthOf(names); j++) {
	if(j<27) {
		y = 1;
	} else {
		y = 0;
	}
	open(path3+names[j]+".jpg");
	rename("Original");
	open(path1+names[j]+"Myelin.jpg");
	rename("Myelin");
	run("Make Binary");
	run("Duplicate...", " ");
	rename("Myelin Outline");
	run("Outline");
	run("Dilate");
	open(path1+names[j]+"Mye1.jpg");
	rename("Mye1");
	run("Make Binary");
	run("Duplicate...", " ");
	rename("Mye1 Outline");
	run("Outline");
	run("Dilate");
	if(y==1) {
		open(path1+names[j]+"Mye2.jpg");
		rename("Mye2");
		run("Make Binary");
		run("Duplicate...", " ");
		rename("Mye2 Outline");
		run("Outline");
		run("Dilate");
		imageCalculator("Subtract create", "Original", "Mye2 Outline");
		rename("Mye2 Out");
		close("Mye2 Outline");
		selectWindow("Mye2 Out");
		saveAs("Jpeg", path1+names[j]+"Mye2 Outline.jpg");
	}

	imageCalculator("Subtract create", "Original", "Myelin Outline");
	rename("Myelin Out");
	close("Myelin Outline");
	imageCalculator("Subtract create", "Original", "Mye1 Outline");
	rename("Mye1 Out");
	close("Mye1 Outline");
	selectWindow("Myelin Out");
	saveAs("Jpeg", path1+names[j]+"Myelin Outline.jpg");
	selectWindow("Mye1 Out");
	saveAs("Jpeg", path1+names[j]+"Mye1 Outline.jpg");
	run("Close All");	
}
