names = newArray("Nerve 1 bottom half","Nerve 1 top half","Nerve 2 grey","Nerve 3 bottom half copy","Nerve 3 bottom half","Nerve 3 top half","Nerve 4 bottom","Nerve 4 bottom middle","Nerve 4 top middle","Nerve 4 top","Nerve 5 left","Nerve 5 middle","Nerve 5 right","Nerve 6 bottom","Nerve 6 top","Nerve 7 bottom","Nerve 7 middle","Nerve 7 top","Nerve 8 bottom","Nerve 8 middle","Nerve 8 top","Nerve 9 bottom","Nerve 9 top","Nerve 10 bottom","Nerve 10 top","Nerve 11 bottom","Nerve 11 top","Nerve 12 bottom","Nerve 12 middle","Nerve 12 top","Nerve 13 bottom","Nerve 13 middle","Nerve 13 top","Nerve 14 bottom","Nerve 14 middle","Nerve 14 top","Nerve 15 bottom","Nerve 15 middle","Nerve 15 top","Nerve 16 left","Nerve 16 left copy","Nerve 16 middle","Nerve 16 right","Nerve 16 right copy","Nerve 17 bottom","Nerve 17 bottom copy","Nerve 17 middle copy","Nerve 17 middle","Nerve 17 top","Nerve 18 middle","Nerve 18 top","Nerve 18 bottom","Nerve 19 left","Nerve 19 middle","Nerve 19 middle lt copy","Nerve 19 right","Nerve 19 right lt copy","Nerve 20 middle","Nerve 20 top","Nerve 20 bottom","Nerve 21 middle","Nerve 21 top","Nerve 21 bottom","Nerve 22 left left","Nerve 22 left","Nerve 22 middle","Nerve 22 right","Nerve 24 bottom","Nerve 24 middle","Nerve 24 top","Nerve 25 bottom", "Nerve 25 middle", "Nerve 25 top","Nerve 26 bottom","Nerve 26 bottom middle","Nerve 26 top middle","Nerve 26 top");
path1 = "J:\\Research\\Histology-Imaging\\Images\\Fibers\\";
for(a=0; a<lengthOf(names); a++) {
	open(path1+names[a]+"Axon0.jpg");
	run("Make Binary");
	run("Analyze Particles...", "clear");
	n = nResults;
	print(names[a]+"-"+n);
	close(names[a]+"Axon0.jpg");
}