function action(input, output, filename) {
	open(input + filename);
	run("Split Channels");
	setAutoThreshold("Default dark");
	//run("Threshold...");
	setThreshold(14000, 1000000000000000000000000000000.0000);
	setOption("BlackBackground", false);
	run("Convert to Mask");
	run("Analyze Particles...", "size=0.50-7.00 show=Overlay add");
	saveAs("Tiff", output + filename);
	close();
	//selectWindow("C1-" + filename);
	roiManager("Measure");

	saveAs("Results",  output + filename + ".csv");
	roiManager("reset");	
	run("Close"); 
	
		
		
		close();
}

input = "D:\\210701\\SIP\\";

output = "D:\\210701\\SIP\\data\\";

setBatchMode(true);
list = getFileList(input);
for (i = 0; i < list.length; i++)
		action(input, output, list[i]);
setBatchMode(false);