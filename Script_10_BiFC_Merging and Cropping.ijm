//Dynamicaly define input files, output folders and file type
#@ File (label = "DAPI file", style = "open") image_DAPI
#@ File (label = "YFP file", style = "open") image_YFP
#@ File (label = "Output directory", style = "directory") output
#@ String (label = "File type", value = ".tif") suffix

//Define pre- and suffixe for processed pictures
#@ String (label = "Merge suffix", value = "_merge") output_sufMerge
#@ String (label = "Crop suffix", value = "") output_sufCrop


//Open DAPI-Picture
open(image_DAPI);
SaveDAPI = File.nameWithoutExtension //Stores the Name of the file for saving later
Image_DAPI = getTitle(); //Allows for selection of the opened file
run("Duplicate...", "title=Duplicate");

//Split channels, keep only the blue one
run("Split Channels");
selectImage("Duplicate (green)");
close();
selectImage("Duplicate (red)");
close();


//Open YFP-Picture
open(image_YFP);
SaveYFP = File.nameWithoutExtension //Stores the Name of the file for saving later
Image_YFP = getTitle() //Allows for selection of the opened file
selectImage(Image_YFP);
run("Duplicate...", "title=Duplicate2");

//Split channels, keep only the green one
run("Split Channels");
selectImage("Duplicate2 (blue)");
close();
selectImage("Duplicate2 (red)");
close();


//Merging of blue DAPI and green YFP channels
run("Merge Channels...", "c2=[Duplicate2 (green)] c3=[Duplicate (blue)]");
selectImage("RGB");

//Add selection for cropping
run("Specify...", "width=750 height=750 x=1 y=1");
waitForUser("Position ROI, then hit OK"); 
//You need to select the area to crop manually. Select a well visible/overlayed area for the cropping


//Crop and save Merged image
roiManager("Add");
run("Crop");
saveAs("PNG", output + File.separator + File.nameWithoutExtension + output_sufMerge);

//Crop and save YFP and DAPI images
selectImage(Image_DAPI);
roiManager("Select", 0);
run("Crop");
saveAs("PNG", output + File.separator + SaveDAPI + output_sufCrop);
selectImage(Image_YFP);
roiManager("Select", 0);
run("Crop");
saveAs("PNG", output + File.separator + SaveYFP + output_sufCrop);


//Cleanup for next pixture
roiManager("Select", 0);
roiManager("Deselect");
roiManager("Delete");
run("Close All");
