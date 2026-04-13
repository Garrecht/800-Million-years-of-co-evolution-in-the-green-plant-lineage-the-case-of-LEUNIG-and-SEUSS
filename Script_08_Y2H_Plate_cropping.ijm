//Dynamicaly define input and output folders and file type
#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output
#@ String (label = "File type", value = ".jpg") suffix

//Define pre- and suffixe for processed pictures
#@ String (label = "Output prefix", value = "") output_prefix
#@ String (label = "Output suffix", value = "") output_suffix

// function to scan folders/subfolders/files to find files with correct suffix
processFolder(input);
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix))
			processFile(input, output, list[i]);
	}
}

//Execution of macro for each file
function processFile(input, output, file) {


//image handeling and preparation
open(input + File.separator + file);
original_image = getTitle(); //Get image name for later selection
selectWindow(original_image);
	
//Manual cropping
makeRectangle(1032, 612, 3888, 2580);
waitForUser("Position selection EXACTLY on internal border (edge media - plate wall) of the plates!"); 
run("Crop");
// Sometimes, to reduce noise during the analysis, a slight blur is advised
//run("Gaussian Blur...", "sigma=1");

//Saving the image
filename_pure = output_prefix + File.nameWithoutExtension + output_suffix; //Adds pre- and suffixe to the file name
saving_prefix = output + File.separator + filename_pure; //Defines save location
saveAs("JPEG", saving_prefix);

//Clean-Up for next picture
run("Close All");

}
