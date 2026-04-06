// 1st open images as hyperstacks

// Save name
nameImg = getTitle();

//run("Stack to Hyperstack...", "order=xyctz channels=3 slices=1 frames=1 display=Color");

// Update dimensions
getDimensions(width, height, channels, slices, frames);

// Case: 1 channeñ, 3 slices wrongly codified
if (channels == 1 && slices == 3) {
    run("Stack to Hyperstack...", "order=xyzct channels=3 slices=1 frames=1 display=Composite");
}

// if Z slices, max project
getDimensions(width, height, channels, slices, frames);
if (slices > 1) {
    run("Z Project...", "projection=[Max Intensity]");
    prefixN = "MAX_";
} else {
    prefixN = "";
}

run("Split Channels");

// Function to scale 16-bit to 8-bit without saturation and enhance contrast
function enhanceAndConvert(imTitle) {
    selectImage(imTitle);
    run("Enhance Contrast", "saturated=1");
    run("8-bit");
}

// Apply to 3 channels
enhanceAndConvert("C3-" + prefixN + nameImg); // → R
enhanceAndConvert("C1-" + prefixN + nameImg); // → G
enhanceAndConvert("C2-" + prefixN + nameImg); // → B

// Combine in RGB with order C3-R, C1-G, C2-B
run("Merge Channels...", 
  "c1=[C3-" + prefixN + nameImg + "] " +
  "c2=[C1-" + prefixN + nameImg + "] " +
  "c3=[C2-" + prefixN + nameImg + "] " +
  "create ignore");

// Convert to RGB image (8-bit)
run("RGB Color");

days = "4";
name = "KIKO 777";
datePhoto = "100226";
saveAs("Tiff", "F:/Lab/MouseMusclePOGLUT1/data/rawImages/tifs/p"+days+"/"+ name +" p" +days+" "+datePhoto+".tif");
close();
close();