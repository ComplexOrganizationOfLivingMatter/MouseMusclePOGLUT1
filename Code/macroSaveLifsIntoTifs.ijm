//nameImg = "Carpeta 170322 Para Anlizar WT.lif - WT 51 P20 369 Tiles_Merged";
nameImg = getTitle();
selectImage(nameImg);
run("Z Project...", "projection=[Max Intensity]");
run("Split Channels");
run("Merge Channels...", "c1=[C3-MAX_" + nameImg + "] c2=[C1-MAX_" + nameImg + "] c3=[C2-MAX_" + nameImg + "] create");
run("RGB Color");
days = "20";
name = "WT 4"
saveAs("Tiff", "F:/Lab/MouseMusclePOGLUT1/data/rawImages/tifs/"+days+" days/"+ name +" p" +days+".tif");
