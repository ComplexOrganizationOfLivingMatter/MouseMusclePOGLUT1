// Obtener el nombre de la imagen original
nameImg = getTitle();
selectImage(nameImg);

// Proyección Z (máxima intensidad)
run("Z Project...", "projection=[Max Intensity]");

// Separar canales
run("Split Channels");

// Función para escalar canal 16-bit a 8-bit sin saturar
// Función para mejorar contraste y convertir a 8-bit correctamente
function enhanceAndConvert(imTitle) {
    selectImage(imTitle);
    run("Enhance Contrast", "saturated=1");
    run("8-bit");
}

// Aplicar a los 3 canales
enhanceAndConvert("C3-MAX_" + nameImg); // → R
enhanceAndConvert("C1-MAX_" + nameImg); // → G
enhanceAndConvert("C2-MAX_" + nameImg); // → B

// Combinar en RGB con orden C3-R, C1-G, C2-B
run("Merge Channels...", 
  "c1=[C3-MAX_" + nameImg + "] " +
  "c2=[C1-MAX_" + nameImg + "] " +
  "c3=[C2-MAX_" + nameImg + "] " +
  "create");

// Convertir a imagen RGB (8-bit por canal)
run("RGB Color");

days = "10";
name = "KIKO 345";
datePhoto = "200525";
saveAs("Tiff", "F:/Lab/MouseMusclePOGLUT1/data/rawImages/tifs/p"+days+"/"+ name +" p" +days+" "+datePhoto+".tif");
close();
close();