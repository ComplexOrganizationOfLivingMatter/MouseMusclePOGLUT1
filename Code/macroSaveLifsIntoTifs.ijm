// 1st open images as hyperstacks

// Guardar nombre
nameImg = getTitle();

//run("Stack to Hyperstack...", "order=xyctz channels=3 slices=1 frames=1 display=Color");

// Actualizar dimensiones
getDimensions(width, height, channels, slices, frames);

// Caso: 1 canal, 3 slices mal codificados
if (channels == 1 && slices == 3) {
    run("Stack to Hyperstack...", "order=xyzct channels=3 slices=1 frames=1 display=Composite");
}

// Si hay Z, proyectar
getDimensions(width, height, channels, slices, frames);
if (slices > 1) {
    run("Z Project...", "projection=[Max Intensity]");
    prefixN = "MAX_";
} else {
    prefixN = "";
}

run("Split Channels");

// Función para escalar canal 16-bit a 8-bit sin saturar
// Función para mejorar contraste y convertir a 8-bit correctamente
function enhanceAndConvert(imTitle) {
    selectImage(imTitle);
    run("Enhance Contrast", "saturated=1");
    run("8-bit");
}

// Aplicar a los 3 canales
enhanceAndConvert("C3-" + prefixN + nameImg); // → R
enhanceAndConvert("C1-" + prefixN + nameImg); // → G
enhanceAndConvert("C2-" + prefixN + nameImg); // → B

// Combinar en RGB con orden C3-R, C1-G, C2-B
run("Merge Channels...", 
  "c1=[C3-" + prefixN + nameImg + "] " +
  "c2=[C1-" + prefixN + nameImg + "] " +
  "c3=[C2-" + prefixN + nameImg + "] " +
  "create ignore");

// Convertir a imagen RGB (8-bit por canal)
run("RGB Color");

days = "4";
name = "KIKO 777";
datePhoto = "100226";
saveAs("Tiff", "F:/Lab/MouseMusclePOGLUT1/data/rawImages/tifs/p"+days+"/"+ name +" p" +days+" "+datePhoto+".tif");
close();
close();