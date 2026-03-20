run("Bio-Formats Macro Extensions");
// Paso 1: Obtener archivo .lif manualmente (compatibilidad completa)
dir = getDirectory("Selecciona el directorio del archivo .lif");
name = getString("Nombre exacto del archivo (.lif):", "archivo.lif");
path = dir + name;

// Paso 2: Preparar Bio-Formats
Ext.setId(path);
Ext.getSeriesCount(seriesCount);

print("🔍 Se encontraron " + seriesCount + " series en el archivo:\n");

// Paso 3: Inicializar arrays
sizeArray = newArray(seriesCount);
indexArray = newArray(seriesCount);
labelArray = newArray(seriesCount);

// Paso 4: Recorrer y guardar información de cada serie
for (i = 0; i < seriesCount; i++) {
    Ext.setSeries(i);
    Ext.getSizeX(x);
    Ext.getSizeY(y);
    Ext.getSizeZ(z);
    Ext.getSizeC(c);
    Ext.getSizeT(t);
    pixels = x * y;
    sizeArray[i] = pixels;
    indexArray[i] = i;
    labelArray[i] = "Serie " + i + ": " + x + "x" + y + " Z=" + z + " C=" + c + " T=" + t + " → " + pixels + " píxeles";
}
Ext.close();

// Paso 5: Ordenar las series por tamaño XY (burbuja)
for (i = 0; i < seriesCount-1; i++) {
    for (j = i+1; j < seriesCount; j++) {
        if (sizeArray[j] > sizeArray[i]) {
            // Intercambiar tamaño
            temp = sizeArray[i]; sizeArray[i] = sizeArray[j]; sizeArray[j] = temp;
            // Intercambiar índice
            temp = indexArray[i]; indexArray[i] = indexArray[j]; indexArray[j] = temp;
            // Intercambiar label
            temp = labelArray[i]; labelArray[i] = labelArray[j]; labelArray[j] = temp;
        }
    }
}

// Paso 6: Mostrar resultados ordenados
print("📋 Series ordenadas por resolución XY (de mayor a menor):");
for (i = 0; i < seriesCount; i++) {
    print((i+1) + ". " + labelArray[i]);
}
