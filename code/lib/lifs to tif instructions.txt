1. run mergedSeriesDetection.ij to detect the series where the whole image within the .lif is reconstructed.

2. Open the .lif manually in ImageJ specifying the specific series to load (always include the Number + 1). Load only the largest series, avoiding the individual tiles.

3. Once opened. Run macroSaveLifsIntoTifs.ij

4. Created the histogram matched version of the .tif images