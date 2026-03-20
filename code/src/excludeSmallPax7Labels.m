rootPath = "F:\Lab\MouseMusclePOGLUT1\data\trainingDataset\individual label\Pax7_v2_Mar";
folderAnnotatedImages = "after_curation";
folderFilteredImages = "after_curation_cleanArtifacts";

%% Parameters
micronsPerPixel = 0.5055443;
minArea_um2 = 10; % minimum nuclei size in um^2
minArea_pixels = minArea_um2 / (micronsPerPixel^2);

%% Get all tif images recursively from annotated folder
annotatedPath = fullfile(rootPath, folderAnnotatedImages);
imageFiles = dir(fullfile(annotatedPath, "**", "*.tif"));

%% Process each image
for k = 1:length(imageFiles)
    
    % Read image
    imgPath = fullfile(imageFiles(k).folder, imageFiles(k).name);
    img = imread(imgPath);
    
    % Ensure binary mask (0 or 255)
    imgBinary = img > 0;             % convert to logical
    imgBinary = bwareaopen(imgBinary, ceil(minArea_pixels)); % remove small objects
    imgFiltered = uint8(imgBinary) * 255; % convert back to uint8
    
    % Determine relative path to folderAnnotatedImages
    relPath = strrep(imgPath, annotatedPath, ""); % e.g., /train/mask/xyz.tif
    relPath = fullfile(folderFilteredImages, relPath); % target path in filtered folder
    
    % Make sure the folder exists
    saveDir = fullfile(rootPath,fileparts(relPath));
    if ~exist(saveDir, "dir")
        mkdir(saveDir);
    end
    
    % Save the filtered image
    imwrite(imgFiltered, fullfile(rootPath,relPath));
    fprintf("Processed and saved: %s\n", relPath);
end

disp("All images processed and filtered!");