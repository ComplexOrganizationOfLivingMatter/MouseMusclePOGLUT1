%% Generate training dataset and test dataset

cropSize = 1020;
nPairsToCrop = 50; % 250 for train and 50 for test
dataPath = fullfile('..','..','data');
pathTifs = dir(fullfile(dataPath,'annotatedMuscleMasks','**','*.tif'));
pathFilteredMask = fullfile('..','..','results', "binaryNucleiStemCells_filtered");
pathRaw = fullfile(dataPath, "rawImages","tifsHistMatched");
%% Load images
maskFiles = dir(fullfile(pathFilteredMask, '*.tif'));
rawFiles  = dir(fullfile(pathRaw, '**','*.tif'));

path2SaveDataset = fullfile(dataPath,"trainingDataset","individual label","Pax7_v2","backup");
outputTrainRaw = fullfile(path2SaveDataset,"train","raw");
outputTrainMask = fullfile(path2SaveDataset,"train","mask");
outputTestRaw = fullfile(path2SaveDataset,"test","raw");
outputTestMask = fullfile(path2SaveDataset,"test","mask");

% Create output folders if they don't exist
folders = {outputTrainRaw, outputTrainMask, outputTestRaw, outputTestMask};
for f = 1:numel(folders)
    if ~exist(folders{f}, 'dir')
        mkdir(folders{f});
    end
end



%% Prepare list of matching mask–raw pairs
rawFiles  = dir(fullfile(pathRaw, '**', '*.tif'));
maskFiles = dir(fullfile(pathFilteredMask, '*.tif'));

rawNames = {rawFiles.name};
maskNames = {maskFiles.name};
selectedNames = {pathTifs.name};

% Keep only files that exist both as filtered mask and as raw image
validPairs = {};
for i = 1:numel(selectedNames)
    name = selectedNames{i};
    if any(strcmp(maskNames, name)) && any(strcmp(rawNames, name))
        validPairs{end+1} = name; %#ok<SAGROW>
    end
end

fprintf('Found %d valid image pairs.\n', numel(validPairs));
if isempty(validPairs)
    error('No matching mask–raw image pairs found.');
end

%% Random cropping and saving
nCropsSaved = 0;
maxAttempts = nPairsToCrop * 10;
attempt = 0;
rng('shuffle');

cropHistory = containers.Map(); % keys = filenames, values = Nx4 arrays of [x y cropSize cropSize]

while nCropsSaved < nPairsToCrop && attempt < maxAttempts
    attempt = attempt + 1;

    % Randomly pick a valid pair
    idx = randi(numel(validPairs));
    fileName = validPairs{idx};

    % Load mask and raw
    maskPath = fullfile(pathFilteredMask, fileName);
    mask = imread(maskPath);

    rawIdx = find(strcmp(rawNames, fileName), 1);
    rawPath = fullfile(rawFiles(rawIdx).folder, rawFiles(rawIdx).name);
    raw = imread(rawPath);

    % Check image size
    [h, w, ~] = size(mask);
    if h < cropSize || w < cropSize
        continue;
    end

    % Random crop coordinates
    x = randi([1, w - cropSize + 1]);
    y = randi([1, h - cropSize + 1]);
    newRect = [x, y, cropSize, cropSize];

    % Check overlap with previous crops of the same image
    if isKey(cropHistory, fileName)
        prevRects = cropHistory(fileName);
        if any(rectint(prevRects, newRect) > 0)
            continue; % overlap detected, skip
        end
    else
        prevRects = [];
    end

    % Crop the images
    maskCrop = mask(y:y+cropSize-1, x:x+cropSize-1);
    rawCrop = raw(y:y+cropSize-1, x:x+cropSize-1, :);

    % Skip empty mask crops
    if all(maskCrop(:) == 0)
        continue;
    end

    % Save crop
    nCropsSaved = nCropsSaved + 1;
    if nCropsSaved <= 250
        outMask = fullfile(outputTrainMask, sprintf('mask_%03d.tif', nCropsSaved));
        outRaw  = fullfile(outputTrainRaw, sprintf('raw_%03d.tif', nCropsSaved));
    else
        outMask = fullfile(outputTestMask, sprintf('mask_%03d.tif', nCropsSaved-250));
        outRaw  = fullfile(outputTestRaw, sprintf('raw_%03d.tif', nCropsSaved-250));
    end

    imwrite(maskCrop, outMask);
    imwrite(rawCrop, outRaw);

    % Update crop history
    cropHistory(fileName) = [prevRects; newRect];

    fprintf('Saved crop %d/%d from %s\n', nCropsSaved, nPairsToCrop, fileName);
end


if nCropsSaved < nPairsToCrop
    warning('Only %d crops were generated (insufficient non-empty regions).', nCropsSaved);
else
    disp('✅ Dataset generation completed successfully.');
end