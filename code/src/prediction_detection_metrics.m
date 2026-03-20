close all
clear all
%% ================= USER INPUT =================
% Ground truth folder containing labeled masks
rootPath = 'F:\Lab\MouseMusclePOGLUT1\results\detectionMetrics';

gtFolder   = fullfile(rootPath, 'imagesTest','GT emi_vicky combined');

% Folder containing predicted masks
% predFolder = fullfile(rootPath, 'imagesTest','modelV1');
namePredImages = 'modelV1';
predFolder = fullfile(rootPath, 'imagesTest',namePredImages);



distThresh = 10;       % Distance threshold (pixels) to consider a centroid matched
outExcel   = fullfile(rootPath, ['detection_evaluation_' namePredImages '_autobinary_160326.xlsx']); % Output Excel file
imgExt     = '*.tif'; % Image file extension
binaryThr  = 250;     % Threshold to binarize predicted masks

%%Min nuclei size
micronsPerPixel = 0.5055443;
minArea_um2 = 10; % minimum nuclei size in um^2
minArea_pixels = minArea_um2 / (micronsPerPixel^2);

%% ==============================================

% List all ground truth files
gtFiles = dir(fullfile(gtFolder, imgExt));

nImages = length(gtFiles); % Total number of images

% Initialize arrays to store per-image metrics
imageNames = string({gtFiles.name})'; % Image file names
TPs = zeros(nImages,1);              % True positives
FPs = zeros(nImages,1);              % False positives
FNs = zeros(nImages,1);              % False negatives
GTs = zeros(nImages,1);              % Number of ground-truth objects
PRs = zeros(nImages,1);              % Number of predicted objects
Precisions = zeros(nImages,1);       % Precision per image
Recalls = zeros(nImages,1);          % Recall per image
F1s = zeros(nImages,1);              % F1-score per image

%% ================= MAIN LOOP =================
for i = 1:nImages
    
    % --- Load images ---
    gt  = imread(fullfile(gtFolder, gtFiles(i).name)); % Ground truth
    pr  = imread(fullfile(predFolder, gtFiles(i).name)); % Prediction
    if ~isa(pr,'uint8')
        BW = uint8(pr*255); % Scale prediction to 0-255 if needed
    else
        BW=pr;
    end
    % --- Binarize images ---
    gt = gt > 0; % Ground truth binary mask
    % Predicted mask: threshold + fill holes + morphological close
    % pr = imclose(bwmorph(BW >= binaryThr,'fill'), strel('disk',1));
    pr = imclose(bwmorph(imbinarize(BW),'fill'), strel('disk',1));


    % Remove small nuclei
    pr = bwareaopen(pr, ceil(minArea_pixels)); % remove small objects

    % imwrite(uint8(pr*255),fullfile(rootPath, 'imagesTest',gtFiles(i).name))

    % --- Label connected components ---
    gtL = bwlabel(gt); % GT labeled objects
    prL = bwlabel(pr); % Predicted labeled objects

    % --- Extract centroids ---
    gtProps = regionprops(gtL, 'Centroid'); % GT centroids
    prProps = regionprops(prL, 'Centroid'); % Predicted centroids

    nGT = numel(gtProps); % Number of GT objects
    nPR = numel(prProps); % Number of predicted objects

    GTs(i) = nGT; % Store number of GT objects
    PRs(i) = nPR; % Store number of predicted objects

    % --- Edge cases: zero objects ---
    if nGT == 0 && nPR == 0
        % Nothing in GT or prediction → perfect
        TPs(i) = 0; FPs(i) = 0; FNs(i) = 0;
        Precisions(i) = 1; Recalls(i) = 1; F1s(i) = 1;
        continue
    end

    if nGT == 0
        % No GT objects → all predictions are false positives
        TPs(i) = 0; FPs(i) = nPR; FNs(i) = 0;
        Precisions(i) = 0; Recalls(i) = 1; F1s(i) = 0;
        continue
    end

    if nPR == 0
        % No predicted objects → all GT are false negatives
        TPs(i) = 0; FPs(i) = 0; FNs(i) = nGT;
        Precisions(i) = 1; Recalls(i) = 0; F1s(i) = 0;
        continue
    end

    % --- Build centroid arrays ---
    gtC = cat(1, gtProps.Centroid); % nGT x 2 array
    prC = cat(1, prProps.Centroid); % nPR x 2 array

    % --- Compute distance matrix between GT and predicted centroids ---
    D = pdist2(gtC, prC); % nGT x nPR

    matchedGT = false(nGT,1); % Track which GT are matched
    matchedPR = false(nPR,1); % Track which predictions are matched

    % --- Greedy one-to-one matching ---
    TP_count = 0;
    while true
        [minD, idx] = min(D(:)); % Find closest pair
        if minD > distThresh
            break % No pair within threshold
        end
        [g,p] = ind2sub(size(D), idx); % Get GT index g and predicted index p

        % Mark as matched
        TP_count = TP_count + 1;
        matchedGT(g) = true;
        matchedPR(p) = true;

        % Remove matched row and column from distance matrix
        D(g,:) = inf;
        D(:,p) = inf;
    end

    % --- Store per-image results ---
    TPs(i) = TP_count;
    FNs(i) = sum(~matchedGT);
    FPs(i) = sum(~matchedPR);

    % Compute Precision / Recall / F1 per image
    Precisions(i) = TPs(i) / (TPs(i) + FPs(i) + eps);
    Recalls(i) = TPs(i) / (TPs(i) + FNs(i) + eps);
    F1s(i) = 2 * Precisions(i) * Recalls(i) / (Precisions(i) + Recalls(i) + eps);
end

%% ================= CREATE TABLE =================
timestamp = datetime('now');

% Per-image table
T = table(imageNames, TPs, FPs, FNs, GTs, PRs, Precisions, Recalls, F1s, ...
    repmat(distThresh,nImages,1), ...
    'VariableNames', {'ImageName','TP','FP','FN','Total_GT','Total_Pred','Precision','Recall','F1','DistanceThreshold_px'});

% --- Add a final row with averages across all images ---
avgT = table( ...
    "Average", ...                  % ImageName
    mean(TPs), ...                  % Mean TP
    mean(FPs), ...                  % Mean FP
    mean(FNs), ...                  % Mean FN
    mean(GTs), ...                  % Mean Total_GT
    mean(PRs), ...                  % Mean Total_Pred
    mean(Precisions), ...           % Mean Precision
    mean(Recalls), ...              % Mean Recall
    mean(F1s), ...                  % Mean F1
    distThresh, ...                 % DistanceThreshold_px
    'VariableNames', T.Properties.VariableNames);

T = [T; avgT]; % Append average row

%% ================= SAVE TO EXCEL =================
if isfile(outExcel)
    % Read existing table
    Told = readtable(outExcel);

    % Ensure matching types for safe concatenation
    if ~isa(Told.ImageName,'string'), Told.ImageName = string(Told.ImageName); end
    numericCols = {'TP','FP','FN','Total_GT','Total_Pred','Precision','Recall','F1','DistanceThreshold_px'};
    for k = 1:length(numericCols)
        col = numericCols{k};
        if ~isa(Told.(col),'double'), Told.(col) = double(Told.(col)); end
    end

    % Concatenate old and new tables
    T = [Told; T];
end

% Write table to Excel
writetable(T, outExcel);

%% ================= REPORT =================
fprintf('\n=== Detection results ===\n');
fprintf('Number of images: %d\n', nImages);
fprintf('Saved per-image metrics and average to: %s\n', outExcel);
