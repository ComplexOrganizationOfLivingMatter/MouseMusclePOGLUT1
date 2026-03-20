%% Histogram Matching function

%% as reference image we will consider KIKO_1_7um
dirTifs = fullfile('..','data','rawImages','tifs');
daysFolder = 'p4';
imgReference = imread(fullfile(dirTifs,'Initial dataset','KIKO 1 p20.tif'));

listImages=dir(fullfile(dirTifs,daysFolder,'*26.tif'));

% Display the file names
for i = 1:length(listImages)
    disp(fullfile(listImages(i).folder, listImages(i).name));

    img2HM = imread(fullfile(listImages(i).folder, listImages(i).name));
    matchedImage = imhistmatch(img2HM, imgReference);
    targetDir = fullfile(dirTifs,'..','tifsHistMatched',daysFolder);
    if ~exist(targetDir,'dir')
        mkdir(targetDir)
    end
    imshow([img2HM,matchedImage])
    close;
    imwrite(matchedImage,fullfile(targetDir,listImages(i).name))
end

