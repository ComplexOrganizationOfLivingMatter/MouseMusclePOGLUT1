%measure number of cells in whole biopsy

clear all;

addpath(genpath('lib'))

pathProjectionsImages = dir(fullfile('..','data','maxProjectedImages','*.tif'));
pathMasks = dir(fullfile('..','data','matlabMask','*.mat'));

boneZone = 6;

rowTable = table();

for nImg = 1:size(pathProjectionsImages,1)

    imageName = strrep(pathProjectionsImages(nImg).name,'.tif','');
    disp(imageName)
    
    infoImg = imfinfo(fullfile(pathProjectionsImages(nImg).folder,[imageName '.tif']));
    
    %select mask delimiting the muscle zones to analyse
    load(fullfile(pathMasks(nImg).folder,[imageName '.mat']),'labelsWholeImg');
    labelRegions = relabelCategoriesToImg(labelsWholeImg);
    muscleZone = labelRegions>0 & labelRegions~=boneZone;
    
    areaMuscle = sum(muscleZone(:)==1);
    areaMuscleMicrometers = areaMuscle/(infoImg.XResolution*infoImg.YResolution);
    
    nucleiPax7Img = imread(fullfile('..','results','binaryNucleiStemCells',[imageName '.tif']));
    nucleiPax7Img(muscleZone==0)=0;
    labelNucleiPax7Img = bwlabel(nucleiPax7Img);
    areaNuclei = regionprops(labelNucleiPax7Img,'Area');
    totalNumberNuclei = sum(vertcat(areaNuclei.Area)>0);
    areaPerNuclei = areaMuscleMicrometers/totalNumberNuclei;
    
    rowTable = [rowTable; array2table([areaMuscleMicrometers,totalNumberNuclei,areaPerNuclei],'RowNames',{imageName}, 'VariableNames',{'totalArea','totalNumberOfNuclei','totalArea/totalNuclei'})];
end

writetable(rowTable,fullfile('..','results', ['totalNucleiPax7_' date '.xls']))