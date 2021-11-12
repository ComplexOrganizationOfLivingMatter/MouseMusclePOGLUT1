%select image to analyse
clear all;

addpath(genpath('lib'))

pathProjectionsImages = dir(fullfile('..','data','maxProjectedImages','*.tif'));
pathMasks = dir(fullfile('..','data','matlabMask','*.mat'));
pathPredictions = dir(fullfile('..','results','predictedImages','*.tif'));

boneZone = 6;

for nImg = 1:size(pathProjectionsImages,1)

    imageName = strrep(pathProjectionsImages(nImg).name,'.tif','');
    
    disp(imageName)
    
    infoImg = imfinfo(fullfile(pathProjectionsImages(nImg).folder,[imageName '.tif']));

    imgProjection = imread(fullfile(pathProjectionsImages(nImg).folder,[imageName '.tif']));
    %select mask delimiting the muscle zones to analyse
    load(fullfile(pathMasks(nImg).folder,[imageName '.mat']),'labels');
    maskZones = relabelCategoriesToImg(labels);

    imgPrediction = imread(fullfile(pathPredictions(nImg).folder,[imageName '.tif']));
    
    %binarize with autothreshold
    BW = imbinarize(imgPrediction);
    
    %make compact the detected objects   
    BW_dilated = imdilate(BW,strel('disk',1));
    BW_erode = imerode(BW_dilated,strel('disk',1)); 
    finalNuclei = BW_erode;
    
    taggedNuclei = bwlabel(BW_erode,4);
    areaNuclei = regionprops(taggedNuclei,'Area');  

    %limit a minimum threshold size of 5um^2
    thresholdNucleiArea = 5;
    areaMicrometers = vertcat(areaNuclei.Area)/(infoImg.XResolution*infoImg.YResolution);
    
    labels2del = find(areaMicrometers < thresholdNucleiArea);
    
    finalNuclei(ismember(taggedNuclei,labels2del))=0;
    
    
    clearvars BW BW_dilated BW_erode taggedNuclei imgPrediction
    
    path2saveBinaryImage = fullfile('..','results','binaryNucleiStemCells');
    if ~exist(fullfile(path2saveBinaryImage,[imageName '.tif']),'file')
        imwrite(finalNuclei,fullfile(path2saveBinaryImage,[imageName '.tif']))
    end
    
    %%calculate area per muscle zone
    areaZones= regionprops(maskZones,'Area');
    areaZonesMicrometers=vertcat(areaZones.Area)/(infoImg.XResolution*infoImg.YResolution);
    
    %count number of nuclei per zone
    centroidsNuclei = regionprops(finalNuclei,'Centroid');
    centroidsNuclei = round(vertcat(centroidsNuclei.Centroid));
    
    nLabels = maskZones(sub2ind(size(maskZones),centroidsNuclei(:,2),centroidsNuclei(:,1)));
    nNucleiPerZone = arrayfun(@(x) sum(ismember(nLabels,x)),[1:length(areaZonesMicrometers)]);
    
    %generate voronoi diagram
    finalNucleiLabel = bwlabel(finalNuclei);
    maskTotalVoronoi = zeros(size(finalNuclei));
    maxAreaVoronoiCellsPerZone = zeros(1,length(areaZonesMicrometers))';
    minAreaVoronoiCellsPerZone = zeros(1,length(areaZonesMicrometers))';
    stdAreaVoronoiCellsPerZone = zeros(1,length(areaZonesMicrometers))';
    for nZone =  1:length(areaZonesMicrometers)
        muscleBinaryZone = maskZones==nZone;
        if (nZone ~= boneZone) && areaZonesMicrometers(nZone)>0
            auxImg = VoronoizateCells(muscleBinaryZone,finalNucleiLabel.*muscleBinaryZone);
            voronoiArea = regionprops(auxImg,'Area');
            areaVoronoiCells = vertcat(voronoiArea.Area);
            areaVoronoiCells(areaVoronoiCells==0)=[];
            
            maxAreaVoronoiCellsPerZone(nZone)=max(areaVoronoiCells/(infoImg.XResolution*infoImg.YResolution));
            minAreaVoronoiCellsPerZone(nZone)=min(areaVoronoiCells/(infoImg.XResolution*infoImg.YResolution));
            stdAreaVoronoiCellsPerZone(nZone)=std(areaVoronoiCells/(infoImg.XResolution*infoImg.YResolution));
            
            maskTotalVoronoi(auxImg>0)=auxImg(auxImg>0);
        end
    end
    
    clearvars finalNucleiLabel muscleBinaryZone auxImg
    
    path2saveVoronoiImage = fullfile('..','results','voronoiImages');

    imwrite(uint16(maskTotalVoronoi),fullfile(path2saveVoronoiImage, [imageName '.tif']))
    
    T = array2table([areaZonesMicrometers,nNucleiPerZone',[areaZonesMicrometers./nNucleiPerZone'],stdAreaVoronoiCellsPerZone,maxAreaVoronoiCellsPerZone,minAreaVoronoiCellsPerZone],'VariableNames',{'muscleGroupArea','numberOfStemCellNuclei','areaPerNuclei','stdAreaVoronoiCell','maxAreaVoronoiCell','minAreaVoronoiCell'});
    writetable(T,fullfile('..','results',[imageName '.xls']));
 
    
end