%select image to analyse
clear all;

addpath(genpath('lib'))

pathProjectionsImages = dir(fullfile('..','data','rawImages','tifs','**','*.tif'));
pathMasks = dir(fullfile('..','data','annotatedMuscleMasks','**','*.tif'));
pathPredictions = dir(fullfile('..','data','predictions','poglut1_pax7_2D_20240729-151810_PREDICT_OverHistMatched','**','per_image','*.tif'));

boneZone = 6;

%11L is labelled as 111. Regions with a label larger than maxMuscleZone
%will be reassignated to 15.
maxMuscleZone = 14;

numPropertiesPerImage = 2;
propertiesNames = {'totalNuclei', 'totalNucleiByAreaEstimation'};
imagesNames ={pathProjectionsImages(:).name};
T_GeneralResults = array2table(zeros(size(pathProjectionsImages,1), numPropertiesPerImage),'VariableNames',propertiesNames,'RowNames',imagesNames);

for nImg = 1:size(pathProjectionsImages,1)
    
    infoImg = imfinfo(fullfile(pathProjectionsImages(nImg).folder,pathProjectionsImages(nImg).name));

    imgProjection = imread(fullfile(pathProjectionsImages(nImg).folder,pathProjectionsImages(nImg).name));
    %select mask delimiting the muscle zones to analyse
        %%%% load(fullfile(pathMasks(nImg).folder,[imageName '.mat']),'labels');
        %%%% maskZones = relabelCategoriesToImg(labels);
    
    imgPrediction = imread(fullfile(pathPredictions(1).folder,pathProjectionsImages(nImg).name));
    
    %binarize with autothreshold
    BW = imbinarize(imgPrediction);
    
    %make compact the detected objects   
    BW_dilated = imdilate(BW,strel('disk',1));
    BW_erode = imerode(BW_dilated,strel('disk',1)); 
    finalNuclei = BW_erode;
    
    taggedNuclei = bwlabel(BW_erode,4);
    areaNuclei = regionprops(taggedNuclei,'Area');  

    %limit a minimum threshold size of 5um^2. If area>50um2 then, try
    %watershed.
    minNucleiArea = 5;
    maxNucleiArea = 70;

    areaMicrometers = vertcat(areaNuclei.Area)/(infoImg.XResolution*infoImg.YResolution);
    
    % %Watershed large nuclei to try to separate them
    % labels2Watershed = find(areaMicrometers >= maxNucleiArea);
    % maskWS = ismember(taggedNuclei,labels2Watershed);
    % D = bwdist(~maskWS);
    % L = watershed(-D);
    % L(~maskWS)=0;

    labels2del = find(areaMicrometers < minNucleiArea);
    
    finalNuclei(ismember(taggedNuclei,labels2del))=0;
    
    
    clearvars BW BW_dilated BW_erode taggedNuclei imgPrediction
    
    path2saveBinaryImage = fullfile('..','results','binaryNucleiStemCells');
    if ~exist(path2saveBinaryImage,'dir'), mkdir(path2saveBinaryImage); end

    if ~exist(fullfile(path2saveBinaryImage,pathProjectionsImages(nImg).name),'file')
        imwrite(uint8(finalNuclei)*255,fullfile(path2saveBinaryImage,pathProjectionsImages(nImg).name))
    end
    
    
    %count number of nuclei per zone
    propsNuclei = regionprops(finalNuclei,'Centroid','Area');
    centroidsNuclei = round(vertcat(propsNuclei.Centroid));

    areaFinalNuclei = round(vertcat(propsNuclei.Area));
    areaMicrometersFinal = vertcat(areaFinalNuclei)/(infoImg.XResolution*infoImg.YResolution);

    totalNucleiWholeImg = size(centroidsNuclei,1); 

    %count nuclei divided by max area
    nucleiByArea = arrayfun(@(x) ceil(x/maxNucleiArea), areaMicrometersFinal);
    totalNucleiByArea = sum(nucleiByArea);


    try 
        maskZones = imread(fullfile(pathMasks(1).folder,pathProjectionsImages(nImg).name));
        %
        maskZones(maskZones>maxMuscleZone) = maxMuscleZone+1;
        %%calculate area per muscle zone
        areaZones= regionprops(maskZones,'Area');
        areaZonesMicrometers=vertcat(areaZones.Area)/(infoImg.XResolution*infoImg.YResolution);
        
        %count number of nuclei per zone
        nLabels = maskZones(sub2ind(size(maskZones),centroidsNuclei(:,2),centroidsNuclei(:,1)));
        nNucleiPerZone = arrayfun(@(x) sum(ismember(nLabels,x)),[1:length(areaZonesMicrometers)]);
        nNucleiPerZoneByNucleiArea = arrayfun(@(x) sum(ismember(nLabels,x).*nucleiByArea),[1:length(areaZonesMicrometers)]);


        % %generate voronoi diagram
        % finalNucleiLabel = bwlabel(finalNuclei);
        % maskTotalVoronoi = zeros(size(finalNuclei));
        % maxAreaVoronoiCellsPerZone = zeros(1,length(areaZonesMicrometers))';
        % minAreaVoronoiCellsPerZone = zeros(1,length(areaZonesMicrometers))';
        % stdAreaVoronoiCellsPerZone = zeros(1,length(areaZonesMicrometers))';
        % for nZone =  1:length(areaZonesMicrometers)
        %     muscleBinaryZone = maskZones==nZone;
        %     if (nZone ~= boneZone) && areaZonesMicrometers(nZone)>0
        %         auxImg = VoronoizateCells(muscleBinaryZone,finalNucleiLabel.*muscleBinaryZone);
        %         voronoiArea = regionprops(auxImg,'Area');
        %         areaVoronoiCells = vertcat(voronoiArea.Area);
        %         areaVoronoiCells(areaVoronoiCells==0)=[];
        % 
        %         maxAreaVoronoiCellsPerZone(nZone)=max(areaVoronoiCells/(infoImg.XResolution*infoImg.YResolution));
        %         minAreaVoronoiCellsPerZone(nZone)=min(areaVoronoiCells/(infoImg.XResolution*infoImg.YResolution));
        %         stdAreaVoronoiCellsPerZone(nZone)=std(areaVoronoiCells/(infoImg.XResolution*infoImg.YResolution));
        % 
        %         maskTotalVoronoi(auxImg>0)=auxImg(auxImg>0);
        %     end
        % end
        
        clearvars finalNucleiLabel muscleBinaryZone auxImg
        
        % path2saveVoronoiImage = fullfile('..','results','voronoiImages');
        % if ~exist(path2saveVoronoiImage,'dir'), mkdir(path2saveVoronoiImage); end
        % 
        % imwrite(uint16(maskTotalVoronoi),fullfile(path2saveVoronoiImage, pathProjectionsImages(nImg).name))
        T = array2table([[1:length(areaZonesMicrometers)]',areaZonesMicrometers,nNucleiPerZone',nNucleiPerZoneByNucleiArea',areaZonesMicrometers./nNucleiPerZone',areaZonesMicrometers./nNucleiPerZoneByNucleiArea'],'VariableNames',{'muscleGroups','muscleGroupArea','numberOfStemCellNuclei','numberOfStemCellNucleiByAreaEstimation','areaPerNuclei','areaPerNucleiByAreaEstimation'});
        % T = array2table([areaZonesMicrometers,nNucleiPerZone',nNucleiPerZoneByNucleiArea',areaZonesMicrometers./nNucleiPerZone',areaZonesMicrometers./nNucleiPerZoneByNucleiArea',stdAreaVoronoiCellsPerZone,maxAreaVoronoiCellsPerZone,minAreaVoronoiCellsPerZone],'VariableNames',{'muscleGroupArea','numberOfStemCellNuclei','numberOfStemCellNucleiByAreaEstimation','areaPerNuclei','areaPerNucleiByAreaEstimation','stdAreaVoronoiCell','maxAreaVoronoiCell','minAreaVoronoiCell'});
        writetable(T,fullfile('..','results','numberOfPax7PerMusclePerImage.xlsx'),'sheet',strrep(pathProjectionsImages(nImg).name,'.tif',''));

    catch
        disp([pathProjectionsImages(nImg).name '-- Does not have a paired painted muscle mask'])
    end
        
  
    T_GeneralResults.totalNuclei(nImg)=totalNucleiWholeImg;
    T_GeneralResults.totalNucleiByAreaEstimation(nImg)=totalNucleiByArea;
    
end

writetable(T_GeneralResults,fullfile('..','results','numberOfPax7PerImage.xlsx'),'WriteRowNames', true);