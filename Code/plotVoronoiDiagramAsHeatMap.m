clear all;

addpath(genpath('lib'))

pathProjectionsImages = dir(fullfile('..','data','maxProjectedImages','*.tif'));
pathMasks = dir(fullfile('..','data','matlabMask','*.mat'));
pathVoronoi = dir(fullfile('..','results','voronoiImages','*.tif'));

boneZone = 6;

cmp = colormap('hot');
cmpHotInverse =  cmp(end:-1:1,:);

for nImg = 2%1:size(pathProjectionsImages,1)

    imageName = strrep(pathProjectionsImages(nImg).name,'.tif','');
    disp(imageName)
    
    infoImg = imfinfo(fullfile(pathProjectionsImages(nImg).folder,[imageName '.tif']));
    
    %select mask delimiting the muscle zones to analyse
    load(fullfile(pathMasks(nImg).folder,[imageName '.mat']),'labels');
    maskZones = relabelCategoriesToImg(labels);
    perimZones = imdilate(bwperim(maskZones>0),strel('disk',5));
    
    imgVoronoi = imread(fullfile(pathVoronoi(nImg).folder,[imageName '.tif']));
   
    areaVor = regionprops(imgVoronoi,'Area');
    areaVor = vertcat(areaVor.Area);
    
    maxTotalArea = max(areaVor);
    minTotalArea = min(areaVor);
    
    %normalize
    normAreaTotal = (areaVor-minTotalArea)/(maxTotalArea-minTotalArea);
    labelsVor = find(areaVor>0);
    maskAreaVorNormTotal = zeros(size(maskZones));
    for nLab = labelsVor'
        maskAreaVorNormTotal(imgVoronoi==nLab)= normAreaTotal(nLab);
    end     
    combiImage = maskAreaVorNormTotal;
    combiImage(perimZones)=1;
    fig=figure;
    imshow(combiImage);
    colormap(cmpHotInverse)
    colorbar()
    exportgraphics(gca,fullfile('..','results','voronoiImages','heatMaps','areaNormalizeAllImage',[imageName '.pdf']),'Resolution',600)
    
    clearvars combiImage maskAreaVorNormTotal
    close all

    %normalize imagePerRegion
    maskAreaVorPerZone = zeros(size(maskZones));
    for nZone =  1:max(unique(maskZones(:)))
        muscleBinaryZone = maskZones==nZone;
        if (nZone ~=boneZone) && (length(unique(muscleBinaryZone(:)))>1)
            auxVorImage = zeros(size(maskZones));
            auxVorImage(muscleBinaryZone>0) = imgVoronoi(muscleBinaryZone>0);
            areaVor = regionprops(auxVorImage,'Area');
            areaVor = vertcat(areaVor.Area);
            
            labelsVor = find(areaVor>0);
            
            maxTotalArea = max(areaVor);
            minTotalArea = min(areaVor);
            %normalize
            normAreaZone = (areaVor-minTotalArea)/(maxTotalArea-minTotalArea);
            
            for nLab = labelsVor'
                maskAreaVorPerZone(imgVoronoi==nLab)= normAreaZone(nLab);
            end           
        end
    end
    

    combiImage = maskAreaVorPerZone;
    combiImage(perimZones)=1;
    fig2 = figure;
    imshow(combiImage);
    colormap(cmpHotInverse)
    colorbar()
    exportgraphics(gca,fullfile('..','results','voronoiImages','heatMaps','areaNormalizePerMuscleGroup',[imageName '.pdf']),'Resolution',600)

    clearvars combiImage maskAreaVorPerZone muscleBinaryZone auxVorImage perimZones
    close all

    
end