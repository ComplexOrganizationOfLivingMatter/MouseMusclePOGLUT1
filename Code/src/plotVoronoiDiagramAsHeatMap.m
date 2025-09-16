clear all;

addpath(genpath('lib'))

pathMasks = fullfile('..','..','data','annotatedMuscleMasks');
pathVoronoi = fullfile('..','..','results','voronoiImages');

boneZone = 6;

cmp = colormap('hot');
cmpHotInverse =  cmp(end:-1:1,:);

%%iterate across all the WT to get their average, min and max "cell size"

ages={'p4','p7', 'p10','p15','p20','p35'};
areaVorPerAge=cell(size(ages));
for pAge = 1:length(ages)
    pathVImages = dir(fullfile(pathVoronoi,ages{pAge},'*.tif'));
    for nImg = 1:size(pathVImages,1)
        if ~contains(pathVImages(nImg).name, 'KI')
            imgVoronoi = imread(fullfile(pathVImages(nImg).folder,pathVImages(nImg).name));
            areaVor = regionprops(imgVoronoi,'Area');
            areaVor = vertcat(areaVor.Area);
            areaVorPerAge{pAge} = [areaVorPerAge{pAge};areaVor(areaVor>0)];

        end
    end
end

meanAreaAge = cellfun(@mean, areaVorPerAge);
minAreaAge = cellfun(@min, areaVorPerAge);
maxAreaAge = cellfun(@max, areaVorPerAge);

for pAge = 1:length(ages)
    pathVImages = dir(fullfile(pathVoronoi,ages{pAge},'*.tif'));
    for nImg = 1:size(pathVImages,1)
        try
            imageName = strrep(pathVImages(nImg).name,'.tif','');
            disp(imageName)
                
            imgVoronoi = imread(fullfile(pathVImages(nImg).folder,pathVImages(nImg).name));

            %select mask delimiting the muscle zones to analyse
            maskZones=imread(fullfile(pathMasks,ages{pAge},pathVImages(nImg).name));
            maskZones = imresize(maskZones,size(imgVoronoi),'nearest');

            perimZones = imdilate(bwperim(maskZones>0),strel('disk',5));
            
           
            areaVorInit = regionprops(imgVoronoi,'Area');
            areaVorInit = vertcat(areaVorInit.Area);
            areaVor=areaVorInit(areaVorInit>0);
        
            averageWTarea = meanAreaAge(pAge);
            normArea = log10(areaVor/averageWTarea);
            % figure;hist(normArea)
        
        
            maxTotalArea = 1;%max(areaVor);
            minTotalArea = -1;%min(areaVor);
            normArea(normArea>maxTotalArea)=maxTotalArea;
            normArea(normArea<minTotalArea)=minTotalArea;
            
            %normalize
            % normAreaTotal = (areaVor-minTotalArea)/(maxTotalArea-minTotalArea);
            labelsVor = find(areaVorInit>0);
            maskAreaVorNormTotal = zeros(size(maskZones));
            for nLab = 1:length(labelsVor)
                maskAreaVorNormTotal(imgVoronoi==labelsVor(nLab))= normArea(nLab);
            end     
            combiImage = maskAreaVorNormTotal;
            combiImage(perimZones)=1;
            combiImage(maskZones==0)=-1;
            combiImage(maskZones==6)=-1;
            fig=figure;
            imshow(combiImage);
            colormap(cmpHotInverse)
            caxis([-1 1]); % Define el rango de valores a mostrar en el colormap
            colorbar('Ticks',[-1:0.25:1])
            hcb=colorbar;
            hcb.Label.String = 'log 10 (area\_i / area\_mean\_WT\_age)';  % Añadir el tag o etiqueta
            
        
            folder2save = fullfile('..','..','results','voronoiImages','heatMaps','areaNormalizedByAvgWT_Age');
            if ~exist(folder2save,'dir')
                mkdir(folder2save)
                % mkdir(fullfile('..','..','results','voronoiImages','heatMaps','areaNormalizePerMuscleGroup'))
            end
            exportgraphics(gca,fullfile(folder2save,ages{pAge},[imageName '.pdf']),'Resolution',600)
            
            clearvars combiImage maskAreaVorNormTotal
            close all
        catch
            disp(['error in ' imageName])
        end

    % %normalize imagePerRegion
    % maskAreaVorPerZone = zeros(size(maskZones));
    % for nZone =  1:max(unique(maskZones(:)))
    %     muscleBinaryZone = maskZones==nZone;
    %     if (nZone ~=boneZone) && (length(unique(muscleBinaryZone(:)))>1)
    %         auxVorImage = zeros(size(maskZones));
    %         auxVorImage(muscleBinaryZone>0) = imgVoronoi(muscleBinaryZone>0);
    %         areaVor = regionprops(auxVorImage,'Area');
    %         areaVor = vertcat(areaVor.Area);
    % 
    %         labelsVor = find(areaVor>0);
    % 
    %         maxTotalArea = max(areaVor);
    %         minTotalArea = min(areaVor);
    %         %normalize
    %         normAreaZone = (areaVor-minTotalArea)/(maxTotalArea-minTotalArea);
    % 
    %         for nLab = labelsVor'
    %             maskAreaVorPerZone(imgVoronoi==nLab)= normAreaZone(nLab);
    %         end           
    %     end
    % end
    % 
    % 
    % combiImage = maskAreaVorPerZone;
    % combiImage(perimZones)=1;
    % fig2 = figure;
    % imshow(combiImage);
    % colormap(cmpHotInverse)
    % colorbar()
    % exportgraphics(gca,fullfile('..','..','results','voronoiImages','heatMaps','areaNormalizePerMuscleGroup',[imageName '.pdf']),'Resolution',600)
    % 
    % clearvars combiImage maskAreaVorPerZone muscleBinaryZone auxVorImage perimZones
    % close all

    
    end
end