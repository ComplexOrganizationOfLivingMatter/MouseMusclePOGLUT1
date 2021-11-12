function [voronoiImg] = VoronoizateCells(binaryMask,imgCells)

    voronoiImg=imgCells.*cast(binaryMask,class(imgCells));

    perimCells=bwperim(voronoiImg>0);
    
    %Get bounded valid pixels
    idsToFill = find(binaryMask==1 & imgCells==0);
    [row, col] = ind2sub(size(binaryMask),idsToFill);
    labelPerId = zeros(size(idsToFill));
    
    idsPerim = find(perimCells==1);
    [rowPer, colPer] = ind2sub(size(binaryMask),idsPerim);
    labelsPerimIds = voronoiImg(perimCells);
    
    %From valid pixels get closest seed (add this value)
    %tic
%   disp('generating 2D Voronoi')

    parfor nId = 1:length(idsToFill)
        distCoord = pdist2([col(nId),row(nId)],[colPer,rowPer]);
        [~,idSeedMin]=min(distCoord);
        labelPerId(nId) = labelsPerimIds(idSeedMin);
    end

    voronoiImg(idsToFill)=labelPerId;
end

