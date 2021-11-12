%select image to analyse
clear all;

addpath(genpath('lib'))

pathImages = dir(fullfile('..','data','trainingDataSet_MouseMusclePoglut','muscleFibresOutlines','y_CellPose','*.png'));

for nImg = 18:size(pathImages,1)
    
    imageName = strrep(pathImages(nImg).name,'.png','');
    imgCells = imread(fullfile(pathImages(nImg).folder,pathImages(nImg).name));
    imgCells = imgCells(:,:,1);
    uniqueCells = unique(imgCells(:));
    uniqueCells(uniqueCells==0)=[];
    imgPerims = false(size(imgCells));
    for nCell = uniqueCells'
        imgAux = bwperim(imgCells==nCell);
        imgPerims(imgAux>0)=1;    
    end
    dilatedPerims = imdilate(imgPerims,strel('disk',1));
    
%     imwrite(dilatedPerims,fullfile(pathImages(nImg).folder,'outlines',[imageName '.tif']))
    imgColor = imread(fullfile('..','data','trainingDataSet_MouseMusclePoglut','muscleFibresOutlines','x',[imageName '.tif']));
    imgMembranes = imgColor(:,:,2);
    
    volumeSegmenter(cat(3,imgMembranes,imgMembranes),cat(3,dilatedPerims,dilatedPerims))
    pause()
    labels = bwareaopen(labels(:,:,1),100);
    labels = 1-bwareaopen(1-labels,100);
    labels(1:2,:)=0;
    labels(end-1:end,:)=0;
    labels(:,1:2)=0;
    labels(:,end-1:end)=0;
    imwrite(labels(:,:,1),fullfile('..','data','trainingDataSet_MouseMusclePoglut','muscleFibresOutlines','y_CellPose','outlines',[imageName '.tif']))
end