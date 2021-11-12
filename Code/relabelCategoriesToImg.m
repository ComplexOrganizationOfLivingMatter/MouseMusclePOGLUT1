function imgLabel = relabelCategoriesToImg(labels)

    allCategories = categories(labels);
    numCategories = cellfun(@(x) str2double(strrep(x,'Label','')), allCategories);

    imgUint8 = uint8(labels);
    indValues = unique(imgUint8(:));
    indValues(indValues==0)=[];

    auxImg = imgUint8;

    for nCateg = 1:length(indValues)
       imgUint8(auxImg==indValues(nCateg))=numCategories(nCateg);    
    end

    %delete empty slice
    imgLabel=imgUint8(:,:,1);

end