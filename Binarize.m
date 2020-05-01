function [binCilia] = Binarize(ciliaI, PARAMS)

ciliaLogical = imbinarize(ciliaI);
%[~,threshold] = edge(ciliaLogical,'sobel');
%BWs = edge(ciliaLogical,'sobel',threshold * PARAMS.fudgeFactor);  
se90 = strel('line',3,90);
se0 = strel('line',3,0);
BWsdil = imdilate(ciliaLogical,[se90 se0]);
BWdfill = imfill(BWsdil,'holes');
BWnobord = imclearborder(BWdfill,4);
seD = strel('diamond', PARAMS.cleaning);
BWfinal = imerode(BWnobord,seD);
BWfinal = imerode(BWfinal,seD);
binCilia = bwareaopen(BWfinal, PARAMS.minArea);

end
