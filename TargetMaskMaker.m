function [targetMask] = TargetMaskMaker(n, ciliaMask, targetI, PARAMS)

    currentCiliaMask = cell2mat(ciliaMask(n));
    regionMean = cell2mat(struct2cell(regionprops(currentCiliaMask, targetI, 'MeanIntensity')));
    regionMax = cell2mat(struct2cell(regionprops(currentCiliaMask, targetI, 'MaxIntensity')));
    targetI(~currentCiliaMask) = 0;
    imshow(targetI);
    targetMask = targetI > 2*regionMean;
    se = strel('diamond', 3);
    targetMask = imdilate(targetMask, se);
    targetMask = imfill(targetMask, 'holes');
    targetMask = bwareafilt(targetMask,1);
    imshow(labeloverlay(targetI, targetMask)); 
end