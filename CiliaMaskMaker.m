function [ciliaMask, allCilia] = CiliaMaskMaker(PARAMS, binCilia, ciliaI0)

[L, num] = bwlabel(binCilia);

for k = 1 : num
    
        singleCilia = ismember(L, k);
        trueCilia(k) = 1;
        angle = cell2mat(struct2cell(regionprops(singleCilia, 'Orientation')));
        se = strel('line', PARAMS.labelSize, angle);
        se2 = strel('line', PARAMS.labelSize/2, angle + 90);
        ciliaMask1 = imdilate(singleCilia, se2);
        ciliaMask0 = imdilate(ciliaMask1, se);
        imshow(labeloverlay(ciliaI0, ciliaMask0));
        ciliaMask{k} = ciliaMask0;
        allCilia = sum(trueCilia);
        disp(allCilia);
        
end
end
