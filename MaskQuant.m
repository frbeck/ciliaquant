function [nResults] = MaskQuant(n, ciliaMask, targetMask, targetI)

currentCiliaMask = cell2mat(ciliaMask(n));
inCiliaMean = cell2mat(struct2cell(regionprops(currentCiliaMask, targetI, 'MeanIntensity')));
inCiliaMax = cell2mat(struct2cell(regionprops(currentCiliaMask, targetI, 'MaxIntensity')));
aroundCiliaMean = mean(cell2mat(struct2cell(regionprops(targetMask, targetI, 'MeanIntensity'))));
aroundCiliaMax = mean(cell2mat(struct2cell(regionprops(targetMask, targetI, 'MaxIntensity'))));
nResults = [inCiliaMean inCiliaMax aroundCiliaMean aroundCiliaMax];

end