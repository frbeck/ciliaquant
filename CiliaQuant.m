filename = 'n3f_mcherry_arl13b_clr_dox_noshh';
addpath('/Users/francisbeckert/Desktop/matlab folder/bfmatlab');
data = bfopen(filename);
PARAMS.totalFields = 2;
PARAMS.cleaning = 1;
PARAMS.fudgeFactor = 2;
PARAMS.minArea = 200;
PARAMS.labelSize = 100;
PARAMS.dapiNum=0; 
PARAMS.ciliaNum=2; 
PARAMS.targetNum=1; 
PARAMS.channelNum=2; 
PARAMS.max = false; 
iResults = [];
results = [];

for i = 1:PARAMS.totalFields
    
    PARAMS.numFields = i;
    
    [fields,zstacks, marker_images]=FileReader(data, PARAMS);

    ciliaI0 = cell2mat(marker_images(PARAMS.ciliaNum, 1));
    
    imshow(ciliaI0);
    
    targetI = cell2mat(marker_images(PARAMS.targetNum, 1));

    [ciliaI] = artifactCleaner(ciliaI0);

    [binCilia] = Binarize(ciliaI, PARAMS);

    [ciliaMask, allCilia] = CiliaMaskMaker(PARAMS, binCilia, ciliaI0);
    
    for n = 1:allCilia
        
        [targetMask] = TargetMaskMaker(n, ciliaMask, targetI, PARAMS);
        
        [nResults] = MaskQuant(n, ciliaMask, targetMask, targetI);
        
        iResults = vertcat(iResults, nResults);
    end
end

titles = {'Mean in cilia' 'Max in cilia' 'Mean around cilia' 'Max around cilia'};

output = [titles; num2cell(iResults)];

outputFilename = strcat(filename, '_results.csv');

writetable(cell2table(output(2:end,:),'VariableNames',output(1,:)), outputFilename);

close all
