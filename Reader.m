%function [f] = Reader(filename, thresh, i, n)
filename = '2_10_2020_tomD4_MBCD_extract_dox_shh_100uL_dish2_well4_1_MMStack_Pos0.ome.tif';

f = [];
g = [];
i = [];
h = [];

    
iInit = imread(filename, 400);
brightI = iInit*10;
figure ;
imshow(brightI);
roi = drawfreehand('Label', 'ROI');
h.FaceAlpha = 1;
h.FaceSelectable = false; 
info = imfinfo(filename);
numframe = length(info);
L = createMask(roi);

for j = 1 : numframe
    I = imread(filename, j);
    cellF = regionprops(L, I, 'MeanIntensity');
    f = [f, cellF];
end

F = f';

figure ;
imshow(brightI);
roi = drawfreehand('Label', 'ROI');
h.FaceAlpha = 1;
h.FaceSelectable = false; 
info = imfinfo(filename);
numframe = length(info);
L2 = createMask(roi);

for j = 1 : numframe
    I = imread(filename, j);
    cellF = regionprops(L2, I, 'MeanIntensity');
    g = [g, cellF];
end

G = g';

for j = 1 : numframe
    fnum = cell2mat(struct2cell(f(j)));
    gnum = cell2mat(struct2cell(g(j)));
    r = gnum./fnum;
    i = [i, r];
end

i = i./i(1);
i = i';
