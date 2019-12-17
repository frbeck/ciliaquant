%This file contains the function that takes in a single field from a lif
%file, finds the cilia in the field based on user-defined
%parameters/thresholds, and quantifies the signal of the
%protein-of-interest, such as SMO or PTCH1
function [true_cilia_postn,cilia_channel,poi_channel,poi_intensities,true_cilia_indices,possb_cilia_info,cilia_masks]=quant_cilia_mask_for_field(field,PARAMS)


cilia_channel=field{PARAMS.ciliaNum,1}; %get max-z projection values for the cilia staining

poi_channel=field{PARAMS.targetNum,1}; %get max-z projection values for the protein-of-interest (POI) staining

%get rid of salt and pepper noise in the first channel and then create a
%black and white mask by thresholding
M=(medfilt2(cilia_channel,[3 3])>PARAMS.Threshold); %identify pixels with real signal and store them in M

%save the cilia mask as tif file for troubleshooting
imwrite(M, [PARAMS.output_name '_field' num2str(PARAMS.currentField,'%d') '_cilia_mask.tif']);

%identify and number each object that could be a cilium in the mask.
%Objects are identified based on adjacent pixels that show some signal in
%the mask.
% if (PARAMS.dilate > 0)
%     SE = strel('disk', PARAMS.dilate);
%     Mdil = imdilate(M,SE);
    
L=bwlabel(M);
cilia_masks = L;
%Ldil=bwlabel(Mdil);

%define the features/thresholds of objects that are cilia-like
MinArea=PARAMS.MinArea;
MaxArea=PARAMS.MaxArea;
%eccentricity represents an object's deviation from being a circle
MaxEccentricity=.9999;
MinEccentricity=.7;
%how filled in is the object?
MinSolidity=.6;
%how bright is the brightest pixel?
minMax=PARAMS.MinMax;
%calculate background by taking an average over 100x100 square
b=imfilter(poi_channel,ones(100)/10000,'symmetric'); %used later to reduce noise from the Smo staining to eliminate pixel values which are unrepresentative of their surroundings.

%calculate properties for each object

possb_cilia_info=regionprops(L,'Eccentricity','Area','Solidity','PixelIdxList', 'BoundingBox', 'MajorAxisLength');
%possb_cilia_info_dil=regionprops(Ldil,'Eccentricity','Area','Solidity','PixelIdxList', 'BoundingBox', 'MajorAxisLength');


%Go through the objects and determine if each is a cilium based on its
%properties. If the object is a cilium, calculate the mean fluorescence
%from the protein-of-interest channel and subtract the background.
true_cilia_postn=L; %matrix that tracks the positions of the true cilia
poi_intensities = []; %vector to hold the intensities of the POI for each determined cilium
for i=1:length(possb_cilia_info)
    possb_cilia_info(i).positiveCilia = 1; %initialize the positiveCilia feature for the object. This meant to track if the object is a true cilium 
    maxI=max(cilia_channel(possb_cilia_info(i).PixelIdxList)); %get pixel with max intensity from cilia staining for this possible "cilia"
    possb_cilia_info(i).maxI = maxI;
    %first check if any other objects are nearby to weed out any
    %"cilium-like" objects created from noisy acetylated-tubulin staining.
    boundingBox = possb_cilia_info(i).BoundingBox;
    %calculate the x coordinate for the upper left corner of expanded bounding box
    upperLeftX = boundingBox(1) - PARAMS.perimeterThres;
    if upperLeftX < 1
        upperLeftX = 1;
    end
    %calculate the y coordinate for the upper left corner of expanded bounding box
    upperLeftY = boundingBox(2) - PARAMS.perimeterThres;
    if upperLeftY < 1
        upperLeftY = 1;
    end
    %calculate the new width of the bounding box
    matSize = size(true_cilia_postn);
    newWidth = boundingBox(3) + (2*PARAMS.perimeterThres);
    if (upperLeftX + newWidth) > matSize(1)
        newWidth = matSize(1) - upperLeftX;
    end
    %calculate the new length of the bounding box
    newLength = boundingBox(4)+ (2*PARAMS.perimeterThres);
    if (upperLeftY + newLength) > matSize(2)
        newLength = matSize(2) - upperLeftY;
    end
    %determine the values needed for the expanded bounding box
    exp_startx = floor(upperLeftX);
    exp_starty = floor(upperLeftY);
    exp_endx = floor(upperLeftX+newWidth);
    exp_endy = floor(upperLeftY+newLength);
    expandedBox = true_cilia_postn(exp_starty:exp_endy, exp_startx:exp_endx);
    elements = numel(unique(expandedBox)); %get the number of unique values in this expanded box
    %imshow(expandedBox);
%     if elements > 2 %if there is more than 2 unique values in the expanded box, throw out this object because that means there are 2 or more objects in this box
%         true_cilia_postn(possb_cilia_info(i).PixelIdxList)=0; %change this "object's" mapping number to 0
%         possb_cilia_info(i).positiveCilia = 0; %set cilia marker to 0/FALSE
%     %check if this "cilium" posseses the feature/meets the thresholds for
    %being labeled a "cilium"
    if maxI<minMax | possb_cilia_info(i).Area<MinArea | possb_cilia_info(i).Area>MaxArea | possb_cilia_info(i).Eccentricity>MaxEccentricity | possb_cilia_info(i).Eccentricity<MinEccentricity | possb_cilia_info(i).Solidity<MinSolidity
        true_cilia_postn(possb_cilia_info(i).PixelIdxList)=0; %if not, change this "cilia's" mapping number to 0
        possb_cilia_info(i).positiveCilia = 0; %set cilia marker to 0/FALSE
    else
        %this is a true cilium, so get mean POI staining intensity and
        %subtract the mean background
        %poi_intensities(i)=mean(poi_channel(possb_cilia_info_dil(i).PixelIdxList))-mean(b(possb_cilia_info_dil(i).PixelIdxList));
        poi_intensities(i)=mean(poi_channel(possb_cilia_info(i).PixelIdxList))-mean(b(possb_cilia_info(i).PixelIdxList));
    end
end

%find all of the non-zero elements in poi_intensities, which means find the
%indices for the true-cilium objects
true_cilia_indices=find(poi_intensities); 
%trim poi_intensities so it only holds the mean intensities for the
%true-cilium objects
poi_intensities=poi_intensities(true_cilia_indices); 

display(['done quantifying cilia for field_' num2str(PARAMS.currentField,'%d')]);
end 

