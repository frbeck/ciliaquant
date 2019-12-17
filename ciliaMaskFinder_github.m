%PARAMETERS that can be altered based on the imaging
PARAMS.Threshold=170;%This is the threshold for considering a pixel as signal or noise in the cilia staining channel
PARAMS.MinMax=50; %This is the threshold for the minimum maximum intensity value for a potential "cilium" object to be considered as a true cilium
PARAMS.MinArea=100;%This is the threshold for the minimum area for a potential-cilium" object to be considered as a true cilium
PARAMS.MaxArea=50000;%Specify the order of the images in the Leica (.lif) file
PARAMS.dapiNum=1; %to specify the position of the DAPI channel
PARAMS.ciliaNum=3; %to specify the position of the cilia-marker channel
PARAMS.targetNum=2; %to specify the position of the channel for the POI
PARAMS.channelNum=3; %to specify the number of channels in the image
PARAMS.perimeterThres = 0; % Distance to search around each potential-cilium object for other objects
PARAMS.dilate = 0; %radius of SE used to dilate the cilia mask
PARAMS.max = false; %Whether to output normalized poi as max of z-plane values

PARAMS.numFields = 3; %to specify the number of fields captured in this Leica file
PARAMS.numMarkers= 3; %to specify the number of markers searched for; usually just 3 because one staining each for DAPI, cilia, POI
%specify the leica file to analyze
leicaFile = 'smonullmef_smoWT_siAP2m1_shh_evc2stain.lif';
%generate output file name
PARAMS.output_name=leicaFile(1:end-4);
%get the max-Z projections for each marker
[fields,zstacks] = get_max_z_projections_github(leicaFile, PARAMS);

%Begin quantifying cilia and POI in cilia

for i = 1:PARAMS.numFields %go through each field to quantify POI in each cilium
    PARAMS.currentField=i;
    
    series = zstacks{i,1};
    
    z = size(series);
    z = z(1);
    
    true_cilia_area = []; %vector to hold the area of each of cilia
    true_cilia_length = []; %vector to hold the length of each of cilia
    true_poi_cilia = []; %cell to hold ...
    true_poi_sub = []; %cell to hold ...
    true_cilia_intensity = []; %cell to hold ...
    normalized_poi_intensity = []; %cell to hold...
    
    %call to function to actually quantify
    [true_cilia_postn,cilia_channel,poi_channel,poi_intensities,true_cilia_indices,possb_cilia_info,cilia_masks]=quant_cilia_mask_for_field(fields{i,1}, PARAMS);
    %call to function to output snapshots of each chosen cilium and the POI
    %signal in each cilium
    ciliaWindow = snapShots(possb_cilia_info, true_cilia_postn,poi_channel,poi_intensities, PARAMS);

    
    
    %aggregate data from different fields
    for j = 1:numel(true_cilia_indices)
        index = true_cilia_indices(j);
        true_cilia_area(j) = possb_cilia_info(index).Area;
        true_cilia_length(j) = possb_cilia_info(index).MajorAxisLength;
        cilia_mask = cilia_masks == index;
        se = strel('disk',0);
        cilia_mask = imdilate(cilia_mask,se);
        se = strel('disk', 12);
        dilated_mask = imdilate(cilia_mask, se);
        sub_mask = dilated_mask - cilia_mask;
        
        cilia_int = zeros((z / PARAMS.channelNum),1);
        
        poi_cilia_int = zeros((z / PARAMS.channelNum),1);
        poi_sub_int = zeros((z / PARAMS.channelNum),1);
        
        l = 1;
        
        for k = 1:z
            if mod(k - PARAMS.ciliaNum * 1, PARAMS.channelNum) > 0
                continue
            end
            
            cilia_index = k;
            poi_index = k + (PARAMS.targetNum - PARAMS.ciliaNum);
            
            true_cilia_intensity(j,l) = mean(series{cilia_index,1}(find(cilia_mask > 0)));
            
            true_poi_cilia(j,l) = mean(series{poi_index,1}(find(cilia_mask > 0)));
            
            true_poi_sub(j,l) =  mean(series{poi_index,1}(find(sub_mask > 0)));
            
            normalized_poi_intensity(j,l) = true_poi_cilia(j,l) - true_poi_sub(j,l);
            
            l = l + 1;
        end
    end
    
    %disp(['Quantified ',num2str(length(true_poi_auc)), ' cilia']);
    %output statistics on each identified cilia
    
    if PARAMS.max == true
        normalized_poi_intensity= max(normalized_poi_intensity.');
    end
    
    csvwrite(strcat(PARAMS.output_name,'_field', num2str(i),'_z_stats.csv'), [true_poi_cilia; true_poi_sub; true_cilia_intensity]);
    csvwrite(strcat(PARAMS.output_name,'_field', num2str(i),'_cilia_stats.csv'), [true_cilia_area; true_cilia_length]);
    csvwrite(strcat(PARAMS.output_name,'_field',num2str(i),'_POI_backgrCorrect.csv'), [normalized_poi_intensity]);
    
end
