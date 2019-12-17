%This function takes in the path to the Leica file that holds the imaging
%data and information for the order of the channels and number of fields,
%separates the fields, does max-z projections for each channel for each field, and returns
%a cell-array that contains the max-z projections for each field.
function [fields,zstacks]=get_max_z_projections_github(leicaFile, PARAMS)

%add path to bio-formats folder
addpath('/Users/francisbeckert/Desktop/matlab folder/bfmatlab');
leicaFile = 'HAEvc2_cltcclone5a_shh.lif';
%read in the Leica imaging file
data = bfopen(leicaFile);
PARAMS.numFields = 3;
fields = {PARAMS.numFields,15}; %create a cell array to hold the max-z-projections from each field
zstacks = {};
%get all of the max-z-projections for each channel for each field
for(f = 1:PARAMS.numFields)
    
    fields{f,2} = ['field_' num2str(f,'%d')]; %put the description of the field in the cell array
    
    series = data{f,1}; %this gets the planes for field f, which includes the cilia imaging, target protein imaging, and dapi signal
    
    zstacks{f,1} = series;
    
    plane1 = series{1,1}; %get the first plane from the series
    plane_size = size(plane1); %get the size of the first plane
    
    %create a cell array to hold the max-z-projections. First column holds
    %the max-z-projections, second column holds a description of the
    %channel, third column holds the colormap for the max-z-projections
    PARAMS.numMarkers= 2;
    PARAMS.dapiNum=0; %to specify the position of the DAPI channel
    PARAMS.ciliaNum=2; %to specify the position of the cilia-marker channel
    PARAMS.targetNum=1; %to specify the position of the channel for the POI
    PARAMS.channelNum=2; %to specify the number of channels in the image
    marker_images = cell(PARAMS.numMarkers, 3);
    for(i = 1:PARAMS.numMarkers) %initialize the matrices for each marker's max-z-projection
        marker_images{i,1} = zeros(plane_size(1),'uint8');
        if(i == PARAMS.dapiNum) %also add in the name of the channels to the second column of the cell array
            marker_images{i,2} = 'DAPI_Channel';
        elseif(i == PARAMS.ciliaNum)
            marker_images{i,2} = 'Cilia_Channel';
        elseif(i == PARAMS.targetNum)
            marker_images{i,2} = 'Desired_Target_Protein_Channel';
        else
            marker_images{i,2} = 'Other_Target_Protein_Channel';
        end
        %get the colormap for this marker/channel
        marker_images{i,3} = data{f,3}{f,i};
        
    end

    %go through the series data and generate a max-z-projection for each marker
    seriesSize = size(series);
    for( i = 1:seriesSize(1) )
        marker_image_index = mod(i, PARAMS.numMarkers); %get the index for the matrix that will hold the z-projection values for this marker
        if(marker_image_index == 0) % whenever the result is zero, we want the index to be equal to PARAMS.numMarkers
            marker_image_index = PARAMS.numMarkers;
        end 
        current = marker_images{marker_image_index, 1}; %get the actual matrix for the current marker
        plane = series{i, 1}; %get the plane that we want to compare to our current 2-D array
        marker_images{marker_image_index, 1} = max(current, plane); %get the max values for each pixel between the two matrixes
    end
    disp(['done with max-z-projections for field_' num2str(f,'%d')] );

    %want to print out the individual images into a tif file if possible
    for(i =1:PARAMS.numMarkers)
        colorMap = marker_images{i,3};
        PARAMS.output_name=leicaFile(1:end-4);
        filename = [PARAMS.output_name '_field' num2str(f,'%d') '_' marker_images{i,2} '.tif'];
        imwrite(marker_images{i,1}, colorMap, filename);
    end
    fields{f,1} = marker_images; %put max-z-projections for field f in column one
    
end

disp('Done with max-z-projections for all fields');

end