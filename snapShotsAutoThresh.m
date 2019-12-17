%This function outputs two tif files: one contains snapshots for each of
%the true-cilia determined in the "quant_cilia_for_field.m" function and
%the other contains snapshots of what the signal from the POI channel
%looks like in each of the true-cilia. These files are for the user to
%verify that the pipeline is really identifying cilia and quantifying the
%correct POI signal.
function ciliaWindow = snapShots(R, L1,J1,V1, PARAMS)
ciliaWindow = [];
numShots = numel(R); %get number of possible cilium-objects
ciliaLogical = [R.positiveCilia];
% out of R, find positive cilia 
allCilia = sum(ciliaLogical); %get the number of true-cilia
if allCilia == 0 %check to prevent error if no cilia are found
    disp(['Zero cilia in field_' num2str(PARAMS.currentField, '%d')]);
    saveas(figure(1), [PARAMS.output_name '_field' num2str(PARAMS.currentField,'%d') '_cilia_chosen.tif']);
    saveas(figure(2), [PARAMS.output_name '_field' num2str(PARAMS.currentField,'%d') '_target_protein_.tif']);
    return;
end
%create vectors to hold protein data
POIarea = [];
global_density = [];
regional_mean_density = [];
IntDensity = [];

positiveCiliaCounter = 0;
%initialize counter
ciliaWindow(allCilia).box = [];
for M = 1 : numShots
    if R(M).positiveCilia == 1 
        positiveCiliaCounter = positiveCiliaCounter + 1; %add to counter
        %get the coordinates of the snapshot
        startx = floor(R(M).BoundingBox(1));
        starty = floor(R(M).BoundingBox(2));
        endx = startx + R(M).BoundingBox(3);
        endy = starty + R(M).BoundingBox(4);
        %a couple of checks to make sure pictures print out ok
        if startx == 0
            startx = 1;
        end
        if starty == 0
            starty = 1;
        end
        ciliaWindow(positiveCiliaCounter).box = (J1(starty:endy, startx:endx));
        fig4=figure(4);
        imshow(ciliaWindow(positiveCiliaCounter).box);
        fig1=figure(1);
        subplot(10, ceil(allCilia/10), positiveCiliaCounter)
        imshow(L1(starty:endy, startx:endx))
        fig2=figure(2);
        subplot(10, ceil(allCilia/10), positiveCiliaCounter)
        imshow(J1(starty:endy, startx:endx));
        fig3=figure(3);
        imshow(J1(starty:endy, startx:endx));
%        saveas(fig2, [PARAMS.output_name '_fieldN_single_protein_.tif']);
        %imwrite(J1(starty:endy, startx:endx),'myF.tif','tif');
        %Create an image(protein cut-out) "I" of the protein in the same box as the
        %identified cilia
        I = J1(starty:endy, startx:endx);%imread('smonullmef_smoWT_siAP2m1_shh_evc2stain_fieldN_single_protein_.tif');
        set(0,'defaultfigurecolor',[0 0 0])
        %Create binarized images of the various protein cut-outs with
        %different binarization threshholds to remove different levels of
        %noise
        BW = imbinarize(I,'global');
        imshow(BW)
        BW2 = bwperim(BW,18);
        BW3 = imbinarize(I,0.8);
        BW4 = imbinarize(I,0.6);
        %Create a matrix with the binarized image expressed as 1/0
        [a,b]= find(BW==1);
        %disp([a,b]);

        s = regionprops(BW4,I,{'Centroid','WeightedCentroid'});
        imshow(I)
        title('Weighted (red) and Unweighted (blue) Centroids'); 
        hold on
        numObj = numel(s);
        for k = 1 : numObj
            plot(s(k).WeightedCentroid(1), s(k).WeightedCentroid(2), 'r*')
            plot(s(k).Centroid(1), s(k).Centroid(2), 'bo')
        end
        hold off

        s = regionprops(BW4,I,{'Centroid','PixelValues','BoundingBox'});
        imshow(I)
        title('Mean by Region')
        hold on
        for k = 1 : numObj
            s(k).mean = mean(double(s(k).PixelValues));
            text(s(k).Centroid(1),s(k).Centroid(2), ...
                sprintf('%2.1f', s(k).mean), ...
                'EdgeColor','b','Color','r');
        end
        hold off
        

        pxarray = bwareaopen(BW,1);
        total = sum(pxarray,'all');
        disp("area:");
        disp(total);
        POIarea = [POIarea,total];
        
        AreaMeasured = boundary([a,b]);
        density = mean(I(:));
        disp("global density:")
        disp(density);
        global_density = [global_density,density];

        RegionalMean = mean([s.mean],'all');
        disp("regional mean:")
        disp(RegionalMean);
        regional_mean_density = [regional_mean_density,RegionalMean];
        
        ID = (RegionalMean*total);
        
        IntDensity = [IntDensity, ID];

        %figure
        %bar(1:numObj,[s.mean])
        %xlabel('Region Label Number')
        %ylabel('Standard Deviation')

        %figure
        %imshowpair(I,BW2,'montage');
        
    end
end

%save the figures
saveas(fig1, [PARAMS.output_name '_field' num2str(PARAMS.currentField,'%d') '_cilia_chosen.tif']);
saveas(fig2, [PARAMS.output_name '_field' num2str(PARAMS.currentField,'%d') '_target_protein_.tif']);
saveas(fig4, [PARAMS.output_name '_field' num2str(PARAMS.currentField,'%d') '_single_protein2_.tif']);
csvwrite(strcat(PARAMS.output_name,'_field', num2str(PARAMS.currentField,'%d'),'_protein_stats.csv'), [POIarea; global_density; regional_mean_density; IntDensity]);

close(fig1);
close(fig2);

end