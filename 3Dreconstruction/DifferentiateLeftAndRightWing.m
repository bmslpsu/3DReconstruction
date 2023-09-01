function [oneWingFlag]=DifferentiateLeftAndRightWing(init,frames)
% add a merged wing flag when it is difficult to sep the two wings
%seperates the two wings from the single wing reconstruction

%% NOTES:
%% directories
DLt_Coef = init.data.DLT;
Rec_Cord_UI = load(init.paths.volume, 'Rec_Cord');
%% Video readers
Vreader1 = VideoReader(init.paths.vid{1});
Vreader2 = VideoReader(init.paths.vid{2}); 
Vreader3 = VideoReader(init.paths.vid{3}); 
%% figures
figure_show=figure;
%% start clustering loop
oneWingFlag=zeros(frames(end),1); %Updated: serves as an error flag. This prompts a full manual reconstruction
for ii=frames
    %% load the wing and body reconstruction
    % Load the body reconstruction data
    disp(ii)
    reconpath = fullfile(init.folders.reconstruction, ['frame_' num2str(ii) '.mat']);
    load(reconpath, 'body_xyz');
    ReconstructionDataBody = body_xyz;
    load(reconpath, 'WingsCoordinates');
    load(reconpath, 'x_axis');
    %% implement k-mean clustering for the wings (also uses another algorithm to check it)
    [cluster1,~]=kmeans(double(WingsCoordinates),2); %kmeans produces 2 clusters
    cluster1_len=length(unique(cluster1));
    cluster2=clusterdata(double(WingsCoordinates),'maxclust',2); %this one produces at most 2
    cluster2_len=length(unique(cluster2));
    
    %% Place each wing in a certain variable without knowing which wing it is
    Bodycom=mean(ReconstructionDataBody); %find the center of mass of the body
    wing1_cluster=WingsCoordinates(cluster2==1,:);
    wing2_cluster=WingsCoordinates(cluster2==2,:);
    if length(wing1_cluster)<100 || length(wing2_cluster)<100 || cluster2_len ~=2
        oneWingFlag(ii)=1; %sets the flag to 1
        frame1 = read(Vreader1, ii);
        frame2 = read(Vreader2, ii);
        frame3 = read(Vreader3, ii);
        [wing1_cluster,wing2_cluster]=CorrectWingSegmentation(WingsCoordinates,DLt_Coef,frame1,frame2,frame3);
    end
    wing1_COM=mean(wing1_cluster);
    wing2_COM=mean(wing2_cluster);
    %% determine left or right wing
    vec_Wing1=wing1_COM-Bodycom; %vector from Body com to each wing point
    vec_Wing2=wing2_COM-Bodycom;
    dot_wing1=dot(vec_Wing1,x_axis); %dot product with the x-axis
    dot_wing2=dot(vec_Wing2,x_axis);
    
    %% Sep the wings into left and right and add to mat file
    if dot_wing1>0 &&dot_wing2<0 %this means both wings are easy to split
        wingLeft=wing1_cluster;
        wingRight= wing2_cluster;
        
    elseif dot_wing1<0 &&dot_wing2>0 %easy to split but wings switched in clustering (happens...)
        wingLeft=wing2_cluster;
        wingRight=wing1_cluster;
    else %will probably happen for damaged wings
        disp('this shouldnt happen')
        oneWingFlag(ii)=1; %sets the flag to 1
    end
    save(reconpath,'wingLeft','-append')
    save(reconpath,'wingRight','-append')
    %% scatter plot for verification purposes
    scatter3(wingLeft(:,1),wingLeft(:,2),wingLeft(:,3),'b')
    hold on
    scatter3(wingRight(:,1),wingRight(:,2),wingRight(:,3),'r')
    scatter3(ReconstructionDataBody(:,1),ReconstructionDataBody(:,2),ReconstructionDataBody(:,3),'g')
    hold off
    drawnow
end
close(figure_show)
end
%% Functions---------------
function [wing_coordinates1,wing_coordinates2]=CorrectWingSegmentation(wing_reconstruction,DLt_Coef,frame1,frame2,frame3)
%% a simple user interface that allows the user to refine the reconstruction of
%the wing based on various metrics that indicate a faulty wing
%reconstruction
disp('select an roi that has wing 1')
[~,~,~,wing_coordinates1]=UserInputReconstructionV3(wing_reconstruction,DLt_Coef,frame1,frame2,frame3);
disp('select an roi that has wing 2')
[~,~,~,wing_coordinates2]=UserInputReconstructionV3(wing_reconstruction,DLt_Coef,frame1,frame2,frame3);
end