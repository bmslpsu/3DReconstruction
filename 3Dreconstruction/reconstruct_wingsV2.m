function [] = reconstruct_wingsV2(init, frames)
%% A newer faster version of the 3D reconstruction algorithm
% the reconstruction is coordinate based (Older version contain a voxel based one)
% more user interface to differentiate and crop wings
% different clustering algorithms to better sep the wings
% loop for error detection in case of large reconstructions

maskdir = init.folders.mask;
crop_reigon = init.data.crop; 
img_sz = init.data.image_size{1};
DLt_Coef = init.data.DLT;
Rec_Cord_UI = load(init.paths.volume, 'Rec_Cord');
Rec_Cord_UI = Rec_Cord_UI.Rec_Cord;

mask_path = maskReader_image(maskdir);
image_type = '.png';
MOV = cell(3,1);

% start reconstruction
disp('Generating reconstruction using specific masks')
n_frame = length(frames);
for n = 1:n_frame
    if ~mod(frames(n),10) || (frames(n) == frames(1)) || (frames(n) == frames(end))
        disp(frames(n))
    end
    
    % Load the wing masks for each view
    fI = num2str(frames(n));
    for f = 1:3
        Icrop = imread([mask_path(f).wing fI image_type]); % read mask image
        MOV{f} = pad_image(Icrop, crop_reigon{f}, img_sz); % pad so mask matches full image size
        Icrop_body = imread([mask_path(f).body fI image_type]); % read mask image
        MOV_body{f} = pad_image(Icrop_body, crop_reigon{f}, img_sz); % pad so mask matches full image size
    end
    
    % Load the body reconstruction data
    reconpath = fullfile(init.folders.reconstruction, ['frame_' fI '.mat']);
    load(reconpath, 'body_xyz');
    ReconstructionData = body_xyz;
    mov1_wing = MOV{1};
    mov2_wing = MOV{2};
    mov3_wing = MOV{3};  
    flag_eps = 0;
    eps = 1000; % initial minimum mask size to extract information
    while 1 % loop twice.  increasing eps to remove a view in the next iteration
        % initialization of loop variabls
        volume_cord = Rec_Cord_UI; %generates a volume grid
        
        % determine which views to use for the reconstruction
        view_keep = double( [sum(mov1_wing,'all') >= eps ,  sum(mov2_wing,'all') >= eps , sum(mov3_wing,'all') >= eps] );
        Pixel_Count = 3; % determines in how many images the pixel should be active
        % get the object mask. Overlay body if image doesnt have info
        mov1_wing = MOV{1}+MOV_body{1}*(1-view_keep(1));
        mov2_wing = MOV{2}+MOV_body{2}*(1-view_keep(2));
        mov3_wing = MOV{3}+MOV_body{3}*(1-view_keep(3));
        % generate the image that has the wing and body mask together
        % get the pixel value of each voxel
        [uv1] = round(dlt_inverse(DLt_Coef(:,1),volume_cord));
        [uv2] = round(dlt_inverse(DLt_Coef(:,2),volume_cord));
        [uv3] = round(dlt_inverse(DLt_Coef(:,3),volume_cord));
        
        % remove pixels that exceed the range of the image (active pixels
        % not detected here)
        index_remove1=uv1(:,1)<1 | uv1(:,1)>img_sz(2) | uv1(:,2)<1 | uv1(:,2)>img_sz(1);
        index_remove1=find(index_remove1>0);
        index_remove2=uv2(:,1)<1 | uv2(:,1)>img_sz(2) | uv2(:,2)<1 | uv2(:,2)>img_sz(1);
        index_remove2=find(index_remove2>0);
        index_remove3=uv3(:,1)<1 | uv3(:,1)>img_sz(2) | uv3(:,2)<1 | uv3(:,2)>img_sz(1);
        index_remove3=find(index_remove3>0);
        index_remove=unique([index_remove1 ; index_remove2 ; index_remove3]);
        uv1(index_remove,:)=[];
        uv2(index_remove,:)=[];
        uv3(index_remove,:)=[];
        volume_cord(index_remove,:)=[];
%         voxel_cord(index_remove,:)=[]; %till here both voxel and volume coordinates are the same

        % check if each pixel is 1 or 0
        ind1 = sub2ind(size(mov1_wing),uv1(:,2),uv1(:,1));
        ind2 = sub2ind(size(mov2_wing),uv2(:,2),uv2(:,1));
        ind3 = sub2ind(size(mov3_wing),uv3(:,2),uv3(:,1));
        pixel_value1=mov1_wing(ind1);
        pixel_value2=mov2_wing(ind2);
        pixel_value3=mov3_wing(ind3);
        
        %find the pixels that aren't part of the fly and wings
        SumPixels=pixel_value1+pixel_value2+pixel_value3;
        Index_keep=find(SumPixels>=Pixel_Count); %everything less than 3 is removed since it shouldnt exist in at least on image
        Rec_Cord=volume_cord(Index_keep,:);%remove the inactive voxels from coordinate frame
        %% Clustering and checking for reconstruction
        if ~isempty(Rec_Cord) %if the volume is empty then clustering won't work  
            try %possible error here if there are no points. Rare but happens
            % implement 3D volume filtering to remove reconstruction noise
            % clustering and cluster check
            ClusteringComponents2 = clusterdata(Rec_Cord,'Maxclust',2);
            Clustervalues2=unique(ClusteringComponents2);
            removed_index=[];
            for i=1:length(Clustervalues2)
                if sum(double(ClusteringComponents2==Clustervalues2(i)))<150
                    Rec_Cord(ClusteringComponents2==Clustervalues2(i),:)=[]; %gets rid of those points
                    ClusteringComponents2(ClusteringComponents2==Clustervalues2(i),:)=[];
                    removed_index=[removed_index i];
                end
            end
            Clustervalues2(removed_index)=[]; %need to check this line
            NumberOfClusters2=length(Clustervalues2);
            catch
                NumberOfClusters2=0;
            end
            % sep wings based on clustering
            if NumberOfClusters2==2
                wing1_coordinates=Rec_Cord(ClusteringComponents2==Clustervalues2(1),:);
                wing2_coordinates=Rec_Cord(ClusteringComponents2==Clustervalues2(2),:);
                
            elseif NumberOfClusters2>2 && length(Rec_Cord)>4500 %usually occurs due to tether occlusion
                warning('Too many clusters detected')
                disp('Select one wing and then the other')
                [~,~,~,wing1_coordinates]=UserInputReconstructionV2(Rec_Cord,DLt_Coef,mov1_wing,mov2_wing,mov3_wing);
                [~,~,~,wing2_coordinates]=UserInputReconstructionV2(Rec_Cord,DLt_Coef,mov1_wing,mov2_wing,mov3_wing);
                
            elseif NumberOfClusters2==1
                wing1_coordinates=Rec_Cord;
                wing2_coordinates=[0 0 0];
            else
                wing1_coordinates=[0 0 0];
                wing2_coordinates=[0 0 0];
            end
            
            % calculate wing dimensions and determine if reconstruction is the right size
            flag_size1=FindWingSizeFlag(wing1_coordinates,ReconstructionData);
            flag_size2=FindWingSizeFlag(wing2_coordinates,ReconstructionData);
        else
            flag_size1=1;
            flag_size2=1;
            wing1_coordinates=[0 0 0];
            wing2_coordinates=[0 0 0];
        end
     	% modify eps to remove the view with the least information
        if (flag_size1 || flag_size2) && flag_eps==0
            eps=min([sum(mov1_wing,'all') sum(mov2_wing,'all') sum(mov3_wing,'all')])+1; %makes eps larger than the smallest mask
            flag_eps=1;
            Rec_coordinate_wing1_Old=wing1_coordinates; %keeps the older images in case on of them is a good reconstruction to be used later
            Rec_coordinate_wing2_Old=wing2_coordinates;
            flag_size1_old=flag_size1;
            flag_size2_old=flag_size2;
        elseif flag_eps==1 && flag_size1 || flag_size2
            flag_eps=0;
            break
        else
            flag_eps=0;
            break
        end
    end
    % interesting issue solved here. wing ordering doesn't always remain the same after reomving a view
    % this might have to become a function in case there is an issue with
    % sep the bodies with different reconstructions
    % there is an issue where changes in size lead to change in order of
    % extracted objects. So function needs to be made
    if exist('Rec_coordinate_wing1_Old','var')
        
        [wing1_coordinates, wing2_coordinates]=CheckWhichWingV2(wing1_coordinates, wing2_coordinates, Rec_coordinate_wing1_Old, Rec_coordinate_wing2_Old);
       if length(wing1_coordinates)<=500
           flag_size1_old=0;
       end
       if length(wing2_coordinates)<=500
           flag_size2_old=0;
       end
        if ~flag_size1_old
            wing1_coordinates=Rec_coordinate_wing1_Old;
        end
        
        if ~flag_size2_old
            wing2_coordinates=Rec_coordinate_wing2_Old;
        end
    end

    % place in struct
    Rec_Cord_forPlot= [wing1_coordinates; wing2_coordinates]; %variable only used for ploting the scatter
    scatter3(Rec_Cord_forPlot(:,1),Rec_Cord_forPlot(:,2),Rec_Cord_forPlot(:,3),'r')
    hold on
    scatter3(ReconstructionData(:,1),ReconstructionData(:,2),...
        ReconstructionData(:,3))
    hold off
    drawnow
    WingsCoordinates=[wing1_coordinates; wing2_coordinates];
    save(reconpath,'WingsCoordinates','-append')
    clear Rec_Voxel_wing1_Old Rec_Voxel_wing2_Old Rec_Voxel3 Rec_Voxel_wing1 Rec_Voxel_wing2
    clear volume_cord voxel_cord wing2_coordinates wing1_coordinates
    clear Rec_coordinate_wing1_Old Rec_coordinate_wing2_Old
end
end

%% FUNCTIONS
function flag_size = FindWingSizeFlag(wing1_coordinates, ReconstructionData)
flag_size = 1;
if length(wing1_coordinates) <= 1200
    flag_size = 1;
else
    body_midpt=mean(ReconstructionData);
    wing_body_vec=wing1_coordinates-body_midpt;
    [~,wing_tipI]=max(vecnorm(wing_body_vec,2,2));
    [~,wing_bottomI]=min(vecnorm(wing_body_vec,2,2));
    span_vec=wing1_coordinates(wing_tipI,:)-wing1_coordinates(wing_bottomI,:);
    span_dist=norm(span_vec);
    
    wing_mid=mean(wing1_coordinates);
    vec_allPoints=wing1_coordinates-wing_mid; % vectors from midpt to all points
    dotProduct_all=dot(vec_allPoints,repmat(span_vec, length(vec_allPoints),1),2);
    [~,min_dotProductI]=min(abs(dotProduct_all));
    Normal_vectors=vec_allPoints(min_dotProductI,:);
    wing_chord=min(vecnorm(Normal_vectors,2,2));
%     scatter3(wing1_coordinates(:,1),wing1_coordinates(:,2),wing1_coordinates(:,3))
%     hold on
%     scatter3(wing_mid(:,1),wing_mid(:,2),wing_mid(:,3))
%     scatter3(wing1_coordinates(index_max,1),wing1_coordinates(index_max,2),wing1_coordinates(index_max,3))
    
    if span_dist<0.002/1.4 || wing_chord<0.00035/1.4
        flag_size=1;
    else
        flag_size=0;
    end
end
end