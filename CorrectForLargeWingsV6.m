function []=CorrectForLargeWingsV6(init,frames,flag_stroke_reversal,voxel_size)
%% function corrects for large wings and can be used to correct for wings that weren't properly reconstructed
%% Initialize the directories of the folders
DLt_Coef = init.data.DLT;
%% Video readers
Vreader1 = VideoReader(init.paths.vid{1});
Vreader2 = VideoReader(init.paths.vid{2});
Vreader3 = VideoReader(init.paths.vid{3});

f1=figure;
%% size correction
for ii=frames
    disp(['Showing frame: ' num2str(ii)])
    %% read the frames
    image1 = read(Vreader1, ii);
    image2 = read(Vreader2, ii);
    image3 = read(Vreader3, ii);
    %% read the wing 3D coordinates
    reconpath = fullfile(init.folders.reconstruction, ['frame_' num2str(ii) '.mat']);
    load(reconpath, 'wingLeft');
    load(reconpath, 'wingRight');
    load(reconpath, 'body_xyz');
    load(reconpath, 'error_flag'); %determines if user will conduct wing reconstruction
    %% redo wing 1 and wing 2. When reconstruction cannot be recovered
    StillError='y';  %% added by mls frame 2687
    if error_flag || (length(wingLeft)<200 || length(wingRight)<200) %|| ii<=frames(1) %last part of statement refines the reconstruction at the first 3 frames
        disp('Select one wing then the next')
        %user reconstruction for horrible wings
        wing_1=UserWingGeneration(voxel_size,init,DLt_Coef,ii);
        wing_2=UserWingGeneration(voxel_size,init,DLt_Coef,ii);
        %the next lines correct wing orientation
        flag_stroke_reversal(ii)=1; %forces the code to refine the fit
        figure_GUI=figure;
        scatter3(wing_1(:,1),wing_1(:,2),wing_1(:,3),'r') %assume wing 1 is right
        hold on
        scatter3(wing_2(:,1),wing_2(:,2),wing_2(:,3),'b')
        scatter3(body_xyz(:,1),body_xyz(:,2),body_xyz(:,3),'g')
        x_switch=input('Flip wings? Blue should be left(y/n): ','s');
        StillError=input('Still error ?(y/n): ','s');  %% added by mls
        switch x_switch %correct which wing is right
            case 'n'
                wingRight=wing_1;
                wingLeft=wing_2;
            case 'y' %flips the wings
                wingRight=wing_2;
                wingLeft=wing_1;
        end
        close(figure_GUI)
    end
    %% correct wing 1 and 2
    wing_coordinates1=CorrectByUserV2(wingRight,body_xyz,DLt_Coef,image1,image2,image3,flag_stroke_reversal(ii),StillError);
    wing_coordinates2=CorrectByUserV2(wingLeft,body_xyz,DLt_Coef,image1,image2,image3,flag_stroke_reversal(ii),StillError);
    %% plot the data
    subplot(2,1,1)
    scatter3(wingRight(:,1),wingRight(:,2),wingRight(:,3))
    hold on
    scatter3(wingLeft(:,1),wingLeft(:,2),wingLeft(:,3))
    scatter3(body_xyz(:,1),body_xyz(:,2),body_xyz(:,3),'g')
    hold off
    subplot(2,1,2)
    if ~isempty(wing_coordinates1)
        scatter3(wing_coordinates1(:,1),wing_coordinates1(:,2),wing_coordinates1(:,3))
    else
        scatter3(wingRight(:,1),wingRight(:,2),wingRight(:,3),'r')
    end
    hold on
    if ~isempty(wing_coordinates2)
        scatter3(wing_coordinates2(:,1),wing_coordinates2(:,2),wing_coordinates2(:,3))
    else
        scatter3(wingLeft(:,1),wingLeft(:,2),wingLeft(:,3),'r')
    end
    scatter3(body_xyz(:,1),body_xyz(:,2),body_xyz(:,3),'g')
    %% correct the structure data
    if ~isempty(wing_coordinates1)
        wing1_corrected=wing_coordinates1;
    else
        wing1_corrected=wingRight;
    end
    if ~isempty(wing_coordinates2)
        wing2_corrected=wing_coordinates2;
    else
        wing2_corrected=wingLeft;
    end
    %% save and misc
    save(reconpath,'wing1_corrected','-append') %wing 1 is the right wing
    save(reconpath,'wing2_corrected','-append') %wing 2 is the left wing
    hold off
    clear wing_coordinates2 wing_coordinates1
end
close(f1)
end
%% FUNCTIONS ---------------------------------------------------------------------
function wing_corrected=UserWingGeneration(voxel_size,init,DLt_Coef,ii)
Vread = vidReader(init.paths.vid);
voxel_size_reduced = voxel_size*10;
space_range_initial = [-0.004 0.004]; % the initial range
volume_cord = CreateVoxelSpave(space_range_initial, voxel_size);
volume_cord_red = CreateVoxelSpave(space_range_initial, voxel_size_reduced);

uv = cell(3,1);
Mask_rect = cell(3,1);
index_remove = false;
for n = 1:3
    uv{n} = round(dlt_inverse(DLt_Coef(:,n),volume_cord)); % inverse DLT normal voxel resolution
    uv_red = round(dlt_inverse(DLt_Coef(:,n), volume_cord_red)); % inverse DLT reduced
    image_Original = read(Vread{n}, ii); % load first image
    Mask_rect{n} = FindMask(image_Original, uv_red); % find the ROI of the fly in each image
    image_res = size(Mask_rect{n}); % ROI size
    index_remove = index_remove | ...
        ( (uv{n}(:,1) < 1) | (uv{n}(:,1) > image_res(2)) | (uv{n}(:,2) < 1) | (uv{n}(:,2) > image_res(1)) );
end

volume_cord(index_remove,:) = [];
pixel_value = cell(3,1);
for n = 1:3
    uv{n}(index_remove,:) = [];
    ind = sub2ind(size(Mask_rect{n}), uv{n}(:,2), uv{n}(:,1));
    pixel_value{n} = Mask_rect{n}(ind);
end

% Find the pixels that aren't part of the fly and wings
SumPixels = sum(cat(2,pixel_value{:}),2);
Index_keep = SumPixels == 3; % everything less than 3 is removed since it shouldn't exist in at least one image
wing_corrected = volume_cord(Index_keep,:); % remove the false voxels from coordinate frame
end

function wing_coordinates=CorrectByUserV2(wing_coordinates,body_xyz,DLt_Coef,image1,image2,image3,FlagStrokeReversal,StillError)
%% a simple user interface that allows the user to refine the reconstruction of
%the wing based on various metrics that indicate a faulty wing
%reconstruction
%% code to determine possible issues in wing reconstruction and fix it
if ((FlagStrokeReversal && length(wing_coordinates)>100) || length(wing_coordinates)>200) && (StillError~='n')
    try
        disp('select an roi that has the wing leading edge and tip')
        [~,~,~,wing_coordinatesLeadingEdge1]=UserInputReconstructionV3(wing_coordinates,DLt_Coef,image1,image2,image3);
        disp('select an roi that has the wing trailing edge')
        [~,~,~,wing_coordinatesTrailingEdge1]=UserInputReconstructionV3(wing_coordinates,DLt_Coef,image1,image2,image3);
        min_length_sample=min(length(wing_coordinatesLeadingEdge1),length(wing_coordinatesTrailingEdge1));
        wing_coordinatesLeadingEdge1=datasample(wing_coordinatesLeadingEdge1,3*min_length_sample,1);
        wing_coordinatesTrailingEdge1=datasample(wing_coordinatesTrailingEdge1,3*min_length_sample,1);
        
        PointsToFit=[wing_coordinatesLeadingEdge1; wing_coordinatesTrailingEdge1];
        %fit a plane to the selected points
        fitobject = fit([PointsToFit(:,1),...
            PointsToFit(:,2)],PointsToFit(:,3),...
            'poly11');
        coeff=coeffvalues(fitobject);
        
        planeC0eff=[coeff(1) coeff(2) coeff(3) -1]; %[D A B C] Ax+By+Cz+D=0
        for i=1:length(wing_coordinates)
            DistancesToPlane(i)=abs(planeC0eff(2)*wing_coordinates(i,1)+planeC0eff(3)*wing_coordinates(i,2)...
                +planeC0eff(4)*wing_coordinates(i,3)+planeC0eff(1))/...
                sqrt(planeC0eff(2)^2+planeC0eff(3)^2+planeC0eff(4)^2);
        end
        wing_coordinates(DistancesToPlane>1.0e-04,:)=[];
    catch
    end
else
end
end

function Mask_rect = FindMask(image_Original, uv)
f1 = figure;
imshow(image_Original)
hold on
scatter(uv(:,1), uv(:,2))
Mask_rect = roipoly;
close(f1)
end

function voxel_cor = CreateVoxelSpave(space_range_initial, voxel_size)
x = space_range_initial(1):voxel_size:space_range_initial(2);
y = space_range_initial(1):voxel_size:space_range_initial(2);
z = space_range_initial(1):voxel_size:space_range_initial(2);
[xx,yy,zz] = meshgrid(x,y,z);
voxel_cor = [xx(:),yy(:),zz(:)]; % generates a volume grid
end