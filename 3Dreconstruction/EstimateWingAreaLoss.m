function [Area_ratio]=EstimateWingAreaLoss(init,frames,voxel_size, damage_flag)
%% This code estimates the area change or area ratio of the damaged wing following wing damage

% damage_flag is set as 1 for damaged wing flies and 0 for intact wing
% flies


%% Initialize the directories of the folders
DLt_Coef = init.data.DLT;
%% Video readers
Vreader1 = VideoReader(init.paths.vid{1});
Vreader2 = VideoReader(init.paths.vid{2});
Vreader3 = VideoReader(init.paths.vid{3});

f1=figure;
sI=frames(1);
eI=frames(end);
%% load frame until user is happy with which frames to use for reconstruction

frame_measure=[];
for ii=[sI eI]
    frame_start=ii; % the frame a which the code starts checking
    while 1
        %% read the frames
        image1 = read(Vreader1, frame_start);
        image2 = read(Vreader2, frame_start);
        image3 = read(Vreader3, frame_start);
        
        montage({image1,image2,image3}) % show the images
        disp('Please select an option')
        disp('Advance one frame: a')
        disp('Go back one frame: b')
        disp('Jump 10 frames forward: 10a')
        disp('Jump 10 frames back : 10b')
        decision=input('Accept: y',"s");
        
        switch decision
            case 'a'
                frame_start=frame_start+1;
            case 'b'
                frame_start=frame_start-1;
            case '10a'
                frame_start=frame_start+10;
            case '10b'
                frame_start=frame_start-10;
            case 'y' %selects the frame and ends the while loop
                frame_measure=[frame_start frame_measure];
                break
            otherwise
                disp('Input not valid')
                why
        end
        if frame_start<sI % makes sure the user does not go beyond the range of the data
            frame_start=sI;
            disp('Index beyond lower range of data, frame index reset to start index')
        elseif frame_start>eI
            frame_start=eI;
            disp('Index beyond upper range of data, frame index reset to end index')
        end
    end
end

%% begin manual reconstruction and correction by the user

wing_1=UserWingGeneration(voxel_size,init,DLt_Coef,frame_measure(1));
wing_2=UserWingGeneration(voxel_size,init,DLt_Coef,frame_measure(2));
%% correct the two wings so that area estimations become easier
wing_coordinates1=CorrectByUserV2(wing_1,DLt_Coef,image1,image2,image3,1);
wing_coordinates2=CorrectByUserV2(wing_2,DLt_Coef,image1,image2,image3,1);
%% project wing unto the plane that it was fit to. Basically turns it into a thin plate
end

%% FUNCTIONs-------------------------------
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

function wing_coordinates=CorrectByUserV2(wing_coordinates,DLt_Coef,image1,image2,image3,FlagStrokeReversal)
%% a simple user interface that allows the user to refine the reconstruction of
%the wing based on various metrics that indicate a faulty wing
%reconstruction
%% code to determine possible issues in wing reconstruction and fix it
if (FlagStrokeReversal && length(wing_coordinates)>1200) || length(wing_coordinates)>5000
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
        wing_coordinates(DistancesToPlane>1.2e-04,:)=[];
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