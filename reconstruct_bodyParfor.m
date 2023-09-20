function [body_angle] = reconstruct_bodyParfor(init, frames,plot_flag)
%% reconstruct_body: reconstructs the fly body in the 3D plan and uses it to determine the heading of the fly
% Faster reconstruction using a predefined smaller volume
% Uses PCA for finding the fly heading
% OUPUT:
% ReconstructionData_space: The reconstruction data in metric units
% ReconstructionData_Voxel: The reconstruction in voxel coordinates
% load and analyze video
% determines changes in heading and the orthogonal system for each
% reconstruction
if plot_flag
    f1_body_plot=figure;
end
maskdir = init.folders.mask;
crop_reigon = init.data.crop; 
img_sz = init.data.image_size;
DLt_Coef = init.data.DLT;
load(init.paths.volume, 'Rec_Cord')
load(init.paths.z_axis, 'z_axis')

recondir = fullfile(init.folders.root, 'reconstruction');
mkdir(recondir)

mask_path = maskReader_image(maskdir);
image_type = '.png';
MOV = cell(3,1);

tic
sI = frames(1);
n_frame = length(frames); %length of the frame vector to iterate each frame
body_angle=zeros(1,length(frames));

%% check if exist files
rec_frame_files = dir(fullfile(recondir, '*.mat'));  % list all .mat file
fileNamesList = {};
for i = 1:length(rec_frame_files)
    fileName = rec_frame_files(i).name;
    fileNamesList{end+1} = fileName;
end
ask_user = 0;
for ii = 2:n_frame
    fI = num2str(frames(ii));
    if ~ismember(['frame_' fI '.mat'], fileNamesList)
        only_calculate_body_angle = 'n';
    else
        ask_user = 1;
    end
end
%% ask if file exists
if ask_user
    while true
        fprintf('Do you only need to calculate body angle?');
        only_calculate_body_angle = input('[y/n]:', 's');
        if only_calculate_body_angle == 'n' || only_calculate_body_angle == 'N'|| only_calculate_body_angle == 'y'|| only_calculate_body_angle == 'Y'
            break;
        else
            warning(['The input is neither ', 'y ', 'nor ', 'n ', 'please re-enter.\n']);
        end
    end
end

%% initialize first frame as the baseline
recon_initial=IntializeHeadingAndAngle(frames, mask_path, image_type, crop_reigon, img_sz,...
    DLt_Coef, Rec_Cord, z_axis, recondir, only_calculate_body_angle);
body_angle(1)=0; %initialize the first body angle

%% continue to rest of frame
parfor ii = 2:n_frame
    frame = frames(ii);
    % Load the tether masks for each view
    fI = num2str(frames(ii));
    Icrop = imread([mask_path(1).body fI image_type]); % read mask image
    MOV1 = pad_image(Icrop, crop_reigon{1}, img_sz{1}); % pad so mask matches full image size
    Icrop = imread([mask_path(2).body fI image_type]); % read mask image
    MOV2 = pad_image(Icrop, crop_reigon{2}, img_sz{2});
    Icrop = imread([mask_path(3).body fI image_type]); % read mask image
    MOV3 = pad_image(Icrop, crop_reigon{3}, img_sz{3});
    % Recontruct the body from all 3 views
    body_xyz = reconstruct_from_imagesParFor(MOV1,MOV2,MOV3, DLt_Coef, Rec_Cord, false);
     
    % fit 3D line to fly
    body_vector = pca(body_xyz); % sed to find the fly heading
    body_vector = body_vector(:,1); % take the strongest component
    body_vector = body_vector/norm(body_vector);
    if  dot(body_vector, z_axis) < 0 %fly heading.z-axis should always be positive
        body_vector =- body_vector;
    end
    % y-axis and x-axis calculation + save axes into a structure
    x_axis = cross(body_vector, z_axis);
    x_axis = x_axis / norm(x_axis);
    y_axis = cross(z_axis,x_axis);

    %%use the vec multiplication method. Gives the wrapped angle
    % det_vec=recon_initial.x_axis(1)*x_axis(2)+x_axis(1)*recon_initial.x_axis(2);
    % dot_vec=dot(recon_initial.x_axis(1:2), x_axis(1:2));
    % body_angle(ii)=atan2(det_vec,dot_vec);

    %% another solution
    cross_prod=cross(recon_initial.x_axis,x_axis);
    body_angle(ii)=atan2(dot(cross_prod,z_axis), dot(recon_initial.x_axis,x_axis));
 	%% Save
    savepath = fullfile(recondir, ['frame_' fI '.mat']);  
    if only_calculate_body_angle == 'n' || only_calculate_body_angle == 'N'
        parsave(savepath, frame,body_xyz,body_vector,x_axis,y_axis,z_axis,body_angle(ii))
    end
end

if plot_flag
    close(f1_body_plot)
end
toc
angle_dir=fullfile(init.folders.root, 'body_angle_all');
mkdir(angle_dir)
angle_savepath=fullfile(angle_dir, 'BodyAngle.mat');
if exist(angle_savepath, 'file')
    while true
        fprintf('The file: %s exists, do you want to overwrite? ', angle_savepath);
        ask_overwrite = input('[y/n]:', 's');
        if ask_overwrite == 'n' || ask_overwrite == 'N'
            original_sI = input('The original sI (enter an integer): ');
            if frames(1) < original_sI
                warning('The first frame: %d is smaller than the original sI, please check.\n', frames(1));
                return
            end
            original_eI = input('The original eI (enter an integer): ');
            bodyangle_file = load(angle_savepath);
            if frames(1)>= original_sI && frames(end) <= original_eI % incase there may be more functions, we check the range again.
                for ii = 2:length(frames)
                    bodyangle_file.body_angle(frames(ii)-original_sI+1) = body_angle(ii);
                    body_angle = bodyangle_file.body_angle;
                    save(angle_savepath, 'body_angle');
                end
            else
                warning('The first frame: %d and the last frame: %d are not within the original sI and eI, please check.\n', frames(1), frames(end));
                return
            end
            break
        elseif ask_overwrite == 'y' || ask_overwrite == 'Y'
            fprintf('overwrite: %s\n', angle_savepath);
            save(angle_savepath, 'body_angle');
            save(fullfile(angle_dir, 'angle_frames.mat'), 'frames')
            break
        else
            warning(['The input is neither ', 'y ', 'nor ', 'n ', 'please re-enter.\n']);
        end
    end
else
    save(angle_savepath, 'body_angle');
    save(fullfile(angle_dir, 'angle_frames.mat'), 'frames')
end
end
%% FUNCTIONS-----------------
function parsave(savepath, frame,body_xyz,body_vector,x_axis,y_axis,z_axis,body_angle)
  save(savepath, 'frame','body_xyz','body_vector','x_axis','y_axis','z_axis','body_angle' )
end

function recon_initial=IntializeHeadingAndAngle(frames,mask_path,image_type,crop_reigon,img_sz,...
    DLt_Coef, Rec_Cord, z_axis, recondir, only_calculate_body_angle)

recon_initial.frame = frames(1);

fI = num2str(frames(1));
for f = 1:3
    Icrop = imread([mask_path(f).body fI image_type]); % read mask image
    MOV{f} = pad_image(Icrop, crop_reigon{f}, img_sz{f}); % pad so mask matches full image size
end

% Recontruct the body from all 3 views
coordinate_fly = reconstruct_from_images(MOV, DLt_Coef, Rec_Cord, false);
recon_initial.body_xyz = coordinate_fly;

% fit 3D line to fly
body_vector = pca(coordinate_fly); % sed to find the fly heading
body_vector = body_vector(:,1); % take the strongest component
body_vector = body_vector/norm(body_vector); %normalize the body vector
recon_initial.body_vector = body_vector; % checks if the body vector has flipped from one frame to the next
%% User interface to determine if heading of fit vector is correct
recon_initial.body_vector=recon_initial.body_vector*sign(dot(z_axis,body_vector)); %makes sure the body vector is in the +ve z-axis
f1_body=figure;
scatter3(coordinate_fly(:,1),coordinate_fly(:,2),coordinate_fly(:,3))
hold on
body_center=mean(coordinate_fly);
body_vector_red=recon_initial.body_vector/200;
plot3([body_center(1) body_center(1)+body_vector_red(1)],...
    [body_center(2) body_center(2)+body_vector_red(2)]...
    ,[body_center(3) body_center(3)+body_vector_red(3)],'k')
flip_axis=input('Flip body-axis y/n? ','s');
switch flip_axis
    case 'y'
        recon_initial.body_vector=-recon_initial.body_vector;
    case 'n'
        %Do nothing
end

close(f1_body)
clear body_center body_vector_red
%% rest of the axes
% y-axis and x-axis calculation + save axes into a structure
x_axis = cross(recon_initial.body_vector, z_axis);
x_axis = x_axis / norm(x_axis);
y_axis = cross(z_axis,x_axis);
recon_initial.x_axis = x_axis;
recon_initial.y_axis = y_axis;
recon_initial.z_axis = z_axis;

det_vec=recon_initial.x_axis(1)*x_axis(2)+x_axis(1)*recon_initial.x_axis(2);
dot_vec=dot(recon_initial.x_axis(1:2), x_axis(1:2));
body_angle=atan2(det_vec,dot_vec);
recon_initial.body_angle = body_angle;

%% Save
if only_calculate_body_angle == 'n' || only_calculate_body_angle == 'N'
    savepath = fullfile(recondir, ['frame_' fI '.mat']);
    save(savepath, '-struct', 'recon_initial')
end
end