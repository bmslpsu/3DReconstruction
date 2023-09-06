%% run this code for a new video than run_main for hull reconstruction
%% Improved code for laser
%% Modifications:
%Saves all body angles to one vector in a sep folder
%Uses parallel processing for image analysis
%Uses par processing for reconstruction
%Detects flapping frequency from 2D image
%Estimates area ratio change in damaged-wing flies

%% Code starts here
clc
clear
close all

%% intialize folder paths for functions (maybe later)

%% Load the DLT coeffs 7 backgrounds, set video directories
% root = 'S:\Public\Wael\Chapter 1\DamagedWingExp\AnalyzedData\';
root = 'S:\Public\Lingsheng\Laser_RigidExp';
fly = 1000; %change fly number to match the number of fly you are analyzing 

%% initialize directory
% hint:copy the root before run this line
% choose videos from cam 1 to cam 3 from the root and fly number above. Than, choose background in the same order of the videos.
% drag from the center of the body to the end of the pin along the pin 3
% hit space when done. 3 times.
% crop the images and leave enough space for all yaw angles. double click
% the box when done. 3 times.
[~] = initialize(root,fly);

%% tracking start and end index
sI = 20; %start frame
eI = 8020; %end frame
laser_activation_frame=8000; % the exact frame in which the wing got hit by laser
init = load([root '\fly_' num2str(fly) '\init\init.mat']); %load intialization
load([root '\fly_' num2str(fly) '\init\init.mat']);
voxel_size = 5e-5; % the resolution of the reconstruction (voxel size)
% %% Alter the root directory match the current machine one (yay)
% init=AlterRootAndPaths(init,root,fly); %only for when the project was made on a different machine

%% Body masks (used to create body 3D reconstruction)
clc
tic
parfor n = 1:3 %parfor enables faster mask creation
    if n == 3
        bflag = true;
    else
        bflag = false;
    end
    make_body_mask_image(false, paths.vid{n}, sI, eI, bflag, data.tether{n}, ...
        data.background{n}, data.crop{n}, folders.mask, n); %function for body mask creation
end
toc

%% Wing masks
close all ; clc
tic
parfor n = 1:3
    make_wing_mask_image(true, paths.vid{n}, sI, eI, data.background{n}, data.crop{n}, folders.mask, n);
end
toc

%% Correct tether occlusion for wing masks
close all ; clc
tic
parfor n = 1:3
    correct_tether_occlusion_imageV2(init, true,n,sI,eI);
end
toc
%% find stroke reversal and flapping frequency before and after wing damage
%This is done using the 2d wing motion rather than 3d
clc; close all;
tic
inst_freq = find_stroke_reversalParFor(init,sI,eI,laser_activation_frame);
toc
%% 3D part and analysis of the code
%% --------------------------------------
%% Initialize voxel range and reconstruction resolution (reduces computational requirement)
close all ; clc
[~] = find_voxel_range(init, voxel_size);

%% 3D reconstruction portion of the code
%% consider turning this into a sep function
%find axis of rotation based on tether location and motion (might be best to use the entire trial for this one)
close all ; clc
tic
[~] = find_rotation_axis(init, sI, eI, true);
toc

%% Find the body 3D reconstruction and heading
close all ; clc
frames = sI:eI;
[body_angle] = reconstruct_bodyParfor(init, frames,true);
%% Wing 3D reconstruction
close all ; clc
frames = sI:sI+eI;
reconstruct_wingsV2(init, frames);

%% plotting function (Ben's work)
% clf
% clc
% frame_range = [1 1000];
% frame_step = 1;
% vidpath = 'E:\EXPERIMENTS\3DMAGNO\Experiment_SS\data\fly_3\test2.mp4';
% make_3d_display(init.paths.vid, init.data.background, frame_range, frame_step, [], 10)