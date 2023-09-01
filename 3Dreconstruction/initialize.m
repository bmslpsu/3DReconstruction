function [data, paths, folders, names] = initialize(root, fly)
%% initialize_data_paths: sets directories, loads DLT & background data, lets user set crop & tether points
fly = num2str(fly);
folders.root = fullfile(root, ['fly_' fly]);
folders.init = fullfile(folders.root, 'init'); % data directory
folders.mask = fullfile(folders.root, 'mask'); % mask directory
folders.background = fullfile(folders.root, 'background'); % background directory
folders.DLT = fullfile(folders.root, 'DLT'); % DLT directory
folders.reconstruction = fullfile(folders.root, 'reconstruction'); % reconstruction directory
folders.z_axis = fullfile(folders.reconstruction, 'z_axis'); % z-axis directory
folders.volume = fullfile(folders.reconstruction, 'volume'); % volume directory

names.DLT = 'DLT.csv';
paths.DLT = fullfile(folders.DLT, names.DLT);
data.DLT = csvread(paths.DLT); % read the DLT coeffs from the file

names.z_axis = 'z_axis.mat';
paths.z_axis = fullfile(folders.z_axis, names.z_axis);

names.volume = 'volume.mat';
paths.volume = fullfile(folders.volume, names.volume);

mkdir(folders.init)

% Select videos
fly = num2str(fly);
for n = 1:3 % each camera view
    [names.vid{n,1}, folders.vid{n,1}] = uigetfile(fullfile(root, ['fly_' fly '\vid\*.avi']), ...
        ['Select video from camera #' num2str(n)]);
    paths.vid{n,1} = fullfile(folders.vid{n}, names.vid{n}); 
end

% Select and load the backgrounds
se = strel('squar',4);
for n = 1:3 % each camera view
    [names.background{n,1}] = uigetfile({'*.tif','files'}, ...
        ['Select background for camera #' num2str(n)], folders.background, 'MultiSelect','off');
    paths.background{n,1} = fullfile(folders.background, names.background{n});
    
    try
        I = imread(paths.background{n});
        I = imcomplement(I);
        I = imdilate(I, se);
    catch
        I = [];
        warning('No background images detected. Analysis can still proceed')
    end
    data.background{n,1} = I;
    data.image_size{n,1} = size(I);
end

% Find initial tether points
n_vid = 3;
data.tether = cell(n_vid,1);
for n = 1:n_vid
    data.tether{n} = set_tether_points(paths.vid{n});
end

% Set crop reigons
for n = 1:3 % each camera view
    Vreader = VideoReader(paths.vid{n});
    frame = read(Vreader, 1);
    [~, rectout] = imcrop(frame);
    xy = rectout(3:4);
    xy = 2*floor(xy ./ 2) + 1;
    rectout(3:4) = xy;
    try
        crop_reigon = round(rectout);
    catch
        warning('not cropping images')
        crop_reigon = [];
    end
    data.crop{n,1} = crop_reigon;
    %Icrop = imcrop(frame, crop_reigon);
    %[Iout] = pad_image(Icrop, crop_reigon, frame);
end
close

paths.init = fullfile(folders.init, 'init.mat');
save(paths.init, 'data', 'folders', 'paths', 'names', '-v7.3')

end