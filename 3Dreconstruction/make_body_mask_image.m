function [] = make_body_mask_image(showplot, root_vid, sI, eI, bottomFlag, ...
    tether_pts, background, crop_reigon, savedir, label)
%% Generates the tether and body mask while subtracting the background
% INPUTS:
%   showplot    : generate plots to track what is going on when NOT running in batch mode
%   root_vid    : video path
%   sI          : start frame for analysis
%   eI          : end frame for analysis
%   bottomFlag  : false side views, true for bottom view
%   tether_pts  : location of the tether start and end points
%   background  : background image for each frame
%   crop_reigon : crop area
%   savedir     : directory to save masks
%   label       : trailing numeric label for mask subfolder (corresponds to camera view)
%
% OUTPUT:
%   -
%

Vreader = VideoReader(root_vid); % root video path
frame = read(Vreader, 1); % 1st frame
% dim = size(frame);
frame_window = 14; % # of frames to take before and after current frame for buffer
thresh = round( 0.1*(2*frame_window + 1) ); % any pixel that retains less than this much white pixels is considered part of the fly
if frame_window >= sI
    disp('Not enough frames to create initial mask.')
    disp('Please change start index of number of frames to use.')
    return
end
frames = sI:eI;
n_frame = length(frames);

label = char(string(label));
% maskpath = fullfile(savedir, ['mask_frames_' label '.mat']);

% Check directories
dircheck = exist(savedir, 'dir');
if ~dircheck
    disp('Making mask directory')
    mkdir(savedir)
end

%% Initialize image buffer
window = (sI - frame_window) : (sI + frame_window);
n_window = length(window);
[y_image , x_image] = size(frame); % the size of the image in pixels
buffer = false(y_image, x_image, n_window);
for n = 1:n_window % loop for each frame
    frame = read(Vreader, window(n));
    bw = process_image(frame, background); % clean up bw image
    buffer(:,:,n) = bw;
    %montage({frame, bw})
    %pause(0.2)
end

%% Generate Masks
disp('Generating masks for body')
if showplot
    fig = figure;
    set(fig, 'Units', 'inches', 'Position', [2 2 8 6])
    movegui(fig, 'center')
end

bodydir = fullfile(savedir, ['body_' label]);
tetherdir = fullfile(savedir, ['tether_' label]);

mkdir(bodydir)
mkdir(tetherdir)

bodypath = fullfile(bodydir, 'frame_');
tetherpath = fullfile(tetherdir, 'frame_');

se = strel('square',3);
image_type = 'png';
for n = 1:n_frame % each frame
    if ~mod(frames(n),10) || (frames(n) == sI) || (frames(n) == eI)
        disp(frames(n))
    end
    bw_mask = uint8(sum(buffer,3));
    mask_body_tether = bw_mask <= thresh; % create the combined body & tether mask by removing wing motion
    [mask_body, mask_tether] = seperate_body_tether(mask_body_tether, tether_pts, bottomFlag); % seperate body & tether masks
    
    % Clean masks
    mask_body = process_mask(mask_body, se);
    mask_tether = process_mask(mask_tether, []);
    
    if showplot
        disp_frame = read(Vreader, frames(n));
        if n ~= 1
            subplot(2,2,1) ; set(H.raw, 'CData', disp_frame)
            subplot(2,2,2) ; set(H.bw, 'CData', double(bw_mask))
            subplot(2,2,3) ; set(H.body, 'CData', mask_body)
            subplot(2,2,4) ; set(H.tether, 'CData', mask_tether)
        else
            subplot(2,2,1) ; H.raw = imshow(disp_frame);
            subplot(2,2,2) ; H.bw = imshow(double(bw_mask));
            subplot(2,2,3) ; H.body = imshow(mask_body);
            subplot(2,2,4) ; H.tether = imshow(mask_tether);
        end
        drawnow
    end
    
    mask_body_crop = imcrop(mask_body, crop_reigon);
    mask_tether_crop = imcrop(mask_tether, crop_reigon);
    
    if ~all( size(mask_body_crop) == size(mask_tether_crop) )
        warning('why?')
    end
    
    % Save the current mask for body & tether
    imwrite(mask_body_crop, [bodypath num2str(frames(n)) '.' image_type], image_type)
    imwrite(mask_tether_crop, [tetherpath num2str(frames(n)) '.' image_type], image_type)

    % Get the new frame for the buffer
    new_frame_idx = frames(n) + frame_window + 1; % new frame index
    new_frame = read(Vreader, new_frame_idx); % new frame
    buffer = circshift(buffer, -1, 3);  % shift the buffer back by 1 frame
    buffer(:,:,end) = process_image(new_frame, background); % replace last image in buffer with new frame
end
end

function [frame] = process_image(frame, background)
% process_image: clean binarized image 
%   
if length(background)<20
   background=0; 
end
frame = frame + background; % subtract background
frame = medfilt2(frame); % median filter to remove noise
frame = imadjust(frame); % adjust contrast
frame = imbinarize(frame); % binarize the image
frame = imcomplement(frame); % invert
%frame = bwmorph(frame, 'clean');
frame = bwareaopen(frame, 300); % remove small objects
frame = imfill(frame, 'holes'); % fill holes
frame = imcomplement(frame);  % invert
end

function [frame] = process_mask(frame, se)
% process_mask: clean mask image 
%  
frame = logical(frame);
frame = bwareafilt(frame, 2);
if ~isempty(se)
    frame = imdilate(frame, se);
end
frame = bwmorph(frame, 'thicken');
end
