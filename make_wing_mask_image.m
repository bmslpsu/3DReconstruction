function [] = make_wing_mask_image(showplot, root_vid, sI, eI, background, crop_reigon, maskpath, label)
%% Finds the wing mask by extracting the body and tethers from the image

disp('Creating wing mask')

frames = sI:eI;
n_frame = length(frames);

Vread = VideoReader(root_vid); % reads raw video

label = char(string(label));

bodypath = fullfile(maskpath, ['body_' label], 'frame_');
tetherpath = fullfile(maskpath, ['tether_' label], 'frame_');

wingdir = fullfile(maskpath, ['wing_' label]);
mkdir(wingdir)
wingpath = fullfile(wingdir, 'frame_');

fig_thresh = figure;
se = strel('rectangle',[10 10]);
image_type = 'png';
for n = 1:n_frame
    if ~mod(frames(n),10) || (frames(n) == sI) || (frames(n) == eI)
        disp(frames(n))
    end
    fI = num2str(frames(n));
    
    frame = read( Vread, frames(n));
    
    mask_body = imread([bodypath fI '.' image_type]);
    mask_tether = imread([tetherpath fI '.' image_type]);
    if length(background)<10
        background=0;
    end
    frame = frame + background; % subtract the background
    
    if ~isempty(crop_reigon)
        frame = imcrop(frame, crop_reigon);
    end
    
    frame = medfilt2(frame); % median filter to remove noise from background subtraction and tether removal
    
    frame(mask_body) = 255; % removes the body and tethers from the system
    frame(mask_tether) = 255;
    
    thresh = 2.5*graythresh(frame); % is selected randomly. Might not be even needed
    if thresh > 1 % might have to be manually tuned. Sometimes the automatic thresholding is just too high and converts certain part of the background to foreground
        thresh = 0.75;
    end
        
    image_wing = imbinarize(frame, thresh); % compliments image then thresholds
    image_wing = imcomplement(image_wing);
    image_wing = bwareaopen(image_wing, 200); % remove pixels smaller than a certain size
    image_wing = imfill(image_wing, 'holes'); % closes the holes
    image_wing = bwareafilt(image_wing, 4); % the filtered image of the wing (keep largest 4 objects. 4 in case wing is obstructed by tether)
    image_wing = imfill(image_wing, 'holes'); % closes the holes
    %image_wing = imclose(image_wing, se);
    
    imwrite(image_wing, [wingpath num2str(frames(n)) '.' image_type], image_type)

    if showplot
        if n ~=1
            subplot(2,1,1) ; set(H.raw, 'CData', frame)
            subplot(2,1,2) ; set(H.mask, 'CData', image_wing)
        else
            subplot(2,1,1) ; H.raw = imshow(frame);
            subplot(2,1,2) ; H.mask = imshow(image_wing);
        end
        drawnow
    end
end
close(fig_thresh)
end