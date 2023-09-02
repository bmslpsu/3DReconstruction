function [inst_freq] = find_stroke_reversalParFor(init,sI,eI, laser_activation_frame)
%% Locate the approximate location where the wing reversal occurs
maskdir = init.folders.mask;
%% finds the location of the wing mask in the image
parfor n = sI:eI
    % Extact the images
    maskdir_wing1=fullfile(maskdir, ['wing_' num2str(1)]);
    maskdir_wing_frame1=fullfile(maskdir_wing1,['frame_',num2str(n) '.png']);
    wing_frame1 = imread(maskdir_wing_frame1);
    [Wing, ~] = find(wing_frame1); %[row, col]=find(wingmask) returns the y-cor of the 1 values in the mask
    if isempty(Wing) %if no wing is detected set this index as naa
        COM(n) = nan;
    else
        COM(n) = mean(Wing); %if there is a wing add center to array
    end
end
%% We start putting data in COM at index 15. However that means 1 to 14 are zero or nan so we can't use those
%% Interpolate nan values (just in case), filter, & smooth the data
t_wing = sI:eI; %vec of time that starts at sI 
COM=COM(sI:eI); %COM is cropped so that elements before sI are gone

cutOffLow = 500; % cut-off frequency
cutOffHigh = 100; % cut-off frequency
[b,a] = butter(3, cutOffLow/(8000/2));
[bH,aH] = butter(3, cutOffHigh/(8000/2),'high');
% Interpolate missing points
COM_int = fillmissing(COM, 'spline', 'SamplePoints', t_wing);
COM_int=dtrend(COM_int,3);
% Filter the data
COM_int = filtfilt(b, a, COM_int);
COM_int = filtfilt(bH, aH, COM_int);
%% Locate the stroke reversal areas
del_time=1/8000; %frame time difference
[~,locs]=findpeaks(COM_int,'MinPeakDistance',15); %index of the stroke reversal
inst_freq=8000./gradient(locs); %gradient finds diff between reversals in index. 8000/num of frames=freq of stroke
time_locs=t_wing(locs); %frame location of each peak/reversal
plot(time_locs/8000,inst_freq);
%% save location 
freq_dir=fullfile(init.folders.root, 'wing frequency');
mkdir(freq_dir)
angle_savepath=fullfile(freq_dir, 'WingFreq.mat');
save(angle_savepath, 'inst_freq','time_locs','laser_activation_frame','COM_int') %saves filtered time data for more work
