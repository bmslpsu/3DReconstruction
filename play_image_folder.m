% Get list of all BMP files in this directory
% DIR returns as a structure array.  You will need to use () and . to get
% the file names.
clear
close all
clc

root = 'E:\EXPERIMENTS\3DMAGNO\Experiment_SS\data\fly_3\mask\body_3';
imagefiles = dir([root '\*.png']);
T = struct2table(imagefiles);
[~,sortI] = natsortfiles(T.name);
T = T(sortI,:);
n_file = length(imagefiles);    % Number of files found
frame = imread(fullfile(T.folder{1}, T.name{1}));


idx = nan(n_file,1);
for n = 1:n_file
    fparts = strsplit(T.name{n}, {'_', '.'});
    idx(n) = str2double(fparts{2});
end

dx_frame = diff(idx);

gap = find(dx_frame > 1);

test = idx(gap)

%%
H = imshow(frame);
for n = 1:n_file
    if (n == 1) || (n == n_file) ||~mod(n,100)
        disp(n)
    end
    fpath = fullfile(T.folder{n}, T.name{n});
    frame = imread(fpath);
    set(H, 'CData', frame)
    drawnow
end