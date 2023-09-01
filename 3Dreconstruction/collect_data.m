% Get list of all BMP files in this directory
% DIR returns as a structure array.  You will need to use () and . to get
% the file names.
clear
close all
clc

root = 'E:\EXPERIMENTS\3DMAGNO\Experiment_SS\data\fly_3\reconstruction';
files = dir([root '\*.mat']);
T = struct2table(files);
[~,sortI] = natsortfiles(T.name);
T = T(sortI,:);
n_file = length(files);

idx = nan(n_file,1);
all.frames = nan(n_file,1);
all.body_angle = nan(n_file,1);
for n = 1:n_file
    fname = fullfile(T.folder{n}, T.name{n});
    data = load(fname, 'body_angle', 'frame');
    all.frames(n) = data.frame;
    all.body_angle(n) = data.body_angle;
end

%%
collected_dir = fullfile(root, 'collected');
mkdir(collected_dir)
collected_fpath = fullfile(collected_dir, 'data.mat');
save(collected_fpath, 'all', '-v7.3')
