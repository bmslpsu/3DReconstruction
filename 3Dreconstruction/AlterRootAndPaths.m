function init=AlterRootAndPaths(init,root_new,fly_number)

%% change the paths
init.paths.DLT=fullfile(root_new, ['fly_' num2str(fly_number)],'DLT','DLT.csv');
init.paths.z_axis=fullfile(root_new, ['fly_' num2str(fly_number)],'reconstruction','z_axis','z_axis.mat');
init.paths.volume=fullfile(root_new, ['fly_' num2str(fly_number)],'reconstruction','volume','volume.mat');



temp=strsplit(init.paths.vid{1},'\');
init.paths.vid{1}=fullfile(root_new,['fly_' num2str(fly_number)], 'Videos',temp{end-1},temp{end});
temp=strsplit(init.paths.vid{2},'\');
init.paths.vid{2}=fullfile(root_new,['fly_' num2str(fly_number)], 'Videos',temp{end-1},temp{end});
temp=strsplit(init.paths.vid{3},'\');
init.paths.vid{3}=fullfile(root_new,['fly_' num2str(fly_number)], 'Videos',temp{end-1},temp{end});

temp=strsplit(init.paths.background{1},'\');
init.paths.background{1}=fullfile(root_new,['fly_' num2str(fly_number)], 'background',temp{end});
temp=strsplit(init.paths.background{2},'\');
init.paths.background{2}=fullfile(root_new,['fly_' num2str(fly_number)], 'background',temp{end});
temp=strsplit(init.paths.background{3},'\');
init.paths.background{3}=fullfile(root_new,['fly_' num2str(fly_number)], 'background',temp{end});

init.paths.init=fullfile(root_new,'init','init.mat');

%% change folder location
init.folders.root=fullfile(root_new,['fly_' num2str(fly_number)]);
init.folders.init=fullfile(root_new,['fly_' num2str(fly_number)],'init');
init.folders.init=fullfile(root_new,['fly_' num2str(fly_number)],'init');
init.folders.mask=fullfile(root_new,['fly_' num2str(fly_number)],'mask');
init.folders.background=fullfile(root_new,['fly_' num2str(fly_number)],'background');
init.folders.DLT=fullfile(root_new,['fly_' num2str(fly_number)],'DLT');
init.folders.reconstruction=fullfile(root_new,['fly_' num2str(fly_number)],'reconstruction');
init.folders.z_axis=fullfile(root_new,['fly_' num2str(fly_number)],'z_axis');
init.folders.volume=fullfile(root_new,['fly_' num2str(fly_number)],'volume');

temp=strsplit(init.folders.vid{1},'\');
init.folders.vid{1}=fullfile(root_new,['fly_' num2str(fly_number)],'Videos',temp{end-1});
temp=strsplit(init.folders.vid{2},'\');
init.folders.vid{2}=fullfile(root_new,['fly_' num2str(fly_number)],'Videos',temp{end-1});
temp=strsplit(init.folders.vid{3},'\');
init.folders.vid{3}=fullfile(root_new,['fly_' num2str(fly_number)],'Videos',temp{end-1});