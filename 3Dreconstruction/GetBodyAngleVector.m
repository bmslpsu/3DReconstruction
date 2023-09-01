function body_angle_all=GetBodyAngleVector(init,frames)
body_angle_all=zeros(1,length(frames));

for ii=frames
    reconpath = fullfile(init.folders.reconstruction, ['frame_' num2str(ii) '.mat']);
    load(reconpath, 'body_angle');
    body_angle_all(ii)=body_angle;
end

% angle_dir=fullfile(init.folders.root, 'BodyAngle');
% mkdir(angle_dir)
% angle_savepath=fullfile(angle_dir, 'BodyAngle.mat');
% save(angle_savepath, '-struct', 'angle')