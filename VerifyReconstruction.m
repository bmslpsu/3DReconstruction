function []= VerifyReconstruction(init,sI,eI)
figure_check=figure;
for ii=sI:eI
    %% load the wing and body hulls
    disp(['frame: ' num2str(ii)])
    %% load 3D data
    reconpath = fullfile(init.folders.reconstruction, ['frame_' num2str(ii) '.mat']);
    load(reconpath, 'wing1_corrected');
    load(reconpath, 'wing2_corrected');
    load(reconpath, 'body_xyz');

    scatter3(body_xyz(:,1),body_xyz(:,2),body_xyz(:,3),'g')
    hold on
    scatter3(wing1_corrected(:,1),wing1_corrected(:,2),wing1_corrected(:,3),'b')
    scatter3(wing2_corrected(:,1),wing2_corrected(:,2),wing2_corrected(:,3),'r')
    pause()
    hold off
end