function []=Rotate3DVoxelsV3(init,sI,eI,plot_data)
%% due to calibration, the axes of the camera system is not always aligned with the global system
%this code checks and corrects for this sort of condition
%% INPUT:
%ReconstructionData: the body reconstruction structure with all the data
%wing_reconstruction_refined: the wing reconstruction structure (voxel and
%volume)
%% OUTPUT:
% Corrected body and wing in voxel and volume space

f1=figure;
for ii=sI:eI
    disp(['frame: ' num2str(ii)])
        %% axese and rotation matrix
    reconpath = fullfile(init.folders.reconstruction, ['frame_' num2str(ii) '.mat']);
    load(reconpath, 'x_axis');
    load(reconpath, 'y_axis');
    load(reconpath, 'z_axis');
    rotm=[-x_axis;  y_axis; z_axis];
    %% load the body and wing data
    load(reconpath, 'wing2_corrected');
    load(reconpath, 'wing1_corrected');
    load(reconpath, 'body_xyz');
    
    % Body rotation
    image_voxels_transformed=RotateHull(body_xyz,rotm);
    
    %wing rotation
    wing_voxels_transfromedRight=RotateHull(wing1_corrected,rotm); %extract the 3D image
    wing_voxels_transfromedLeft=RotateHull(wing2_corrected,rotm); %extract the 3D image
    %% reset fly position
    index_offset=-0; %the offset. If 0 then does nothing
    image_voxels_transformed(:,1)=image_voxels_transformed(:,1)-index_offset;
    image_voxels_transformed(:,2)=image_voxels_transformed(:,2)-index_offset;
    image_voxels_transformed(:,3)=image_voxels_transformed(:,3)-index_offset;
    
    wing_voxels_transfromedRight(:,1)=wing_voxels_transfromedRight(:,1)-index_offset;
    wing_voxels_transfromedRight(:,2)=wing_voxels_transfromedRight(:,2)-index_offset;
    wing_voxels_transfromedRight(:,3)=wing_voxels_transfromedRight(:,3)-index_offset;
    
    wing_voxels_transfromedLeft(:,1)=wing_voxels_transfromedLeft(:,1)-index_offset;
    wing_voxels_transfromedLeft(:,2)=wing_voxels_transfromedLeft(:,2)-index_offset;
    wing_voxels_transfromedLeft(:,3)=wing_voxels_transfromedLeft(:,3)-index_offset;
    
    %% plots for verification
    if plot_data
        scatter3(image_voxels_transformed(:,1),image_voxels_transformed(:,2),image_voxels_transformed(:,3))
        hold on
        scatter3(wing_voxels_transfromedRight(:,1),wing_voxels_transfromedRight(:,2),wing_voxels_transfromedRight(:,3))
        scatter3(wing_voxels_transfromedLeft(:,1),wing_voxels_transfromedLeft(:,2),wing_voxels_transfromedLeft(:,3))
        %plot the old data
        %     scatter3(ReconstructionData(ii).Voxel_Fly(:,1),ReconstructionData(ii).Voxel_Fly(:,2),ReconstructionData(ii).Voxel_Fly(:,3))
        %     scatter3(wing_reconstruction_refined(ii).WingVolume(:,1),wing_reconstruction_refined(ii).WingVolume(:,2),wing_reconstruction_refined(ii).WingVolume(:,3))
        center_old=mean(body_xyz);
        plot3([center_old(1) center_old(1)+.1*x_axis(1)],...
            [center_old(2) center_old(2)+.1*x_axis(2)],...
            [center_old(3) center_old(3)+.1*x_axis(3)],'r')
        plot3([center_old(1) center_old(1)+.1*y_axis(1)],...
            [center_old(2) center_old(2)+.1*y_axis(2)],...
            [center_old(3) center_old(3)+.1*y_axis(3)],'g')
        plot3([center_old(1) center_old(1)+.1*z_axis(1)],...
            [center_old(2) center_old(2)+.1*z_axis(2)],...
            [center_old(3) center_old(3)+.1*z_axis(3)],'k')
        drawnow
    end
    
    save(reconpath,'image_voxels_transformed','-append')
    save(reconpath,'wing_voxels_transfromedRight','-append')
    save(reconpath,'wing_voxels_transfromedLeft','-append')
end
close(f1)
end
function image_voxels_transformed=RotateHull(image_voxels,rotm)
image_voxels_transformed=zeros(size(image_voxels));
for jj=1:size(image_voxels,1)
    image_voxels_transformed(jj,:)=rotm*double(image_voxels(jj,:)');
end
% image_voxels_transformed=[round(image_voxels_transformed); ceil(image_voxels_transformed); floor(image_voxels_transformed)];
image_voxels_transformed=unique(image_voxels_transformed,'row');
end