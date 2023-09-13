function []=ConvertToVoxelSpacev2(init,sI,eI,plot_data,voxel_size)

offset=300; %how much we offset the fly 
f1=figure;
for ii=sI:eI
    disp(['frame: ' num2str(ii)])
    reconpath = fullfile(init.folders.reconstruction, ['frame_' num2str(ii) '.mat']);
    load(reconpath, 'wing_voxels_transfromedRight');
    load(reconpath, 'wing_voxels_transfromedLeft');
    load(reconpath, 'image_voxels_transformed');
    
    %% coverts the wing into voxel space
    wingRV=round(wing_voxels_transfromedRight/voxel_size+offset);
    wingLV=round(wing_voxels_transfromedLeft/voxel_size+offset);
    BodyRecV=round(image_voxels_transformed/voxel_size+offset);
    %% save data
    save(reconpath,'wingRV','-append')
    save(reconpath,'wingLV','-append')
    save(reconpath,'BodyRecV','-append')
    if plot_data
        scatter3(BodyRecV(:,1),BodyRecV(:,2),BodyRecV(:,3),'g')
        hold on
        scatter3(wingRV(:,1),wingRV(:,2),wingRV(:,3),'r')
        scatter3(wingLV(:,1),wingLV(:,2),wingLV(:,3),'b')
        drawnow
    end
end
close(f1)

