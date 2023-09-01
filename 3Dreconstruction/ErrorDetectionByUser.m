function []=ErrorDetectionByUser(init,frames)

f1=figure;
for ii=frames
    %% load the wing and body reconstruction
    % Load the body reconstruction data
    reconpath = fullfile(init.folders.reconstruction, ['frame_' num2str(ii) '.mat']);
    load(reconpath, 'wingLeft');
    load(reconpath, 'wingRight');
    load(reconpath, 'body_xyz');
    load(reconpath, 'error_flag');
    if error_flag==0
        %% plot the data
        scatter3(wingRight(:,1),wingRight(:,2),wingRight(:,3),'r')
        hold on
        scatter3(wingLeft(:,1),wingLeft(:,2),wingLeft(:,3),'b')
        scatter3(body_xyz(:,1),body_xyz(:,2),body_xyz(:,3),'g')
        
        %% user input
        input_error=input(['Frame ' num2str(ii) 'Error/Imporovement in wing reconstruction y/n? '],'s');
        switch input_error
            case 'y'
                error_flag=1;
            case 'n'
                %does nothing
        end
    end
    clf(f1)
    save(reconpath,'error_flag','-append')
end