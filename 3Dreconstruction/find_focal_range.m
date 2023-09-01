function [pixel_fly,coordinate_fly,voxel_fly] = find_focal_range(range,delta_range,mov1,mov2,mov3, DLt_Coef)
%% Output
% pixel_fly: pixels in each image from each camera that contain the fly
% coordinate_fly: coordinate of each reconstructed pixel in 3D

% Initialize pixels and location
pixel_fly = [];
coordinate_fly = [];
voxel_fly = [];
dim = size(mov1);

% inverse to find associated uv
k = 1;
win = range(1):delta_range:range(2);
n_win = length(win);
coordinate_fly = zeros(10,3);
for x = win
    for y = win
        for z = win
            % Convert cooridnates to 
            [uv1] = dlt_inverse(DLt_Coef(:,1),[x y z]);
            [uv2] = dlt_inverse(DLt_Coef(:,2),[x y z]);
            [uv3] = dlt_inverse(DLt_Coef(:,3),[x y z]);
            uv1 = round(uv1); % round to the nearest pixel value
            uv2 = round(uv2);
            uv3 = round(uv3);
            
            % find the pixels within range
            if (uv1(1) > 0) && (uv1(2) > 0) && (uv1(1) <= dim(2)) && (uv1(2) <= dim(1))
                if (uv2(1) > 0) && (uv2(2) > 0) && (uv2(1) <= dim(2)) && (uv2(2) <= dim(1))
                    if (uv3(1) > 0) && (uv3(2) > 0) && (uv3(1) <= size(mov1,2)) && (uv3(2) <= dim(1))
                        if mov1(uv1(2),uv1(1)) && mov2(uv2(2),uv2(1)) && mov3(uv3(2),uv3(1))
                            %pixel_fly(k,:) = [uv1 uv2 uv3];
                            coordinate_fly(k,:) = [x y z];
                            %voxel_fly(k,:) = ([x y z]-range(1))/delta_range+1; % this seems to work but not sure about indexing
                            k = k + 1;
%                             subplot(3,1,1)
%                             hold on
%                             imshow(mov1.cdata)
%                             plot(uv1(1),uv1(2),'*')
%                             subplot(3,1,2)
%                             hold on
%                             imshow(mov2.cdata)
%                             plot(uv2(1),uv2(2),'*')
%                             subplot(3,1,3)
%                             hold on
%                             imshow(mov3.cdata)
%                             plot(uv3(1),uv3(2),'*')
%                             clf
                        end
                    end
                end
            end
        end
    end
    
end

try 
    scatter3(coordinate_fly(:,1), coordinate_fly(:,2), coordinate_fly(:,3),1)
    xlim([-6 6]*10^-3)
    ylim([-6 6]*10^-3)
    zlim([-6 6]*10^-3)
catch
    error('Couldnt plot the 3D scatter')
    dbstop if error
end
voxel_fly=uint8(voxel_fly);

% test points taken from original DLT
test=0;
if test ==1
    xyz=[2 -.75 -2]*10^-3;
    [uv_1] = dlt_inverse(DLt_Coef(:,1),xyz);
    [uv_2] = dlt_inverse(DLt_Coef(:,2),xyz);
    [uv_3] = dlt_inverse(DLt_Coef(:,3),xyz);
    uv_1=round(uv_1);
    uv_2=round(uv_2);
    uv_3=round(uv_3);
    figure
    subplot(3,1,1)
    hold on
    imshow(mov1.cdata)
    plot(uv_1(1),uv_1(2),'*')
    subplot(3,1,2)
    hold on
    imshow(mov2.cdata)
    plot(uv_2(1),uv_2(2),'*')
    subplot(3,1,3)
    hold on
    imshow(mov3.cdata)
    plot(uv_3(1),uv_3(2),'*')
end

end