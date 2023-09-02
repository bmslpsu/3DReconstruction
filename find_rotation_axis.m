function z_axis = find_rotation_axis(init, sI, eI, showplot)
%% find_rotation_axis:
% Updated to take modified voxel space for faster reconstruction
% Finding the z-axis is now done using PCA
% 
% Finds the axis of rotation and the z-axis of the fly from the
% magnetic tether.
%
% INPUTS:
%   Mask1,2,3: The mask structures than contains the tether and body mask
%   image_resolution: the size of the video images
%   DLt_Coef: the DLT coefficients arranged as [Cam1 Cam2 Cam3]
%   start_index & end_index: range of videos to be analyzed
%   Rec_Cord: The volume in which the fly and tether are present. This is
%   determined from a previous function
%
% OUTPUT:
%   z_axis: The normalized vector of the z_axis or yaw axis of the fly
%

maskdir = init.folders.mask;
crop_reigon = init.data.crop; 
img_sz = init.data.image_size;
DLt_Coef = init.data.DLT;
load(init.paths.volume, 'Rec_Cord')

mask_path = maskReader_image(maskdir);

% Start the 3D reconstruction
disp('Generating z-axis coordinates')
j = 1;
frame_incriment = (eI - sI) / 10;
frame_incriment = round(frame_incriment);
MOV = cell(3,1);
frames = sI:frame_incriment:eI;
n_frame = length(frames);
Tether_vec = zeros(n_frame,3);
for n = frames
    % Load the tether masks for each view
    for f = 1:3
        Icrop = imread([mask_path(f).tether num2str(n) '.png']); % read mask image
        MOV{f} = pad_image(Icrop, crop_reigon{f}, img_sz{f}); % pad so mask matches full image size
    end

    % Recontruct the tether from all 3 views
    recon = reconstruct_from_images(MOV, DLt_Coef, Rec_Cord, false);

%     % Fits using linear regressions
%     [x_reg, y_reg, z_reg, p0, d] = LinearRegression_3D(recon); % used to fit  the tether to a st line
%     TetherData(j).x_reg = x_reg;
%     TetherData(j).y_reg = y_reg;
%     TetherData(j).z_reg = z_reg; % coordinates of the line for a certain range of points
%     TetherData(j).d = d;
%     TetherData(j).p0 = p0;`

    % Fits 3D line using PCA
    PCA_Tether = pca(recon); % used to find the fly heading
    Tether_vec(j,:) = PCA_Tether(:,1); % take the strongest component
    if j > 1 % makes sure all the vectors for the z axis are along the same direction
        Tether_vec(j,:) = sign(dot(Tether_vec(j,:), Tether_vec(j-1,:)))*Tether_vec(j,:);
    end
    j = j + 1;
    
    if showplot
        scatter3(recon(:,1), recon(:,2), recon(:,3), 1)
        hold on
        xlim(2*[-2 2]*10^-3)
        ylim(2*[-2 2]*10^-3)
        zlim(2*[-2 2]*10^-3)
        drawnow
    end
end

% Estimate the yaw axis for multiple frames and normalize the axis to a unit vector
z_axis = sum(Tether_vec, 1);
z_axis = z_axis / norm(z_axis);

%% correct the axis orientation
MOV_body = cell(3,1);
for f = 1:3
    Icrop_body = imread([mask_path(f).body num2str(sI) '.png']); % read mask image
    MOV_body{f} = pad_image(Icrop_body, crop_reigon{f}, img_sz{f}); % pad so mask matches full image size
end

% Recontruct the tether from all 3 views
recon = reconstruct_from_images(MOV_body, DLt_Coef, Rec_Cord, false);
fly_center=mean(recon);
scatter3(recon(:,1),recon(:,2),recon(:,3))
hold on
plot3([fly_center(1), fly_center(1)+z_axis(1)], [fly_center(2), fly_center(2)+z_axis(2)], [fly_center(3), fly_center(3)+z_axis(3)], 'k', 'LineWidth', 2)

flip_axis=input('Flip z-axis y/n? ','s');
switch flip_axis
    case 'y'
        z_axis=-z_axis;
    case 'n'
        %Do nothing
end
%% plots and save

mkdir(init.folders.z_axis)
save(init.paths.z_axis, 'z_axis')
disp('Done')

end