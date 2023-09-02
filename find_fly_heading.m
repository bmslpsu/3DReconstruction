function [ReconstructionData_space, angle] = find_fly_heading(maskdir, DLt_Coef, crop_reigon, img_sz, z_axis, range, volume_cord)
%% find_fly_heading: reconstructs the fly body in the 3D plan and uses it to determine the heading of the fly
% Uses PCA for finding the fly heading
%
% OUPUT:
%   ReconstructionData_space: The reconstruction data in metric units
%   ReconstructionData_Voxel: The reconstruction in voxel coordinates

% Load and analyze video
% determines changes in heading and the orthogonal system for each reconstruction
tic

frames = 1:4000;
% Generate the reconstruction
% ReconstructionData_space = ReconstructionFunctionV2(start_index,end_index,DLt_Coef,Mask1,Mask2,Mask3,...
%     range,delta_range,image_resolution,volume_cord);

ReconstructionData_space = reconstruction_body(maskdir, DLt_Coef, crop_reigon, img_sz, volume_cord, frames);

% Place reconstruction in structure and find heading of the fly
% Fig_body = figure;
n_frame = length(frames);
angle = zeros(n_frame,1);
sI = 1;
fig = figure;
ax = subplot(1,1,1) ; hold on
for n = frames
    % Place data in structure
    coordinate_fly = ReconstructionData_space(n).CoordinateReconstruction;
    
    % Fit 3D line to fly
    Body_vec = pca(coordinate_fly); % used to find the fly heading
    Body_vec = Body_vec(:,1); % take the strongest component
    Body_vec = -Body_vec/norm(Body_vec);
    ReconstructionData_space(n).Body_vec = Body_vec; % checks if the body vector has flipped from one frame to the next
    if (n > 1) && ( dot(ReconstructionData_space(n).Body_vec,ReconstructionData_space(n-1).Body_vec) < 0 )
        ReconstructionData_space(n).Body_vec = -ReconstructionData_space(n).Body_vec;
    end
    
    % y-axis and x-axis calculation & save axes into a structure
    x_axis = cross(Body_vec,z_axis);
    x_axis = x_axis/norm(x_axis);
    y_axis = cross(z_axis,x_axis);
    ReconstructionData_space(n).xAxis = x_axis;
    ReconstructionData_space(n).yAxis = y_axis;
    ReconstructionData_space(n).zAxis = z_axis;
    
    % Plot axes
    hold on 
    scatter3(coordinate_fly(:,1), coordinate_fly(:,2), coordinate_fly(:,3))

    xlim(range)
    ylim(range)
    zlim(range)
    mean_fly=mean(coordinate_fly);
    plot3([mean_fly(1) mean_fly(1)+x_axis(1)],[mean_fly(2) mean_fly(2)+x_axis(2)],[mean_fly(3) mean_fly(3)+x_axis(3)],'r')
    plot3([mean_fly(1) mean_fly(1)+y_axis(1)],[mean_fly(2) mean_fly(2)+y_axis(2)],[mean_fly(3) mean_fly(3)+y_axis(3)],'b')
    plot3([mean_fly(1) mean_fly(1)+z_axis(1)],[mean_fly(2) mean_fly(2)+z_axis(2)],[mean_fly(3) mean_fly(3)+z_axis(3)],'k')
    hold off
    
    % Calculate angle of the fly, this starts with the first frame and moves forward
    %angle(start_index)=0;
    if n > 0
       angle(n)= atan2(norm(cross(ReconstructionData_space(sI).xAxis,ReconstructionData_space(n).xAxis)),...
           dot(ReconstructionData_space(sI).xAxis, ReconstructionData_space(n).xAxis));
    end
end
toc 