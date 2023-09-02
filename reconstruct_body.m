function [recon, angle] = reconstruct_body(init, frames,plot_flag)
%% reconstruct_body: reconstructs the fly body in the 3D plan and uses it to determine the heading of the fly
% Faster reconstruction using a predefined smaller volume
% Uses PCA for finding the fly heading
% OUPUT:
% ReconstructionData_space: The reconstruction data in metric units
% ReconstructionData_Voxel: The reconstruction in voxel coordinates
% load and analyze video
% determines changes in heading and the orthogonal system for each
% reconstruction
if plot_flag
    f1_body_plot=figure;
end
maskdir = init.folders.mask;
crop_reigon = init.data.crop; 
img_sz = init.data.image_size;
DLt_Coef = init.data.DLT;
load(init.paths.volume, 'Rec_Cord')
load(init.paths.z_axis, 'z_axis')

recondir = fullfile(init.folders.root, 'reconstruction');
mkdir(recondir)

mask_path = maskReader_image(maskdir);
image_type = '.png';
MOV = cell(3,1);

tic
sI = frames(1);
n_frame = length(frames);
body_angle_all=zeros(1,length(frames));
for n = 1:n_frame
    if ~mod(frames(n),10) || (frames(n) == frames(1)) || (frames(n) == frames(end))
        disp(frames(n))
    end
    recon.frame = frames(n);
    
    % Load the tether masks for each view
    fI = num2str(frames(n));
    for f = 1:3
        Icrop = imread([mask_path(f).body fI image_type]); % read mask image
        MOV{f} = pad_image(Icrop, crop_reigon{f}, img_sz{f}); % pad so mask matches full image size
    end

    % Recontruct the body from all 3 views
    coordinate_fly = reconstruct_from_images(MOV, DLt_Coef, Rec_Cord, false);
    recon.body_xyz = coordinate_fly;
     
    % fit 3D line to fly
    body_vector = pca(coordinate_fly); % sed to find the fly heading
    body_vector = body_vector(:,1); % take the strongest component
    body_vector = body_vector/norm(body_vector);
    recon.body_vector = body_vector; % checks if the body vector has flipped from one frame to the next
    if (frames(n) > sI) && ( dot(recon.body_vector, recon_prev.body_vector) < 0 )
        recon.body_vector = -recon.body_vector;
    end
    if frames(n)==sI %determines if the first body vector is in the right direction (from abdomen to head)
        recon.body_vector=recon.body_vector*sign(dot(z_axis,body_vector)); %makes sure the body vector is in the +ve z-axis
        f1_body=figure;
        scatter3(coordinate_fly(:,1),coordinate_fly(:,2),coordinate_fly(:,3))
        hold on
        body_center=mean(coordinate_fly);
        body_vector_red=recon.body_vector/200;
        plot3([body_center(1) body_center(1)+body_vector_red(1)],...
            [body_center(2) body_center(2)+body_vector_red(2)]...
            ,[body_center(3) body_center(3)+body_vector_red(3)],'k')
        flip_axis=input('Flip body-axis y/n? ','s');
        switch flip_axis
            case 'y'
                recon.body_vector=-recon.body_vector;
            case 'n'
                %Do nothing
        end
        close(f1_body)
        clear body_center body_vector_red
    end
    % y-axis and x-axis calculation + save axes into a structure
    x_axis = cross(recon.body_vector, z_axis);
    x_axis = x_axis / norm(x_axis);
    y_axis = cross(z_axis,x_axis);
    recon.x_axis = x_axis;
    recon.y_axis = y_axis;
    recon.z_axis = z_axis;
    
    % Calculate angle of the fly
    if frames(n) > sI
       angle = atan2(norm(cross(recon_init.x_axis, recon.x_axis)),...
           dot(recon_init.x_axis, recon.x_axis));
    else %initialize zero position of the fly as the zero yaw angle
        angle = 0;
        recon_init = recon;
    end
    recon.body_angle = angle;
    body_angle_all(n)=angle;
    recon_prev = recon;
    %% Plot
    if plot_flag
        
        scatter3(coordinate_fly(:,1),coordinate_fly(:,2),coordinate_fly(:,3))
        hold on
        body_center=mean(coordinate_fly);
        body_vector_red=recon.body_vector/200;
        x_axis_red=x_axis/200;
        y_axis_red=y_axis/200;
        z_axis_red=z_axis/200;
        plot3([body_center(1) body_center(1)+body_vector_red(1)],...
            [body_center(2) body_center(2)+body_vector_red(2)]...
            ,[body_center(3) body_center(3)+body_vector_red(3)],'k')
        plot3([body_center(1) body_center(1)+x_axis_red(1)],...
            [body_center(2) body_center(2)+x_axis_red(2)]...
            ,[body_center(3) body_center(3)+x_axis_red(3)],'r')
        plot3([body_center(1) body_center(1)+y_axis_red(1)],...
            [body_center(2) body_center(2)+y_axis_red(2)]...
            ,[body_center(3) body_center(3)+y_axis_red(3)],'b')
        plot3([body_center(1) body_center(1)+z_axis_red(1)],...
            [body_center(2) body_center(2)+z_axis_red(2)]...
            ,[body_center(3) body_center(3)+z_axis_red(3)],'g')  
    end
 	%% Save
    savepath = fullfile(recondir, ['frame_' fI '.mat']);
    save(savepath, '-struct', 'recon')
end
if plot_flag
    close(f1_body_plot)
end
toc
angle_dir=fullfile(init.folders.root, 'body_angle_all');
mkdir(angle_dir)
angle_savepath=fullfile(angle_dir, 'BodyAngle.mat');
save(angle_savepath, '-struct', 'body_angle_all')
end