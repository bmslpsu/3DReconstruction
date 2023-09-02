function Rec_Cord = find_voxel_range(init, voxel_size)
%% find_voxel_range: user interface function to help reduce the voxel space and reduce the time needed to
% one of the main problems when doing 3D work is that the 3D voxels require
% a long time for full reconstruction. This function reduces that volume
% using the user data

Vread = vidReader(init.paths.vid);
DLt_Coef = init.data.DLT;

voxel_size_reduced = voxel_size*10;
space_range_initial = [-0.01 0.01]; % the initial range
volume_cord = CreateVoxelSpave(space_range_initial, voxel_size);
volume_cord_red = CreateVoxelSpave(space_range_initial, voxel_size_reduced);

uv = cell(3,1);
Mask_rect = cell(3,1);
index_remove = false;
for n = 1:3
    uv{n} = round(dlt_inverse(DLt_Coef(:,n),volume_cord)); % inverse DLT normal voxel resolution
    uv_red = round(dlt_inverse(DLt_Coef(:,n), volume_cord_red)); % inverse DLT reduced
    image_Original = read(Vread{n}, 1); % load first image
    Mask_rect{n} = FindMask(image_Original, uv_red); % find the ROI of the fly in each image
    image_res = size(Mask_rect{n}); % ROI size
    index_remove = index_remove | ...
        ( (uv{n}(:,1) < 1) | (uv{n}(:,1) > image_res(2)) | (uv{n}(:,2) < 1) | (uv{n}(:,2) > image_res(1)) );
end

volume_cord(index_remove,:) = [];
pixel_value = cell(3,1);
for n = 1:3
   uv{n}(index_remove,:) = [];
   ind = sub2ind(size(Mask_rect{n}), uv{n}(:,2), uv{n}(:,1));
   pixel_value{n} = Mask_rect{n}(ind);
end

% Find the pixels that aren't part of the fly and wings
SumPixels = sum(cat(2,pixel_value{:}),2);
Index_keep = SumPixels == 3; % everything less than 3 is removed since it shouldn't exist in at least one image
Rec_Cord = volume_cord(Index_keep,:); % remove the false voxels from coordinate frame

mkdir(init.folders.volume)
save(init.paths.volume, 'Rec_Cord')

end

%% Functions
function Mask_rect = FindMask(image_Original, uv)
f1 = figure;
imshow(image_Original)
hold on
scatter(uv(:,1), uv(:,2))
Mask_rect = roipoly;
close(f1)
end

function voxel_cor = CreateVoxelSpave(space_range_initial, voxel_size)
x = space_range_initial(1):voxel_size:space_range_initial(2);
y = space_range_initial(1):voxel_size:space_range_initial(2);
z = space_range_initial(1):voxel_size:space_range_initial(2);
[xx,yy,zz] = meshgrid(x,y,z);
voxel_cor = [xx(:),yy(:),zz(:)]; % generates a volume grid
end
