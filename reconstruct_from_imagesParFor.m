function [recon] = reconstruct_from_imagesParFor(MOV1,MOV2,MOV3, DLt_Coef, volume_cord, showplot)
%% reconstruct_from_images: get the 3D pixel cooridnates from masks in each camera view
% Updated to 
% 
% INPUTS:
%   MOV         : cell array with each cell containing a 2D mask frame
%   DLt_Coef    : DLT coefficents for each frame
%   volume_cord	: volume in which the fly and tether are present
%   showplot    : (boolean) show debug plot
%
% OUTPUT:
%   recon       : reconstruction volume coordinates (x,y,z)
%
MOV = cell(3,1);
MOV{1}=MOV1;
MOV{2}=MOV2;
MOV{3}=MOV3;
if nargin < 6
    showplot = true;
end

dim = cellfun(@(x) size(x), MOV, 'UniformOutput', false);
dim = cat(3, dim{:});

if ~all(dim)
    error('masks for each view must be the same size')
else
   img_sz = dim(:,:,1); 
end

% Create the xy coorindtaes from the given volume
n_view = length(MOV);
uv = cell(n_view,1);
index_remove = false;
for f = 1:3
    uv{f} = round(dlt_inverse(DLt_Coef(:,f), volume_cord)); % inverse DLT normal voxel resolution
    index_remove = index_remove | ...
        ( (uv{f}(:,1) < 1) | (uv{f}(:,1) > img_sz(2)) | (uv{f}(:,2) < 1) | (uv{f}(:,2) > img_sz(1)) );
end

% Get the pixel value for each valid xy coorindate
volume_cord(index_remove,:) = [];
pixel_value = cell(3,1);
for f = 1:3
   uv{f}(index_remove,:) = [];
   ind = sub2ind(img_sz, uv{f}(:,2), uv{f}(:,1));
   pixel_value{f} = MOV{f}(ind);
end

% Find the pixels that apear in every view
all_pixel = cat(2, pixel_value{:}); % pixel values from each view
all_on_pixel = all(all_pixel, 2); % pixel values that apear in every view
% all_on_pixel = sum(all_pixel,2) >= 2;
uv_on = cellfun(@(x) x(all_on_pixel,:), uv, 'UniformOutput', false); % xy pixel cooridinates in each view
recon = volume_cord(all_on_pixel,:); % volume xyz cooridinates that appear in every view

if showplot
   for f = 1: n_view
      subplot(1,3,f)
      imshow(MOV{f}) ; hold on ; axis on
      plot(uv{f}(:,1), uv{f}(:,2), '.b', 'MarkerSize', 1)
      plot(uv_on{f}(:,1), uv_on{f}(:,2), '.r', 'MarkerSize', 5);
   end
end

end