function [img_body_dilated, img_tether_dilated] = seperate_body_tether(Image_bw, pts_teth, bottomFlag)
%% Seperates body & tether masks
%
% INPUT:
%   Image_bw: binarized image or body and tether mask
%   bottomFlag: is this the bottom view
%   pts_teth: location of tether
%
% OUTPUT
%   img_body: body image or body mask
%   img_tether: tether image or mask
%

% Extract body and tether mask
y = round(pts_teth(:,2)); % y and x location of tether points
if ~bottomFlag % side view
    img_body = Image_bw;
    img_tether = Image_bw;
    img_tether(y(1):end,:) = 0;
    img_body(1:y(1),:) = 0;
    
    % remove the noise thats below the fly
    CC = bwconncomp(img_body);
    lengthP = zeros(CC.NumObjects,1);
    for i = 1:CC.NumObjects
       lengthP(i) = length(CC.PixelIdxList{i});
    end
    [~,index_max_pixelSize] = max(lengthP);
    img_body = zeros(size(img_body));
    img_body(CC.PixelIdxList{index_max_pixelSize}) = 1;

elseif bottomFlag % bottom view which needs special methods
    Image_bw(y(2):end,:) = 0; % removes the brass rod(not particularly needed now)
    Image_bw = bwareafilt(Image_bw,3); % removes the noise and small objects
    Image_bw_noBody = Image_bw;
    Image_bw_noBody(1:y(1),:) = 0;
%     Vector_sum = sum(Image_bw_noBody,2);
%     index_sepPoint = find(Vector_sum > 10, 1, 'first'); % finds the first point at which the body enters the tether mask
%     if ~isempty(index_sepPoint)
%         Image_bw_noBody(1:index_sepPoint,:) = 0; % removes the body
%     end
    img_tether = Image_bw_noBody;
    img_body = Image_bw-img_tether;
end
se = strel('square',3);
img_tether_dilated = imdilate(img_tether,se);
img_body_dilated = imdilate(img_body,se);

end