function [Iout] = pad_image(Icrop, crop_reigon, I)
% pad_image: takes a cropped image, the crop reactangualr reigon, & 
% matrix of the orginal image size (I) and pads the cropped image with 
% 0's to match the size of the orginal image.
%
if all(size(I) > 1)
    dim = size(I); % matrix with orginal image size
else
    dim = I; % orginal image size directly
end
topR = crop_reigon(2) - 1; % top padding
leftR = crop_reigon(1) - 1; % left padding
bottomR = dim(1) - crop_reigon(4) - crop_reigon(2); % bottom padding
rightR = dim(2) - crop_reigon(3) - crop_reigon(1); % right padding

Iout = padarray(Icrop, [topR leftR], 0, 'pre'); % add padding to top & left
Iout = padarray(Iout, [bottomR rightR], 0, 'post'); % add padding to bottom & right
end