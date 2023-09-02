function [H] = imshow3Dview(MOV, H)
%% imshow3D: dispaly 3 images on x-y-z planes
MOV = cellfun(@(x) uint8(flipud(x)), MOV, 'UniformOutput', false);
scl = 0.005;
if isempty(H)
    H(1) = surface(scl*[-1 1; -1 1], scl*[1 1; 1 1],   scl*[-1 -1; 1 1], ...
        'FaceColor', 'texturemap', 'CData', MOV{1});
    H(2) = surface(scl*[1 1; 1 1],   scl*[-1 1; -1 1], scl*[-1 -1; 1 1], ...
        'FaceColor', 'texturemap', 'CData', MOV{2});
    H(3) = surface(scl*[-1 1; -1 1], scl*[-1 -1; 1 1], scl*[-1 -1; -1 -1], ...
        'FaceColor', 'texturemap', 'CData', MOV{3});

    view(3);
    axis image
    axis off
    grid on
    colormap(gray)
else
    set(H(1), 'CData', MOV{1})
    set(H(2), 'CData', MOV{2})
    set(H(3), 'CData', MOV{3})
end

drawnow

end