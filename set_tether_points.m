function [pt_sep] = set_tether_points(root_vid)
Vreader = VideoReader(root_vid);
I = read(Vreader, 1);
fig = figure;
imshow(I)
title('Select the start and end point of the tether: Fly to top. Click space when done.')
roi = drawline(gca, 'Color', 'r', 'LineWidth', 0.5);
pause
pt_sep = roi.Position;
close(fig)
end