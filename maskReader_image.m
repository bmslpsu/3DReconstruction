function [paths] = maskReader_image(maskdir)
% maskReader_image: sets up mask paths for body, tether, & wings
%   

for n = 1:3
    label = num2str(n);
    paths(n).body   = fullfile(maskdir, ['body_' label], 'frame_');
    paths(n).tether	= fullfile(maskdir, ['tether_' label], 'frame_');
    paths(n).wing   = fullfile(maskdir, ['wing_' label], 'frame_');
end

end

