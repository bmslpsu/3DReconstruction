function [Vreader] = vidReader(vidpaths)
% vidReader: sets up paths to read raw videos
%   

Vreader = cell(3,1);
for n = 1:3
    Vreader{n} = VideoReader(vidpaths{n});
end

end

