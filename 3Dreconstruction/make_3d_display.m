function [] = make_3d_display(vidpths, background, recon, frame_range, frame_step, vidpath, slow_scale)
%% make_3d_display: 

Vread = vidReader(vidpths);
% n_frame = Vread(1).body.NumFrames;
frames = frame_range(1):frame_step:frame_range(2);
n_frame = length(frames);

if ~isempty(vidpath)
    export = true;
    fps = 8000 / frame_step;
    if ~isempty(slow_scale)
        vidfps = fps / slow_scale;
    end
    Vwrite = VideoWriter(vidpath,'MPEG-4');
    Vwrite.FrameRate = vidfps;
    Vwrite.Quality = 100;
    open(Vwrite)
else
    export = false;
end

fig = figure (1) ; clf
set(fig, 'Color', 'w')
fig_sz = fig.Position;
MOV = cell(3,1);
for n = 1:n_frame
    disp(frames(n))
    for f = 1:3
        I = read(Vread{f}, frames(n)) + background{f};
        MOV{f} = flipud(I);
    end
    
    if n == 1
        H(1) = surface([-1 1; -1 1], [1 1; 1 1], [-1 -1; 1 1], ...
            'FaceColor', 'texturemap', 'CData', MOV{1} );
        H(2) = surface([1 1; 1 1], [-1 1; -1 1], [-1 -1; 1 1], ...
            'FaceColor', 'texturemap', 'CData', MOV{2} );
        H(3) = surface([-1 1; -1 1], [-1 -1; 1 1], [-1 -1; -1 -1], ...
            'FaceColor', 'texturemap', 'CData', MOV{3} );

        view(3);
        axis image
        axis off
        colormap(gray)
    else
        set(H(1), 'CData', MOV{1})
        set(H(2), 'CData', MOV{2})
        set(H(3), 'CData', MOV{3})
    end
    
    drawnow
    pause(0.001)
    
    
    if export
        fig_frame = getframe(fig);
        fig_frame.cdata = fig_frame.cdata(...
            round(0.1*fig_sz(2)):end-round(0.1*fig_sz(2)), ...
            round(0.12*fig_sz(1)):end-round(0.07*fig_sz(1)), :);
        writeVideo(Vwrite, fig_frame);
    end
    
end
if export
    close(Vwrite)
end

end