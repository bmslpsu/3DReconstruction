function [OcclusionFlag] = correct_tether_occlusion_imageV2(init, showplot,label,sI,eI)
%% Detect the location and type of wing occlusions & remake masks
se_correction = strel('disk',150);
se = strel('square',2);


fig_correction = figure;

maskdir = init.folders.mask;
label = char(string(label));
maskdir_wing=fullfile(maskdir, ['wing_' label]);
maskdir_tether=fullfile(maskdir, ['tether_' label]);
n_frame=eI;
OcclusionFlag = uint8(zeros(n_frame,1));
for n = sI:n_frame
    if ~mod(n,10) || (n == sI) || (n == n_frame)
        disp(n)
    end
    
    %% Wing and tether mask directory for each frame
    maskdir_wing_frame=fullfile(maskdir_wing,['frame_',num2str(n) '.png']);
    maskdir_tether_frame=fullfile(maskdir_tether,['frame_',num2str(n) '.png']);
    %% Read tether & wing masks
    mask_wing = imread(maskdir_wing_frame);
    mask_tether = imread(maskdir_tether_frame);
    mask_tether_morp = bwmorph(mask_tether, 'thicken');
    
	% Determine occlusion range for automatic corrections
    mask_WingTetherInter = mask_tether_morp & mask_wing; % remove the wing and tether but keeping the interesection range
    mask_WingTetherInter = imdilate(mask_WingTetherInter,se);
    CC_WTmask = bwconncomp(mask_WingTetherInter);
    
    if CC_WTmask.NumObjects >= 2
        OcclusionFlag(n) = 2; % automated fix
    elseif CC_WTmask.NumObjects == 1
        OcclusionFlag(n) = 1; % wing is behind a portion of the tether (hard to fix)
    else
        OcclusionFlag(n) = 0; % no need for fix
    end
    
    % Correction statements based on occlusion flag
    if OcclusionFlag(n) == 2
        % performs a closing to reduce tether occlusion
        bw = imclose(mask_WingTetherInter, se_correction);        
        mask_wing_correct = mask_wing | bw;
        imwrite(mask_wing_correct,maskdir_wing_frame);
    else
        mask_wing_correct = mask_wing; % does nothing
    end
    
    if showplot
        if n ~= sI
            subplot(1,2,1) ; set(H.raw, 'CData', mask_wing);
            subplot(1,2,2) ; set(H.corrct, 'CData', mask_wing_correct);
        else
            subplot(1,2,1) ; H.raw = imshow(mask_wing);
            subplot(1,2,2) ; H.corrct = imshow(mask_wing_correct);
        end
        drawnow
    end    
end
try
    close(fig_correction)
catch
end
end
