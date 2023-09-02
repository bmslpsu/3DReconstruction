function Flag_reversal = find_stroke_reversal(init,sI,eI, showplot)
%% Locate the approximate location where the wing reversal occurs
maskdir = init.folders.mask;


COM = cell(1,2);
COM{1} = zeros(eI,1);
COM{2} = zeros(eI,1);
Wing = cell(2,1);
for n = sI:eI
    if ~mod(n,10) || (n == 1) || (n == eI)
        disp(n)
    end
    
    % Extact the images
    for f = 1:2
        maskdir_wing1=fullfile(maskdir, ['wing_' num2str(f)]);
        maskdir_wing_frame1=fullfile(maskdir_wing1,['frame_',num2str(n) '.png']);
        wing_frame1 = imread(maskdir_wing_frame1);
        [Wing{f}, ~] = find(wing_frame1);
        
        if isempty(Wing{f})
            COM{f}(n) = nan;
        else
            COM{f}(n) = mean(Wing{f});
        end
    end
end

% Interpolate nan values (just in case), filter, & smooth the data
t_wing = 1:eI;
cutOff = 700; % cut-off frequency
[b,a] = butter(3, cutOff/(8000/2));
COM_int = COM;
COM_vel = COM;
for d = 1:2
    % Interpolate
    COM_int{d} = fillmissing(COM{d}, 'spline', 'SamplePoints', t_wing);
    
    % Filter
    COM_int{d} = filtfilt(b, a, COM_int{d});
    
    % Velocity
    %COM_vel{d} = abs(diff(COM_int{d}));
    COM_vel{d} = abs(central_diff(COM_int{d}));
end

% Locate the areas of the stroke reversal
error_eps = 1.2;
index_rev = (COM_vel{1} <= error_eps) | (COM_vel{2} <= error_eps);
%index_rev = [index_rev ; index_rev(end)];

% Create the flag for the reversal
index_offset = 1;
Flag_reversal = false(eI,1);
for n = 2:eI-index_offset
    if sum(index_rev(n:n+index_offset)) > 0
        Flag_reversal(n-1:n) = true;
    else
       %Flag_reversal(n) = 0;
   end
end

for n = 1:eI-4
    if index_rev(n)
        Flag_reversal(n:n+4) = true;
    end
end

% Debug plot
if showplot
    f1 = figure;
    subplot(2,2,1) ; hold on ; title('Position of wing mask COM in pixels')
       plot(COM_int{1})
       plot(COM_int{2})
   
    subplot(2,2,3) ; hold on ; title('Velocity of wing mask COM in pixels')
       plot(COM_vel{1})
       plot(COM_vel{2})

	subplot(2,2,[2 4]) ; hold on ; title('Stoke reversal boolean')
    imagesc(Flag_reversal)
    axis tight
end

end
