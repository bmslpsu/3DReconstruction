function [FilteredAngleR, FilteredAngleL]=FilterEulerAngles(EulerAnglesR, EulerAnglesL,cutFreq,SamplingFreq,init)
%This function filters and corrects the Euler angles of both wings

%% start by correcting the stroke angle (might only need to correct this one)
%% Right wing
for ii=2:length(EulerAnglesR)
    if abs(EulerAnglesR(ii,1)-EulerAnglesR(ii-1,1))>2.5
        EulerAnglesR(ii,1)=EulerAnglesR(ii,1)-2*pi*sign(EulerAnglesR(ii,1)-EulerAnglesR(ii-1,1));
        
    end
    
%     if abs(EulerAnglesR(ii,2)-EulerAnglesR(ii-1,2))>pi/2
%         EulerAnglesR(ii,2)=EulerAnglesR(ii,2)*sign(EulerAnglesR(ii-1,2));
%         
%     end
end

%% Left wing
for ii=2:length(EulerAnglesL)
    if abs(EulerAnglesL(ii,1)-EulerAnglesL(ii-1,1))>2.5
        EulerAnglesL(ii,1)=EulerAnglesL(ii,1)-2*pi*sign(EulerAnglesL(ii,1)-EulerAnglesL(ii-1,1));
        
    end
    
%     if abs(EulerAnglesL(ii,2)-EulerAnglesL(ii-1,2))>pi/2
%         EulerAnglesL(ii,2)=EulerAnglesL(ii,2)-2*pi*sign(EulerAnglesL(ii,2)-EulerAnglesL(ii-1,2));
%         
%     end
end

%% angle filtering
FilteredAngleR=FilterAngles(EulerAnglesR,cutFreq,SamplingFreq);
FilteredAngleL=FilterAngles(EulerAnglesL,cutFreq,SamplingFreq);
FilteredAngleL(:,1)=FilteredAngleL(:,1);
%% change signs and offset

FilteredAngleR(:,1)=FilteredAngleR(:,1)+pi;
FilteredAngleL(:,1)=-FilteredAngleL(:,1)+0;
FilteredAngleR(:,2)=-FilteredAngleR(:,2);
FilteredAngleL(:,2)=-FilteredAngleL(:,2);
FilteredAngleR=FilteredAngleR*180/pi;
FilteredAngleL=FilteredAngleL*180/pi;
%% save data
dir_angles=fullfile(init.folders.root,'angles');
mkdir(dir_angles)

save(fullfile(dir_angles,'wing_data.mat'),'FilteredAngleR','FilteredAngleL') %wing 1 is the right wing

%% plot the angles
figure
subplot(3,1,1)
plot(FilteredAngleR(:,1),'b')
hold on
plot(FilteredAngleL(:,1),'r')
title('Stroke angle')
legend('Right wing','Left wing')
subplot(3,1,2)
hold on
plot(FilteredAngleR(:,2),'b')
% plot(-EulerAnglesR(:,2)*180/pi,'b--')
plot(FilteredAngleL(:,2),'r')
% plot(-EulerAnglesL(:,2)*180/pi,'r--')
title('Rotation angle')
subplot(3,1,3)
plot(FilteredAngleR(:,3),'b')
hold on
plot(FilteredAngleL(:,3),'r')
title('Deviation angle')
xlabel('Frame')
end
%% FUNCTIONS
function FilteredAngle=FilterAngles(EulerAngle,cutFreq,SamplingFreq)
%% start with a hampel filter (using def settings)
FilteredAngle(:,1)=hampel(EulerAngle(:,1));
FilteredAngle(:,2)=hampel(EulerAngle(:,2));
FilteredAngle(:,3)=hampel(EulerAngle(:,3));

%% low-pass filtering (variable cut-off freq)
[b_stroke, a_stroke]=butter(4,cutFreq/(SamplingFreq/2));
[b_rotation, a_rotation]=butter(4,1.8*cutFreq/(SamplingFreq/2)); %rotation angle is faster hence a larger cut-off freq
[b_dev, a_dev]=butter(4,1*cutFreq/(SamplingFreq/2)); %smaller cut-off to reduce noise
%% filter data 
FilteredAngle(:,1)=filtfilt(b_stroke,a_stroke,FilteredAngle(:,1));
FilteredAngle(:,2)=filtfilt(b_rotation,a_rotation,FilteredAngle(:,2));
FilteredAngle(:,3)=filtfilt(b_dev,a_dev,FilteredAngle(:,3));
%% spline interpolation?
end