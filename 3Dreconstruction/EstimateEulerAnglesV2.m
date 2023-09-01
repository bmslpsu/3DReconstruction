function [phi1, theta1, eta1, phi2, theta2, eta2]=EstimateEulerAnglesV2(init,sI,eI)

hull_analysis_dir=fullfile(init.folders.root,"hull_analysis");

for ii=sI:eI
    hull_analysis_path = fullfile(hull_analysis_dir, ['frame_' num2str(ii) '.mat']);
    load(hull_analysis_path, 'chord1AltHat');
    load(hull_analysis_path, 'chord1Hat');
    load(hull_analysis_path, 'chord2AltHat');
    load(hull_analysis_path, 'chord2Hat');
    load(hull_analysis_path, 'span1Hat');
    load(hull_analysis_path, 'span2Hat');
    load(hull_analysis_path, 'angle_rotation1_final_unit');
    load(hull_analysis_path, 'angle_rotation2_final_unit');
    
    %% wing 1
    rotAngle=45;
    RotationMatrixX=rotx(rotAngle);
    spanHatRightRotated=RotationMatrixX*span1Hat;
    phi1(ii)   = atan2(spanHatRightRotated(2),spanHatRightRotated(1)); %stroke
    theta1(ii) = asin(spanHatRightRotated(3)); %deviation
    
    %also finds the angle
    eta1(ii) = angle_rotation1_final_unit*180/pi;
 
 %% wing 2
    RotationMatrixX=rotx(rotAngle);
    spanHatLeftRotated=RotationMatrixX*span2Hat;
    phi2(ii)   = atan2(spanHatLeftRotated(2),spanHatLeftRotated(1)); %stroke
    theta2(ii) = asin(spanHatLeftRotated(3)); %deviation
    
    eta2(ii) = angle_rotation2_final_unit*180/pi;

end
%% extra code
phi1(1:sI-1)=[];
theta1(1:sI-1)=[];
eta1(1:sI-1)=[];
phi2(1:sI-1)=[];
theta2(1:sI-1)=[];
eta2(1:sI-1)=[];
eta1=eta1*pi/180;
eta2=eta2*pi/180;
%% figures
figure
plot(theta1*180/pi)
hold on
plot(theta2*180/pi)
figure
plot(phi1)
hold on
plot(phi2)
figure
plot(eta1*180/pi)
hold on
plot(eta2*180/pi)