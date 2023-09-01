function [phi1, theta1, eta1, phi2, theta2, eta2, eta11]=EstimateEulerAngles(init,sI,eI)

hull_analysis_dir=fullfile(init.folders.root,"hull_analysis");

for ii=sI:eI
    hull_analysis_path = fullfile(hull_analysis_dir, ['frame_' num2str(ii) '.mat']);
    load(hull_analysis_path, 'chord1AltHat');
    load(hull_analysis_path, 'chord1Hat');
    load(hull_analysis_path, 'chord2AltHat');
    load(hull_analysis_path, 'chord2Hat');
    load(hull_analysis_path, 'span1Hat');
    load(hull_analysis_path, 'span2Hat');
    
    %% wing 1
    rotAngle=45;
    RotationMatrixX=rotx(rotAngle);
    spanHatRightRotated=RotationMatrixX*span1Hat;
    phi1(ii)   = atan2(spanHatRightRotated(2),spanHatRightRotated(1)); %stroke
    theta1(ii) = asin(spanHatRightRotated(3)); %deviation
    
    
    
    r_rot = vrrotvec([1,0,0],span1Hat);
    m_rot = vrrotvec2mat(r_rot);
    y_wing=[0;1;0]*sign(chord1Hat(2));
    z_wing1=cross(span1Hat,y_wing);
    z_wing1=z_wing1/norm(z_wing1);
    wingNormVecWael = cross(span1Hat,chord1Hat);
    %project z axis on the span-normal plane
    
    z_wing_proj=z_wing1-dot(z_wing1,span1Hat)*span1Hat;
    z_wing_proj=z_wing_proj/norm(z_wing_proj);
    %project onto xy plant
    %also finds the angle
    eta1(ii) = sign(wingNormVecWael(2))*atan2(norm(cross(z_wing_proj,wingNormVecWael)),dot(z_wing_proj,wingNormVecWael));
 
    %% estimates of eta1 sam's code
    phi1Hat = [-sin(phi1(ii)),cos(phi1(ii)),0]; %a vector thats perpendicular to span vector
    phi1Proj = [span1Hat(1),span1Hat(2),0]; %projection of wing span on xy plane
    phi1Proj = phi1Proj / norm(phi1Proj); %normalize it
    theta1Hat = -phi1Proj*sin(theta1(ii)) + [0,0,1]*cos(theta1(ii)); %vector of the deviation angle motion
    eta11(ii) = atan2( abs(dot(chord1Hat,theta1Hat)), abs(dot(chord1Hat,phi1Hat)));
 %% wing 2
    RotationMatrixX=rotx(rotAngle);
    spanHatLeftRotated=RotationMatrixX*span2Hat;
    phi2(ii)   = atan2(spanHatLeftRotated(2),spanHatLeftRotated(1)); %stroke
    theta2(ii) = asin(spanHatLeftRotated(3)); %deviation
    
    
    
    r_rot = vrrotvec([1,0,0],span2Hat);
    m_rot = vrrotvec2mat(r_rot);
    y_wing2=[0;1;0]*sign(chord2Hat(2));
    z_wing2=cross(span1Hat,y_wing2);
    z_wing2=z_wing2/norm(z_wing2);
    wingNormVecWael2 = cross(span2Hat,chord2Hat);
    %project z axis on the span-normal plane
    
    z_wing_proj2=z_wing2-dot(z_wing2,span2Hat)*span2Hat;
    z_wing_proj2=z_wing_proj2/norm(z_wing_proj2);
    %project onto xy plant
    %also finds the angle
    eta2(ii) = sign(wingNormVecWael2(2))*atan2(norm(cross(z_wing_proj2,wingNormVecWael2)),dot(z_wing_proj2,wingNormVecWael2));

end
%% extra code
phi1(1:sI-1)=[];
theta1(1:sI-1)=[];
eta1(1:sI-1)=[];
phi2(1:sI-1)=[];
theta2(1:sI-1)=[];
eta2(1:sI-1)=[];
eta11(1:sI-1)=[];
figure
plot(theta1*180/pi)
hold on
plot(theta2*180/pi)
figure
plot(phi1)
hold on
plot(phi2)
figure
plot(eta1)
hold on
plot(eta2)
figure
plot(eta11)
hold on

