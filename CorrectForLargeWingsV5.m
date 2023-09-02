function wing_reconstruction=CorrectForLargeWingsV5(init,frames,flag_stroke_reversal)

%% Initialize the directories of the folders
DLt_Coef = init.data.DLT;
%% Video readers
Vreader1 = VideoReader(init.paths.vid{1});
Vreader2 = VideoReader(init.paths.vid{2});
Vreader3 = VideoReader(init.paths.vid{3});

f1=figure;
%% size correction
for ii=frames
    disp(['Showing frame: ' num2str(ii)])
    %% read the frames
    image1 = read(Vreader1, ii);
    image2 = read(Vreader2, ii);
    image3 = read(Vreader3, ii);
    %% read the wing 3D coordinates
    reconpath = fullfile(init.folders.reconstruction, ['frame_' num2str(ii) '.mat']);
    load(reconpath, 'wingLeft');
    load(reconpath, 'wingRight');
    load(reconpath, 'body_xyz');
    %% correct wing 1 and 2
    wing_coordinates1=CorrectByUserV2(wingRight,DLt_Coef,image1,image2,image3,flag_stroke_reversal(ii));
    wing_coordinates2=CorrectByUserV2(wingLeft,DLt_Coef,image1,image2,image3,flag_stroke_reversal(ii));
    %% plot the data
    subplot(2,1,1)
    scatter3(wingRight(:,1),wingRight(:,2),wingRight(:,3))
    hold on
    scatter3(wingLeft(:,1),wingLeft(:,2),wingLeft(:,3))
    scatter3(body_xyz(:,1),body_xyz(:,2),body_xyz(:,3),'g')
    subplot(2,1,2)
    if ~isempty(wing_coordinates1)
        scatter3(wing_coordinates1(:,1),wing_coordinates1(:,2),wing_coordinates1(:,3))
    else
        scatter3(wingRight(:,1),wingRight(:,2),wingRight(:,3),'r')
    end
    hold on
    if ~isempty(wing_coordinates2)
        scatter3(wing_coordinates2(:,1),wing_coordinates2(:,2),wing_coordinates2(:,3))
    else
        scatter3(wingLeft(:,1),wingLeft(:,2),wingLeft(:,3),'r')
    end
    scatter3(body_xyz(:,1),body_xyz(:,2),body_xyz(:,3),'g')
    %% correct the structure data
    if ~isempty(wing_coordinates1)
        wing1_corrected=wing_coordinates1;
    else
        wing1_corrected=wingRight;
    end
    if ~isempty(wing_coordinates2)
        wing2_corrected=wing_coordinates2;
    else
        wing2_corrected=wingLeft;
    end
    wing_reconstruction(ii).WingsCoordinates=[wing_coordinates1;wing_coordinates2];
    %% save and misc
    save(reconpath,'wing1_corrected','-append')
    save(reconpath,'wing2_corrected','-append')
    hold off
    clear wing_coordinates2 wing_coordinates1
end
close(f1)
end

function wing_coordinates=CorrectByUserV2(wing_coordinates,DLt_Coef,image1,image2,image3,FlagStrokeReversal)
%% a simple user interface that allows the user to refine the reconstruction of
%the wing based on various metrics that indicate a faulty wing
%reconstruction
%% code to determine possible issues in wing reconstruction and fix it
if FlagStrokeReversal && length(wing_coordinates)>3500
    
    disp('select an roi that has the wing leading edge and tip')
    [~,~,~,wing_coordinatesLeadingEdge1]=UserInputReconstructionV3(wing_coordinates,DLt_Coef,image1,image2,image3);
    disp('select an roi that has the wing trailing edge')
    [~,~,~,wing_coordinatesTrailingEdge1]=UserInputReconstructionV3(wing_coordinates,DLt_Coef,image1,image2,image3);
    PointsToFot=[wing_coordinatesLeadingEdge1; wing_coordinatesTrailingEdge1];
    %fit a plane to the selected points
    fitobject = fit([PointsToFot(:,1),...
        PointsToFot(:,2)],PointsToFot(:,3),...
        'poly11');
    coeff=coeffvalues(fitobject);
    
    planeC0eff=[coeff(1) coeff(2) coeff(3) -1]; %[D A B C] Ax+By+Cz+D
    for i=1:length(wing_coordinates)
        DistancesToPlane(i)=abs(planeC0eff(2)*wing_coordinates(i,1)+planeC0eff(3)*wing_coordinates(i,2)...
            +planeC0eff(4)*wing_coordinates(i,3)+planeC0eff(1))/...
            sqrt(planeC0eff(2)^2+planeC0eff(3)^2+planeC0eff(4)^2);
    end
    wing_coordinates(DistancesToPlane>mean(DistancesToPlane)/2,:)=[];
else
end
end

