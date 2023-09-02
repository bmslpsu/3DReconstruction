function []=CorrectChordAutomatic(init,sI,eI)
%automatic corrector to keep things consistant

%% start analysis
hull_analysis_dir=fullfile(init.folders.root,"hull_analysis");
for ii=sI:eI
    hull_analysis_path = fullfile(hull_analysis_dir, ['frame_' num2str(ii) '.mat']);
    load(hull_analysis_path, 'chord1AltHat');
    load(hull_analysis_path, 'chord1Hat');
    load(hull_analysis_path, 'chord2AltHat');
    load(hull_analysis_path, 'chord2Hat');
    
    if ii==sI %initialize chord1 and chord2
        chord1=chord1Hat;
        chord2=chord2Hat;
        
    elseif ii>sI
        %% correct chord 1
        dot_c_1=dot(chord1,chord1Hat);
        dot_cAlt_1=dot(chord1,chord1AltHat);
        if abs(dot_c_1)>=dot_cAlt_1
            chord1_new=chord1Hat;
        elseif abs(dot_c_1)<dot_cAlt_1
            chord1_new=chord1Hat;
        end
        
        if dot(chord1_new,chord1)<0
            chord1_new=-chord1_new;
        end
        
        %% correct chord 2
        dot_c_2=dot(chord2,chord2Hat);
        dot_cAlt_2=dot(chord2,chord2AltHat);
        if abs(dot_c_2)>=dot_cAlt_2
            chord2_new=chord2Hat;
        elseif abs(dot_c_2)<dot_cAlt_2
            chord2_new=chord2Hat;
        end
        
        if dot(chord2_new,chord2)<0
            chord2_new=-chord2_new;
        end
        chord1=chord1_new;
        chord2=chord2_new;
    else
        disp('Error in code')
        break
    end
    
    %% load 3D data
    reconpath = fullfile(init.folders.reconstruction, ['frame_' num2str(ii) '.mat']);
    load(reconpath, 'wingRV');
    load(reconpath, 'wingLV');
    load(reconpath, 'BodyRecV');
    wing1Coords=wingRV;
    wing2Coords=wingLV;
    bodyCoords=BodyRecV;
    %% intialize figure
    if ii==sI
        h=figure('position',[ 424   447   900   500]) ;
        %az = -10 ;% 265 ;
        
        % calculate axis for dispaly
        axMat = zeros(3,6) ;
        axMat(1,:) = [min(bodyCoords(:,1)) max(bodyCoords(:,1)) min(bodyCoords(:,2)) max(bodyCoords(:,2)) min(bodyCoords(:,3)) max(bodyCoords(:,3))  ] ;
        axMat(2,:) = [min(wing1Coords(:,1)) max(wing1Coords(:,1)) min(wing1Coords(:,2)) max(wing1Coords(:,2)) min(wing1Coords(:,3)) max(wing1Coords(:,3))  ] ;
        axMat(3,:) = [min(wing2Coords(:,1)) max(wing2Coords(:,1)) min(wing2Coords(:,2)) max(wing2Coords(:,2)) min(wing2Coords(:,3)) max(wing2Coords(:,3))  ] ;
        ax        =  [min(axMat(:,1))/1.2 1.2*max(axMat(:,2)) min(axMat(:,3))/1.2 1.2*max(axMat(:,4)) min(axMat(:,5))/1.2 1.2*max(axMat(:,6)) ] ;
        clear axMat ;
        
        colmap = [0 0.8 0 ; 1 0 0 ; 0 0 1 ];
    end
    %% plot the fly with the chord
    % plot clustering results
    figure(h);
    colmap = [0 0.8 0 ; 1 0 0 ; 0 0 1 ];
    clf;
    hold on
    
   scatter3(bodyCoords(:,1),bodyCoords(:,2),bodyCoords(:,3),'g') ;
    scatter3(wing1Coords(:,1),wing1Coords(:,2),wing1Coords(:,3),'r') ;
    scatter3(wing2Coords(:,1),wing2Coords(:,2),wing2Coords(:,3),'b') ;
    % plot centroids of three clusters
    newCentroids=[mean(bodyCoords); mean(wing1Coords); mean(wing2Coords)];
    plot3(newCentroids(:,1), newCentroids(:,2), newCentroids(:,3), ...
        'ko','markerfacecolor','k','markersize',12) ;
    
    A = 60 * (45/45) ;
    B = 30 * (45/45) ;

    % plot right wing chord
    xvec = [0 B*chord1(1)] + newCentroids(2,1) ;
    yvec = [0 B*chord1(2)] + newCentroids(2,2) ;
    zvec = [0 B*chord1(3)] + newCentroids(2,3) ;
    plot3(xvec, yvec, zvec, 'ko-','linewidth',3,'markersize',8, 'markerfacecolor',colmap(2,:)) ;
     
    % plot left wing chord
    xvec = [0 B*chord2(1)] + newCentroids(3,1) ;
    yvec = [0 B*chord2(2)] + newCentroids(3,2) ;
    zvec = [0 B*chord2(3)] + newCentroids(3,3) ;
    plot3(xvec, yvec, zvec, 'ko-','linewidth',3,'markersize',8, 'markerfacecolor',colmap(3,:)) ;
    
    
    hold off ;
    axis equal ;
    grid on ;
    box on ;
    
    %view(92,8) ;
    view([159 10]);
    
    axis(ax) ;
    
    refresh ;
end