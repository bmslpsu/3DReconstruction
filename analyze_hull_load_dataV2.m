function []=analyze_hull_load_dataV2(init,sI,eI,plot_data,correction_flag)
%% Wael's version of the hull analysis code. Lots of body analysis is removed since it's not needed
% This function saves the chord and span vectors which will be corrected
% by another code
% Suggested work:
% - Clear variables that aren't used
% -Delete unused variables

%% Generate directory for hull analysis
hull_analysis_dir=fullfile(init.folders.root,"hull_analysis");
mkdir(hull_analysis_dir);

%% Initial defintions
chordFraction     = 0.33 ; % fraction of the chord voxels used to find chord
delta             = 2.0;  % strip width used in finding the wing chord

wingTipVelocityThreshold = 5 ;  % voxels/frame

wingLength = 45 ; %(wing1Length + wing2Length)/2 * params.pixPerCM / 232
LL = .55 ; %0.66
LL = wingLength * LL ; % used with farthestPoint

Nimages           = eI ;
%% Initialize the variables to be analyzed
% INITIALIZE COORDINATE VARIABLES. KEEP SAME NAMES AS IN ORIGINAL VERSION
% OF THE ALGORITHM (UNLINE THAT VERSION, HERE WE PRE-ALLOCATE).
% IN THE MEANTIME, UNITS ARE PIXELS. NEED TO CONVERT TO REAL-WORLD UNITS AT
% SOME POINT (MULTIPLY VOXELSIZE)
% -----------------------------------------------------------------------
xs = zeros(Nimages,1); % body center of mass
ys = zeros(Nimages,1);
zs = zeros(Nimages,1);

x1s = zeros(Nimages,1); % wing 1 center of mass
y1s = zeros(Nimages,1);
z1s = zeros(Nimages,1);

phi1s = zeros(Nimages,1); % right wing Euler angles
theta1s = zeros(Nimages,1);
eta1s = zeros(Nimages,1);

x2s = zeros(Nimages,1); % wing 2 center of mass
y2s = zeros(Nimages,1);
z2s = zeros(Nimages,1);

phi2s   = zeros(Nimages,1); % left wing Euler angles
theta2s = zeros(Nimages,1);
eta2s   = zeros(Nimages,1);

rightWingTips = zeros(Nimages,3) ;
leftWingTips  = zeros(Nimages,3) ;

% Memory allocations for unit vectors:
% ------------------------------------
rollHats    = zeros(Nimages,3) ;

psiHats     = zeros(Nimages,3) ;
normRolls   = zeros(Nimages,3) ; %#ok<NASGU> % currently allocated just to avoid bugs in gui

diag11Right = zeros(Nimages,1);
diag12Right = zeros(Nimages,1);
diag21Left = zeros(Nimages,1);
diag22Left = zeros(Nimages,1);
% more memory allocations
rightChordTopProjections = zeros(Nimages,1) ;
leftChordTopProjections  = zeros(Nimages,1) ;
%% begin hull analysis
for ii=sI:eI
    %% load the wing and body hulls
    disp(['frame: ' num2str(ii)])
    %% load 3D data
    reconpath = fullfile(init.folders.reconstruction, ['frame_' num2str(ii) '.mat']);
    load(reconpath, 'wingRV');
    load(reconpath, 'wingLV');
    load(reconpath, 'BodyRecV');
    wing1Coords=wingRV;
    wing2Coords=wingLV;
    bodyCoords=BodyRecV;
    %% Plotting flag and basic plotting
    %% plotting and plot flag
    if plot_data && ii==sI
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
    %% initialize centroid
    newCentroids = zeros(3,3) ; % (x,y,z) for 3 clusters (rows)
    % calculate centroids
    newCentroids(1,:) = mean(bodyCoords) ;
    cBody = newCentroids(1,:);
    xs(ii)   = newCentroids(1,1) ;
    ys(ii)   = newCentroids(1,2) ;
    zs(ii)   = newCentroids(1,3) ;
    %% Start wing 1
    if ~isempty(wing1Coords)
        wing1LargestCC_initial = findLargestHullCC (wing1Coords); %WS: might want to get rid of this one
        
        % from the largest CC, find the farpoint and then select the voxels
        % that are close to the farpoint.
        [farPoint1, ~, list1] = farthestPointWaelV1(wing1LargestCC_initial, cBody, LL) ;
        
        % find the largest connected component within the voxels in list1.
        % list1 are the indices of the voxles in wing1LargestCC_initial whose distance
        % from "farPoint" is smaller than LL.
        %
        % largest component is calculated AGAIN, to account for a case where cutting
        % the LL distance left us with two connected components
        wing1LargestCC = findLargestHullCC (wing1LargestCC_initial(list1,:));
        
        clear wing1LargestCC_initial
        
        % make sure the far point is included in the connected components
        % if not, take another one
        if (~ismember(farPoint1, wing1LargestCC,'rows'))
            %keyboard
            
            tmp = setdiff( wing1Coords(list1,:), wing1LargestCC, 'rows') ;
            wing1LargestCC = findLargestHullCC(tmp) ;
            clear tmp
        end
        newCentroids(2,:) = ...
            calcCentroidFromTopView(wing1LargestCC, farPoint1, wingLength/3) ;
        
    end
    
    %% Start wing 2
    if ~isempty(wing2Coords)
        % first find the largest connected component of the wing voxels
        % (first use of this function)
        wing2LargestCC_initial = findLargestHullCC (wing2Coords);
        
        % from the largest CC, find the farpoint and then select the voxels
        % that are close to the farpoint.
        [farPoint2, ~, list2] = farthestPointWaelV1(wing2LargestCC_initial, cBody, LL) ;
        
        % find the largest connected component within the voxels in list2.
        % list2 are the indices of the voxles in wing2LargestCC_initial whose distance
        % from "farPoint" is smaller than LL.
        %
        % largest component is calculated AGAIN, to account for a case where cutting
        % the LL distance left us with two connected components
        wing2LargestCC = findLargestHullCC (wing2LargestCC_initial(list2,:));
        
        clear wing2LargestCC_initial
        % make sure the far point is included in the connected components
        % if not, take another one
        if (~ismember(farPoint2, wing2LargestCC,'rows'))
            %keyboard
            tmp = setdiff( wing2Coords(list2,:), wing2LargestCC, 'rows') ;
            wing2LargestCC = findLargestHullCC(tmp) ;
            clear tmp
        end
        
        newCentroids(3,:) = ...
            calcCentroidFromTopView(wing2LargestCC, farPoint2, wingLength/3) ;
        
        % combine three hulls into "coords" and specify bodyRows, wing1Rows and
        % wing2Rows
        n1 = size(bodyCoords,1) ;
        n2 = size(wing1Coords,1) ;
        n3 = size(wing2Coords,1) ;
        
        % note there might be a duplicity between "coords", "wingXLargestCC"
        coords = [ bodyCoords ; wing1Coords ; wing2Coords ] ;
        nn = n1+n2+n3 ;
        % use logical indexing (should be faster than numeric)
        bodyRows  = false(nn,1) ; bodyRows(1:n1)                  = true ; % bodyRows  = 1:n1 ;
        wing1Rows = false(nn,1) ; wing1Rows((n1+1):(n1+n2))       = true ;
        wing2Rows = false(nn,1) ; wing2Rows((n1+n2+1):(n1+n2+n3)) = true ;
        
        clear n1 n2 n3 nn
    end
    %% define center of mass for each body
    cRight=newCentroids(2,:);
    cLeft=newCentroids(3,:);
    %% Find wing 1 tip,span, chord
    farPointR    = farthestPointWaelV1(wing1LargestCC, cBody, LL) ;
    spanHatRight = farPointR - cRight ; %vector of the span of right wing
    spanHatRight = spanHatRight' ;
    spanHatRight = spanHatRight / norm(spanHatRight) ; %normalize span
    
    % makes wing span point outward
    if dot(cRight-cBody,spanHatRight) < 0
        spanHatRight = -spanHatRight;
    end
    
    % FIND WING TIP
    % -------------
    
    rightWingTip = findWingTip(wing1LargestCC, spanHatRight', cRight);
    if (isnan(rightWingTip(1)))
        rightWingTip = farPointR ;
    end
    
    % recalculate span vector based on the refined wing tip
    spanHatRight = rightWingTip - cRight ;
    spanHatRight = spanHatRight' ;
    spanHatRight = spanHatRight / norm(spanHatRight) ;
    
    span1Hat     = spanHatRight ; % get rid of this duplicity
    
    
    % store right wing center of mass (centroid)
    x1s(ii) = cRight(1) ; % newCentroids(rightWingInd,1) ;
    y1s(ii) = cRight(2) ; % newCentroids(rightWingInd,2) ;
    z1s(ii) = cRight(3) ; % newCentroids(rightWingInd,3) ;
    
    % right wing phi and theta.
    % NOTE: here phi is with respect to the x axis in LAB FRAME OF REF.
    phi1   = atan2(spanHatRight(2),spanHatRight(1));
    theta1 = asin(spanHatRight(3));
    
    phi1Proj = [span1Hat(1),span1Hat(2),0];
    phi1Proj = phi1Proj / norm(phi1Proj);
    phi1Hat = [-sin(phi1),cos(phi1),0];
    theta1Hat = -phi1Proj*sin(theta1) + [0,0,1]*cos(theta1);
    
    try
        
        rightWingVoxels = double(wing1LargestCC) ; %converts voxels from int to double
        Nvox            = size(rightWingVoxels,1) ;
        
        mat1 = rightWingVoxels - repmat(cRight, Nvox,1) ;
        mat2 = repmat(spanHatRight', Nvox, 1) ;
        
        distFromMidSpan = abs(sum(mat1.*mat2,2) ) ;
        clear mat1 mat2
        
        chordRowsInd = find(distFromMidSpan<delta) ;
        
        if (isempty(chordRowsInd))
            % first try a larger delta
            chordRowsInd = find(distFromMidSpan<3*delta) ;
            % check if still empty
            if (isempty(chordRowsInd))
                error('hullAnalysis:Chord','Bad clustering - empty right chord') ;
            end
        end
        clear distFromMidSpan
        % among the chordRowsInd voxels, find the pair that is the most distant
        % apart. this pair defines the direction of the chord
        % original code calculated a distance matrix between all pairs. here we
        % do it a bit more efficiently by excluding most of the voxels before
        % calculating the distance matrix
        
        % caculate the distance^2 of each voxel from the wing centroid
        Nvox    = length(chordRowsInd);
        mat1    = (rightWingVoxels(chordRowsInd,:) - repmat(cRight, Nvox,1)).^2 ;
        distVec = (sum(mat1,2)).^0.5 ;
        clear mat1
        
        % select only the top quarter of the voxels, i.e the most distant from
        % wing centroid
        
        [~, sortedInd] = sort(distVec,'descend') ;
        selectedInd    = chordRowsInd(sortedInd(1:ceil(Nvox*chordFraction))) ;
        %selectedIndRight = selectedInd ;
        clear distVec
        
        % find the most distant pair
        
        distMat = squareform (pdist (rightWingVoxels(selectedInd,:))) ;
        
        [maxRowVec Irow] = max(distMat,[],1) ;
        [~, Icol] = max(maxRowVec) ;
        Irow = Irow(Icol) ;
        
        if (distMat(Irow, Icol)~=max(distMat(:)))
            disp('Error with finding max. plz check.') ;
            disp('problem 6?') ;
            %keyboard ;
        end
        
        
        % the following indices give voxel coordinate:
        vox1IndRight = selectedInd(Irow) ; % voxel position is rightWingVoxels(vox1IndRight,:)
        vox2IndRight = selectedInd(Icol) ; % voxel position is rightWingVoxels(vox2IndRight,:)
        
        chord1Hat = rightWingVoxels(vox1IndRight,:)' - rightWingVoxels(vox2IndRight,:)' ;
        
        % force chord to be vertical to the span vector
        chord1Hat = chord1Hat - span1Hat * dot(span1Hat, chord1Hat) ;
        
        diag11 = norm(chord1Hat) ; % will be used later
        
        chord1Hat = chord1Hat / norm(chord1Hat) ;
        
        
        % FIND THE SECOND DIAGONAL OF THE RIGHT WING PARALLELOGRAM (ADDED FEB 2012)
        % ------------------------------------------------------------------------
        
        % see http://mathworld.wolfram.com/Point-PlaneDistance.html
        mat1    = rightWingVoxels(chordRowsInd,:) - repmat(cRight, Nvox,1) ;
        % mat1 contains the vectors connecting the voxels
        % in chordRowsInd to the wing centoid.
        
        % calc the vector normal to the span and chord
        wingNormVec = cross(span1Hat, chord1Hat) ;
        
        % calc the (signed) distance from each point in mat1 to the wing plane
        mat2 = repmat(wingNormVec', Nvox, 1) ;
        distVec = sum(mat1.*mat2,2) ;
        
        % find the largest positive and largest negative distances, which
        % correspond to the farthest voxels on each side of the plane
        
        [maxval, indmax] = max(distVec) ; % index into rightWingVoxels
        [minval, indmin] = min(distVec) ; % index into rightWingVoxels
        
        if (maxval<=0 || minval>=0)
            disp('error in finding alternative chord vector for right wing') ;
            %keyboard ;
        end
        
        indmin = chordRowsInd(indmin) ;
        indmax = chordRowsInd(indmax) ;
        chord1AltHat = rightWingVoxels(indmax,:)' - rightWingVoxels(indmin,:)' ;
        
        % force alternative chord to be vertical to the span vector
        chord1AltHat = chord1AltHat - span1Hat * dot(span1Hat, chord1AltHat) ;
        
        diag12 = norm(chord1AltHat) ; % will be used later
        % now normalize
        chord1AltHat = chord1AltHat / norm(chord1AltHat) ;
        
        % choose the sign of the alternative chord vector such that it is
        % has positive overlap with the "main" chord vector
        
        %if ( dot(chord1AltHat, chord1Hat) < 0 )
        %    chord1AltHat = - chord1AltHat ;
        %end
        
        
        diag11Right(ii) = diag11 ;
        diag12Right(ii) = diag12 ;
        
        
        % if one of the diagonals is siginficantly longer, choose the longer
        % one and do not proceed to the velocity criterion below
        diagSwapFlag     = false ;
        velocitySwapFlag = false ;
        
        % find wingtip 'velocity' with respect to the body
        if(ii>sI)
            % previous version calculate wing centroid velocity:
            vWingCM = [  (x1s(ii) - xs(ii)) - (x1s(ii-1)-xs(ii-1)) ; ...
                (y1s(ii) - ys(ii)) - (y1s(ii-1)-ys(ii-1))  ; ...
                (z1s(ii) - zs(ii)) - (z1s(ii-1)-zs(ii-1)) ] ;
            
            vWing =  ( rightWingTip - [xs(ii) ys(ii) zs(ii)] ) - ...
                ( rightWingTips(ii-1,:) - [xs(ii-1) ys(ii-1) zs(ii-1)]) ;
            
            % keep only the component perpendicular to the span vector
            vWing = vWing - span1Hat' * dot(span1Hat, vWing) ;
            
            
            nrm = norm(vWing) ;
            
        else
            nrm = 0 ;
        end
        
        if (ii>sI)
            if (nrm~=0)
                vWing = vWing / nrm ;
                dot1 = dot(chord1Hat, vWing) ;
                dot2 = dot(chord1AltHat, vWing) ;
                
                if (dot2>dot1) % swap
                    velocitySwapFlag = true ;
                end
                if (dot1<0 && dot2<0 && nrm>=wingTipVelocityThreshold && ~velocitySwapFlag)
                    chord1Hat = - chord1Hat ;
                end
            end
        end
        
        swapFlag = (velocitySwapFlag && nrm>=wingTipVelocityThreshold) || ... % believe velocity if |v|>2
            (diagSwapFlag && nrm<wingTipVelocityThreshold) ;
        
        % probably need a smarter way to "weigh" the two types of swaps
        % see, e.g. frame 14 where the wings should have been swapped
        if (swapFlag)
            tmp = chord1Hat ;
            chord1Hat = chord1AltHat ;
            chord1AltHat = tmp ;
            
            tmp = diag11Right(ii) ;
            diag11Right(ii) = diag12Right(ii) ;
            diag12Right(ii) = tmp ;
            
            clear tmp ;
        end
        
        % ----------
        
        
        % makes chord look good on right wing during downstroke
        %chord1Hat = abs(dot(chord1Hat,phi1Hat))*phi1Hat + abs(dot(chord1Hat,theta1Hat))*theta1Hat;
        
        chord1Hat2 = -abs(dot(chord1Hat,phi1Hat))*phi1Hat + abs(dot(chord1Hat,theta1Hat))*theta1Hat;
  
        
        
    catch %#ok<CTCH>
        disp('error occured while trying to find right chord. taking values from the previous frame.')
        % this ignores the case where the error occures on the first frame...
        vox1IndRight = [] ;
        vox2IndRight = [] ;
    end
    
    [p1R p2R rightChordTopProj] = chordTopViewProjection(spanHatRight, cRight, wing1LargestCC) ;
    
    %% Wing 2 tip, span, and chord
    % use a new "definition" of the span vector, as the unit vector from
    % the wing's center of mass to its farthest point.
    farPointL = farthestPointWaelV1(wing2LargestCC, cBody, LL) ;
    spanHatLeft = farPointL - cLeft ;
    spanHatLeft = spanHatLeft' ;
    spanHatLeft = spanHatLeft / norm(spanHatLeft) ;
    
    % makes wing span point outward
    if dot(cLeft-cBody,spanHatLeft) < 0
        spanHatLeft = -spanHatLeft;
    end
    
    % FIND WING TIP
    % -------------
    
    leftWingTip = findWingTip(wing2LargestCC, spanHatLeft', cLeft);
    if (isnan(leftWingTip(1)))
        leftWingTip = farPointL ;
    end
    
    % recalculate span vector based on the refined wing tip
    spanHatLeft = leftWingTip - cLeft ;
    spanHatLeft = spanHatLeft' ;
    spanHatLeft = spanHatLeft / norm(spanHatLeft) ;
    
    span2Hat    = spanHatLeft ; % get rid of this duplicity
    
    % store right wing center of mass (centroid)
    x2s(ii) = cLeft(1) ; % newCentroids(rightWingInd,1) ;
    y2s(ii) = cLeft(2) ; % newCentroids(rightWingInd,2) ;
    z2s(ii) = cLeft(3) ; % newCentroids(rightWingInd,3) ;
    
    % right wing phi and theta.
    % NOTE: here phi is with respect to the x axis in LAB FRAME OF REF.
    phi2   = atan2(spanHatLeft(2),spanHatLeft(1));
    theta2 = asin(spanHatLeft(3));
    
    phi2Proj  = [span2Hat(1),span2Hat(2),0];
    phi2Proj  = phi2Proj / norm(phi2Proj);
    phi2Hat   = [-sin(phi2),cos(phi2),0];
    theta2Hat = -phi2Proj*sin(theta2) + [0,0,1]*cos(theta2);
    try
        %leftWingVoxels = coords(leftWingRows,:) ;
        leftWingVoxels  = double(wing2LargestCC) ;
        
        Nvox            = size(leftWingVoxels,1) ;
        
        mat1 = leftWingVoxels - repmat(cLeft, Nvox,1) ;
        mat2 = repmat(spanHatLeft', Nvox, 1) ;
        
        distFromMidSpan = abs(sum(mat1.*mat2,2) ) ;
        clear mat1 mat2
        
        chordRowsInd = find(distFromMidSpan<delta) ;
        
        if (isempty(chordRowsInd))
            % first try a larger delta
            chordRowsInd = find(distFromMidSpan<3*delta) ;
            % check if still empty
            if (isempty(chordRowsInd))
                error('hullAnalysis:Chord','Bad clustering - empty left chord') ;
            end
        end
        
        clear distFromMidSpan
        
        % among the chordRowsInd voxels, find the pair that is the most distant
        % apart. this pair defines the direction of the chord
        % original code calculated a distance matric between all pairs. here we
        % do it a bit more efficiently by excluding most of the voxels before
        % calculating the distance matrix
        
        % caculate the distance^2 of each voxel from the wing centroid
        Nvox    = length(chordRowsInd);
        mat1    = (leftWingVoxels(chordRowsInd,:) - repmat(cLeft, Nvox,1)).^2 ;
        distVec = (sum(mat1,2)).^0.5 ;
        clear mat1
        
        % select only the top quarter of the voxels, i.e the most distant from
        % wing centroid
        
        [~, sortedInd] = sort(distVec,'descend') ;
        selectedInd    = chordRowsInd(sortedInd(1:ceil(Nvox*chordFraction))) ;
        clear distVec
        
        %selectedIndLeft = selectedInd ;
        % find the most distant pair
        
        distMat = squareform (pdist (leftWingVoxels(selectedInd,:))) ;
        
        [maxRowVec Irow] = max(distMat,[],1) ;
        [~, Icol] = max(maxRowVec) ;
        Irow = Irow(Icol) ;
        
        
        if (distMat(Irow, Icol)~=max(distMat(:)))
            disp('Error with finding max. plz check.') ;
            %disp('problem 7?') ;
            %keyboard ;
        end
        
        % the following indices give voxel coordinate:
        vox1IndLeft = selectedInd(Irow) ; % voxel position is leftWingVoxels(vox1IndLeft,:)
        vox2IndLeft = selectedInd(Icol) ; % voxel position is leftWingVoxels(vox2IndLeft,:)
        
        chord2Hat = leftWingVoxels(vox1IndLeft,:)' - leftWingVoxels(vox2IndLeft,:)' ;
        
        p1 = leftWingVoxels(vox1IndLeft,:)' ;
        p2 = leftWingVoxels(vox2IndLeft,:)' ;
        
        % force chord vector to be perpendicular to the span vector
        chord2Hat = chord2Hat - span2Hat * dot(span2Hat, chord2Hat) ;
        
        diag21 = norm(chord2Hat) ; % will be used later
        
        chord2Hat = chord2Hat / norm(chord2Hat) ;
        
        if (chord2Hat(3)<0)
            chord2Hat = - chord2Hat ;
        end
        
        % makes chord look good on left wing during upstroke
        
        %disp ('check this vector manipulation');
        % chord2Hat = abs(dot(chord2Hat,phi2Hat))*phi2Hat + abs(dot(chord2Hat,theta2Hat))*theta2Hat;
        
        % FIND THE SECOND DIAGONAL OF THE LEFT WING PARALLELOGRAM (ADDED FEB 2012)
        % ------------------------------------------------------------------------
        
        % see http://mathworld.wolfram.com/Point-PlaneDistance.html
        mat1    = leftWingVoxels(chordRowsInd,:) - repmat(cLeft, Nvox,1) ;
        % mat1 contains the vectors connecting the voxels
        % in chordRowsInd to the wing centoid.
        
        % calc the vector normal to the span and chord
        wingNormVec = cross(span2Hat, chord2Hat) ;
        
        % calc the (signed) distance from each point in mat1 to the wing plane
        mat2 = repmat(wingNormVec', Nvox, 1) ;
        distVec = sum(mat1.*mat2,2) ;
        
        % find the largest positive and largest negative distances, which
        % correspond to the farthest voxels on each side of the plane
        
        [maxval, indmax] = max(distVec) ; % index into rightWingVoxels
        [minval, indmin] = min(distVec) ; % index into rightWingVoxels
        
        if (maxval<=0 || minval>=0)
            disp('error in finding alternative chord vector for left wing') ;
            %keyboard ;
        end
        
        indmin = chordRowsInd(indmin) ;
        indmax = chordRowsInd(indmax) ;
        chord2AltHat = leftWingVoxels(indmax,:)' - leftWingVoxels(indmin,:)' ;
        
        p3 = leftWingVoxels(indmax,:)' ;
        p4 = leftWingVoxels(indmin,:)' ;
        
        % force alternative chord to be vertical to the span vector
        chord2AltHat = chord2AltHat - span2Hat * dot(span2Hat, chord2AltHat) ;
        
        diag22 = norm(chord2AltHat) ; % will be used later
        
        % now normalize
        chord2AltHat = chord2AltHat / norm(chord2AltHat) ;
        
        % choose the sign of the alt vector
        %if ( dot(chord2AltHat, chord2Hat) < 0 )
        %    chord2AltHat = - chord2AltHat ;
        %end
        
        if (chord2AltHat(3)<0)
            chord2AltHat = - chord2AltHat ;
        end
        
        diag21Left(ii) = diag21 ;
        diag22Left(ii) = diag22 ;
        
        % if one of the diagonals is siginficantly longer, choose the longer
        % one and do not proceed to the velocity criterion below
        diagSwapFlag     = false ;
        velocitySwapFlag = false ;
        zSwapFlag        = false ;
        
        % for the z-swap - decide to swap if alt-chord has larger z range
        
        if ( max([ p3(3), p4(3)]) > max([ p1(3), p2(3)]) && ...
                min([ p3(3), p4(3)]) < min([ p1(3), p2(3)]) )
            zSwapFlag = true ;
        end
        
        if (diag22/diag21 >= 1.3)
            diagSwapFlag = true ;
        end
        
        % find wingtip 'velocity' with respect to the body
        if(ii>1)
            % previous version calculate wing centroid velocity:
            vWingCM = [  (x2s(ii) - xs(ii)) - (x2s(ii-1)-xs(ii-1)) ; ...
                (y2s(ii) - ys(ii)) - (y2s(ii-1)-ys(ii-1))  ; ...
                (z2s(ii) - zs(ii)) - (z2s(ii-1)-zs(ii-1)) ] ;
            
            vWing =  ( leftWingTip - [xs(ii) ys(ii) zs(ii)] ) - ...
                ( leftWingTips(ii-1,:) - [xs(ii-1) ys(ii-1) zs(ii-1)]) ;
            
            % keep only the component perpendicular to the span vector
            vWing = vWing - span2Hat' * dot(span2Hat, vWing) ;
            
            nrm = norm(vWing) ;
        else
            nrm = 0 ;
        end
        
        if (ii>1)
            if (nrm~=0)
                vWing = vWing / nrm ;
                dot1 = dot(chord2Hat, vWing) ;
                dot2 = dot(chord2AltHat, vWing) ;
                
                if (dot2>dot1) % swap
                    velocitySwapFlag = true ;
                end
                if (dot1<0 && dot2<0 && nrm>=wingTipVelocityThreshold && ~velocitySwapFlag)
                    chord2Hat = - chord2Hat ;
                end
            end
        end
        
        swapFlag = (velocitySwapFlag && nrm>=wingTipVelocityThreshold) || ... % believe velocity if |v|>2
            (diagSwapFlag && nrm<wingTipVelocityThreshold) ;
        
        %swapFlag = zSwapFlag || (diagSwapFlag && nrm<wingTipVelocityThreshold) ;
        
        % probably need a smarter way to "weigh" the two types of swaps
        if (swapFlag)
            tmp = chord2Hat ;
            chord2Hat = chord2AltHat ;
            chord2AltHat = tmp ;
            
            tmp = diag21Left(ii) ;
            diag21Left(ii) = diag22Left(ii) ;
            diag22Left(ii) = tmp ;
            clear tmp ;
        end
        
        
        [p1L, p2L, leftChordTopProj] = chordTopViewProjection(spanHatLeft, cLeft, wing2LargestCC) ;
        
        % ----------
        
        chord2Hat2 = -abs(dot(chord2Hat,phi2Hat))*phi2Hat + abs(dot(chord2Hat,theta2Hat))*theta2Hat;
        % in the first quadrant only
        %         eta2 = atan2( abs(dot(chord2Hat,theta2Hat)), abs(dot(chord2Hat,phi2Hat)));
        
        r_rot = vrrotvec([1,0,0],span2Hat);
        m_rot = vrrotvec2mat(r_rot);
        y_wing2=[0;1;0];
        z_wing2=cross(y_wing2,span2Hat);
        z_wing2=z_wing2/norm(z_wing2);
        wingNormVecWael2 = cross(chord2Hat,span2Hat);
        %project z axis on the span-normal plane
        
        z_wing_proj2=z_wing2-dot(z_wing2,span2Hat)*span2Hat;
        z_wing_proj2=z_wing_proj2/norm(z_wing_proj2);
        
    catch %#ok<CTCH>
        disp('error occured while trying to find left chord. taking values from the previous frame.')
        % this ignores the case where the error occures on the first frame...
        vox1IndLeft = [] ;
        vox2IndLeft = [] ;
    end
    
    
    %% Plot the data
    if (plot_data) %&& (ind > 604)
        % plot clustering results
        figure(h);
        colmap = [0 0.8 0 ; 1 0 0 ; 0 0 1 ];
        clf;
        hold on
        % color of inner side of the wings
        plot3(coords(bodyRows,1),coords(bodyRows,2),coords(bodyRows,3),'o','markersize',2,'color',colmap(1,:)) ;
        plot3(coords(wing1Rows,1),coords(wing1Rows,2),coords(wing1Rows,3),'o','markersize',2,'color',colmap(2,:)) ;
        plot3(coords(wing2Rows,1),coords(wing2Rows,2),coords(wing2Rows,3),'o','markersize',2,'color',colmap(3,:)) ;
        
        % plot centroids of three clusters
        plot3(newCentroids(:,1), newCentroids(:,2), newCentroids(:,3), ...
            'ko','markerfacecolor','k','markersize',12) ;
        
        try
            % plot wing tips
            plot3(rightWingTip(1), rightWingTip(2), rightWingTip(3), ...
                'ks','markerfacecolor','k','markersize',10) ;
        catch %#ok<CTCH>
            disp('problem 8?') ;
            %keyboard ;
        end
        
        %{
            % plot chord voxels
            plot3( leftWingVoxels(selectedIndLeft,1),leftWingVoxels(selectedIndLeft,2), ...
                leftWingVoxels(selectedIndLeft,3),'ko') ;
            plot3( rightWingVoxels(selectedIndRight,1),rightWingVoxels(selectedIndRight,2), ...
                rightWingVoxels(selectedIndRight,3),'ko') ;
        %}
        
        A = 60 * (wingLength/45) ;
        B = 30 * (wingLength/45) ;
        % plot body axis
        %         xvec = [0 A*rollHat(1)] + newCentroids(1,1) ;
        %         yvec = [0 A*rollHat(2)] + newCentroids(1,2) ;
        %         zvec = [0 A*rollHat(3)] + newCentroids(1,3) ;
        %         plot3(xvec, yvec, zvec, 'ko-','linewidth',3,'markersize',8,'markerfacecolor',colmap(bodyInd,:)) ;
        
        % plot right wing span
        xvec = [0 A*span1Hat(1)] + newCentroids(2,1) ;
        yvec = [0 A*span1Hat(2)] + newCentroids(2,2) ;
        zvec = [0 A*span1Hat(3)] + newCentroids(2,3) ;
        plot3(xvec, yvec, zvec, 'ko-','linewidth',3,'markersize',8, 'markerfacecolor',colmap(2,:)) ;
        
        % plot right wing chord
        xvec = [0 B*chord1Hat(1)] + newCentroids(2,1) ;
        yvec = [0 B*chord1Hat(2)] + newCentroids(2,2) ;
        zvec = [0 B*chord1Hat(3)] + newCentroids(2,3) ;
        plot3(xvec, yvec, zvec, 'ko-','linewidth',3,'markersize',8, 'markerfacecolor',colmap(2,:)) ;
        
        % plot right wing alternative chord
        xvec = [0 B*chord1AltHat(1)] + newCentroids(2,1) ;
        yvec = [0 B*chord1AltHat(2)] + newCentroids(2,2) ;
        zvec = [0 B*chord1AltHat(3)] + newCentroids(2,3) ;
        plot3(xvec, yvec, zvec, 'ks--','linewidth',2,'markersize',8, 'markerfacecolor',colmap(2,:)) ;
        
        % plot left wing span
        xvec = [0 A*span2Hat(1)] + newCentroids(3,1) ;
        yvec = [0 A*span2Hat(2)] + newCentroids(3,2) ;
        zvec = [0 A*span2Hat(3)] + newCentroids(3,3) ;
        plot3(xvec, yvec, zvec, 'ko-','linewidth',3,'markersize',8, 'markerfacecolor',colmap(3,:)) ;
        
        % plot left wing chord
        xvec = [0 B*chord2Hat(1)] + newCentroids(3,1) ;
        yvec = [0 B*chord2Hat(2)] + newCentroids(3,2) ;
        zvec = [0 B*chord2Hat(3)] + newCentroids(3,3) ;
        plot3(xvec, yvec, zvec, 'ko-','linewidth',3,'markersize',8, 'markerfacecolor',colmap(3,:)) ;
        
        % plot left wing alternative chord
        xvec = [0 B*chord2AltHat(1)] + newCentroids(3,1) ;
        yvec = [0 B*chord2AltHat(2)] + newCentroids(3,2) ;
        zvec = [0 B*chord2AltHat(3)] + newCentroids(3,3) ;
        plot3(xvec, yvec, zvec, 'ks--','linewidth',2,'markersize',8, 'markerfacecolor',colmap(3,:)) ;
        
        plot3( leftWingVoxels(vox1IndLeft,1), leftWingVoxels(vox1IndLeft,2), ...
            leftWingVoxels(vox1IndLeft,3), 'k^','markerfacecolor','y') ;
        plot3( leftWingVoxels(vox2IndLeft,1), leftWingVoxels(vox2IndLeft,2), ...
            leftWingVoxels(vox2IndLeft,3), 'k^','markerfacecolor','y') ;
        
        plot3( rightWingVoxels(vox1IndRight,1), rightWingVoxels(vox1IndRight,2), ...
            rightWingVoxels(vox1IndRight,3), 'k^','markerfacecolor','y') ;
        plot3( rightWingVoxels(vox2IndRight,1), rightWingVoxels(vox2IndRight,2), ...
            rightWingVoxels(vox2IndRight,3), 'k^','markerfacecolor','y') ;
        
        % plot farthest point for both wings
        
        %[p, ~, list] = farthestPoint(rightWingVoxels, cBody, LL) ;
        farPointRight = farPoint1 ; %WS: has to be added because we do the clustering before the analysis
        farPointLeft  = farPoint2 ;
        plot3(farPointRight(1), farPointRight(2), farPointRight(3),'ks','markersize',14,'markerfacecolor','y','linewidth',2) ;
        %plot3(rightWingVoxels(list,1), rightWingVoxels(list,2), rightWingVoxels(list,3),'o','markersize',2,'color',[0.75 0 0]) ;
        % color of outside of the wing
        plot3(wing1LargestCC(:,1), wing1LargestCC(:,2), wing1LargestCC(:,3),'o','markersize',2,'color',[0.6 0 0]) ;
        
        %[p, ~, list] = farthestPoint(leftWingVoxels, cBody, LL) ;
        plot3(farPointLeft(1), farPointLeft(2), farPointLeft(3),'ks','markersize',14,'markerfacecolor','y','linewidth',2) ;
        %plot3(leftWingVoxels(list,1), leftWingVoxels(list,2), leftWingVoxels(list,3),'o','markersize',2,'color',[0 0 0.75]) ;
        % color of outside of the wing
        plot3(wing2LargestCC(:,1), wing2LargestCC(:,2), wing2LargestCC(:,3),'o','markersize',2,'color',[0 0 0.6]) ;
        
        xlabel('x') ; ylabel('y'); zlabel('z') ;
        
        hold off ;
        axis equal ;
        grid on ;
        box on ;
        
        %view(92,8) ;
        view([159 10]);
        
        axis(ax) ;
        
        refresh ;
        
    end
    %% GUI corrections for chord
    if ~correction_flag && ii==sI %only corrected for first wing
        disp('Correct left')
        [chord1Hat,chord1AltHat]=CorrectWingChord(chord1Hat,chord1AltHat);
        disp('Correct right')
        [chord2Hat,chord2AltHat]=CorrectWingChord(chord2Hat,chord2AltHat);
    elseif correction_flag %corrected for all wings (not recommeded)
        disp('Correct left')
        [chord1Hat,chord1AltHat]=CorrectWingChord(chord1Hat,chord1AltHat);
        disp('Correct right')
        [chord2Hat,chord2AltHat]=CorrectWingChord(chord2Hat,chord2AltHat);
    end
    %% angle corrections
    
    %% estimate the location of the new orthogonal system
    rotAngle=45;
    rotation_hat=rotx(rotAngle);
    z_hat=rotation_hat*[0; 0; 1];
    y_hat=rotation_hat*[0; 1; 0];
    %% rotation angle 2
    [vec_min2,R_min2, angle_min2]=MinimizeXComponent(span2Hat,y_hat,-pi:0.0001:0);
    [vec_align2, R_aligned2, angle_aligned2]=AlignVectors(vec_min2,[1 0 0]',z_hat);
    chord2_transformed=R_aligned2*R_min2*[chord2Hat ;0];
    chord2_transformed=chord2_transformed(1:end-1);
    [vec_min2_rot,R_min, angle_min]=MinimizeXComponent_SmallAngle(chord2_transformed,vec_align2);
    angle_rotation2_final(ii)=-angle_min;
    angle_rotation2_final_unit=-angle_min;
    %% rotation angle 1
    [vec_min1,R_min1, angle_min1]=MinimizeXComponent(span1Hat,y_hat,0:0.0001:pi); %motion of the wing is in the opp direction of the right wrt yhat
    [vec_align1, R_aligned1, angle_aligned1]=AlignVectors(vec_min1,[1 0 0]',z_hat);
    chord1_transformed=R_aligned1*R_min1*[chord1Hat ;0];
    chord1_transformed=chord1_transformed(1:end-1);
    [vec_min1_rot,~, angle_min1]=MinimizeXComponent_SmallAngle(chord1_transformed,vec_align1);
    angle_rotation1_final(ii)=angle_min1;
    angle_rotation1_final_unit=angle_min1;
        %% rotation angle 1
    y_wing1=[0;1;0];
    z_wing1=cross(span1Hat,y_wing1);
    z_wing1=z_wing1/norm(z_wing1);
    wingNormVecWael1 = cross(span1Hat,chord1Hat);
    %project z axis on the span-normal plane
    
    z_wing_proj=z_wing1-dot(z_wing1,span1Hat)*span1Hat;
    z_wing_proj=z_wing_proj/norm(z_wing_proj);
    
    phi1Proj = [span1Hat(1),span1Hat(2),0];
    phi1Proj = phi1Proj / norm(phi1Proj);
    phi1Hat = [-sin(phi1),cos(phi1),0];
    theta1Hat = -phi1Proj*sin(theta1) + [0,0,1]*cos(theta1);
    eta1 = atan2( abs(dot(chord1Hat,theta1Hat)), abs(dot(chord1Hat,phi1Hat)));
    eta1_vec(ii)=eta1;
    rotation_angle_wael1=EstimateWingRotation(span1Hat,chord1Hat,y_wing1);
    %% rotation angle 2
    
    y_wing2=[0;1;0];
    z_wing2=cross(y_wing2,span2Hat);
    z_wing2=z_wing2/norm(z_wing2);
    wingNormVecWael2 = cross(chord2Hat,span2Hat);
    %project z axis on the span-normal plane
    
    z_wing_proj2=z_wing2-dot(z_wing2,span2Hat)*span2Hat;
    z_wing_proj2=z_wing_proj2/norm(z_wing_proj2);

    rotation_angle_wael2=-EstimateWingRotation(span2Hat,chord2Hat,y_wing2);
    rotation_angle_vec2(ii)=rotation_angle_wael2;
    phi2Proj  = [span2Hat(1),span2Hat(2),0];
    phi2Proj  = phi2Proj / norm(phi2Proj);
    phi2Hat   = [-sin(phi2),cos(phi2),0];
    theta2Hat = -phi2Proj*sin(theta2) + [0,0,1]*cos(theta2);
    eta2 = atan2( abs(dot(chord2Hat,theta2Hat)), abs(dot(chord2Hat,phi2Hat)));
    eta2_vec(ii)=eta2;
    %% second plot 
    if (plot_data)
        % plot clustering results
        figure(h);
        colmap = [0 0.8 0 ; 1 0 0 ; 0 0 1 ];
        clf;
        hold on
        
        plot3(coords(bodyRows,1),coords(bodyRows,2),coords(bodyRows,3),'o','markersize',2,'color',colmap(1,:)) ;
        plot3(coords(wing1Rows,1),coords(wing1Rows,2),coords(wing1Rows,3),'o','markersize',2,'color',colmap(2,:)) ;
        plot3(coords(wing2Rows,1),coords(wing2Rows,2),coords(wing2Rows,3),'o','markersize',2,'color',colmap(3,:)) ;
        % plot centroids of three clusters
        plot3(newCentroids(:,1), newCentroids(:,2), newCentroids(:,3), ...
            'ko','markerfacecolor','k','markersize',12) ;
        
        try
            % plot wing tips
            plot3(rightWingTip(1), rightWingTip(2), rightWingTip(3), ...
                'ks','markerfacecolor','k','markersize',10) ;
        catch %#ok<CTCH>
            disp('problem 8?') ;
            %keyboard ;
        end
        
        A = 60 * (wingLength/45) ;
        B = 30 * (wingLength/45) ;
        
        % plot right wing span
        xvec = [0 A*span1Hat(1)] + newCentroids(2,1) ;
        yvec = [0 A*span1Hat(2)] + newCentroids(2,2) ;
        zvec = [0 A*span1Hat(3)] + newCentroids(2,3) ;
        plot3(xvec, yvec, zvec, 'ko-','linewidth',3,'markersize',8, 'markerfacecolor',colmap(2,:)) ;
        
        % plot right wing chord
        xvec = [0 B*chord1Hat(1)] + newCentroids(2,1) ;
        yvec = [0 B*chord1Hat(2)] + newCentroids(2,2) ;
        zvec = [0 B*chord1Hat(3)] + newCentroids(2,3) ;
        plot3(xvec, yvec, zvec, 'ko-','linewidth',3,'markersize',8, 'markerfacecolor',colmap(2,:)) ;
        
        % plot right wing alternative chord
        xvec = [0 B*chord1AltHat(1)] + newCentroids(2,1) ;
        yvec = [0 B*chord1AltHat(2)] + newCentroids(2,2) ;
        zvec = [0 B*chord1AltHat(3)] + newCentroids(2,3) ;
        plot3(xvec, yvec, zvec, 'ks--','linewidth',2,'markersize',8, 'markerfacecolor',colmap(2,:)) ;
        
        % plot left wing span
        xvec = [0 A*span2Hat(1)] + newCentroids(3,1) ;
        yvec = [0 A*span2Hat(2)] + newCentroids(3,2) ;
        zvec = [0 A*span2Hat(3)] + newCentroids(3,3) ;
        plot3(xvec, yvec, zvec, 'ko-','linewidth',3,'markersize',8, 'markerfacecolor',colmap(3,:)) ;
        
        % plot left wing chord
        xvec = [0 B*chord2Hat(1)] + newCentroids(3,1) ;
        yvec = [0 B*chord2Hat(2)] + newCentroids(3,2) ;
        zvec = [0 B*chord2Hat(3)] + newCentroids(3,3) ;
        plot3(xvec, yvec, zvec, 'ko-','linewidth',3,'markersize',8, 'markerfacecolor',colmap(3,:)) ;
        
        % plot left wing alternative chord
        xvec = [0 B*chord2AltHat(1)] + newCentroids(3,1) ;
        yvec = [0 B*chord2AltHat(2)] + newCentroids(3,2) ;
        zvec = [0 B*chord2AltHat(3)] + newCentroids(3,3) ;
        % plot left wing y-axis
        xvec = [0 B*y_wing1(1)] + newCentroids(2,1) ;
        yvec = [0 B*y_wing1(2)] + newCentroids(2,2) ;
        zvec = [0 B*y_wing1(3)] + newCentroids(2,3) ;
        plot3(xvec, yvec, zvec, 'b-','linewidth',2,'markersize',8, 'markerfacecolor',colmap(2,:))
        % plot left  wing z-axis
        xvec = [0 B*z_wing1(1)] + newCentroids(2,1) ;
        yvec = [0 B*z_wing1(2)] + newCentroids(2,2) ;
        zvec = [0 B*z_wing1(3)] + newCentroids(2,3) ;
        plot3(xvec, yvec, zvec, 'r-','linewidth',2,'markersize',8, 'markerfacecolor',colmap(2,:));
        % plot left wing normal vector
        xvec = [0 B*wingNormVecWael1(1)] + newCentroids(2,1) ;
        yvec = [0 B*wingNormVecWael1(2)] + newCentroids(2,2) ;
        zvec = [0 B*wingNormVecWael1(3)] + newCentroids(2,3) ;
        plot3(xvec, yvec, zvec, 'y-','linewidth',2.5,'markersize',8, 'markerfacecolor',colmap(2,:));
        plot3(xvec, yvec, zvec, 'ks--','linewidth',2,'markersize',8, 'markerfacecolor',colmap(2,:)) ;
         % plot right wing y-axis
        xvec = [0 B*y_wing2(1)] + newCentroids(3,1) ;
        yvec = [0 B*y_wing2(2)] + newCentroids(3,2) ;
        zvec = [0 B*y_wing2(3)] + newCentroids(3,3) ;
        plot3(xvec, yvec, zvec, 'b-','linewidth',2,'markersize',8, 'markerfacecolor',colmap(2,:))
        % plot right wing z-axis
        xvec = [0 B*z_wing2(1)] + newCentroids(3,1) ;
        yvec = [0 B*z_wing2(2)] + newCentroids(3,2) ;
        zvec = [0 B*z_wing2(3)] + newCentroids(3,3) ;
        plot3(xvec, yvec, zvec, 'r-','linewidth',2,'markersize',8, 'markerfacecolor',colmap(2,:));
        % plot right wing normal vector
        xvec = [0 B*wingNormVecWael2(1)] + newCentroids(3,1) ;
        yvec = [0 B*wingNormVecWael2(2)] + newCentroids(3,2) ;
        zvec = [0 B*wingNormVecWael2(3)] + newCentroids(3,3) ;
        plot3(xvec, yvec, zvec, 'y-','linewidth',2.5,'markersize',8, 'markerfacecolor',colmap(2,:));
        
        
        plot3( leftWingVoxels(vox1IndLeft,1), leftWingVoxels(vox1IndLeft,2), ...
            leftWingVoxels(vox1IndLeft,3), 'k^','markerfacecolor','y') ;
        plot3( leftWingVoxels(vox2IndLeft,1), leftWingVoxels(vox2IndLeft,2), ...
            leftWingVoxels(vox2IndLeft,3), 'k^','markerfacecolor','y') ;
        
        plot3( rightWingVoxels(vox1IndRight,1), rightWingVoxels(vox1IndRight,2), ...
            rightWingVoxels(vox1IndRight,3), 'k^','markerfacecolor','y') ;
        plot3( rightWingVoxels(vox2IndRight,1), rightWingVoxels(vox2IndRight,2), ...
            rightWingVoxels(vox2IndRight,3), 'k^','markerfacecolor','y') ;
        
        % plot farthest point for both wings
        
        %[p, ~, list] = farthestPoint(rightWingVoxels, cBody, LL) ;
        farPointRight = farPoint1 ; %WS: has to be added because we do the clustering before the analysis
        farPointLeft  = farPoint2 ;
        plot3(farPointRight(1), farPointRight(2), farPointRight(3),'ks','markersize',14,'markerfacecolor','y','linewidth',2) ;
        %plot3(rightWingVoxels(list,1), rightWingVoxels(list,2), rightWingVoxels(list,3),'o','markersize',2,'color',[0.75 0 0]) ;
        plot3(wing1LargestCC(:,1), wing1LargestCC(:,2), wing1LargestCC(:,3),'o','markersize',2,'color',[0.6 0 0]) ;
        
        %[p, ~, list] = farthestPoint(leftWingVoxels, cBody, LL) ;
        plot3(farPointLeft(1), farPointLeft(2), farPointLeft(3),'ks','markersize',14,'markerfacecolor','y','linewidth',2) ;
        %plot3(leftWingVoxels(list,1), leftWingVoxels(list,2), leftWingVoxels(list,3),'o','markersize',2,'color',[0 0 0.75]) ;
        plot3(wing2LargestCC(:,1), wing2LargestCC(:,2), wing2LargestCC(:,3),'o','markersize',2,'color',[0 0 0.6]) ;
        
        xlabel('x') ; ylabel('y'); zlabel('z') ;
        
        hold off ;
        axis equal ;
        grid on ;
        box on ;
        
        %view(92,8) ;
        view([159 10]);
        
        axis(ax) ;
        
        refresh ;
        
    end
    %% save data
    hull_analysis_save = fullfile(init.folders.root, 'hull_analysis', ['frame_' num2str(ii) '.mat']);
    save(hull_analysis_save,'chord2Hat','chord2AltHat','chord1Hat','chord1AltHat',...
        'span1Hat','span2Hat','rotation_angle_wael2','rotation_angle_wael1','eta1','eta2',...
        'angle_rotation1_final_unit','angle_rotation2_final_unit','farPointLeft','farPointRight')
end


end
%% ------------------------
% FUNCTIONS
function [p, dst, list] = farthestPointWaelV1(coords, p0, L)
% finds the point in coords which is farthest from p0
% returns the point coordinate p and its distance from p0 in dst. Furtherst
% point is detected by finding an n number of points that are the furthest
% from the body and taking the average. This should help account for
% potential noise in reconstruction
% if there are several points with the same distance, return only one.
% also finds the indices of voxels in coords that whose distance from p is
% smaller or equal to L

% find p
Nvox = size(coords,1) ;
mat1 = double(coords) - repmat(p0, Nvox,1) ;
dst2vec  = sum (mat1 .* mat1, 2) ;
[dst, ind] = maxk(dst2vec,8) ;
p = coords(ind,:) ;
p=mean(p);
% modify p so it belongs to the wing
mat2 = double(coords) - repmat(p, Nvox,1);
dst2vec_P  = sum (mat2 .* mat2, 2);
[dst_P, ind_P] = min(dst2vec_P);
p = coords(ind_P,:) ;
% find distance from new p to center of body
mat3 = double(coords) - repmat(p, Nvox,1) ; %matrix that contains the point coordinates
dst2vec  = sum (mat3 .* mat3, 2) ;
[dst, ind] = max(dst2vec) ;
% find list
mat1 = double(coords) - repmat(p, Nvox, 1) ;
dst2vec  = sum (mat1 .* mat1, 2) ;
list = (dst2vec <= L^2) ;

end

function [c xyunique ] = calcCentroidFromTopView(coords, farPoint, L)
% find the centroid of the voxel set "coords" in two steps:
% 0. include only the points whose distnace to the farPoint is more than L
%    L is typically wingLength/3
% 1. consider only the unique (x,y) pairs, ignoring z coordinate. Find their
%    center of mass on the xy plane (xc, yc)
% 2. consider the voxels (x,y,z) for which (x-xc)^2+(y-yc)^2<=delta2.
%    i.e. take only the distance in the xy plane, ignoring their z.
%    Find the center of mass along z of this set only.
%
% retrun the center of wing as well as the unique coordinates in the xy
% plane. this will be used in the function xx to calculate the projection
% of the chord in the xy plane as seen from the top view.

c = zeros(1,3) ;
delta2 = 4 ;

% step 0
% ------
Nvox = size(coords,1) ;
dst = myNorm ( coords - repmat(farPoint,Nvox,1) ) ;
ind0 = dst > L ;

% step 1
% ------
xyunique = unique(coords(ind0,1:2),'rows') ;
cxy      = mean(xyunique) ;
c(1:2)   = cxy ;

% step 2
% ------
coords2 = coords(ind0,:) ;

ind = (coords2(:,1)-c(1)).^2 + (coords2(:,2)-c(2)).^2 <= delta2 ;

c(3) = mean(coords2(ind,3)) ;

end

function tip = findWingTip(coords, spanVec, wingCM)

% find distances of coords from wingCM
Nvox = size(coords,1) ;

fraction = 0.05 ;

% 1st screening
% -------------
% consider only voxels in a cone of 30 degrees around the line starting
% from wingCM along the direction of spanVec
mat1 = coords - repmat(wingCM, Nvox, 1) ;
mat1 = mat1 ./ repmat( myNorm(mat1), 1,3) ;
mat2 = repmat(spanVec, Nvox, 1) ;
dotprod = dot(mat1, mat2, 2) ;

ind1 = dotprod > 0.5 ; % cos(30 deg) = 0.5

% 2nd screening
% -------------
% take the farthest "fraction" of the voxels and calculate their mean position

coords = coords(ind1,:) ;
Nvox   = size(coords,1) ;

dst  = myNorm ( coords - repmat(wingCM,Nvox,1) ) ;

[~, sortedInd] = sort(dst,'descend') ;
Nvox1 = ceil(Nvox*fraction) ;
selectedInd    = sortedInd(1:Nvox1) ;
if (isempty(selectedInd))
    tip = [NaN NaN NaN] ;
else
    if (numel(selectedInd)==1)
        tip = coords(selectedInd,:) ;
    else
        tip = mean(coords(selectedInd,:)) ;
    end
end

end

function [p1 p2 d12] = chordTopViewProjection(spanVec, wingCM, coords)

%[p1 p2 d12] = chordTopViewProjection(spanHatRight, cRight, rightWingLargestCC) ;


% p1 and p2 are the extreme points along the "chord" and d12 is the
% distance between them

plotFlag = false ; % for debugging purposes

delta = 2 ; % width of strip

rcm = double(wingCM(1:2)) ; % work in 2D

sxy = [spanVec(1) spanVec(2) ];   % 2D span vector
sxy = sxy / norm(sxy) ;           % normalize to be a unit vector in xy plane

cxy = [ -spanVec(2) spanVec(1)] ; % rotate spanVec by 90 degrees ("chord" vector)
cxy = cxy / norm(cxy) ;           % normalize to be a unit vector in xy plane

xyunique = unique(coords(:,1:2),'rows') ; % top projection

N = size(xyunique,1) ;

% relative vector w.r.t. wing cm
xyunique = double(xyunique) ;
xyunique(:,1) = xyunique(:,1) - rcm(1) ; % use this shape since it is insensitive
xyunique(:,2) = xyunique(:,2) - rcm(2) ; % to the dimension of rcm (2,1) or (1,2)

% calc the component of each pixel along the sxy direction (sxy is a unit
% vector, so that's fine)
dotprod = dot ( xyunique, repmat(sxy,N,1) ,2) ;

ind = abs(dotprod)<=delta ; % logical indices of pixels inside the delta-strip

clear dotprod ;
% now do the same trick for the pixels in the strip - calculate their
% componenet along the cxy vector

N2 = sum(ind) ; % how many pixels in the strip

if (N2==0)
    disp('no pixels in the delta-strip. assign zero to "chord" xy projection') ;
    p1 = rcm ;
    p2 = rcm ;
    d12 = 0 ;
    disp('problem 9?') ;
    % keyboard ;
    return ;
end

dotprod = dot(xyunique(ind,:), repmat(cxy, N2,1), 2) ;

a1 = max(dotprod) ;
a2 = min(dotprod) ;

if (a1<0 || a2>0)
    disp('problem 10?') ;
    %keyboard ;
end

p1 = rcm + a1 * cxy ;
p2 = rcm + a2 * cxy ;
d12 = norm(p1 - p2) ;

if (plotFlag)
    Lspan = 20 ;
    Lchord = 10 ;
    
    figure(71)
    plot(xyunique(:,1), xyunique(:,2),'b.') ;
    hold on ;
    
    plot(xyunique(ind,1), xyunique(ind,2),'ko','markersize',10) ;
    plot([0 Lspan*sxy(1)], [0 Lspan*sxy(2)],'k-','linewidth',4) ;
    plot([0 Lchord*cxy(1)], [0 Lchord*cxy(2)],'b-','linewidth',4) ;
    
    plot(a1*cxy(1), a1*cxy(2),'k^','markerfacecolor','y','markersize',10) ;
    plot(a2*cxy(1), a2*cxy(2),'k^','markerfacecolor','r','markersize',10) ;
    
    
    axis equal ; box on ; grid on ;
    title(['d_1_2 = ' num2str(d12)]) ;
    hold off ;
end

end