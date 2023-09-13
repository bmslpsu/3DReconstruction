
%% Improved code for laser
%% Modifications:
%Saves all body angles to one vector in a sep folder
%Uses parallel processing for image analysis
%Uses par processing for reconstruction
%Detects flapping frequency from 2D image
%Estimates area ratio change in damaged-wing flies

%% Code starts here
clc
clear
close all

%% Set mask data directories
% root = 'S:\Public\Wael\Chapter 1\DamagedWingExp\AnalyzedData\';
root = 'S:\Public\Lingsheng\Laser_RigidExp';
fly = 1000; %change fly number to match the number of fly you are analyzing 

%% initialized directory and voxel_size
init = load([root '\fly_' num2str(fly) '\init\init.mat']); %load intialization
load([root '\fly_' num2str(fly) '\init\init.mat']);
voxel_size = 5e-5; % the resolution of the reconstruction (voxel size)

%% Split left and right wing
[oneWingFlag,error_Flag] = DifferentiateLeftAndRightWingV3(init,21:22);
%% UI error detection
clc
ErrorDetectionByUser(init,21:22); 
%% find stroke reversal
flag_stroke_reversal=find_stroke_reversal(init,21,22,false);
%% UI correction for those pesky large wings (important at stroke reversals)
clc; close all
errorlabel = 6744; 
CorrectForLargeWingsV6(init,21:22,zeros(length(flag_stroke_reversal),1),voxel_size); 

%% Be familiar with the code until this section

%% Extract info from code (PART 2)
%% verify reconstruction and correction (just for a visual check, does nothing to change data)
VerifyReconstruction(init,21,21)
%% align the reconstruction with a global reference frame
Rotate3DVoxelsV3(init,21,21,true); 
%% convert the hull coordinates to voxel coordinates 
clc; close all
ConvertToVoxelSpacev2(init,21,21,true,voxel_size);
%% Analyze the hull
clc; close all
analyze_hull_load_dataV2(init,21,21,true,true);
%% Correct chordwise and span vectors (not working)
% clc; close all;
% CorrectChordAutomatic(init,sI,eI)
%% estimates Euler angles
clc
[phi1, theta1, eta1, phi2, theta2, eta2]=EstimateEulerAnglesV2(init,1301,1350);
%% peform angle corrections and filtering
cutFreq=400;
SamplingFreq=8000;
[EulerAnglesR_Filt, EulerAnglesL_Filt]=FilterEulerAngles([phi1', eta1', theta1'], [phi2', eta2', theta2'],cutFreq,SamplingFreq,init);
