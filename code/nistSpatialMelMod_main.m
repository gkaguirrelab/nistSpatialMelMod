%
% nistSpatialMelMod_main
%
%	Create a 2D matrix of spectra that are designed to create spatially
%	varying contrast of the melanopsin photopigment while silencing the
%	cones.
%
%	The spectra are defined using the routine:
%   	nistSpatialMelMod_makeModulationPrimaries 
%   which is derived from the receptorIsolateDemo of the SilentSubstitution
%   toolbox.
%

%% Housekeeping
clear all
close all


%% Hardcoded parameters
% Define the spatial properties of the display and observer
displayPixelResolution = [200 200];
dispaySizeMeters = displayPixelResolution / 500;
observerDistanceMeters = 1.0;

% Define the spatial properties of the stimulis
radiusInnerEdgeAnnulusDeg = 2.5;
radiusOuterEdgeAnnulusDeg = 10;
widthHalfCosineSmoothDeg = 2;
gratingSpatialFreqHz = .2;
gratingSpatialPhaseDeg = 0;
gratingOrientationDeg = 45;

% Define the temporal properties of the stimulis
gratingContrastModulateTemporalFreqHz = 0.1;
gratingFramesPerSec = 10;

% Define the spectral properties of the stimulus
maxContrast = 0.3;
targetedReceptor = 4; % this is melanopsin

% Define the path to the calibration file
codeBaseDir = tbLocateProject('nistSpatialMelMod','verbose',false);
calibrationFileName = fullfile(codeBaseDir, 'demoOneLightCalFile', 'OneLightDemoCal.mat');


%% Setup
% Derive some values from these inputs
displaySizeDeg = rad2deg(atan((dispaySizeMeters/2)/ observerDistanceMeters)*2);
pixelsPerDeg = displayPixelResolution./displaySizeDeg;
if min(abs(diff(pixelsPerDeg))./pixelsPerDeg) > 0.05
    error('The display pixels depart from square by more than 5%');
else
    pixelsPerDeg=mean(pixelsPerDeg);
end


%% Obtain the primaries and associated variables
% The default settings of the called function provide melanopsin-directed
% stimulation while silencing the cones, but ignoring the rods and the
% penumbral cones.
[ backgroundPrimaries, maxPositivePrimaries, B_primary, ambientSpd, T_receptors ] = ...
    nistSpatialMelMod_makeModulationPrimaries( calibrationFileName, 'desiredContrast',maxContrast);

% Obtain the SPDs
maxPositiveModulationSPD = B_primary*(maxPositivePrimaries - backgroundPrimaries);
backgroundSPD = B_primary*backgroundPrimaries;

    
%% Define the temporal domain of the stimulus
% Setting this up for a contrast modulation. This could be readily adapted
% for spatially drifting gratings, rotating gratings, etc.
temporalSupport = 0 : 1000 / gratingFramesPerSec : (1/gratingContrastModulateTemporalFreqHz)*1000 - 1000 / gratingFramesPerSec;
temporalSupportRads = (temporalSupport / ((1/gratingContrastModulateTemporalFreqHz)*1000)) * 2 * pi;
temporalContrastModulation = sin(temporalSupportRads);

%% Create the spatial stimulus matrix
% Obtain the distance from fixation (assumed to be in the center of the
% display) for each point on the display
[Xm, Ym] = meshgrid(1:displayPixelResolution(1), 1:displayPixelResolution(2));
distanceFromFixationDeg = sqrt((Xm - displayPixelResolution(1)/2).^2 + (Ym - displayPixelResolution(2)/2).^2) / pixelsPerDeg;

% Define an anonymous function that calls out to a routine that implements
% the half-cosine weighting of the stimulus inner and outer edges, and the
% zeroing of contrast beyond the annular stimulus broders
returnMaskWeight = @(x) nistSpatialMelMod_createAnnulus(x,radiusInnerEdgeAnnulusDeg, radiusOuterEdgeAnnulusDeg, widthHalfCosineSmoothDeg);

% Create the spatial weighting mask
spatialWeightingMask = arrayfun(returnMaskWeight, distanceFromFixationDeg);

% Loop over the temporal support and obtain the SPD at each pixel
for tt = 1:length(temporalSupport)
    
    % Create 2D sinusoidal grating. This remains in the loop fo allow for 
    % future modifications in which we vary some aspect of the grating over
    % time (e.g., orientation, phase, spatial frequency)
    sinusoidalGrating = maxContrast * temporalContrastModulation(tt) * ...
        nistSpatialMelMod_createGrating( displayPixelResolution, ...
        pixelsPerDeg, ...
        gratingSpatialFreqHz, ...
        gratingOrientationDeg, ...
        gratingSpatialPhaseDeg);
    
    % Compute the called-for contrast ipon the targeted photoreceptor at
    % each point of the display
    thisFrameContrastRelativeToMax = ...
        maxContrast * temporalContrastModulation(tt) * sinusoidalGrating .*spatialWeightingMask';
        
    % reshape the 2D matrix of contrasts into a vector to support matrix
    % multiplication below
    thisFrameContrastRelativeToMaxVec = reshape(thisFrameContrastRelativeToMax,prod(displayPixelResolution),1);
    
    % Obtain the SPD at each pixel, scaled by the called-for contrast
    % relative to max
    thisFrameSPDVec = thisFrameContrastRelativeToMaxVec*maxPositiveModulationSPD';
    
    % Add the background SPD
    thisFrameSPDVec = bsxfun(@plus, thisFrameSPDVec, backgroundSPD');
    
    % As a check, calculate the contrast on the targeted photoreceptor at each pixel
    thisFrameReceptorsVec = thisFrameSPDVec*T_receptors';
    backgroundReceptors = T_receptors*(B_primary*backgroundPrimaries + ambientSpd);
    thisFrameTargetedReceptors = reshape(squeeze(thisFrameReceptorsVec(:,targetedReceptor)),displayPixelResolution(1), displayPixelResolution(2));
    thisFrameTargetedWeberContrast = (thisFrameTargetedReceptors - backgroundReceptors(targetedReceptor)) ./ backgroundReceptors(targetedReceptor);
    
    % Display the resulting contrast image
    if tt==1
        figure
        imshow( thisFrameTargetedWeberContrast, [-maxContrast maxContrast] );
        colormap gray(256);
        axis off; axis image;
    else
        imshow( thisFrameTargetedWeberContrast, [-maxContrast maxContrast] );
    end
    
end
