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
clearvars
close all


% Display, I/O, and flow control parameters
outputMovieFileName = '~/Desktop/nistSpatialMelMod_Movie.mp4';
outputSPDMatrixFileName = '~/Desktop/nistSpatialMelMod_SPDMatrix.mat';
outputSPDWavelengthSupportFileName = '~/Desktop/nistSpatialMelMod_WavelengthSupport.mat';


%% Hardcoded parameters
% Define the spatial properties of the display and observer
displayPixelResolution = [192 256];
dispaySizeMeters = [0.5 0.666];
observerDistanceMeters = 1.0;

% Define the spatial properties of the stimulus
radiusInnerEdgeAnnulusDeg = 2.5;
radiusOuterEdgeAnnulusDeg = 12.5;
widthHalfCosineSmoothDeg = 2;
gratingSpatialFreqHz = .2;
gratingSpatialPhaseDeg = 0;
gratingOrientationDeg = 45;

% Define the temporal properties of the stimulus
gratingContrastModulateTemporalFreqHz = 0.1;
gratingFramesPerSec = 10;

% Define the spectral properties of the stimulus
maxContrast = 0.3;
targetedReceptor = 4; % this is melanopsin

% Define the path to the calibration file
codeBaseDir = tbLocateProject('nistSpatialMelMod','verbose',false);
calibrationFileName = fullfile(codeBaseDir, 'demoOneLightCalFile', 'OneLightDemoCal.mat');

% Pick a pixel to plot the spectrum
pixelCoordsToPlot = [round(displayPixelResolution(2)*0.5) round(displayPixelResolution(1)*0.25)];

%% Setup
% Derive some values from these inputs
displaySizeDeg = rad2deg(atan((dispaySizeMeters/2)/ observerDistanceMeters)*2);
pixelsPerDeg = displayPixelResolution./displaySizeDeg;
if min(abs(diff(pixelsPerDeg))./pixelsPerDeg) > 0.05
    error('The display pixels depart from square by more than 5%');
else
    pixelsPerDeg=mean(pixelsPerDeg);
end

% Get the linear index for the pixel we will later plot
pixelIdxToPlot=sub2ind(displayPixelResolution,pixelCoordsToPlot(2),pixelCoordsToPlot(1));

% Prepare a figure
figHandle=figure;
figPanelA=subplot(1,3,1);
figPanelB=subplot(1,3,2);
figPanelC=subplot(1,3,3);

% Open an avi file for saving if a filename has been given
if ~isempty(outputMovieFileName)
    videoOutHandle = VideoWriter(outputMovieFileName,'MPEG-4');
    videoOutHandle.FrameRate = gratingFramesPerSec;
    videoOutHandle.Quality = 100;
    open(videoOutHandle);
end
    

%% Obtain the primaries and associated variables
% The default settings of the called function provide melanopsin-directed
% stimulation while silencing the cones, but ignoring the rods and the
% penumbral cones.
[ backgroundPrimaries, maxPositivePrimaries, B_primary, ambientSpd, T_receptors, wavelengthSupport ] = ...
    nistSpatialMelMod_makeModulationPrimaries( calibrationFileName, 'desiredContrast',maxContrast);

% Obtain the SPDs
maxPositiveModulationSPD = B_primary*(maxPositivePrimaries - backgroundPrimaries);
backgroundSPD = B_primary*backgroundPrimaries;

    
%% Define the temporal domain of the stimulus
% Setting this up for a contrast modulation. This could be readily adapted
% for spatially drifting gratings, rotating gratings, etc.
temporalSupport = 0 : 1000 / gratingFramesPerSec : (1/gratingContrastModulateTemporalFreqHz)*1000 - 1000 / gratingFramesPerSec;
temporalSupportRads = (temporalSupport / ((1/gratingContrastModulateTemporalFreqHz)*1000)) * 2 * pi;
temporalContrastModulation = cos(temporalSupportRads);

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
    
    % Create 2D sinusoidal grating. This remains in the loop to allow for 
    % future modifications in which we vary some aspect of the grating over
    % time (e.g., orientation, phase, spatial frequency)
    sinusoidalGrating = nistSpatialMelMod_createGrating( ...
        displayPixelResolution, ...
        pixelsPerDeg, ...
        gratingSpatialFreqHz, ...
        gratingOrientationDeg, ...
        gratingSpatialPhaseDeg);
    
    % Compute the called-for contrast upon the targeted photoreceptor at
    % each point of the display
    thisFrameContrastRelativeToMax = ...
        temporalContrastModulation(tt) * sinusoidalGrating .*spatialWeightingMask';
        
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
    
    % At a single pixel, calculate the contrast on each photoreceptor
    thisFramePlotPixelAllReceptorsWeberContrast = (thisFrameReceptorsVec(pixelIdxToPlot,:)-backgroundReceptors')./backgroundReceptors';
    
    % Make and update a figure, which we will save as a movie
    % Panel A -- the spatial contrast image for melanopsin
    subplot(figPanelA);
    imshow( thisFrameTargetedWeberContrast, [-maxContrast maxContrast] );
    colormap gray(256);
    axis off; axis image;
    hold on
    text(pixelCoordsToPlot(1),pixelCoordsToPlot(2),'X','HorizontalAlignment','center','VerticalAlignment','middle')
    title(['Relative Weber contrast [\pm' num2str(maxContrast) ']']);
    hold off
    % Panel B -- The spectrum at a point on the display
    subplot(figPanelB);
    plot(wavelengthSupport,backgroundSPD,'k','LineWidth',2);
    hold on
    plot(wavelengthSupport,thisFrameSPDVec(pixelIdxToPlot,:),'r','LineWidth',2);
    title('Spectrum (red) at X');
    xlim([380 780]);
    ylim([0 0.1]);
    xlabel('Wavelength');
    ylabel('Power');
    pbaspect([1 displayPixelResolution(1)/displayPixelResolution(2) 1]);
    hold off
    % Panel C -- a column plot of contrasts at a point on the display
    subplot(figPanelC);
    bar(thisFramePlotPixelAllReceptorsWeberContrast(1:4));
    hold on
    bar(1,thisFramePlotPixelAllReceptorsWeberContrast(1),'r');
    bar(2,thisFramePlotPixelAllReceptorsWeberContrast(2),'g');
    bar(3,thisFramePlotPixelAllReceptorsWeberContrast(3),'b');
    bar(4,thisFramePlotPixelAllReceptorsWeberContrast(4),'c');
    title('Contrast by receptor at X');
    ylim([-0.5 0.5]);
    xlim([0 5]);
    set(gca,'XTickLabel',{'L','M','S','Mel'})
    xlabel('Photoreceptor class');
    ylabel('Weber contrast');
    pbaspect([1 displayPixelResolution(1)/displayPixelResolution(2) 1]);
    hold off
    drawnow
    
    if ~isempty(outputMovieFileName)
        videoFrame=getframe(figHandle);
        writeVideo(videoOutHandle,videoFrame)
    end
    
end % loop over temporalSupport frames

% Close the video object (if there is one)
if ~isempty(outputMovieFileName)
    close(videoOutHandle);
end

% For the last frame in memory, reshape into a 3D array which has the
% dimensions of pixelsX, pixelsY, and wavelengths. Save out this and its
% corresponding wavelength support
exampleFrameSPD = reshape(thisFrameSPDVec, displayPixelResolution(1), displayPixelResolution(2), length(backgroundSPD));
save(outputSPDMatrixFileName,'exampleFrameSPD');
save(outputSPDWavelengthSupportFileName,'wavelengthSupport');
