function [ grating ] = nistSpatialMelMod_createGrating( displayPixelResolution, pixelsPerDeg, gratingSpatialFreqHz, gratingOrientationDeg, gratingSpatialPhaseDeg)
% function [ grating ] = nistSpatialMelMod_createGrating( displayPixelResolution, pixelsPerDeg, gratingSpatialFreqHz, gratingOrientationDeg, gratingSpatialPhaseDeg)
%
% Creates a 2D sinusoidal grating
%

gratingSpatialFrequencyPixels = (gratingSpatialFreqHz / pixelsPerDeg) * max(displayPixelResolution);

X0 = ((1:max(displayPixelResolution)) / max(displayPixelResolution)) - .5;
[Xm, Ym] = meshgrid(X0, X0);
Xt = Xm * cos(deg2rad(gratingOrientationDeg));                % compute proportion of Xm for given orientation
Yt = Ym * sin(deg2rad(gratingOrientationDeg));                % compute proportion of Ym for given orientation
XYt = Xt + Yt;                      % sum X and Y components
XYf = XYt * gratingSpatialFrequencyPixels * 2*pi;                % convert to radians and scale by frequency
grating = sin( XYf + deg2rad(gratingSpatialPhaseDeg));                   % make 2D sinewave

% Crop the grating back to the original dimensions

grating = grating(1:displayPixelResolution(1),1:displayPixelResolution(2));

end % function
