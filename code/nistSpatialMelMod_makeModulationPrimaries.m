function [ backgroundPrimary, modulationPrimary, B_primary, ambientSpd, T_receptors, wavelengthSupport ] = nistSpatialMelMod_makeModulationPrimaries( calibrationFileName, varargin )
% function [ backgroundPrimary, modulationPrimary, B_primary, ambientSpd, T_receptors ] = nistSpatialMelMod_makeModulationPrimaries( calibrationFileName, varargin )
%
%  Returns primaries and associated variables that may be used to compute
%  spectral modulations for the purpose of silent substitution.
%
%  Much of the code is derived from ReceptorIsolateDemo.m within the
%  SilentSubstitutionToolbox
%
% INPUTS:
%   calibrationFileName: full path to a .mat file that contains a
%       calibration of the display device
%
% OPTIONAL - I/O and flow control
%   verbose - true or false
%
% OPTIONAL - device parameters:
%   primaryHeadRoom -  This parameter enforces a constraint that we don't
%       go right to the edge of the gamut.  The head room parameter is
%       defined in the [0-1] device primary space.  Using a little head
%       room keeps us a bit away from the hard edge of the device.
%   maxPowerDiff - This is a smoothness constraint that determins the
%       maximum absolute value of the change in spectral power between two
%       adjacent wavelength samples. The appropriate value depends on the
%       overall power of the viewed light as well as on the wavelength
%       sampling step.
%
% OPTIONAL - photoreceptor classes and modulation direction
%   photoreceptorClasses - A cell array of photoreceptor classes. This is
%       passed to the function GetHumanPhotoreceptorSS, which should know
%       what to do with it.
%   whichReceptorsToTarget - Vector of indices to photoreceptorClasses that
%       identify photoreceptors to stimulate
%	whichReceptorsToIgnore - Vector of indices to photoreceptorClasses that
%       identify photoreceptors to ignore
%       (All remaining photoreceptors are silenced)
%   desiredContrast - vector of contrast to be sought for each of the 
%       targeted receptor classes. Set to empty to maximize contrast.
%
% OPTIONAL - observer and physiologic parameters:
%   observerAgeInYears - impacts spectral sensitivity through lens density
%   fieldSizeDegrees - impacts spectral sensitivity as a consequence of
%       macular pigment
%   pupilDiameterMm - impacts retinal irradiance and thus bleaching
%   vesselOxyFraction - Assumed fraction of oxygenated hemoglobin within
%       retinal blood vessels, used for calculation of spectral sensitivity
%       of penumbral cones
%   vesselOverallThicknessUm - Retinal blood vessel thickness assumed for
%       calculation of spectral sensitivity of penumbral cones
%   correctBleaching - Correct for cone photopigment bleaching?

%% Parse vargin for options passed here
p = inputParser;

% Required
p.addRequired('calibrationFileName',@ischar);

% Optional - I/O and flow control
p.addParameter('verbose',false,@islogical);

% Optional - device parameters
p.addParameter('primaryHeadRoom',0.02,@isnumeric);
p.addParameter('maxPowerDiff',10^-1.5,@isnumeric);

% Optional - photoreceptor classes and modulation direction
p.addParameter('photoreceptorClasses', ...
    {'LConeTabulatedAbsorbance', ...
    'MConeTabulatedAbsorbance', ...
    'SConeTabulatedAbsorbance', ...
    'Melanopsin', ...
    'Rods', ...
    'LConeTabulatedAbsorbancePenumbral', ...
    'MConeTabulatedAbsorbancePenumbral', ...
    'SConeTabulatedAbsorbancePenumbral'},@iscell);
p.addParameter('whichReceptorsToTarget',[4],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('whichReceptorsToIgnore',[5 6 7 8],@(x)(isempty(x) | isnumeric(x)));
p.addParameter('desiredContrast',[0.4],@(x)(isempty(x) | isnumeric(x)));

% Optional - observer and physiologic parameters
p.addParameter('whichModel','HumanPhotopigments',@ischar);
p.addParameter('observerAgeInYears',32,@isnumeric);
p.addParameter('fieldSizeDegrees',27.5,@isnumeric);
p.addParameter('pupilDiameterMm',4.7,@isnumeric);
p.addParameter('vesselOxyFraction',0.85,@isnumeric);
p.addParameter('vesselOverallThicknessUm',5,@isnumeric);
p.addParameter('correctBleaching',true,@islogical);


%% Parse and check the parameters
p.parse(calibrationFileName, varargin{:});


%% Load the calibration file and extract the descrption of spectral primaries
% LoadCalFile needs the path to the cal file spoon-fed. We do so here.
tmpSplitPath=split(calibrationFileName, filesep);
cal = LoadCalFile(char(tmpSplitPath(end)),[],char(join(tmpSplitPath(1:end-1),filesep)));
S = cal.describe.S;
B_primary = cal.computed.pr650M;
ambientSpd = cal.computed.pr650MeanDark;

% S contains a specification of wavelength support in the form of:
%   [start delta n]. The psychtoolbox function SToWLs expands this to a
%   support vector of wavelengths.
wavelengthSupport = SToWls(S);

% Half on in OneLight primary space
backgroundPrimary = 0.5*ones(size(B_primary,2),1);

% Do not pin any of the primaries
whichPrimariesToPin = [];



% The routines that do these computations are in the
% ContrastSplatter directory of the SilentSubstitutionToolbox. They
% provide pre-defined receptor types and compute spectral
% sensitivities using the routines provided in the Psychtoolbox.
% The routines here, however, also allow computation of fraction
% cone bleached, which may be used to adjust pigment peak optical
% density.  They can also compute photopigment variants corrected
% for filtering by blood vessels.

% Note that we don't typically vary or pass the blood vessel
% parameters but rather simply accept the defaults used by
% GetHumanPhotoreceptorSS.  It's mainly for fun that we show
% how to do this here.

% Define photoreceptor classes that we'll consider.
% ReceptorIsolate has a few more built-ins than these.

% Correct for pigment bleaching if desired.  This is done
% separately for open-field and penumbral cones.  The bleaching
% correction routine only knows about L, M, and S cones.
if p.Results.correctBleaching
    % The hard work is done by routine
    % GetConeFractionBleachedFromSpectrum, which is in the
    % ContrastSplatter directory. We set the verbose flag here.
    [fractionBleachedFromIsom, fractionBleachedFromIsomHemo] = ...
        GetConeFractionBleachedFromSpectrum(S, B_primary*backgroundPrimary + ambientSpd, p.Results.fieldSizeDegrees, p.Results.observerAgeInYears, p.Results.pupilDiameterMm, [], p.Results.verbose);
    
    % Assign the fraction bleached for each photoreceptor class.
    for ii = 1:length(p.Results.photoreceptorClasses)
        switch p.Results.photoreceptorClasses{ii}
            case 'LCone'
                fractionBleached(ii) = fractionBleachedFromIsom(1);
            case 'MCone'
                fractionBleached(ii) = fractionBleachedFromIsom(2);
            case 'SCone'
                fractionBleached(ii) = fractionBleachedFromIsom(3);
            case 'LConeHemo'
                fractionBleached(ii) = fractionBleachedFromIsomHemo(1);
            case 'MConeHemo'
                fractionBleached(ii) = fractionBleachedFromIsomHemo(2);
            case 'SConeHemo'
                fractionBleached(ii) = fractionBleachedFromIsomHemo(3);
            otherwise
                fractionBleached(ii) = 0;
        end
    end
else
    fractionBleached = [];
end

% Make sensitivities.  The wrapper routine is
% GetHumanPhotoreceptorSS, which is in the ContrastSplatter
% directory.  Each row of the matrix T_receptors provides the
% spectral sensitivity of the photoreceptor class in the
% corresponding entry of the cell array photoreceptorClasses.
%
% The last two arguments are the oxygenation fraction and the
% vessel thickness. We set them to be empty here, prompting the
% user to enter these values later.
T_receptors = GetHumanPhotoreceptorSS(S, p.Results.photoreceptorClasses, p.Results.fieldSizeDegrees, p.Results.observerAgeInYears, p.Results.pupilDiameterMm, [], fractionBleached, [], []);


% User chooses whether to maximize contrast in targeted receptor classes or
% or get it as close to a specified value as possible.
%
% If we target, here we specify the same contrast for all targeted classes.
% This is not necessary, they can differ.  It just makes the demo code a
% bit simpler to yoke them since we only have to prompt for one number.


% Nice message for user
if p.Results.verbose
    fprintf('\nGenerating stimuli which isolate receptor classes');
    for i = 1:length(p.Results.whichReceptorsToTarget)
        fprintf('\n  - %s', p.Results.photoreceptorClasses{p.Results.whichReceptorsToTarget(i)});
    end
    fprintf('\nGenerating stimuli which ignore receptor classes');
    if (~length(p.Results.whichReceptorsToIgnore) == 0)
        for i = 1:length(p.Results.whichReceptorsToIgnore)
            fprintf('\n  - %s', p.Results.photoreceptorClasses{p.Results.whichReceptorsToIgnore(i)});
        end
    else
        fprintf('\n  - None');
    end
    fprintf('\nThe remaining classes will be silenced\n');
end

%% Call the optimization routine.
%
% Careful examaination of the arguments will reveal that the initialGuess
% for the primaries is set to the background value for the primaries, so
% that the constraints are all met at the start of the search.  The
% optimization routine is much happier when it is done this way -- bad
% things happen if you start with a guess that violates constraints.
modulationPrimary = ReceptorIsolate(T_receptors,p.Results.whichReceptorsToTarget, p.Results.whichReceptorsToIgnore, [], ...
    B_primary, backgroundPrimary, backgroundPrimary, whichPrimariesToPin,...
    p.Results.primaryHeadRoom, p.Results.maxPowerDiff, p.Results.desiredContrast, ambientSpd);

%% Compute the contrasts that we got.
backgroundReceptors = T_receptors*(B_primary*backgroundPrimary + ambientSpd);
modulationReceptors = T_receptors*B_primary*(modulationPrimary - backgroundPrimary);
contrastReceptors = modulationReceptors ./ backgroundReceptors;

if p.Results.verbose
    for j = 1:size(T_receptors,1)
        fprintf('\t%s: contrast = %0.4f\n',p.Results.photoreceptorClasses{j},contrastReceptors(j));
    end
end


end

