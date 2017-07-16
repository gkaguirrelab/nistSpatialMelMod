
codeBaseDir = tbLocateProject('nistSpatialMelMod','verbose',false);
whichModel = 'HumanPhotopigments';
whichPrimaries = 'OneLight';

        % Get a OneLight calibration file, stored here for demo purposes.
        % Extract the descrption of spectral primaries, which is what we
        % need for this demo.  Gamma correction of primary values to
        % settings would need to be handled to get the device to actually
        % produce the spectrum, but here we are just finding the desired
        % modulation.
        calPath = fullfile(codeBaseDir, 'demoOneLightCalFile');
        cal = LoadCalFile('OneLightDemoCal.mat',[],calPath);
        S = cal.describe.S;
        B_primary = cal.computed.pr650M;
        ambientSpd = cal.computed.pr650MeanDark;
        
        % Half on in OneLight primary space
        backgroundPrimary = 0.5*ones(size(B_primary,2),1);
        
        % Don't pin any primaries.  Do enforce a constraint that we don't
        % go right to the edge of the gamut.  The head room parameter is
        % defined in the [0-1] device primary space.  Using a little head
        % room keeps us a bit away from the hard edge of the device.
        whichPrimariesToPin = [];
        primaryHeadRoom = 0.02;
        
        % Set smoothness constraint value.  This is a magic constant
        % defined by hand, that determins the maximum absolute value of the
        % change in spectral power between two adjacent wavelength samples.
        % Thus it's appropriate value depends on the overall power of the
        % viewed light as well as on the wavelength sampling step.
        maxPowerDiff = 10^-1.5;
        
        
        
        
        
                % The routines that do these computations are in the
        % ContrastSplatter directory of the SilentSubstitutionToolbox. They
        % provide pre-defined receptor types and compute spectral
        % sensitivities using the routines provided in the Psychtoolbox.
        % The routines here, however, also allow computation of fraction
        % cone bleached, which may be used to adjust pigment peak optical
        % density.  They can also compute photopigment variants corrected
        % for filtering by blood vessels.
        
        % Prompt user for key parameters that affect the spectral
        % sensitivities.
        %
        % Note that we don't typically vary or pass the blood vessel
        % parameters but rather simply accept the defaults used by
        % GetHumanPhotoreceptorSS.  It's mainly for fun that we show
        % how to do this here.
        observerAgeInYears = GetWithDefault('\tObserver age in years?', 32);
        fieldSizeDegrees = GetWithDefault('\tField size in degrees?', 27.5);
        pupilDiameterMm = GetWithDefault('\tPupil diameter?', 4.7);
        vesselOxyFraction = GetWithDefault(['\tOxygenation fraction for vessel hemoglobin [typical 0.85]?'], 0.85);
        vesselOverallThicknessUm = GetWithDefault(['\tVessel thickness [typical 5 um]?'], 5);
        correctBleaching = GetWithDefault(['\tCorrect for cone photopigment bleaching [1 = yes, 0 = no]?'],1);
        
        % Define photoreceptor classes that we'll consider.
        % ReceptorIsolate has a few more built-ins than these.
        photoreceptorClasses = {'LConeTabulatedAbsorbance', 'MConeTabulatedAbsorbance', 'SConeTabulatedAbsorbance', 'Melanopsin', 'Rods', ...
            'LConeTabulatedAbsorbancePenumbral', 'MConeTabulatedAbsorbancePenumbral', 'SConeTabulatedAbsorbancePenumbral'};
        
        % Correct for pigment bleaching if desired.  This is done
        % separately for open-field and penumbral cones.  The bleaching
        % correction routine only knows about L, M, and S cones.
        if correctBleaching
            % The hard work is done by routine
            % GetConeFractionBleachedFromSpectrum, which is in the
            % ContrastSplatter directory. We set the verbose flag here.
            verbose = true;
            [fractionBleachedFromIsom, fractionBleachedFromIsomHemo] = ...
                GetConeFractionBleachedFromSpectrum(S, B_primary*backgroundPrimary + ambientSpd, fieldSizeDegrees, observerAgeInYears, pupilDiameterMm, [], verbose);
            
            % Assign the fraction bleached for each photoreceptor class.
            for p = 1:length(photoreceptorClasses)
                switch photoreceptorClasses{p}
                    case 'LCone'
                        fractionBleached(p) = fractionBleachedFromIsom(1);
                    case 'MCone'
                        fractionBleached(p) = fractionBleachedFromIsom(2);
                    case 'SCone'
                        fractionBleached(p) = fractionBleachedFromIsom(3);
                    case 'LConeHemo'
                        fractionBleached(p) = fractionBleachedFromIsomHemo(1);
                    case 'MConeHemo'
                        fractionBleached(p) = fractionBleachedFromIsomHemo(2);
                    case 'SConeHemo'
                        fractionBleached(p) = fractionBleachedFromIsomHemo(3);
                    otherwise
                        fractionBleached(p) = 0;
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
        oxygenationFraction = [];
        vesselThickness = [];
        T_receptors = GetHumanPhotoreceptorSS(S, photoreceptorClasses, fieldSizeDegrees, observerAgeInYears, pupilDiameterMm, [], fractionBleached, oxygenationFraction, vesselThickness);
        
        %% Let user choose a photoreceptor class to target
        fprintf('Available photoreceptor classes to target:\n');
        fprintf('\t[1]  Melanopsin, silence open-field cones; ignore rods and penumbral cones\n');
        fprintf('\t[2]  Melanopsin, silence open-field and penumbral cones; ignore rods)\n');
        fprintf('\t[3]  S cones, silence open-field L and M cones, melanopsin, and prenumbral L and M cones; ignore rods and penumbral S cones\n');
        fprintf('\t[4]  S cones, silence open field L and M cones, ingore all others\n');
        fprintf('\t[5]  Penumbral L and M cones, silence open-field cones, melanopsin, and prenumbral S cones; ignore rods\n');
        fprintf('\t[6]  Rods\n');
        whichDirectionNumber = GetWithDefault('Enter direction',1);
        
        % Depending on which direction is chosen, specify the indices
        % into the rows of T_receptors to define the various classes.
        %  Targeted photoreceptors - contrast maximized or driven to specified target contrast
        %  Ignored photoreceptors - ignored in calculation
        %  Minmized photoreceptors - legacy variable no longer used
        %  All remaining photoreceptors are silenced.
        switch (whichDirectionNumber)
            case 1
                whichDirection = 'MelanopsinDirectedLegacy';
                whichReceptorsToTarget = [4];
                whichReceptorsToIgnore = [5 6 7 8];
                whichReceptorsToMinimize = [];
            case 2
                whichDirection = 'MelanopsinDirected';
                whichReceptorsToTarget = [4];
                whichReceptorsToIgnore = [5];
                whichReceptorsToMinimize = [];
            case 3
                whichDirection = 'SDirected';
                whichReceptorsToTarget = [3];
                whichReceptorsToIgnore = [5 8];
                whichReceptorsToMinimize = [];
            case 4
                whichDirection = 'SDirected';
                whichReceptorsToTarget = [3];
                whichReceptorsToIgnore = [4 5 6 7 8];
                whichReceptorsToMinimize = [];
            case 5
                whichDirection = 'PenumbralLM';
                whichReceptorsToTarget = [6 7];
                whichReceptorsToIgnore = [5];
                whichReceptorsToMinimize = [];
            case 6
                whichDirection = 'Rods';
                whichReceptorsToTarget = [5];
                whichReceptorsToIgnore = [6 7 8];
                whichReceptorsToMinimize = [];
            otherwise
                error('Unknown direction entered');
        end
        
        
        
        
        
% User chooses whether to maximize contrast in targeted receptor classes or
% or get it as close to a specified value as possible.
%
% If we target, here we specify the same contrast for all targeted classes.
% This is not necessary, they can differ.  It just makes the demo code a
% bit simpler to yoke them since we only have to prompt for one number.
maximizeTargetContrast = GetWithDefault('\tMaximize contrast? [1 = yes, 0 = no]', 1);
if maximizeTargetContrast
    desiredContrast = [];
else
    desiredContrast = GetWithDefault('\tDesired contrast (applies to all targeted classes)?', 0.45)*ones(size(whichReceptorsToTarget));
end

% Nice message for user
fprintf('\nGenerating stimuli which isolate receptor classes');
for i = 1:length(whichReceptorsToTarget)
    fprintf('\n  - %s', photoreceptorClasses{whichReceptorsToTarget(i)});
end
fprintf('\nGenerating stimuli which ignore receptor classes');
if (~length(whichReceptorsToIgnore) == 0)
    for i = 1:length(whichReceptorsToIgnore)
        fprintf('\n  - %s', photoreceptorClasses{whichReceptorsToIgnore(i)});
    end
else
    fprintf('\n  - None');
end
fprintf('\nThe remaining classes will be silenced\n');

%% Call the optimization routine.
%
% Careful examaination of the arguments will reveal that the initialGuess
% for the primaries is set to the background value for the primaries, so
% that the constraints are all met at the start of the search.  The
% optimization routine is much happier when it is done this way -- bad
% things happen if you start with a guess that violates constraints.
modulationPrimary = ReceptorIsolate(T_receptors,whichReceptorsToTarget, whichReceptorsToIgnore, whichReceptorsToMinimize, ...
    B_primary, backgroundPrimary, backgroundPrimary, whichPrimariesToPin,...
    primaryHeadRoom, maxPowerDiff, desiredContrast, ambientSpd);

%% Compute the contrasts that we got.
backgroundReceptors = T_receptors*(B_primary*backgroundPrimary + ambientSpd);
modulationReceptors = T_receptors*B_primary*(modulationPrimary - backgroundPrimary);
contrastReceptors = modulationReceptors ./ backgroundReceptors;
for j = 1:size(T_receptors,1)
    fprintf('\t%s: contrast = %0.4f\n',photoreceptorClasses{j},contrastReceptors(j));
end

%% Plots
plotDir = fullfile(codeBaseDir,'outputPlots');
if ~isdir(plotDir)
    mkdir(plotDir);
end
curDir = pwd;
cd(plotDir);

% Photoreceptor sensitivities
theFig1 = figure; clf; hold on
plot(SToWls(S),T_receptors,'LineWidth',2);
xlabel('Wavelength (nm)')
ylabel('Sensitivity');
title('Normalized photoreceptor sensitivities');
saveas(theFig1,sprintf('%s_%s_%s_Sensitivities.pdf',whichModel,whichPrimaries,whichDirectionNumber),'pdf');

% Modulation spectra
theFig2 = figure; hold on
plot(SToWls(S),B_primary*modulationPrimary,'r','LineWidth',2);
plot(SToWls(S),B_primary*backgroundPrimary,'k','LineWidth',2);
title('Modulation spectra');
xlim([380 780]);
xlabel('Wavelength');
ylabel('Power');
pbaspect([1 1 1]);
saveas(theFig2,sprintf('%s_%s_%s_Modulation.pdf',whichModel,whichPrimaries,whichDirectionNumber),'pdf');

% Primaries
theFig3 = figure; hold on
plot(modulationPrimary,'r','LineWidth',2);
plot(backgroundPrimary,'k','LineWidth',2);
title('Primary settings');
xlim([0 length(backgroundPrimary)]);
ylim([0 1]);
xlabel('Primary Number (nominal)');
ylabel('Setting');
saveas(theFig3,sprintf('%s_%s_%s_Primaries.pdf',whichModel,whichPrimaries,whichDirectionNumber),'pdf');


