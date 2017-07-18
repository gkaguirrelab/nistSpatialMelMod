
codeBaseDir = tbLocateProject('nistSpatialMelMod','verbose',false);

calibrationFileName = fullfile(codeBaseDir, 'demoOneLightCalFile', 'OneLightDemoCal.mat');

[ backgroundPrimary, modulationPrimary, B_primary, ambientSpd, T_receptors ] = nistSpatialMelMod_makeModulationPrimaries( calibrationFileName, 'desiredContrast',0.3, 'verbose', true );

figure
plot(backgroundPrimary,'-k');
hold on
plot(modulationPrimary,'-r');

[ backgroundPrimary, modulationPrimary, B_primary, ambientSpd, T_receptors ] = nistSpatialMelMod_makeModulationPrimaries( calibrationFileName, 'desiredContrast',0.15, 'verbose', true );

plot(modulationPrimary,'-.b');


backgroundReceptors = T_receptors*(B_primary*backgroundPrimary + ambientSpd);
modulationReceptors = T_receptors*B_primary*(modulationPrimary - backgroundPrimary);
contrastReceptors = modulationReceptors ./ backgroundReceptors;




% %% Plots
% plotDir = fullfile(codeBaseDir,'outputPlots');
% if ~isdir(plotDir)
%     mkdir(plotDir);
% end
% curDir = pwd;
% cd(plotDir);
% 
% % Photoreceptor sensitivities
% theFig1 = figure; clf; hold on
% plot(SToWls(S),T_receptors,'LineWidth',2);
% xlabel('Wavelength (nm)')
% ylabel('Sensitivity');
% title('Normalized photoreceptor sensitivities');
% %saveas(theFig1,sprintf('%s_%s_%s_Sensitivities.pdf',whichModel,whichPrimaries,whichDirectionNumber),'pdf');
% 
% % Modulation spectra
% theFig2 = figure; hold on
% plot(SToWls(S),B_primary*modulationPrimary,'r','LineWidth',2);
% plot(SToWls(S),B_primary*backgroundPrimary,'k','LineWidth',2);
% title('Modulation spectra');
% xlim([380 780]);
% xlabel('Wavelength');
% ylabel('Power');
% pbaspect([1 1 1]);
% %saveas(theFig2,sprintf('%s_%s_%s_Modulation.pdf',whichModel,whichPrimaries,whichDirectionNumber),'pdf');
% 
% % Primaries
% theFig3 = figure; hold on
% plot(modulationPrimary,'r','LineWidth',2);
% plot(backgroundPrimary,'k','LineWidth',2);
% title('Primary settings');
% xlim([0 length(backgroundPrimary)]);
% ylim([0 1]);
% xlabel('Primary Number (nominal)');
% ylabel('Setting');
% %saveas(theFig3,sprintf('%s_%s_%s_Primaries.pdf',whichModel,whichPrimaries,whichDirectionNumber),'pdf');

