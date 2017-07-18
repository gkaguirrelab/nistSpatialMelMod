
codeBaseDir = tbLocateProject('nistSpatialMelMod','verbose',false);

calibrationFileName = fullfile(codeBaseDir, 'demoOneLightCalFile', 'OneLightDemoCal.mat');

[ backgroundPrimary, modulationPrimary ] = nistSpatialMelMod_makeModulationSPD( calibrationFileName );