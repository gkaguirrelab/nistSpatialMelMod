function nistSpatialMelModLocalHook
% nistSpatialMelModLocalHook - Configure things for working on Joe Rice's
%  hyper-spectral display at NIST.
%
% For use with the ToolboxToolbox.  If you copy this into your
% ToolboxToolbox localToolboxHooks directory (by default,
% ~/localToolboxHooks) and delete "LocalHooksTemplate" from the filename,
% this will get run when you execute tbUse({'OLFlickerSensitivityConfig'}) to set up for
% this project.  You then edit your local copy to match your local machine.
%
% The main thing that this does is define Matlab preferences that specify input and output
% directories.
%
% You will need to edit the project location and i/o directory locations
% to match what is true on your computer.

% Say hello
fprintf('Running nistSpatialMelMod Local Hook\n');

% Find the project directory, add it to the path, save this as a
%  pref, and then make this the current directory
projectDir = fullfile(tbLocateProject('nistSpatialMelMod'));
addpath(genpath(projectDir));
setpref('nistSpatialMelMod', 'projectDir', projectDir);