
function ecogGammaAddPaths
% Set Matlab directory path for ECOG Gamma project
%
%     ecogGammaAddPaths


rootPath = ecogGammaRootPath;
% fprintf('ecogGammaRootPath root directory: %s\n',rootPath)

% Adds the root directory to the user's path
addpath(rootPath);

% Generates a list of the directories below the root path.
addpath(fullfile(rootPath, 'data'))
addpath(genpath(fullfile(rootPath, 'Chronux')))
addpath(fullfile(rootPath, 'scripts_to_make_figures'))
addpath(fullfile(rootPath, 'support_functions'))

return;

