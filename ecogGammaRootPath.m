function rootPath = ecogGammaRootPath
% Return the path to the root ECoG Gamma directory
%
% This function must reside in the directory at the base of the ECoG Gamma
% directory structure.  It is used to determine the location of various
% sub-directories.
% 
% Example:
%   fullfile(ecogGammaRootPath,'data')

rootPath=which('ecogGammaRootPath');

rootPath = fileparts(rootPath);

return