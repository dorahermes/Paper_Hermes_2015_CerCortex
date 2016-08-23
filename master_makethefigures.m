
% download the Chronux toolbox (Mitra and Bokil, 2008) http://www.chronux.org 
% make sure it is in your path (added in ecogGammaAddPaths)

% check whether chronux is in the path:
if exist('mtspecgramc','file')==0
    disp('ERROR: make sure Chronux is in your Matlab path')
end


%% set the root paths and add the correct paths:
ecogGammaAddPaths


%% Make figure 1

make_figure1
% reproduces panels B and C of figure 1 from the paper

%% Make figure 2

make_figure2
% reproduces panels A and B of figure 2 from the paper

%% Make figure 3

make_figure3
% reproduces panels B and C of figure 3 from the paper
% function generates 14 figures, one for each pannel, with narrowband and broadband 
% weights for all electrodes for each stimulus

%% Make figure 4

make_figure4
% reproduces panel B from figure 4 from the paper


%% Make figure 5

make_figure5
% reproduces panels A and B from figure 5 from the paper
