# Paper_Hermes_2015_CerCortex

README

The code and data in this directory are a supplement to:

Hermes, D., Miller, K.J., Wandell, B.A., Winawer, J. (2014). Stimulus dependence of gamma oscillations in human visual cortex. Cerebral Cortex.

All code in this repository is written in MATLAB (Mathworks) and, together with the included data, can be used to reproduce the plots from the publication.

The code reproduces the ECoG data panels in the main text of the paper and one panel from the supplement: Main figures 1, 2, 3, 4, 5 and Figure S2.

Code and data are provided as part of the goal of ensuring that computational methods are reproducible by other researchers. 

******************
** DEPENDENCIES **
******************
Matlab toolboxes needed:
	optimization
	statistics
	signal processing

Other toolboxes needed:
	Chronux (Mitra and Bokil, 2008), freely available at http://www.chronux.org
	
***************
** EXAMPLE ****
***************
Reproduce plots from Figure 1 with the following call from the Matlab command window.
	ecogGammaAddPaths;
	make_figure1;

Reproduce all plots with the following call from the Matlab command window. 
	ecogGammaAddPaths;
	master_makethefigures;
	
***************
** CONTENTS ***
***************

Top level:
- README: This is where you are now ;) 
- ecogGammaAddPaths: function to add paths needed for computations and plots
- ecogGammaRootPath: function to returns path to directory with data and functions
- master_makethefigures: master script that calls the separate functions to make the figures
 
scripts_to_make_figures: Scripts to reproduce the main data plots in the publication
- make_figure1
- make_figure2
- make_figure3
- make_figure4
- make_figure5
 
data: Matlab data files ('*.mat') needed for plots. 
- subj1_cortex.mat		% subject 1 cortex rendering
	% mesh with right hemisphere rendering of subject 1
- subj1_electrodes.mat		% subject 1 electrode positions
	% XYZ coordinates of electrode positions
- subj1_resamp_stats.mat	% subject 1 resampled data from for figure 3
	% 4D matrix: electrode X stimulus type X resampling number X fitted value
	% 1)    111 electrodes
	% 2)    stimulus type 1:8
	%       {'white n','pink n','brown n','gr 4','gr 8','gr 16','gr 32','blank'}
	% 3)    100 resamples
	% 4)    fitted value: exponent, broadband weight, narrowband weight, narrowband frequency 
- subj1data_chan112_fh.mat 	% subject 1 data from for figure 4 and 5
	% data_chan112  % raw data from the V1/V2v channel
	% onset_trial   % stimulus onset in samples
	% offset_trial  % stimulus offset in samples
	% srate         % sampling rate (1000 Hz)
	% stims         % stimulus numbers, one for each onset
	% stimnames     % stimulus names for stimulus 1:2
- subj1data_chan112.mat		% subject 1 data from for figure 1, 2 and 5
- subj2data_chan85.mat		% subject 2 data from for figure 1 and 2 
	% data_chan85   % raw data from the V1/V2v channel
	% onsets        % stimulus onsets in ms
	% srate         % sampling rate (1000 Hz)
	% stims         % stimulus numbers, one for each onset
	% stimnames     % stimulus names for stimuli 1:10

support_functions: Support functions called by the main scripts 
- ctmr_gauss_plot.m
- ecog_baselinesubtract.m
- ecog_notch.m
- ecog_spectra.m
- el_add_sizable.m
- el_add.m
- fit_func3_loglog.m
- fit_gammadata.m
- loc_colormap.mat
- loc_view.m
- tripatch.m
