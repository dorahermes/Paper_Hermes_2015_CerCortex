function ok = make_figure3()
% ok = make_figure3()
%
% Function to reproduce panels B and C of figure 3 from the paper
% 
% Hermes, D., Miller, K.J., Wandell, B.A., Winawer, J. (2014). 
% Stimulus dependence of gamma oscillations in human visual cortex. 
% Cerebral Cortex, 2014
%
% Function generates 14 figures, one for each pannel, with narrowband and 
% broadband weights for all electrodes for each stimulus


%     Copyright (C) 2014  D Hermes
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.


%% RENDER BroadBand and Gaussian weights

%%%% load the cortex rendering
load(['./data/subj1_cortex.mat'])
% mesh with right hemisphere rendering of subject 1

%%%% load the electrode positions
load(['./data/subj1_electrodes.mat'])
% XYZ coordinates of electrode positions

%%%% load the broadband and narrowband weights
load(['./data/subj1_resamp_stats.mat']) 
% 4D matrix: electrode X stimulus type X resampling number X fitted value
% 1)    111 electrodes
% 2)    stimulus type 1:8
%       {'white n','pink n','brown n','gr 4','gr 8','gr 16','gr 32','blank'}
% 3)    100 resamples
% 4)    fitted value: exponent, broadband weight, narrowband weight, narrowband frequency 

cond_labels={'white n','pink n','brown n','gr 4','gr 8','gr 16','gr 32','blank','b white n','plaid'};
cond_plot=[1:7];% plot these conditions

% select a view direction
els=elecmatrix;
v_d=[-60,-30]; % view direction makes sure electrodes are visible on surface
% make electrodes pop out just a little bit to be able to see them all
% WARNING: NO ROTATION AFTER THIS STEP as original location is only
% maintained for this view
a_offset=.1*max(abs(els(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
els=els+repmat(a_offset,size(els,1),1);

%%%%% plot BROADBAND weights
a=resamp_parms(:,:,:,2);% 2:broadband weight
% The broadband weight does not have a good zero, so we calcualte the
% median, and subtract the median of the baseline. This is only for
% plotting purposes, not for the statistics
a=median(a,3); % take the median for all the resamples
for k=1:size(resamp_parms,1)
    a(k,:)=a(k,:)-median(a(k,8)); 
end
% Calculate confidence interval for statistics (significance)
b=resamp_parms(:,:,:,2);% 2:broadband weight
b_up=prctile(b,97.5,3);
b_low=prctile(b,2.5,3);

for k=1:length(cond_plot)
    figure('Color',[1 1 1],'Position',[0 0 400 400])
    ctmr_gauss_plot(cortex,[0 0 0],0)
    w_plot=a(:,cond_plot(k));
    % set non-significant values to zero (overlapping with baseline CI)
    w_plot(b_low(:,cond_plot(k))<b_up(:,8))=0;
    el_add_sizable(els,w_plot,1)
    loc_view(v_d(1),v_d(2))
    title(['BB weight ' cond_labels{k}])
end

%%%%% plot GAUSSIAN weights
a=resamp_parms(:,:,:,3);% 3:gaussian weight
a=median(a,3);
b=resamp_parms(:,:,:,3);% 3:gaussian weight
b_up=prctile(b,97.5,3);
b_low=prctile(b,2.5,3);

for k=1:length(cond_plot)
    figure('Color',[1 1 1],'Position',[0 0 400 400])
    ctmr_gauss_plot(cortex,[0 0 0],0)
    w_plot=a(:,cond_plot(k));
    % set non-significant values to zero (overlapping with baseline CI)
    w_plot(b_low(:,cond_plot(k))<b_up(:,8))=0;
    el_add_sizable(els,w_plot,1.5)
    loc_view(v_d(1),v_d(2))
    title(['Gauss weight ' cond_labels{k}])
end

