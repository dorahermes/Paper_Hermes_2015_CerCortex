function ok = make_figure4()
% ok = make_figure4()
%
% Function to reproduce panels B from figure 4 from the paper

% Hermes, D., Miller, K.J., Wandell, B.A., Winawer, J. (2014). 
% Stimulus dependence of gamma oscillations in human visual cortex. 
% Cerebral Cortex, 2014


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

%%
%%
%% EPOCH AND CALCULATE POWERSPECTRUM
%%
%%
clear all

load('./data/subj1data_chan112_fh','data_chan112','onset_trial','offset_trial','stims','srate','stimnames')
% data_chan112  % raw data from the V1/V2v channel
% onset_trial   % stimulus onset in samples
% offset_trial  % stimulus offset in samples
% srate         % sampling rate
% stims         % stimulus numbers, one for each onset
% stimnames     % stimulus names for stimulus 1:2

data=data_chan112;

%%%% NOTCH FILTER 
% notch filter data at 60, 120 and 180 Hz
data = ecog_notch(data,srate);
    
%%%% MAKE EPOCHS
epoch_l=1.2; % -0.2:1 sec
data_epoch=zeros(size(data,2),length(onset_trial),epoch_l*srate);
data_epoch_off=zeros(size(data,2),length(offset_trial),epoch_l*srate);
for elec=1:size(data,2)
    for l=1:length(onset_trial)
        data_epoch(elec,l,:)=...
            data(onset_trial(l)-.2*srate+1:onset_trial(l)+(epoch_l-.2)*srate,elec)';
    end
    for l=1:length(offset_trial)
        data_epoch_off(elec,l,:)=...
            data(offset_trial(l)-.2*srate+1:offset_trial(l)+(epoch_l-.2)*srate,elec)';
    end

end
t=[1:epoch_l*srate]/srate - 0.2;    

disp('made epochs')

%%
%%
%% Make Figure 4
%%
%%

%%%% DEFINE TIME-FREQUENCY MULTITAPER SETTINGS
movingwin=[.200 .050];
params.pad=-1;
params.tapers=[3 5];
params.fpass=[0 200];
params.Fs=srate;
params.trialave=1;

elec=1; % there is only one electrode

figure('Color',[1 1 1],'Position',[0 0 400 200])

disp('calculating time-frequency plots')

% define baseline from ITI
data2use=squeeze(data_epoch_off(elec,:,t>.25 & t<.5))';
[S1b]=mtspecgramc(data2use,movingwin,params);
S1b=mean(S1b,1); % average over time

%%%%% calculate spectrogram trials
for k=1:2 % faces, houses
    data2use=squeeze(data_epoch(elec,stims==k,:))';
    
    % calculate spectgram
    [S1,t_tf,f]=mtspecgramc(data2use,movingwin,params);
    t_tf=t_tf+t(1);
    % normalize wrt baseline
    S1=S1./repmat(S1b,[size(S1,1) 1]);
    subplot(1,3,k)
    imagesc(t_tf,f,log(S1)',[-3 3])
    axis xy
end

%%%%% calculate spectrogram baseline (ITI)
data2use=squeeze(data_epoch_off(elec,:,:))';
% calculate spectgram
[S1,t_tf,f]=mtspecgramc(data2use,movingwin,params);
t_tf=t_tf+t(1);
% normalize wrt baseline
S1=S1./repmat(S1b,[size(S1,1) 1]);
subplot(1,3,3)
imagesc(t_tf,f,log(S1)',[-3 3])
axis xy

% add titles to the subplots
for k=1:2 
    subplot(1,3,k)
    title(stimnames{k})
end
subplot(1,3,3)
title('ITI')

