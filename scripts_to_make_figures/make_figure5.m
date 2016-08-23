function ok = make_figure5()
% ok = make_figure5()
%
% Function to reproduce panels A and B from figure 5 from the paper

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


%% Figure 5A: gratings and noise individual trials


load('./data/subj1data_chan112','data_chan112','onsets','stims','srate','stimnames')
% data_chan112  % raw data from the V1/V2v channel
% onsets        % stimulus onsets in ms
% srate         % sampling rate
% stims         % stimulus numbers, one for each onset
% stimnames     % stimulus names for stimulus 1:8

data=data_chan112;

%%%% NOTCH FILTER 
% notch filter data at 60, 120 and 180 Hz
data = ecog_notch(data,srate);
    
%%%% MAKE EPOCHS
% choose onsets for epochs
onset_trial=round(onsets);% ms to samples
epoch_l=1.2; % epoch length: -0.2:1 sec
data_epoch=zeros(size(data,2),length(onset_trial),epoch_l*srate);
for elec=1:size(data,2)
    for l=1:length(onset_trial)
        data_epoch(elec,l,:)=...
            data(onset_trial(l)-.2*srate+1:onset_trial(l)+(epoch_l-.2)*srate,elec)';
    end
end
% define t - time vector for each epoch
t=[1:epoch_l*srate]/srate - 0.2;    

% BASELINE CORRECT
data_epoch = ecog_baselinesubtract(data_epoch,t>-.1 & t<0);

% CALCULATE SPECTRA
disp('calculating powerspectra gratings and noise patterns')
% if choosing to regress the erp out, make sure data are baseline corrected first!
fft_w=500; % window width
fft_t=t>0 & t<.5; % time segment for spectrum, stimulus was on for 500 ms
fft_ov=0; % overlap
reg_erp=1; % 1 to regress erp out, 0 not to
[f,data_epoch_spectra,data_epoch_erp] = ...
    ecog_spectra(data_epoch,stims,fft_w,fft_t,fft_ov,srate,1);

% FITTING 
disp('fitting broadband + gaussian to grating/noise trials')
resamp_parms=NaN(size(data_epoch_spectra,2),4);
data_fft=squeeze(data_epoch_spectra(1,:,:));
f_use4fit=[35:57 65:115 126:175 186:200];
f_sel=ismember(f,f_use4fit);

for s=1:size(data_fft,1) % trials    
    data_base=mean(data_fft(stims==8,:),1);
    data_fit=mean(data_fft(s,:),1);

    [out_exp,w_pwr,w_gauss,gauss_f,fit_f2] = ...
        fit_gammadata(f,f_use4fit,data_base,data_fit);

    resamp_parms(s,1)=out_exp;
    resamp_parms(s,2)=w_pwr;
    resamp_parms(s,3)=w_gauss;
    resamp_parms(s,4)=gauss_f;
end

%% plotting single trial data noise and gratings

figure

cond_color={[.6 .6 .6],[.4 .4 .4],[.2 .2 .2],[.2 .2 .2],[.4 .4 .4],[.6 .6 .6],[.8 .8 .8],[.5 .5 .5],[.1 .1 .1],'g'};
cond_symbol={'*','*','*','v','v','v','v','o','s','d'};

count_trials=1;
plotting_order=[3 2 1 4 5 6 7 8];

for s=plotting_order
    subplot(2,1,1),hold on
    plot([count_trials:count_trials+length(find(stims==s))-1],resamp_parms(stims==s,3),cond_symbol{s},'Color',cond_color{s},'LineWidth',1)
    plot([count_trials+5 count_trials+length(find(stims==s))-6],[median(resamp_parms(stims==s,3)) median(resamp_parms(stims==s,3))],...
        'Color',[0 0 0],'LineWidth',1)
    plot([count_trials+5 count_trials+length(find(stims==s))-6],[quantile(resamp_parms(stims==s,3),.25) quantile(resamp_parms(stims==s,3),.25)],...
        '--','Color',[0 0 0],'LineWidth',1)
    plot([count_trials+5 count_trials+length(find(stims==s))-6],[quantile(resamp_parms(stims==s,3),.75) quantile(resamp_parms(stims==s,3),.75)],...
        '--','Color',[0 0 0],'LineWidth',1)
    xlim([-1 length(stims)+1])
    ylim([min(resamp_parms(:,3))-0.1 max(resamp_parms(:,3))+0.1])
    
    subplot(2,1,2),hold on
    plot([count_trials:count_trials+length(find(stims==s))-1],resamp_parms(stims==s,2),cond_symbol{s},'Color',cond_color{s},'LineWidth',1)
    plot([count_trials+5 count_trials+length(find(stims==s))-6],[median(resamp_parms(stims==s,2)) median(resamp_parms(stims==s,2))],...
        'Color',[0 0 0],'LineWidth',1)
    plot([count_trials+5 count_trials+length(find(stims==s))-6],[quantile(resamp_parms(stims==s,2),.25) quantile(resamp_parms(stims==s,2),.25)],...
        '--','Color',[0 0 0],'LineWidth',1)
    plot([count_trials+5 count_trials+length(find(stims==s))-6],[quantile(resamp_parms(stims==s,2),.75) quantile(resamp_parms(stims==s,2),.75)],...
        '--','Color',[0 0 0],'LineWidth',1)
    xlim([-1 length(stims)+1])

    
    count_trials=count_trials+length(find(stims==s));
end

subplot(2,1,1),hold on
ylabel('narrowband weight')

subplot(2,1,2),hold on
ylabel('broadband weight')


%%
%%
%% Figure 5B: faces and houses individual trials
%%
clear all
disp('loading face/house data')
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

% BASELINE CORRECT
data_epoch = ecog_baselinesubtract(data_epoch,t>-.1 & t<0);
data_epoch_off = ecog_baselinesubtract(data_epoch_off,t>-.1 & t<0);

% CALCULATE SPECTRA
disp('calculating face/house spectra')
% if choosing regress erp out, make sure data are baseline corrected first!
fft_w=500; % window width
fft_t=t>0 & t<.5; % time segment for spectrum
fft_ov=0; % overlap
reg_erp=1; % 1 to regress erp out, 0 not to
[f,data_epoch_spectra,data_epoch_erp] = ...
    ecog_spectra(data_epoch,stims,fft_w,fft_t,fft_ov,srate,1);
[f,data_epoch_off_spectra,data_epoch_off_erp] = ...
    ecog_spectra(data_epoch_off,ones(size(data_epoch_off,2),1),fft_w,fft_t,fft_ov,srate,1);

% FITTING:
disp('fitting broadband + gaussian to face/house trials')
% resample 1 electrode all trials
data_fft=squeeze(data_epoch_spectra(1,:,:));
data_fft_base=squeeze(data_epoch_off_spectra(1,:,:));
f_use4fit=[35:57 65:115 126:175 186:200];
f_sel=ismember(f,f_use4fit);

resamp_parms_trials=NaN(size(data_fft,1),4);
resamp_parms_base=NaN(size(data_fft_base,1),4);
%face/house trials
for s=1:size(data_fft,1) % nr trials
    data_base=mean(data_fft_base,1);
    data_fit=data_fft(s,:);
    [out_exp,w_pwr,w_gauss,gauss_f,fit_f2] = ...
        fit_gammadata(f,f_use4fit,data_base,data_fit);
    resamp_parms_trials(s,1)=out_exp;
    resamp_parms_trials(s,2)=w_pwr;
    resamp_parms_trials(s,3)=w_gauss;
    resamp_parms_trials(s,4)=gauss_f;
end
%baseline trials
for s=1:size(data_fft_base,1) % nr trials
    data_base=mean(data_fft_base,1);
    data_fit=data_fft_base(s,:);
    [out_exp,w_pwr,w_gauss,gauss_f,fit_f2] = ...
        fit_gammadata(f,f_use4fit,data_base,data_fit);
    resamp_parms_base(s,1)=out_exp;
    resamp_parms_base(s,2)=w_pwr;
    resamp_parms_base(s,3)=w_gauss;
    resamp_parms_base(s,4)=gauss_f;
end

%% plotting single trial data faces and houses

figure

cond_color={[0 0 0],[.8 .8 .8]};
cond_symbol={'d','s'};
count_trials=1;

for s=1:2
    subplot(2,1,1),hold on
    plot([count_trials:count_trials+length(find(stims==s))-1],resamp_parms_trials(stims==s,3),cond_symbol{s},'Color',cond_color{s},'LineWidth',1)
    plot([count_trials+5 count_trials+length(find(stims==s))-6],[median(resamp_parms_trials(stims==s,3)) median(resamp_parms_trials(stims==s,3))],...
        'Color',[0 0 0],'LineWidth',1)
    plot([count_trials+5 count_trials+length(find(stims==s))-6],[quantile(resamp_parms_trials(stims==s,3),.25) quantile(resamp_parms_trials(stims==s,3),.25)],...
        '--','Color',[0 0 0],'LineWidth',1)
    plot([count_trials+5 count_trials+length(find(stims==s))-6],[quantile(resamp_parms_trials(stims==s,3),.75) quantile(resamp_parms_trials(stims==s,3),.75)],...
        '--','Color',[0 0 0],'LineWidth',1)
    ylim([min(resamp_parms_trials(:,3))-0.1 max(resamp_parms_trials(:,3))+0.1])
    
    subplot(2,1,2),hold on
    plot([count_trials:count_trials+length(find(stims==s))-1],resamp_parms_trials(stims==s,2),cond_symbol{s},'Color',cond_color{s},'LineWidth',1)
    plot([count_trials+5 count_trials+length(find(stims==s))-6],[median(resamp_parms_trials(stims==s,2)) median(resamp_parms_trials(stims==s,2))],...
        'Color',[0 0 0],'LineWidth',1)
    plot([count_trials+5 count_trials+length(find(stims==s))-6],[quantile(resamp_parms_trials(stims==s,2),.25) quantile(resamp_parms_trials(stims==s,2),.25)],...
        '--','Color',[0 0 0],'LineWidth',1)
    plot([count_trials+5 count_trials+length(find(stims==s))-6],[quantile(resamp_parms_trials(stims==s,2),.75) quantile(resamp_parms_trials(stims==s,2),.75)],...
        '--','Color',[0 0 0],'LineWidth',1)

    count_trials=count_trials+length(find(stims==s));
end

subplot(2,1,1),hold on
plot([count_trials:count_trials+size(resamp_parms_base(:,2))-1],resamp_parms_base(:,3),'o','Color',[.5 .5 .5],'LineWidth',1)
plot([count_trials+5 count_trials+size(resamp_parms_base(:,2),1)-6],[median(resamp_parms_base(:,3)) median(resamp_parms_base(:,3))],...
    'Color',[0 0 0],'LineWidth',1)
plot([count_trials+5 count_trials+size(resamp_parms_base(:,2),1)-6],[quantile(resamp_parms_base(:,3),.25) quantile(resamp_parms_base(:,3),.25)],...
    '--','Color',[0 0 0],'LineWidth',1)
plot([count_trials+5 count_trials+size(resamp_parms_base(:,2),1)-6],[quantile(resamp_parms_base(:,3),.75) quantile(resamp_parms_base(:,3),.75)],...
    '--','Color',[0 0 0],'LineWidth',1)
ylim([-0.1 2.3])
xlim([0 count_trials+size(resamp_parms_base(:,2),1)+1])

subplot(2,1,2),hold on
plot([count_trials:count_trials+size(resamp_parms_base(:,2),1)-1],resamp_parms_base(:,2),'o','Color',[.5 .5 .5],'LineWidth',1)
plot([count_trials+5 count_trials+size(resamp_parms_base(:,2),1)-6],[median(resamp_parms_base(:,2)) median(resamp_parms_base(:,2))],...
    'Color',[0 0 0],'LineWidth',1)
plot([count_trials+5 count_trials+size(resamp_parms_base(:,2),1)-6],[quantile(resamp_parms_base(:,2),.25) quantile(resamp_parms_base(:,2),.25)],...
    '--','Color',[0 0 0],'LineWidth',1)
plot([count_trials+5 count_trials+size(resamp_parms_base(:,2),1)-6],[quantile(resamp_parms_base(:,2),.75) quantile(resamp_parms_base(:,2),.75)],...
    '--','Color',[0 0 0],'LineWidth',1)
xlim([0 count_trials+size(resamp_parms_base(:,2),1)+1])


subplot(2,1,1),hold on
ylabel('narrowband weight')

subplot(2,1,2),hold on
ylabel('broadband weight')
