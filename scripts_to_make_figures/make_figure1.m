function ok = make_figure1()
% ok = make_figure1()
%
% Function to reproduce panels B and C of figure 1 from the paper
%
% Hermes, D., Miller, K.J., Wandell, B.A., Winawer, J. (2014). 
% Stimulus dependence of gamma oscillations in human visual cortex. 
% Cerebral Cortex, 2014
%
% The figure produces the mean ECoG time-frequency plots for 7 stimulus
% classes for 2 example electrodes in V1, one in Subject 1 and one in Subject 2.


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


ok = false;

figure('Color',[1 1 1],'Position',[0 0 1200 360])

for subject = 1:2
    
    disp(['plotting data for subject ' int2str(subject)])
    
    %%%% LOAD THE DATA
    if subject==1
        load('subj1data_chan112','data_chan112','onsets','stims','srate','stimnames')
        data=data_chan112;
    elseif subject==2
        load('subj2data_chan85','data_chan85','onsets','stims','srate','stimnames')
        data=data_chan85;
    end
    % data_chanXXX  % raw data from the V1/V2v channel
    % onsets        % stimulus onsets in ms
    % srate         % sampling rate
    % stims         % stimulus numbers, one for each onset
    % stimnames     % stimulus names for stimulus 1:8

    %%%% NOTCH FILTER
    % notch filter data at 60, 120 and 180 Hz
    data = ecog_notch(data,srate);
    
    %%%% MAKE EPOCHS
    % choose onsets for epochs
    onset_trial=round(onsets);%samples
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
    
    disp(['done making epochs subject ' int2str(subject)])
        
    % just for the plot number
    % plotting order, as in paper: grating 4,8,16,32 cycles, plaid, noise patterns
    plot_nr=[8 7 6 1 2 3 4 10 9 5];
    
    nr_cond=max(stims);
    
    %%%% DEFINE TIME-FREQUENCY MULTITAPER SETTINGS
    movingwin=[.200 .05];
    params.pad=-1;
    params.tapers=[3 5];
    params.fpass=[0 200];
    params.Fs=srate; 
    params.trialave=1;
    
    %%%% CALCULATING AND PLOTTING THE SPECTOGRAMS
    subplot(2,1,subject),hold on
    
    disp(['calculating time-frequency plots subject ' int2str(subject)])
    
    % define baseline
    data2use=squeeze(data_epoch(elec,stims==8,t>.25 & t<.5))';
    [S1b,t_tf_b,f_b]=mtspecgramc(data2use,movingwin,params);
    S1b=mean(mean(S1b,3),1);
        
    %%%%% all responses:
    for k=1:nr_cond
        data2use=squeeze(data_epoch(elec,stims==k,:))';
        
        % calculate specgram
        [S1,t_tf,f]=mtspecgramc(data2use,movingwin,params);
        
        % make sure the first time point is correct
        t_tf=t_tf+t(1);
        
        % normalize wrt baseline
        S1=S1./repmat(S1b,[size(S1,1) 1]);
        
        % and plot:
        subplot(2,10,plot_nr(k) + 10 * (subject-1))
        imagesc(t_tf,f,log10(S1)',[-1.3 1.3])
        axis xy
        set(gca,'XTick',[0 .5])
    end
    
    % ad condition labels to each subplot
    cond_label={'white n','pink n','brown n','gr 4','gr 8','gr 16','gr 32','blank','bin white n','plaid'};
    for k=1:nr_cond
        subplot(2,10,plot_nr(k) + 10 * (subject-1))
        title(cond_label{k})
    end
end

subplot(2,10,1)
ylabel('subject 1')
subplot(2,10,11)
ylabel('subject 2')

ok = true;

return
