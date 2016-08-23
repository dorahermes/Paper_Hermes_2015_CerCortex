function ok = make_figure2()
% ok = make_figure2()
%
% Function to reproduce panels A and B of figure 2 from the paper
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

ok=false;
% figure 2:
fH1=figure('Color',[1 1 1],'Position',[0 0 1200 360]);

% figure S2:
fH2=figure('Color',[1 1 1],'Position',[0 0 400 400]);

for subject = 1:2
    
    disp(['plotting data for subject ' int2str(subject)])

if subject==1
    load('subj1data_chan112','data_chan112','onsets','stims','srate','stimnames')
    data=data_chan112;
elseif subject==2
    load('subj2data_chan85','data_chan85','onsets','stims','srate','stimnames')
    data=data_chan85;
end

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

% BASELINE CORRECT
data_epoch = ecog_baselinesubtract(data_epoch,t>-.1 & t<0);

% CALCULATE SPECTRA
% if choosing to regress the erp out, make sure data are baseline corrected first!
fft_w=500; % window width
fft_t=t>0 & t<.5; % time segment for spectrum, stimulus was on for 500 ms
fft_ov=0; % overlap
reg_erp=1; % 1 to regress erp out, 0 not to
[f,data_epoch_spectra,data_epoch_erp] = ...
    ecog_spectra(data_epoch,stims,fft_w,fft_t,fft_ov,srate,1);

disp(['done calculating powerspectrs subject ' int2str(subject)])

%%
%%
%% MAKE FIGURE 2
%%
%%

figure(fH1)

% just for the plot number
% plotting order, as in paper: grating 4,8,16,32 cycles, plaid, noise patterns
plot_nr=[8 7 6 1 2 3 4 10 9 5];

nr_cond=max(stims); % subj1: 8, subj2: 10

elec=1; % there is just one electrode

data_fft=squeeze(data_epoch_spectra(elec,:,:));
stims_use=stims;

f_use4fit=[35:57 65:115 126:175 186:200];
f_sel=ismember(f,f_use4fit);

cond_label={'white n','pink n','brown n','gr 4','gr 8','gr 16','gr 32','blank','b white n','plaid'};
cond_color={[.2 .5 1],[0 .2 1],[0 0 .8],[.5 0 0],[1 0 0],[1 .5 .5],[1 .8 .8],[0 0 0],[.5 .8 1],'g'};
nr_resamp=100;
resamp_parms=NaN(nr_cond,nr_resamp,4);

for s=1:nr_cond  
    subplot(2,10,plot_nr(s) + 10 * (subject-1)),hold on
    
    % data to use for fit
    data_base=mean(data_fft(stims_use==8,:),1); % spectrum baseline
    data_fit=mean(data_fft(stims_use==s,:),1);  % spectrum for stimulus S
    
    % calculate fit based on the average
    [out_exp,w_pwr,w_gauss,gauss_f,fit_f2] = ...
        fit_gammadata(f,f_use4fit,data_base,data_fit);
    % plot fit 
    plot(log10(f),fit_f2,'Color',[.5 .5 .5],'LineWidth',1)  
    % plot baseline data
    plot(log10(f(f_sel)),log10(data_base(f_sel)),'k--','LineWidth',1)      
    % plot stimulus data
    plot(log10(f(f_sel)),log10(data_fit(f_sel)),'-','LineWidth',1,'Color',cond_color{s})       

    % some plot settings:
    xlim([log10(35) log10(200)])
    set(gca,'XTick',[log10(25) log10(50) log10(100) log10(200)])
    set(gca,'XTickLabel',{'25','50','100','200'})
    ylim([-2 2.2])
    title(cond_label{s})
    
    % do the bootstrap
    disp(['bootstrap stim ' int2str(s)])
    for k=1:nr_resamp
        data_fit=data_fft(stims_use==s,:);
        
        tr_boot=randsample(size(data_fit,1),size(data_fit,1),'true');
        data_fit=mean(data_fit(tr_boot,:),1);
        [out_exp,w_pwr,w_gauss,gauss_f,fit_f2] = ...
            fit_gammadata(f,f_use4fit,data_base,data_fit);
        resamp_parms(s,k,1)=out_exp;    % exponent
        resamp_parms(s,k,2)=w_pwr;      % broadband weight 
        resamp_parms(s,k,3)=w_gauss;    % gaussian weight
        resamp_parms(s,k,4)=gauss_f;    % gaussian frequency, note that this is biased since it is on top of the powerlaw
        % WARNING: do not use the Gaussian frequency for any conclusions, it is fitted on top of a powerlaw, and an increase in power will bias it towards lower frequencies 
    end
end

% resamp_parms is el X bootstrap X fitted value

%%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%%
%%%%% and make supplemental figure 2 %%%%%
%%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%% %%%%%
figure(fH2)

subplot(2,2,subject),hold on
for k=1:nr_cond
    plot([plot_nr(k) plot_nr(k)],[prctile(squeeze(resamp_parms(k,:,3)),2.5)...%lower
        prctile(squeeze(resamp_parms(k,:,3)),97.5)],'k');%upper
    plot(plot_nr(k),prctile(squeeze(resamp_parms(k,:,3)),50),'.','Color',cond_color{k},...
        'MarkerSize',20)    
end
title('narrowband weight')
xlim([0.5 10.5])
ylim([-0.07 1.8])
set(gca,'XTick',[1:10],'XTickLabel',[])

subplot(2,2,subject+2),hold on
for k=1:nr_cond
    plot([plot_nr(k) plot_nr(k)],[prctile(squeeze(resamp_parms(k,:,2)),2.5)...%lower
        prctile(squeeze(resamp_parms(k,:,2)),97.5)],'k');%upper
    plot(plot_nr(k),prctile(squeeze(resamp_parms(k,:,2)),50),'.','Color',cond_color{k},...
        'MarkerSize',20)
end
title('broadband weight')
xlim([0.5 10.5])
set(gca,'XTick',[1:10],'XTickLabel',[])

end

figure(fH1)
subplot(2,10,1)
ylabel('subject 1')
subplot(2,10,11)
ylabel('subject 2')

figure(fH2)
subplot(2,2,1)
ylabel('subject 1')
subplot(2,2,3)
ylabel('subject 2')

ok=true;
return


