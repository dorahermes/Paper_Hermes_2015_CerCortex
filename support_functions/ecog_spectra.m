function [f,data_epoch_spectra,data_epoch]=...
    ecog_spectra(data_epoch,stims,fft_w,fft_t,fft_ov,srate,reg_erp)

% [f,data_epoch_spectra,data_epoch]=...
%     ecog_spectra(data_epoch,stims,fft_w,fft_t,fft_ov,srate,reg_erp)
% data_epoch will have erp regressed out if option to choose so
% stims is a vector for each epoch, needed for condition-specific ERP
% calculation
% if choosing regress erp out, make sure data are baseline corrected first!


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



% regress erp out
if reg_erp==1 
    for k=1:size(data_epoch,1)%channels
        disp(['el ' int2str(k)])
        for m=1:size(data_epoch,2)%epochs
            x=squeeze(data_epoch(k,m,:));
            % regress ERP out
            s=stims(m);
            av_erp=squeeze(mean(data_epoch(k,stims==s,:),2));
            [~,~,reg_R] = regress(x,av_erp);
            data_epoch(k,m,:)=reg_R;
        end
    end
end

% calculate spectra
[~,f]=pwelch(squeeze(data_epoch(1,1,fft_t)),fft_w,fft_ov,srate,srate);
data_epoch_spectra=zeros(size(data_epoch,1),size(data_epoch,2),length(f));
clear Pxx

% calculate powerspectra
for k=1:size(data_epoch,1)%channels
    disp(['el ' int2str(k)])
    for m=1:size(data_epoch,2)%epochs
        x=squeeze(data_epoch(k,m,:));
        [Pxx,f]=pwelch(x(fft_t),fft_w,fft_ov,srate,srate);
        data_epoch_spectra(k,m,:)=Pxx;
    end
end