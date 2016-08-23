function data = ecog_notch(data,srate)

% data = ecog_notch(data,srate)
% notch filter data around 60Hz, 120Hz and 180Hz
% data is time X electrodes


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


% 5th order butterworth notch filter, [low/(samplerate/2) high/(samplerate/2)]
[n1_b, n1_a]=butter(5,2*[59 61]/srate,'stop'); %60hz
[n2_b, n2_a]=butter(5,2*[119 121]/srate,'stop'); %120hz
[n3_b, n3_a]=butter(5,2*[179 181]/srate,'stop'); %180hz
disp('notching out 60 120 180')
for elec=1:size(data,2)
    data(:,elec)=filtfilt(n1_b,n1_a,data(:,elec)); %60
    data(:,elec)=filtfilt(n2_b,n2_a,data(:,elec)); %120
    data(:,elec)=filtfilt(n3_b,n3_a,data(:,elec)); %180
end

