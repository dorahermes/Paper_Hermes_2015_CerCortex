function data_epoch = ecog_baselinesubtract(data_epoch,t_base)

% data_epoch = ecog_baselinesubtract(data_epoch,t_base)
% subtracts the mean signal during t_base in each epoch
% data_epoch = electrodes X epoch X t

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

% baseline correct
for k=1:size(data_epoch,1)%channels
    for m=1:size(data_epoch,2)%epochs
        x=squeeze(data_epoch(k,m,:));
        x=x-mean(x(t_base));
        data_epoch(k,m,:)=x;
    end
end
