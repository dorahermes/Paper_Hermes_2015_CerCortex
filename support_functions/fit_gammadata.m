function [out_exp,w_pwr,w_gauss,gauss_f,fit_f2] = ...
    fit_gammadata(f,f_use4fit,data_base,data_fit)

% function fits broadband + gaussian
% [out_exp,w_pwr,w_gauss,gauss_f,fit_f2] = ...
%     fit_gammadata(f,f_use4fit,data_base,data_fit);
%
% input:
% f: frequencies
% f_use4fit=[35:57 65:115 126:175 186:200];
% data_base: used to fit exp: (1/f^exp) - enter power (not log)
% data_fit: used to fit weights and gaussian - enter power (not log)
%
% output (exp weight_pwr weight_gauss gamma_freq fit_f2)

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


f_sel=ismember(f,f_use4fit);
x_base=data_base(f_sel);
x_in=data_fit(f_sel);
f_in=f(f_sel);

% fit exponent to base
p=polyfit(log10(f_in),log10(x_base)',1);
out_exp=-p(1);

my_options=optimset('Display','off','Algorithm','trust-region-reflective');
[x]=lsqnonlin(@(x) fit_func3_loglog(x,log10(x_in),log10(f_in'),out_exp),...
    [0 0 log10(50)],[-Inf 0 log10(35)],[Inf Inf log10(80)],...
    my_options);
w_pwr=x(1);
w_gauss=x(2);
gauss_f=x(3);

% fit to data in log-space
fit_f2=w_pwr-out_exp*log10(f) + ...
    w_gauss*.04*sqrt(2*pi)*normpdf(log10(f),gauss_f,.04);
