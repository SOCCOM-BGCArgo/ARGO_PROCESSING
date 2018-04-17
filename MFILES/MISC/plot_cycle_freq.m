function plot_cycle_freq(floatID,floatType)
% ************************************************************************
% plot_cycle_freq.m
% ************************************************************************
%
% Function to plot cycle frequency for specified float.
% 
%
% INPUTS:  
%          floatID = MBARI float ID
%          floatTYPE = NAVIS or APEX? ('N', or 'A')
%
% OUTPUTS: plots days since cycle 1 against cycle #
%
%
% AUTHOR: Tanya Maurer
%         Monterey Bay Aquarium Research Institute
%         tmaurer@mbari.org
%
% DATE: 03/22/2017
% UPDATES:
% NOTES: 
% ************************************************************************
% ************************************************************************

dirs = SetupDirs_sO2;
DATA = getall_floatdata_sO2(dirs,floatID,floatType);
dayssince = diff(DATA.track(:,1));
dayssince = [0; dayssince];
H = myfig(0.2,0.05,0.8,0.45);
plot(DATA.track(:,2),dayssince,'o-')
set(gca,'fontsize',12)
title(['FLOAT ',floatID,' cycle frequency'],'fontsize',12)
xlabel('Cycle Number','fontsize',12)
ylabel('Days Since Cycle 1','fontsize',12)
axis tight
