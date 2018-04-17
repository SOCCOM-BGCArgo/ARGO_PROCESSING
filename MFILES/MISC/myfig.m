function H = myfig(b,l,w,h)
%----------------------------------------------------
% myfig.m
%----------------------------------------------------
%
% Brings up figure for plot, with specific size/position, and white background.
%
% USE AS: H = myfig(l,b,w,h)
% EXAMPLES: H = myfig(0.25,0.25,0.5,0.5); or for larger:
%           H = myfig(0.08,0.15,0.6,0.7);
% INPUTS: define fractions of the screensize width and height
%         b = bottom location
%         l = left location
%         w = figure width
%         h = figure height
%
% OUTPUTS: H = handle to figure;
%
% AUTHOR: Tanya Maurer
%         MBARI
% DATE: 05/09/13
% UPDATES: 06/28/13, made into a function for ability to control figure
% size.  Position is set.
% ----------------------------------------------------

scrsz = get(0,'ScreenSize');
H = figure('color','white','Position',[scrsz(3)*l scrsz(4)*b scrsz(3)*w scrsz(4)*h]);