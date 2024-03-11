function d = interp_SOLO_ST(PTS,BGC_P,offset,tol)
% PURPOSE:
%   This function interpolates primary T & S onto a secondary BGC pressure
%   axis. T & S are need at the BGC sample pressures to properly calculate
%   oxygen, nitrate, pH or BBP
%
% INPUTS:
%   PTS    - an nx6 matrix of primary axis data:[P, P_QC, T, T_QC, S, S_QC]
%   BGC_P  - secondary BGC pressure axis
%   offset - sensor offset from CTD intake (ie SUNA is attached to the
%            the lower part of float usually with an offset ~ 0.8m
%   tol    - pressure match up tolerance in meters. The BGC pressure sample
%            must be with in X meters of a primary pressure measurement
%            where tol = X
%
% OUTPUTS:
%   d =   a structure with 6 fields, 1 x n, n = size(BGC_P,1)
%   d.BGC_P    = secondary pressure axis
%   d.TEMPi    = Primary interpolated temperatures
%   d.TEMPi_QC = Primary interpolated temperature flags fixed to the
%                higher value. Samples which fail tolerance test also
%                marked bad
%   d.PSALi    = Primary interpolated salinity
%   d.PSALi_QC = Primary interpolated temperature flags fixed to the
%                higher bounding value. Samples which fail tolerance test
%                also marked bad
%
% TM 4/20/22 fixed indexing bug for finding and removing missing T,S values

% ************************************************************************
% TESTING
% P    = 1000:-2:1;
% P_QC = P*0+1;
% T    = linspace(7,23,size(P,2));
% T_QC = T*0+1;
% S    = linspace(35.5,32,size(P,2));
% S_QC = T*0+1;
% PTS  = [P; P_QC; T; T_QC; S; S_QC]';
% 
% BGC_P  = (986.2: -5.2:3)';
% offset = 0.8; % for SUNA
% tol    = 2.0;
% 
% % MAKE SOME BAD PRIMARY DATA
% tbad = PTS(:,1) < 700 & PTS(:,1) > 500;
% PTS(tbad,[4,6]) = 4;


% ************************************************************************
%                         DO SOME PREP WORK
% ************************************************************************
% BUILD EMPTY OUTPUT STRUTURE
% d.BGC_P    = [];
% d.TEMPi    = [];
% d.TEMPi_QC = [];
% d.PSALi    = [];
% d.PSALi_QC = [];
fv.bio = 99999;
fv.qc  = 99;
fill0 =  nan(size(BGC_P,1),1);
% keyboard
d.BGC_P    = fill0;
d.TEMPi    = fill0;
d.TEMPi_QC = fill0;
d.PSALi    = fill0;
d.PSALi_QC = fill0;

try
% INPUT SHAPE & ORDER IS IMPORTANT!! INDICES ARE HARD WIRED
iP = 1; 
iT = 3;
iS = 5;

% USE SENSOR OFFSET TO ADJUST BGC PRES FOR PROPER WATER COLUMN ALIGNMENT
% WITH PRIMARY PRES - % USE THIS P AXIS FOR WORK BUT RETURN ORIGINAL
BGC_PO = BGC_P;
tnan          = isnan(BGC_P);
BGC_PO(~tnan) = BGC_P(~tnan) + offset; 


% FLAG, REMOVE & UNIQUE BAD OR MISSING PRIMARY T & S
%badT     = isnan(PTS(:,iP)) | isnan(PTS(:,iT)) | PTS(:,iT+1) == 4; % REMOVE NaN's & DATA
% badT     = isnan(PTS(:,iP)); % % ONLY REMOVE NaN's
badT     = isnan(PTS(:,iT)); % % ONLY REMOVE NaN's; TM 4/21/22 - this should be temp index (need to check iP too?)
PT       = PTS(~badT,[iP,iT,iT+1]);
[~,ia,~] = unique(PT(:,1));
PT       = PT(ia,:);

%badS     = isnan(PTS(:,iP)) | isnan(PTS(:,iS)) | PTS(:,iS+1) == 4;
% badS     = isnan(PTS(:,iP));
badS     = isnan(PTS(:,iS));
PS       = PTS(~badS,[iP,iS,iS+1]);
[~,ia,~] = unique(PS(:,1));
PS       = PS(ia,:);

% MAKE LOGICAL TOLERANCE ARRAYS
tTOL_T = BGC_PO*0+1; % START WITH ALL GOOD
tTOL_S = tTOL_T;
for ct = 1: size(BGC_P,1)
    tTOL_T(ct) = min(abs(PT(:,1) - BGC_PO(ct))) <= tol;
    tTOL_S(ct) = min(abs(PT(:,1) - BGC_PO(ct))) <= tol;
end

% INTERPOLATE & FILL OUTPUT

d.BGC_P = BGC_P;
d.TEMPi = interp1(PT(:,1),PT(:,2),BGC_PO);
d.PSALi = interp1(PS(:,1),PS(:,2),BGC_PO);

% TAKE WORST QC IN INTERP BOUNDS
TQC1       = interp1(PS(:,1),PS(:,3),BGC_PO,'next');
TQC2       = interp1(PS(:,1),PS(:,3),BGC_PO,'previous');
d.TEMPi_QC = max([TQC1,TQC2],[],2);
d.TEMPi_QC(~tTOL_S) = 4;

SQC1       = interp1(PS(:,1),PS(:,3),BGC_PO,'next');
SQC2       = interp1(PS(:,1),PS(:,3),BGC_PO,'previous');
d.PSALi_QC = max([SQC1,SQC2],[],2);
d.PSALi_QC(~tTOL_S) = 4;
catch
    disp('ERROR IN INTERPOLATING PSAL AND TEMP TO BGC AXES!!!!  PSALi AND TEMPi REMAIN FILL VALUE!!!')
end
%clearvars -except d





