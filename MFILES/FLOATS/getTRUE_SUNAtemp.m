function [spec, CTD_for_no3, CTDmn_for_no3] = getTRUE_SUNAtemp(spec,hiresT,hiresP,hiresS,prevORmean)
%--------------------------------------------------------------------------
% getTRUE_SUNAtemp.m
%
% Function to get a better temp match to suna data, as there is a lag in
% temperature between suna and ctd.  Two options are returned: the CTD T
% record closest in pressure (CTDT_no3) or the mean of the closest and the
% one before it (CTDmnT_no3).  spec.T is replaced with one or the other (as
% defined by user in "prevORmean"), and spec returned as output to the
% function.
%
% inputs:
%  spec = spec from .isus file
%  hiresT  = hi-res temp (already processed, HR.TEMP)
%  hiresP  = hi-res pres (already processed, HR.PRES)
%  hiresS  = hi-res salinity (already processed, HR.PSAL)
%  prevORmean = use previous CTDT record, or mean of 2 previous?
%               0 == previous (CTDT_no3)
%               1 == mean of 2 previous (CTDmnT_no3)
%
% T. Maurer
% MBARI
% 08/22/2018
%--------------------------------------------------------------------------

% % %testing:
% Nmsg = 'Y:\floats\n0949\0949.005.isus';
% spec = parse_NO3msg(Nmsg);
% load X:\ARGO_PROCESSING\DATA\FLOATS\NO_WMO_0949\NO_WMO_0949.005.mat

HRT = hiresT;
HRP = hiresP;
HRS = hiresS;
sT = spec.T;
sP = spec.P;
sS = spec.S;
sP = flipud(sP); %change to ascneding to match HR vector
sT = flipud(sT); %change to ascneding to match HR vector
sS = flipud(sS); %change to ascneding to match HR vector

SP = sP;
ST = sT;
SS = sS;

for i = 1:length(SP)
    tmpP = SP(i);
    [c,inds] = nanmin(abs(HRP-tmpP));
	if c <=2 && i~=length(SP) && inds~=length(HRT)  %This eliminates imposing the n+1 indexing on the coarse resolution data at depth (n+1 would offset the Temp much more than 2m, this fix is not applicable at depth)
        %Temp
        CTD_for_no3(i,1) = HRT(inds+1); %1 record deeper than the matching CTD (isus temp is warmer due to lag)
        CTDmn_for_no3(i,1) = nanmean(HRT(inds+1:inds+2)); %try a mean (isus temp is warmer due to lag)
        %Sal
        CTD_for_no3(i,2) = HRS(inds+1); %1 record deeper than the matching CTD 
        CTDmn_for_no3(i,2) = nanmean(HRS(inds+1:inds+2)); %try a mean 
    else
        %Temp
        CTD_for_no3(i,1) = ST(i); %stick with T in isus file for the deep (coarse res) data
        CTDmn_for_no3(i,1) = ST(i);
        %Sal
        CTD_for_no3(i,2) = SS(i); %stick with S in isus file for the deep (coarse res) data
        CTDmn_for_no3(i,2) = SS(i);
    end
end

if prevORmean == 0
    spec.T = flipud(CTD_for_no3(:,1)); %reassign spec.T, return to descending order
    spec.S = flipud(CTD_for_no3(:,2)); %reassign spec.S, return to descending order
else
    spec.T = flipud(CTDmn_for_no3(:,1)); %reassign spec.T, return to descending order
    spec.S = flipud(CTDmn_for_no3(:,2)); %reassign spec.S, return to descending order
end

end % end function
