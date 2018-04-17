function AIC = calcAIC_sO2Argo(gui,DATA)

% ************************************************************************
% calcAIC_sO2Argo.m
% ************************************************************************
%
% Function to calculate Akaike Information Criteria (AIC) as a metric for
% assessing the appropriate number of breakpoints in the calculation of
% drifts in O2 gains.
%
%
% USE AS: AIC = calcAIC_sO2Argo(gui,DATA)
%
% INPUTS:
%    gui, inputs, and DATA area all structures of handles and inputs to the
%    SAGE_O2 gui.
%
% OUTPUTS:
%
%
% AUTHOR: Tanya Maurer
%         MBARI
%         tmaurer@mbari.org
%
% DATE: 03/14/2017
% UPDATES:
% NOTES: 
% ************************************************************************
%
% ************************************************************************

% Calculate AIC
% if isfield(DATA,'brkRES')
%     SSE = sum(DATA.brkRES.^2);
%     n = length(DATA.brkRES);
%     K = length(gui.tbl.Data(:,1))+1;
%     AIC = log(1/n.*SSE)+2*2*(K)/n;
% else
%     disp('Not enough data, setting AIC to NaN.')
%     AIC = NaN;
% end

% CALCULATE AIC
% modified formula from above.  -tm 4/10/2017
if isfield(DATA,'brkRES')
    SSE = sum(DATA.brkRES.^2);
    n = length(DATA.brkRES);
    m = length(gui.tbl.Data(:,1))-1; %Do not include cycle 1 in # of brkpts
    K = 2.*m+2; 
    if m > ((n/4)-1) % Valid data parameters?  see Jones & Day, 1995
        AIC = NaN;
        disp('n~>>K, cannot calculate AIC.  Setting AIC = NaN')
    else
        AIC = log(SSE./n) + (n+K)./(n-K-2); %from Jones & Day, 1995, later referenced by Owens & Wong, 2009.
    end
else
    disp('Not enough data, setting AIC = NaN.')
    AIC = NaN;
end

end %end function
