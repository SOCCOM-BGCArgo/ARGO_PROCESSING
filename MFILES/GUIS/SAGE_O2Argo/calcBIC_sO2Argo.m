function BIC = calcBIC_sO2Argo(gui,DATA)

% ************************************************************************
% calcBIC_sO2Argo.m
% ************************************************************************
%
% Function to calculate Bayesian Information Criteria (BIC) as a metric for
% assessing the appropriate number of breakpoints in the calculation of
% drifts in O2 gains.
%
%
% USE AS: BIC = calcBIC_sO2(gui,DATA)
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
% DATE: 03/20/2019
% UPDATES:
% NOTES: 
% ************************************************************************
%
% ************************************************************************

% Calculate Bayesian Information Criteria (BIC)

% CALCULATE BIC
errorLim = 0; %Cap on residuals.  useful for noisy data in NO3 and pH series.  Set to 0 for O2.
if isfield(DATA,'brkRES')
    SSE = sum(DATA.brkRES.^2);
    n = length(DATA.brkRES);
    m = length(gui.tbl.Data(:,1))-1; %Do not include cycle 1 in # of brkpts
    K = 2.*m+2; 
    if m > ((n/4)-1) % Valid data parameters?  see Jones & Day, 1995
        BIC = NaN;
        disp('n~>>K, cannot calculate BIC.  Setting BIC = NaN')
    else
        BIC = log(1./n*SSE + errorLim.^2) + K*log(n) / n;
%         AIC = log(SSE./n) + (n+K)./(n-K-2); %from Jones & Day, 1995, later referenced by Owens & Wong, 2009.
    end
else
    disp('Not enough data, setting BIC = NaN.')
    BIC = NaN;
end

end
