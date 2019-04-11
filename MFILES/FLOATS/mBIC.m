function bic = mBIC(residuals, k, errorLim)
% This function calculates the Bayesian Information Criterion 
%
% INPUTS:
%   residuals - an array of data - model residuals 
%   k         - number of fitting parameters
%   errorLim   - cap on residuals
R = nansum((residuals .* residuals)); % NO WEIGHTS!
N = max(size(residuals));
% N should be >> k, if not send a warning
warning_msg = ['Too many fitting parameters (k) for sample size (N)!! ',...
               ' - returning NaN'];
           
% Valid data parameters?           
if k*4 > N-1 
    bic = NaN;
    disp(warning_msg)
    return
end
bic = log(1./N*R + errorLim.^2) + k*log(N) / N;

clearvars -except bic

