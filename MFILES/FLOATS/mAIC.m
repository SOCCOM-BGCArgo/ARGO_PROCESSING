function aic = mAIC(residuals, k, type)
% This function calculates several flavors of the modified least squares
% Akaike Information Criterion defined by the "type" string
%
% INPUTS:
%   residuals - an array of data - model residuals 
%   k         - number of fitting parameters
%   type      - 'standard'
%                   AIC = Nlog(R/N) + 2k
%               'small sample'
%                   AIC = Nlog(R/N) + 2k * N/(N-k-1)
%               'owens'
%                   AIC = Nlog(R/N) + 2k * N/(N-k-2)
%
%       N   = number of residual data points
%       log = natural log
%       R   = weighted residual sum of squares
%               (NOTE: equal weights used here!)
%       k   = number of fitting parameters
%               for 'owens' k = 2m+2 where m is the number of break points
%
% 'small sample' is used when N/k < 40 (crteria also probably for Ownens)
% 'owens' was developed to optimize piece wise linear regression fits to
% residual data. For more details see:
%	Owens, B. W. and A. P. S. Wong, 2009, An improved calibration method
%       for the drift of the conductivity sensor on autonomous CTD
%       profiling floats by ?–S climatology, Deep Sea Res. I, 450-457,
%       doi: 10.1016/j.dsr.2008.09.008.
% The other fromulas were from:
%   Burnham, K.P.,Anderson,D.R.,2002. Model Selection and Multimodel
%       Inference: A Practical Information – Theoretic Approach, second ed.
%       Springer, NewYork, USA, p.488.

% residuals = [-0.0002;-0.0004;0.1203;0.112;-0.3038;-0.0824; ...
%     -0.1169;0.0251;0.2177;0.0245];
%     
% k = 8;
% type = 'owens';



R = sum((residuals .* residuals)); % NO WEIGHTS!
N = max(size(residuals));

% N should be >> k, if not send a warning
warning_msg = ['Too many fitting parameters (k) for sample size (N)!! ',...
               ' - returning NaN'];
           
% Valid data parameters?           
if k*4 > N-1 
    aic = NaN;
    disp(warning_msg)
    return
end
           
switch type
    case 'standard'
        aic = N .* log(R./N) + 2.*k;
    case 'small sample'
        aic = N .* log(R./N) + 2.*k .* N./(N-k-1);
    case 'owens' % really Jones and Day 1995
        %aic = N .* log(R./N) + (N.*(N+k))./(N-k-2);% In paper
        aic = log(R./N) + (N+k)./(N-k-2); % OK to drop N from both terms
end

clearvars -except aic

