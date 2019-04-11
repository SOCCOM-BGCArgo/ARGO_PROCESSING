function T = test0slope(X,Y,alpha)
% Simple function to test if the regression slope from the data pair is
% significantly differnt from zero (Yi = mX + B via Model I regression)
% m = slope = dy/dx & B = Y intercept
%
% INPUTS:
%   X     = independant variable data
%   Y     = dependant variable data
%   alpha = signifigance level, often 0.1 (90%), 0.05 (95%), or 0.01 (99%)
%
% OUPUT: a structur T
%   T.test  = test outcome 1 = different from zero, 0 = not
%   T.X     = independant variable data
%   T.Y     = dependant variable data
%   T.alpha = signifigance level;
%   T.N     = number of samples
%   T.DF    = Degrees of freedom (N-2)
%   T.my    = regression slope, dY/dX;
%   T.by    = regression Y intercept;
%   T.t     = t_stat for slope;
%   T.p     = p value;
%

% TEST DATA
% alpha = 0.05;
% X = data(tZchl,iP);
% Y = CHLA(tZchl);
% X = 1:10;
% Y = (1:10) + (rand(10,1)*10)';
% rand_data = [0.059, 0.681, 0.042, 0.071, 0.521, 0.096, 0.818, 0.817, ...
%     0.722, 0.149];
% Y = Y + rand_data*10;  

% ************************************************************************
% TEST FOR PROPER INPUTS
sX = size(X); sY = size(Y);
if max(sX) ~= max(sY) % CHECK ARRAY SIZES
    fprintf(['INPUT ERROR: X & Y must be the same size: ',...
        'X[%0.0fx%0.0f] & Y[%0.0fx%0.0f]\n'],sX,sY);
    return
end
if sX(1) ~= sY(1) || sX(2) ~= sY(2) % CHECK FOR EQUAL DIMMENSIONS
    fprintf(['INPUT ERROR: X & Y must have the same dimensions: ',...
        'X[%0.0fx%0.0f] & Y[%0.0fx%0.0f]\n'],sX,sY);
    return
end

if sX(2) > sX(1) % IF ROW VECTORS MAKE COLUMN VECTORS
    X = X'; % make column vectors
    Y = Y';
end

tnan = isnan(X) | isnan(Y); % REMOVE NaN LINES IF THEY EXIST
if sum(tnan) > 0
    X(tnan) = [];
    Y(tnan) = [];
    fprintf('%0.0f NaN data rows removed before testing\n',sum(tnan));
end
clear sX sY

% QUICK Y DISTRIBUTION CHECK
Ychk = unique(Y);
if max(size(Ychk)) < 3 % Not enough unique data points
    disp('Number of unique data points < 3 - exiting')
    T = [];
    return
end

% ************************************************************************
% NOW DO SOME PREP WORK

% CREATE ANONAYMOUS FUNCTIONS TO CALCULATE T DISTRIBUTION PROBABILITY FUNCTION
%https://www.mathworks.com/matlabcentral/answers/163578-obtaining-the-p-value-from-the-t-stat
tdist2T = @(t,v) (1-betainc(v/(v+t^2),v/2,0.5)); % 2-tailed t-distribution
tdist1T = @(t,v) 1-(1-tdist2T(t,v))/2;           % 1-tailed t-distribution

null_slope = 0; % SLOPE TO TEST AGAINST

[my,by,ry,smy,sby] = lsqfity(X,Y); % MODEL I LINER REGRESSION

% CALCULATE STANDARD ERROR OF SLOPE
N  = size(X,1); % sample count
DF = N-2; % for Y = mX + b (two constants)
Yi = my* X + by; % Y estimate from regression
Xm = mean(X); % mean of X

top = sum((Y-Yi).^2) ./ DF;  % numerator
bot = sum((X-Xm).^2);        % denominator
SE = sqrt(top) ./ sqrt(bot); % standard error of slope

% CALCULATE T STATISTIC
t_stat = (my - null_slope) ./ SE;

% CALCULATE P VALUE
P2T = 1-tdist2T(t_stat,DF); % two tailed +/- about zero (use this)
P1T = 1-tdist1T(t_stat,DF);

% FILL OUT OUTPUT STRUCTURE
T.test  = 0;
T.X     = X;
T.Y     = Y;
T.alpha = alpha;
T.N     = N;
T.DF    = DF;
T.my    = my;
T.by    = by;
T.t     = t_stat;
T.p     = P2T;

if P2T < alpha
    T.test = 1;
end
    

% plot(X,Y,'bo',X,Yi,'R-')
% title(sprintf('alpha = %0.3f  Tstat = %0.2f  p = %0.4f', alpha, t_stat, P2T))
% hold on
% y_lim =ylim;
% plot(xlim, xlim*0+y_lim(1)+diff(ylim)/2,'k--')
% hold off

clearvars -except T


