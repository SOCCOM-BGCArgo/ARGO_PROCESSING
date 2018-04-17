% lsqfitma.m                                     by:  Edward T Peltzer, MBARI
%                                                revised:  2000 Jan 27.
% 
% M-file to calculate a "MODEL-2" least squares fit.
%
%     The line is fit by MINIMIZING the NORMAL deviates.
%
%     The equation of the line is:     y = mx + b.
%
%     This line is called the MAJOR AXIS.  All points are given EQUAL
%       weight.  The units and range for X and Y must be the same.
%     Equations are from York (1966) Canad. J. Phys. 44: 1079-1086;
%       re-written from Kermack & Haldane (1950) Biometrika 37: 30-41;
%       after a derivation by Pearson (1901) Phil. Mag. V2(6): 559-572.
%
%     Data are input and output as follows:
%
%	    [m,b,r,sm,sb] = lsqfitma(X,Y)
%
%             X    =    x data (vector)
%             Y    =    y data (vector)
%
%             m    =    slope
%             b    =    y-intercept
%             r    =    correlation coefficient
%             sm   =    standard deviation of the slope
%             sb   =    standard deviation of the y-intercept
%
%     Note that the equation passes through the centroid:  (x-mean, y-mean)

function [m,b,r,sm,sb]=lsqfitma(X,Y)

% Determine the size of the vector
 
n = length(X);
 
% Calculate sums and other re-used expressions
 
Sx = sum(X);
Sy = sum(Y);
xbar = Sx/n;
ybar = Sy/n;
U = X - xbar;
V = Y - ybar;
 
Suv = sum(U .* V);
Su2 = sum(U .^2);
Sv2 = sum(V .^2);
 
sigx = sqrt(Su2/(n-1));
sigy = sqrt(Sv2/(n-1));
 
% Calculate m, b, r, sm, and sb
 
m = (Sv2 - Su2 + sqrt(((Sv2 - Su2)^2) + (4 * Suv^2)))/(2 * Suv);
b = ybar - m * xbar;
r = Suv / sqrt(Su2 * Sv2);
 
sm = (m/r) * sqrt((1 - r^2)/n);
sb1 = (sigy - sigx * m)^2;
sb2 = (2 * sigx * sigy) + ((xbar * m * (1 + r))/r^2);
sb = sqrt((sb1 + ((1 - r) * m * sb2))/n);
