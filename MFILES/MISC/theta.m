function ptemp = theta(p,t,s,p0)

% THETA  Computes local potential temperature at reference pressure.
%
%  PTEMP = THETA(P,T,S,P0) is the local potential temperature
%       at reference pressure P0 using Bryden 1973 polynomial for
%       adiabatic lapse rate and Runge-Kutta fourth order integration
%       algorithm.
%
%       Units:
%               Pressure        P, P0   dbar
%               Temperature     T       degC IPTS-78
%               Salinity        S       PSU  PSS-78
%	Defaults:
%		P0		0 dbar
%
%       Checkvalue:
%               THETA(10000,40,40,0) = 36.89072
%
%      P could be a [ M by 1 ] Vector if T,S 2-dimensional, 
%           or a [ 1 x 1 x N ] Vector if T,S 3-dimensional,
%                   (the last Matlab 5.# only)

%       18/02/93, C. Mertens, IfM Kiel, changed to Matlab
% added default p0=0dbar	G.Krahmann, IfM Kiel, Mar 1996


if length(size(t)) == 3

 p_si = size(p);
 if p_si(1:2) == [1 1]

  ok=1;
  try
    [hilf1,hilf2,p]=meshgrid(t(1,:,1),t(:,1,1),p);
  catch
    ok=0;
  end
  if ~ok
   error(lasterr)
  end
    clear hilf1 hilf2
 end

elseif length(size(t)) == 2

 if size(p,2) == 1 
    p = p*ones(1,size(t,2)) ;
 end

 if ( size(p,1) == 1 ) & ( size(p0,1) > 1 )
    p = ones(size(p0,1),1) * p;
 end

end

if nargin<4
  p0=0;
end

p = p/10 ; 
p0 = p0/10 ; 
h = p0 - p ;
x = h.*atg(p,t,s) ;
t = t + 0.5*x ;
q = x ;
p = p + 0.5*h ;
x = h.*atg(p,t,s) ;
t = t + 0.29289322*(x - q) ;
q = 0.58578644*x + 0.121320344*q;
x = h.*atg(p,t,s) ;
t = t + 1.707106781*(x - q) ;
q = 3.414213562*x - 4.121320344*q ;
p = p + 0.5*h ;
x = h.*atg(p,t,s) ;
ptemp = t + (x - 2*q)/6 ;

%**************************************************
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

function a = atg(p,t,s)
%ATG Computes adiabatic temperature gradient (required by THETA).
%       A = ATG(P,T,S)

%       VAX 11/750      1983    J.HOLTORFF
%       18/02/93, C. Mertens, IfM Kiel, changed to Matlab

s = s-35.0 ;

a = (((-2.1687E-13*t + 1.8676E-11).*t - 4.6206E-10).*p ...
   + (( 2.7759E-10*t - 1.1351E-08).*s ...
   + ((-5.4481E-12*t + 8.7330E-10).*t - 6.7795E-08).*t + 1.8741E-06)).*p ...
   +  (-4.2393E-07*t + 1.8932E-05).*s ...
   + (( 6.6228E-09*t - 6.8360E-07).*t + 8.5258E-05).*t + 3.5803E-04 ;

