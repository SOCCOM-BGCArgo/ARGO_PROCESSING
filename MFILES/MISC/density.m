% density.m                                      by:  Edward T Peltzer, MBARI
%                                                  revised:  29 Jan 98.
%
% CALCULATE THE DENSITY OF SEAWATER AT A GIVEN S & T, P = 1 ATM
% Equation of State is from Millero & Poisson (1981) DSR V28: 625-629.
%
% [rho]=density(S,T)
%
% INPUT:       	Salinity (S) in g/kg or pss.
%		Temperature (T) in degrees C.
%
% OUTPUT:	 %% Density [rho] in g/cc: %%
%           Density [rho] in kg/m^3 -jp 12/2009
%			rho = density(S,T).

function [rho]=density(S,T)

% DEFINE CONSTANTS FOR EQUATION OF STATE

  R0=+9.99842594E2;
  R1=+6.793952E-2;
  R2=-9.095290E-3;
  R3=+1.001685E-4;
  R4=-1.120083E-6;
  R5=+6.536332E-9;

  A0=+8.24493E-1;
  A1=-4.0899E-3;
  A2=+7.6438E-5;
  A3=-8.2467E-7;
  A4=+5.3875E-9;

  B0=-5.72466E-3;
  B1=+1.0227E-4;
  B2=-1.6546E-6;

  C=+4.8314E-4;

% CALCULATE RHO

    RHO0=R0+T.*(R1+T.*(R2+T.*(R3+T.*(R4+T.*R5))));

    A=A0+T.*(A1+T.*(A2+T.*(A3+T.*A4)));
    B=B0+T.*(B1+T.*B2);
    RHO=RHO0+S.*(A+B.*sqrt(S)+C.*S); 

% CONVERT KG/M3 TO g/cc

 %rho = RHO ./ 1000; g/cc
 rho = RHO; % kg/m^3
