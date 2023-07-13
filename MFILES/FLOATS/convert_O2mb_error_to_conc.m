function O2error = convert_O2mb_error_to_conc(T,S,O2mb)
% Small helper routie to convert DOXY_ADJUSTED_ERROR from a specification
% of 2mb error in PPOX_DOXY (converted to concentration units using in situ
% profile T and S).
% See Argo Quality Control Document for Dissolved Oxygen: https://archimer.ifremer.fr/doc/00354/46542/82301.pdf
%
% Tanya Maurer
% October, 2021
% MBARI
%--------------------------------------------------------------------------
    TK = T+273.15;
    O2sol = oxy_sol(T,S,0); %flag  = umol/kg, T in deg C for this routine
%     E = exp((0.317.*P)./(8.314.*(TK))); %pretty much = 1! is this needed?
    pH20 = exp(52.57 -(6690.9./TK) - 4.681 .* log(TK));
%     O2 = O2sol.*pO2./E./0.20946./(1013.25-pH20)
    O2tmp = O2sol.*O2mb./0.20946./(1013.25-pH20);
    O2error = O2tmp;
%    END