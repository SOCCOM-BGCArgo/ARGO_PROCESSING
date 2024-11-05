function [OUT] = screen_bad_adjusted(IN)

%--------------------------------------------------------------------------
% Small helper function to perform a final assessment on adjusted data
% fields.  If there is raw data marked bad etc due to sensor malfunction,
% this data should NOT be propagated to Argo Adjusted fields.  Therefore,
% this check looks for bad QC flags in raw data, and then ensures the
% associated adjusted & adjusted-qc fields are filled appropriately (with
% fill value).
%
% Tanya Maurer, MBARI
% Feb 13, 2023
%
% INPUTS: MBARI Argo data structure with subfields containing one or
% numerous BGC parameter names. Depending on the platform type (APEX, NAVIS
% or SOLO), this routine will need to be called varying number of times to
% address all parameters.
%
% OUTPUTS: Same structure as the input structure, but with any bad data
% screened out of the associated adjusted fields.
%
% USE EXAMPLES:
%   [LR] = screen_bad_adjusted(LR); (MBARI APEX example, addresses all parameters on an APEX)
%   [HR] = screen_bad_adjusted(HR); (MBARI NAVIS example, addresses all parameters on a NAVIS)
%   [BGC01] = screen_bad_adjusted(BGC01); (MBARI SOLO example, addresses only DOXY)
%   [BGC04] = screen_bad_adjusted(BGC04); (MBARI SOLO example, addresses all OCR parameters)
%
% We want to screen for all potential BGC adjusted parameters containing bad data, including:
%   1.  DOXY
%   2.  PH_IN_SITU_TOTAL
%   3.  NITRATE
%   4.  CHLA
%   5.  BBP700
%   6.  CDOM
%   7.  DOWNWELLING_PAR
%   8.  DOWN_IRRADIANCE380
%   9.  DOWN_IRRADIANCE412
%   10. DOWN_IRRADIANCE490
%--------------------------------------------------------------------------


potential_params = {'DOXY','PH_IN_SITU_TOTAL','NITRATE','CHLA','BBP700','CDOM','DOWNWELLING_PAR','DOWN_IRRADIANCE380','DOWN_IRRADIANCE412','DOWN_IRRADIANCE490'};
FV = 99999;
FVqc = 99;
for ii = 1:length(potential_params)
    if isfield(IN,potential_params{ii})
        tmpQC = [potential_params{ii},'_QC'];
        tmpA = [potential_params{ii},'_ADJUSTED'];
        tmpAQC = [potential_params{ii},'_ADJUSTED_QC'];
        xx = find(IN.(tmpQC)==4);
        if ~isempty(xx)
            disp('hello')
            IN.(tmpA)(xx) = FV;
            IN.(tmpAQC)(xx) = FVqc;
        end
    end
end
OUT = IN;

%end
%--------------------------------------------------------------------------


