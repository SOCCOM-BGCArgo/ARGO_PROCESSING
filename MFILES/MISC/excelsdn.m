function out = excelsdn(SDN)
%USE:  Simple function to convert Matlab SDN to Excel SDN because I always
%      forget the offset
%
% function out = excelsdn(SDN)
%
% INPUTS:
%     SDN =  an array or matrix of Matlab serial date numbers
%
% OUTPUT:
%     out = an array or matrix of Excel serial date numbers
%
% EXAMPLE:
%     data = excelsdn(735032.48)   


% CHECK INPUTS
if ~isnumeric(SDN)
    error('Input incorrect - should be a numeric array.');
    return
elseif sum(SDN(:) < 0) > 0
    error('Input incorrect - numeric values should be greater than zero.');
    return
end

offset = 693960;
out = SDN - offset;