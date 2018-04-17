function out = sdn2doy(sdn)
% Return day of year given a matlab serial date number array

dvec = datevec(sdn);
out  = sdn - datenum(dvec(:,1),1,1);

clearvars -except out