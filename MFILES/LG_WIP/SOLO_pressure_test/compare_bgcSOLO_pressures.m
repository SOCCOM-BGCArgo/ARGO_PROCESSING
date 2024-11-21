% compare_bgcSOLO_pressures.m
%
%--------------------------------------------------------------------------
% Script to start looking at comparison of 'target pressure' versus
% actual sampling pressure on various sensor axes of the BGC-SOLO.  The
% driving question is - how far off do each of the sensors sample from
% their target pressure, and is this negligible (or not) throughout the
% profile?
%
% Junk script started by T. Maurer 7/15/24 for Logan Grady to take over.
% :-)
%--------------------------------------------------------------------------

filedir = 'C:\Users\lgrady\Documents\ALK\';
alkdata = dir([filedir,'*.alk']);
thefiles = char(alkdata.name);
alkdiffs = nan(500);
for ialk = 1:size(thefiles,1)
    filenm = [filedir,thefiles(ialk,:)];
    disp(['Parsing file ',filenm]);
    Dalk = parse_BGCsio_v2(filenm, 1, 'alk');
    Indhdr = strfind(Dalk.hdr,'TARGET PRESSURE');
    iTARGET = find(not(cellfun('isempty',Indhdr)));
    Indhdr = strfind(Dalk.hdr,'PRESSURE (dbar)');
    iPRES = find(not(cellfun('isempty',Indhdr)));
    Tpres = Dalk.data(:,iTARGET);
    pres = Dalk.data(:,iPRES);
    diffpres = Tpres - pres;
    if ialk == 1 % grab target pressure for first column
        alkdiffs(1:length(Tpres),1) = Tpres;
    end
    alkdiffs(1:length(diffpres),ialk+1) = diffpres;
end

% Clean up array:
alkdiffs(all(isnan(alkdiffs),2),:) = [];    %Remove ALL-NaN rows
alkdiffs(:, all(isnan(alkdiffs),1)) = [];   %Remove ALL-NaN columns

