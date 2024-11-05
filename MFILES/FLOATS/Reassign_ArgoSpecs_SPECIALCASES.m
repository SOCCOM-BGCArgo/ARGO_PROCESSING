function [LR, HR, INFO] = Reassign_ArgoSpecs_SPECIALCASES(LR,HR,INFO,FLOATS)
%
% Small helper function called by Process_APEX float, used to reassign
% <PARAM>_ADJUSTED_ERROR fields for floats with special case issues
% (ie Ross Sea floats qc'd to ~500m).  Error inflations assigned on a float-by-float
% basis, based on regional comparisons & analysis.  Updates are user-defined and held
% in "Define_ArgoSpecs_SPECIALCASES".
%
% Tanya Maurer
% MBARI
% 4/18/23
%
% INPUTS: WMOin    -- WMOid
%         CYCin        -- Cycle being processed
%         LR     -- LR data matrix, as defined up to this point in the processing code
%         HR     -- HR data matrix, as defined up to this point in the processing code
%	      INFO   -- INFO matrix; contains Scientific-calib-comments for Argo files
%         FLOATS -- Input structure outlining what needs modification; with the following fields.
%					These definitions should be stored in "Define_ArgoSpecs_SPECIALCASES.m"
%					.WMO: wmo of float
%					.parameter (ie 'NITRATE'; must follow Argo paramter fields)
%					.error_inflate (ie, 0.5 to add a static 0.5 umol/kg to the nitrate adjusted error)
%					.add_comment (ie 'Ross Sea Float; QC assessed at 600m')
% OUTPUTS:PH_ERRout     -- Updated PH_ADJUSTED_ERROR
%         LR     -- LR data matrix, reflecting updated fields as defined by inputs
%         HR     -- HR data matrix, reflecting updated fields as defined by inputs
%		  INFO   -- INFO data matrix, reflecting updated fields as defined by inputs
%--------------------------------------------------------------------------

%See if float is on the list
WMOin = str2num(INFO.WMO_ID);
for i = 1:length(FLOATS)
    tt = FLOATS{i}.WMO;
    if tt == WMOin
        grabind = i;
        MYFLOAT = FLOATS{grabind};
    else
        continue %do nothing; get out of this routine
    end
end
if ~exist('MYFLOAT','var')
    return
end

%Now check which parameter fields need modifying.  If there are more then one, must loop through.
for ii = 1:length(MYFLOAT.parameter)
    %First check if the current cycle requires modification
    if INFO.cast >= MYFLOAT.cyclestart{ii}
        if ~isempty(strfind(MYFLOAT.parameter{ii},'PH'))
            tmpparam = [MYFLOAT.parameter{ii},'_IN_SITU_TOTAL_ADJUSTED_ERROR']; %PH adj error uses full Argo name in our processing, but the sci-cal-comment only uses 'PH'
        else
            tmpparam = [MYFLOAT.parameter{ii},'_ADJUSTED_ERROR'];
        end
        LR.(tmpparam) = LR.(tmpparam)+MYFLOAT.error_inflate{ii};
        if isfield(HR,tmpparam)
           HR.(tmpparam) = HR.(tmpparam)+MYFLOAT.error_inflate{ii}; %for Navis
        end
        tmpcomment = [MYFLOAT.parameter{ii},'_SCI_CAL_COM'];
        INFO.(tmpcomment) = [INFO.(tmpcomment) MYFLOAT.add_comment{ii}];
    end
end




