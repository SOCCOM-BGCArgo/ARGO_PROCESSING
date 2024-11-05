% Define_ArgoSpecs_SPECIALCASES.m
%
% Script to generate a matfile holding structures that specify special-case handling of Argo errors for specific floats.
% This script gets called in Process_*_Float.m prior to the calling of Reassign_ArgoSpecs_SPECIALCASES.matfile%
% 
% Tanya Maurer
% MBARI
% 4/18/23
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Ross Sea Floats----------------------------------------------------------------------------------------------------------
FLOATS{1}.WMO = 7900823;
FLOATS{1}.cyclestart = {1,1};
FLOATS{1}.parameter = {'NITRATE','PH'};
FLOATS{1}.error_inflate = {0.5,0.01};
FLOATS{1}.add_comment = {'Ross Sea Shelf Float, QC to 600m.',' Ross Sea Shelf, QC to 600m.'};

FLOATS{2}.WMO = 7900824;
FLOATS{2}.cyclestart = {1,1};
FLOATS{2}.parameter = {'NITRATE','PH'};
FLOATS{2}.error_inflate = {0.5,0.015};
FLOATS{2}.add_comment = {'Ross Sea Shelf Float, QC to 250m.',' Ross Sea Shelf, QC to 250m.'};

FLOATS{3}.WMO = 7900825;
FLOATS{3}.cyclestart = {1,1};
FLOATS{3}.parameter = {'NITRATE','PH'};
FLOATS{3}.error_inflate = {0.5,0.01};
FLOATS{3}.add_comment = {'Ross Sea Shelf Float, QC to 400m.',' Ross Sea Shelf, QC to 400m.'};

FLOATS{4}.WMO = 7900826;
FLOATS{4}.cyclestart = {1,1};
FLOATS{4}.parameter = {'NITRATE','PH'};
FLOATS{4}.error_inflate = {0.5,0.01};
FLOATS{4}.add_comment = {'Ross Sea Shelf Float, QC to 500m.',' Ross Sea Shelf, QC to 500m.'};

FLOATS{5}.WMO = 7900827;
FLOATS{5}.cyclestart = {1,1};
FLOATS{5}.parameter = {'NITRATE','PH'};
FLOATS{5}.error_inflate = {0.5,0.01};
FLOATS{5}.add_comment = {'Ross Sea Shelf Float, QC to 500m.',' Ross Sea Shelf, QC to 500m.'};


% % ARMOR3D PSAL proxy --------------------------------------------------------------------------------------------------------
% 
FLOATS{6}.WMO = 5905988;
FLOATS{6}.cyclestart = {99,99,99,99};
FLOATS{6}.parameter = {'NITRATE','PH','DOXY','BBP700'};
FLOATS{6}.error_inflate = {0.5, 0.01, 0.1, 0};
FLOATS{6}.add_comment = {' PSALproxy for BGC https://doi.org/10.48670/moi-00052', ' PSALproxy used for BGC https://doi.org/10.48670/moi-00052',' PSALproxy used for BGC https://doi.org/10.48670/moi-00052',' PSALproxy used for BGC https://doi.org/10.48670/moi-00052'};

FLOATS{7}.WMO = 5904674;
FLOATS{7}.cyclestart = {24,24,24};
FLOATS{7}.parameter = {'NITRATE','PH','DOXY'};
FLOATS{7}.error_inflate = {0.2, 0.002, 0.2};
FLOATS{7}.add_comment = {' PSALproxy for BGC https://doi.org/10.48670/moi-00052', ' PSALproxy used for BGC https://doi.org/10.48670/moi-00052',' PSALproxy used for BGC https://doi.org/10.48670/moi-00052'};

FLOATS{8}.WMO = 5904186;
FLOATS{8}.cyclestart = {121};
FLOATS{8}.parameter = {'DOXY'};
FLOATS{8}.error_inflate = {0.05};
FLOATS{8}.add_comment = {' PSALproxy for BGC https://doi.org/10.48670/moi-00052'};

FLOATS{9}.WMO = 5904660;
FLOATS{9}.cyclestart = {88,88,88,88};
FLOATS{9}.parameter = {'NITRATE','PH','DOXY','BBP700'};
FLOATS{9}.error_inflate = {0.2, 0.002, 0.3, 0};
FLOATS{9}.add_comment = {' PSALproxy for BGC https://doi.org/10.48670/moi-00052', ' PSALproxy used for BGC https://doi.org/10.48670/moi-00052',' PSALproxy used for BGC https://doi.org/10.48670/moi-00052',' PSALproxy used for BGC https://doi.org/10.48670/moi-00052'};

FLOATS{10}.WMO = 5905106;
FLOATS{10}.cyclestart = {101,101,101};
FLOATS{10}.parameter = {'NITRATE','DOXY','BBP700'};
FLOATS{10}.error_inflate = {0.3, 1, 0};
FLOATS{10}.add_comment = {' PSALproxy used for BGC https://doi.org/10.48670/moi-00052',' PSALproxy used for BGC https://doi.org/10.48670/moi-00052',' PSALproxy used for BGC https://doi.org/10.48670/moi-00052'};

FLOATS{11}.WMO = 5905077;
FLOATS{11}.cyclestart = {49,49,49,49};
FLOATS{11}.parameter = {'NITRATE','PH','DOXY','BBP700'};
FLOATS{11}.error_inflate = {0.1, 0.002, 0.2, 0};
FLOATS{11}.add_comment = {' PSALproxy used for BGC https://doi.org/10.48670/moi-00052',' PSALproxy used for BGC https://doi.org/10.48670/moi-00052',' PSALproxy used for BGC https://doi.org/10.48670/moi-00052',' PSALproxy used for BGC https://doi.org/10.48670/moi-00052'};

FLOATS{12}.WMO = 5905108;
FLOATS{12}.cyclestart = {46,46,46,46};
FLOATS{12}.parameter = {'NITRATE','PH','DOXY','BBP700'};
FLOATS{12}.error_inflate = {0.2, 0.002, 0.2, 0};
FLOATS{12}.add_comment = {' PSALproxy used for BGC https://doi.org/10.48670/moi-00052',' PSALproxy used for BGC https://doi.org/10.48670/moi-00052',' PSALproxy used for BGC https://doi.org/10.48670/moi-00052',' PSALproxy used for BGC https://doi.org/10.48670/moi-00052'};

FLOATS{13}.WMO = 5905379;
FLOATS{13}.cyclestart = {1,1,1,1};
FLOATS{13}.parameter = {'NITRATE','PH','DOXY','BBP700'};
FLOATS{13}.error_inflate = {0.2, 0.002, 0.2, 0};
FLOATS{13}.add_comment = {' PSALproxy used for BGC https://doi.org/10.48670/moi-00052',' PSALproxy used for BGC https://doi.org/10.48670/moi-00052',' PSALproxy used for BGC https://doi.org/10.48670/moi-00052',' PSALproxy used for BGC https://doi.org/10.48670/moi-00052'};

% % pH pump offset --------------------------------------------------------------------------------------------------------
f = 13;
for ipump = 1:length(pH_pumpoffset_980_floats)
	FLOATS{f+ipump}.WMO = pH_pumpoffset_980_floats(ipump);
	FLOATS{f+ipump}.cyclestart = {1};
	FLOATS{f+ipump}.parameter = {'PH'};
	FLOATS{f+ipump}.error_inflate = {0.005};
	FLOATS{f+ipump}.add_comment = {' pH pump-offset; QC to 980m.'};
end

