function [PH_ERRout, PH_COMMENTout, N_ERRout, N_COMMENTout] = Reassign_ArgoSpecs_LIReqn8(FLOATIDin,CYCin,PH_ERRin,PH_COMMENTin,N_ERRin,N_COMMENTin)
%
% Small helper function called by Process_APEX float, used to reassign 
% <PARAM>_ADJUSTED_ERROR fields for floats with failed optodes that have
% been QC'd using Carter et al 2017 LI(PH/N)R Eqn 8 (without O2,
% T/S/LOCATION inputs only).  Error inflations assigned on a float-by-float
% basis, based on regional comparison of difference between LIR eqn 7 & 8
% using good float data (see analysis emailed by Tanya on 8/25/21) and DMQC
% assessment using SAGE.
%
% Tanya Maurer
% MBARI
% 8/25/21
%
% INPUTS: FLOATIDin    -- Internal MBARI FloatID
%         CYCin        -- Cycle being processed
%         PH_ERRin     -- Original PH_ADJUSTED_ERROR assigned in Process_APEX_float
%         PH_COMMENTin -- Original SCI_CAL_COMMENT for pH, assigned in Process_APEX_float
%         N_ERRin      -- Original NITRATE_ADJUSTED_ERROR assigned in Process_APEX_float
%         N_COMMENTin  -- Original SCI_CAL_COMMENT for nitrate, assigned in Process_APEX_float
% OUTPUTS:PH_ERRout     -- Updated PH_ADJUSTED_ERROR 
%         PH_COMMENTout -- Updated SCI_CAL_COMMENT for pH
%         N_ERRout      -- Updated NITRATE_ADJUSTED_ERROR
%         N_COMMENTout  -- Updated SCI_CAL_COMMENT for nitrate
%--------------------------------------------------------------------------

floatSTR = regexp(FLOATIDin,'\d*','match');
floatNUM = str2num(floatSTR{:});

%-------------------------------------------------------------------------
% Column headers for variable, "NO_O2_errors", are:
% Col 1: float ID
% Col 2: WMO ID
% Col 3: Cycle to start inflating Error for Nitrate
% Col 4: Magnitude of error inflation (static) for Nitrate (umol/kg)
% Col 5: Cycle to start inflating Error for pH
% Col 6: Magnitude of error inflation (static) for pH
NO_O2_errors(1,:) = [9274, 5904655, 1, 0, 14, 0.01];
NO_O2_errors(3,:) = [9642, 5904685, 40, 0.3, 1, 0];
NO_O2_errors(4,:) = [9659, 5904843, 1, 0.3, 1, 0.015];
NO_O2_errors(5,:) = [12551, 5904858, 1, 0.2, 1, 0.005];
NO_O2_errors(6,:) = [12573, 5904982, 62, 0.2, 62, 0.005];
NO_O2_errors(7,:) = [18082, 5906205, 1, 0.25, 1, 0.005];
NO_O2_errors(8,:) = [19327, 5906313, 4, 1, 4, 0.025];
NO_O2_errors(9,:) = [19412, 5906300, 1, 0, 14, 0.005];
NO_O2_errors(10,:) = [1114, 5906306, 8, 0.5, 1, 0];
NO_O2_errors(11,:) = [19142, 5906342, 1, 0.15, 1, 0.003];
NO_O2_errors(12,:) = [1200, 4903273, 1, 0.1, 1, 0.001]; %5/9/22; small error inflation; algorithm seems to do ok compared to LIR with O2.

%-------------------------------------------------------------------------

flt_ind = find(NO_O2_errors(:,1)==floatNUM);

%NITRATE
if CYCin>=NO_O2_errors(flt_ind,3) 
    %cycle uncertainty is affected; update error & comment fields
    N_ERRout = N_ERRin + NO_O2_errors(flt_ind,4);
    N_COMMENTout = [N_COMMENTin,' LINR ref EQ8 used starting cycle ',num2str(NO_O2_errors(flt_ind,3)),' due to optode failure.'];
else
    N_ERRout = N_ERRin;
    N_COMMENTout = N_COMMENTin;
end

%PH_IN_SITU_TOTAL
if CYCin>=NO_O2_errors(flt_ind,5) 
    %cycle uncertainty is affected; update error & comment fields
    PH_ERRout = PH_ERRin + NO_O2_errors(flt_ind,6);
    PH_COMMENTout = [PH_COMMENTin,' LIPHR ref EQ8 used cycle ',num2str(NO_O2_errors(flt_ind,5)),'-on due to optode failure.'];
else
    PH_ERRout = PH_ERRin;
    PH_COMMENTout = PH_COMMENTin;
end

% END

    
