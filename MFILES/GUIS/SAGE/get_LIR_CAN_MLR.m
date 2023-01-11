function [handles, DATA] = get_LIR_CAN_MLR(dirs,handles,DATA)

% function to calculate reference data for use in SAGE


% ********************************************************************
% GET CANYON NEURAL NETWORK, LINR & LIPHER APROXIMATIONs FOR NO3 AND PH
% NEED QC OXYGEN FOR THIS
% out = CANYON_jp(gtime,lat,lon,pres,temp,psal,doxy,param)
if isfield(handles,'new_qc_data')  %references are only recalculated upon new float selection, or after O2 gain is changed.  Thus, if new_qc_data exists when this funtion is called, O2 gain has been updated.
    O2data = handles.new_qc_data.data(:,DATA.iO);
else
    O2data = handles.qc_data.data(:,DATA.iO);
end

if handles.info.qc_flag == 1
    d = handles.qc_data;
else
    d.hdr = {};
end

% CANYON & CANYON-B ESTIMATES FOR NITRATE AND PH
if handles.info.qc_flag == 1 && ~isempty(DATA.iO)
    %                 set(handles.recumpute_text,'Visible','on')
    %                 set(handles.recumpute_text, ...
    %                     'String','LOADING CANYON NEURAL NETWORK NO3  & PH ESTIMATES ....')
    %                 drawnow
    canyonB_no3 = CANYONB(d.data(:,1),d.data(:,4),d.data(:,3), ...
        d.data(:,6),d.data(:,8),d.data(:,10),O2data,{'NO3'});
    canyonB_ph = CANYONB(d.data(:,1),d.data(:,4),d.data(:,3), ...
        d.data(:,6),d.data(:,8),d.data(:,10),O2data,{'pH'});
    %%%Tanya Maurer, June26,2020.  Tested this adjustment to canyon-b
    %%%briefly on 1 float and was VERY MINIMAL! (ie 0.05 millipH)  This
    %%%was at depth, for 1 float.  We should still apply it, and likely
    %%%add as a different variable to the structure (ie
    %%%canyonB_ph.pHspec).
    % %         % CANYON-B pH is not 'in-line' with pH calculated from Alk and DIC,
    % %         % not pH that was spec measured (LIPHR is consistent with pH spec).
    % %         %  So, for consistency, apply adjustment using eqn 1 in Carter et
    % %         %  al, 2018.
    CBspec = canyonB_ph.pH + (canyonB_ph.pH*0.0404 - 0.3168); % TM 040821, uncommented
    canyonB_ph.pH = CBspec;
    
    
    canyon_no3 = CANYON_jp(dirs,d.data(:,1),d.data(:,4),d.data(:,3), ...
        d.data(:,6),d.data(:,8),d.data(:,10),O2data,'NO3');
    canyon_ph = CANYON_jp(dirs,d.data(:,1),d.data(:,4),d.data(:,3), ...
        d.data(:,6),d.data(:,8),d.data(:,10),O2data,'PH');
    handles.canyon.hdr  = [d.hdr([1,2,6]),'canyon_no3','canyon_ph','canyonB_no3','canyonB_ph'];
    handles.canyon.data = [d.data(:,[1,2,6]), canyon_no3, canyon_ph canyonB_no3.NO3 canyonB_ph.pH];
else
    %                 set(handles.recumpute_text,'Visible','on')
    %                 set(handles.recumpute_text, ...
    %                     'String','NO CANYON NEURAL NETWORK NO3  OR PH ESTIMATES!!')
    %                 drawnow
    handles.canyon =[];
end

% GET CANYON NEURAL NETWORK ESTIMATES
DATA.C = handles.canyon;
if ~isempty(DATA.C)
    DATA.iCN   = find(strcmp('canyon_no3',DATA.C.hdr)     == 1);
    DATA.iCPH  = find(strcmp('canyon_ph',DATA.C.hdr)      == 1);
    DATA.iCBN   = find(strcmp('canyonB_no3',DATA.C.hdr)     == 1);
    DATA.iCBPH  = find(strcmp('canyonB_ph',DATA.C.hdr)      == 1);
else
    DATA.iCN   = [];
    DATA.iCPH  = [];
    DATA.iCBN   = [];
    DATA.iCBPH  = [];
end


% LINR & LIPHR ESTIMATES FOR NITRATE AND PH
if handles.info.qc_flag == 1
    % FIRST, always comput LIR without O2 (T,S,position only)
    LXXX_pos = [d.data(:,3), d.data(:,4), d.data(:,DATA.iZ)]; % lon,lat,Z
    
    MeasIDVec    = [1 7]; % PSAL, TEMP, ,
    Measurements = [d.data(:,10), d.data(:,8)];
    
    Equations    = 8; % S, Theta
    
    [NO3_Est_noO2, N_Uncert_Est, MinUncert_Equ] = LINR(LXXX_pos, ...
        Measurements, MeasIDVec,'Equations', Equations); % update 08/11/17
    
    [PH_Est_noO2, PH_Uncert_Est, MinUncert_Equ] = LIPHR(LXXX_pos, ...
        Measurements, MeasIDVec, ...
        'Equations', Equations,'OAAdjustTF',false); % update 08/11/17
    
    %
    handles.LIRnoO2.hdr  = [d.hdr([1,2,6]),'LIRnoO2_no3','LIRnoO2_ph'];
    handles.LIRnoO2.data = [d.data(:,[1,2,6]), NO3_Est_noO2, PH_Est_noO2];
    % SECOND, if oxygen data, compute LIR eqn 7 with O2 input
    if ~isempty(DATA.iO) %LIR with eqn7 (default) requires oxygen!
        
        %         ptemp    = theta(d.data(:,6), d.data(:,8) ,d.data(:,10),0);
        
        %         MeasIDVec    = [1 2 6]; % PSAL, Pot_TEMP, DOXY_ADJ,
        %         Measurements = [d.data(:,10), ptemp, d.data(:,iO)];
        
        MeasIDVec    = [1 6 7]; % PSAL, DOXY_ADJ, TEMP, ,
        Measurements = [d.data(:,10), O2data, d.data(:,8)];
        
        Equations    = 7; % S, Theta, AOU
        
        [NO3_Est, N_Uncert_Est, MinUncert_Equ] = LINR(LXXX_pos, ...
            Measurements, MeasIDVec,'Equations', Equations); % update 08/11/17
        
        [PH_Est, PH_Uncert_Est, MinUncert_Equ] = LIPHR(LXXX_pos, ...
            Measurements, MeasIDVec, ...
            'Equations', Equations,'OAAdjustTF',false); % update 08/11/17
        
        %
        handles.LIR.hdr  = [d.hdr([1,2,6]),'LIR_no3','LIR_ph'];
        %handles.LIR.data = [d.data(:,[1,2,6]), NO3_Est, PH_Est, N_Uncert_Est, PH_Uncert_Est];
        handles.LIR.data = [d.data(:,[1,2,6]), NO3_Est, PH_Est];
        %LIROUTPUTwUncrt = [d.data(:,[1,2,6]), NO3_Est, PH_Est, N_Uncert_Est, PH_Uncert_Est];
        %save('9631test.mat','LIROUTPUTwUncrt');
        
    else
        handles.LIR =[];
    end
end

clear canyon_no3 canyon_ph LXXX_pos  ptemp  MeasIDVec
clear Measurements Equations NO3_Est PH_Est Uncert_Est MinUncert_Equ

% GET LINR & LIPHER NO3 & PH ESTIMATES
DATA.L = handles.LIR ;
if ~isempty(DATA.L)
    DATA.iLN   = find(strcmp('LIR_no3',DATA.L.hdr)     == 1);
    DATA.iLPH  = find(strcmp('LIR_ph',DATA.L.hdr)      == 1);
    DATA.reftemp = DATA.L.data;
    DATA.reftag = 'LIR'; %will be the default on SELECT FLOAT
else
    DATA.iLN   = [];
    DATA.iLPH  = [];
    DATA.reftag = 'NOQC'; %No QC has been done.
end

DATA.LnoO2 = handles.LIRnoO2 ;
if ~isempty(DATA.LnoO2)
    DATA.iLN_noO2   = find(strcmp('LIRnoO2_no3',DATA.LnoO2.hdr)     == 1);
    DATA.iLPH_noO2  = find(strcmp('LIRnoO2_ph',DATA.LnoO2.hdr)      == 1);
else
    DATA.iLN_noO2   = [];
    DATA.iLPH_noO2  = [];
end

% LOAD MLR coeffs and calculate reference fields

DATA.MLR = LoadGuiMLR_GLT;
%PH
MLR = DATA.MLR.Williams_50Sto80S.PH;
tMLR = isnan(O2data) | isnan(handles.qc_data.data(:,DATA.iS));
potT = theta(handles.qc_data.data(:,DATA.iP), handles.qc_data.data(:,DATA.iT), handles.qc_data.data(:,DATA.iS),0);
sig_theta  = density(handles.qc_data.data(:,DATA.iS), potT)-1000; %density at p =0 t= pot temp
DATA.MLRdata.PH.Williams_50Sto80S = MLR.cC + O2data*MLR.cO + handles.qc_data.data(:,DATA.iS)*MLR.cS + ...
    handles.qc_data.data(:,DATA.iT)*MLR.cT + sig_theta*MLR.cST + handles.qc_data.data(:,DATA.iP)*MLR.cP;
DATA.MLRdata.PH.Williams_50Sto80S(tMLR) = NaN;
MLR = DATA.MLR.Williams_30Sto50S.PH;
tMLR = isnan(O2data) | isnan(handles.qc_data.data(:,DATA.iS));
potT = theta(handles.qc_data.data(:,DATA.iP), handles.qc_data.data(:,DATA.iT), handles.qc_data.data(:,DATA.iS),0);
sig_theta  = density(handles.qc_data.data(:,DATA.iS), potT)-1000; %density at p =0 t= pot temp
DATA.MLRdata.PH.Williams_30Sto50S = MLR.cC + O2data*MLR.cO + handles.qc_data.data(:,DATA.iS)*MLR.cS + ...
    handles.qc_data.data(:,DATA.iT)*MLR.cT + sig_theta*MLR.cST + handles.qc_data.data(:,DATA.iP)*MLR.cP;
DATA.MLRdata.PH.Williams_30Sto50S(tMLR) = NaN;
%NO3
MLR = DATA.MLR.Williams_50Sto80S.NO3;
tMLR = isnan(O2data) | isnan(handles.qc_data.data(:,DATA.iS));
potT = theta(handles.qc_data.data(:,DATA.iP), handles.qc_data.data(:,DATA.iT), handles.qc_data.data(:,DATA.iS),0);
sig_theta  = density(handles.qc_data.data(:,DATA.iS), potT)-1000; %density at p =0 t= pot temp
DATA.MLRdata.NO3.Williams_50Sto80S = MLR.cC + O2data*MLR.cO + handles.qc_data.data(:,DATA.iS)*MLR.cS + ...
    handles.qc_data.data(:,DATA.iT)*MLR.cT + sig_theta*MLR.cST + handles.qc_data.data(:,DATA.iP)*MLR.cP;
DATA.MLRdata.NO3.Williams_50Sto80S(tMLR) = NaN;
MLR = DATA.MLR.Williams_30Sto50S.NO3;
tMLR = isnan(O2data) | isnan(handles.qc_data.data(:,DATA.iS));
potT = theta(handles.qc_data.data(:,DATA.iP), handles.qc_data.data(:,DATA.iT), handles.qc_data.data(:,DATA.iS),0);
sig_theta  = density(handles.qc_data.data(:,DATA.iS), potT)-1000; %density at p =0 t= pot temp
DATA.MLRdata.NO3.Williams_30Sto50S = MLR.cC + O2data*MLR.cO + handles.qc_data.data(:,DATA.iS)*MLR.cS + ...
    handles.qc_data.data(:,DATA.iT)*MLR.cT + sig_theta*MLR.cST + handles.qc_data.data(:,DATA.iP)*MLR.cP;
DATA.MLRdata.NO3.Williams_30Sto50S(tMLR) = NaN;
%S, T, and O2 - all empty
DATA.MLRdata.S.Williams_50Sto80S = handles.qc_data.data(:,DATA.iP) * NaN; %replace with nans
DATA.MLRdata.T.Williams_50Sto80S = handles.qc_data.data(:,DATA.iP) * NaN; %replace with nans
DATA.MLRdata.O2.Williams_50Sto80S = handles.qc_data.data(:,DATA.iP) * NaN; %replace with nans
DATA.MLRdata.S.Williams_30Sto50S = handles.qc_data.data(:,DATA.iP) * NaN; %replace with nans
DATA.MLRdata.T.Williams_30Sto50S = handles.qc_data.data(:,DATA.iP) * NaN; %replace with nans
DATA.MLRdata.O2.Williams_30Sto50S = handles.qc_data.data(:,DATA.iP) * NaN; %replace with nans
