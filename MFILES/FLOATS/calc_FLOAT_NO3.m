function NO3 = calc_FLOAT_NO3(spec, ncal,Pcor_flag)
% 
% PURPOSE:
%   Calculate nitrate concentration for each sample in a profile(msg file)
%   for APEX or NAVIS floats. This function calculates nitrate
%   concentration from data supplied as a structure retuned from the
%   function parse_NO3msg.m which includes all needed data.
%   See parse_NO3msg.m for details on output structure
%
% USAGE:
%	NO3 = calc_FLOAT_NO3(spec, ncal,Pcor_flag)
%
% INPUTS:
%	spec      =  a structure built by parse_NO3msg.m
%   ncal      =  a structure built by get_float_cals.m or parseNO3cal.m with NO3 cal data
%   Pcor_flag = 1 to pressure correct ESW, otherwise don't
%
% OUTPUTS:
%   NO3 =   a matrix of nitrate concentrations (umol/L) and other info
%       [SDN, DarkCur,Pres,Temp,Sal,NO3,BL_B,BL_M,RMS ERROR,WL~240,ABS~240]

% HISTORY
%   04/06/2017 Added check for Intensites <= dark current and set to NaN

% *************************************************************************
% TESTING
% cal_file  = ['C:\Users\jplant\Documents\MATLAB\ARGO_PROCESSING\' ...
%                'DATA\CAL\NO3_CAL\6400STNP.cal'];
% ncal      = parseNO3cal(cal_file);
% spec      = parse_NO3msg('C:\temp\6400.054.isus');
% Pcor_flag = 1;

% ************************************************************************
% CONSTANTS
% ************************************************************************
fig_flag = 0; % 1 = show nitrate fit
WL_offset = ncal.WL_offset;  % Optical wavelength offset

% ************************************************************************
% CALCULATE RAW ABSORBANCE 
% CALIBRATION DATA MEASURED FROM ~190 to ~400 nm IN THE LAB. FLOAT CAN'T
% STORE ALL WAVELENGTHS SPECTROPHOTOMETER MEASURES SO ONLY A SUBSET IS
% RETURNED. 
% ************************************************************************
d = spec; clear spec

% MAKE SOME CALIBRATION ADJUSTMENTS BASED ON DATA FLAGS IN CAL STRUCTURE
% 0 based count - pixel registration to wavelength is off by one for some
% instruments. if 217 was pixel 30 it is now pix 31. DOUBLE CHECK THIS!!!
if ncal.pixel_base == 0 
    d.spectra_pix_range = d.spectra_pix_range +1;
    d.pix_fit_win = d.pix_fit_win + 1;
    disp('Pixel registration offset by +1')
end
% CHOOSE REGULAR or SEAWATER DARK CURRENT FOR CALCULATION
DC = d.DC; % default
if ncal.DC_flag == 0
    DC = d.SWDC;
end

% GET PIXEL FIT WINDOW FOR SUNA FLOATS BASED ON WL FIT RANGE
% SUNA FLOATS DON'T HAVE THIS AS A HEADER LINE
if isnan(sum(d.pix_fit_win)) && isfield(d,'WL_fit_win');
    pfit_low = find(ncal.WL >= d.WL_fit_win(1),1,'first');
    pfit_hi  = find(ncal.WL <= d.WL_fit_win(2),1,'last');
    d.pix_fit_win = [pfit_low pfit_hi];
    clear pfit_low pfit_hi
end

%  DO NOT USE PIXEL FIT RANGE IN CAL FILE! USE 1st >= 217 and
%  LAST <= 240 unless cal file has updated fit window

% ADJUST FIT WINDOW IF IT HAS BEEN UPDATED IN THE CAL STRUCTURE
% OTHERWISE USE DEFAULT WINDOW >=217 & <=240
if ~isempty(ncal.min_fit_WL)
    d.pix_fit_win(1) = find(ncal.WL >= ncal.min_fit_WL,1,'first');
else
    d.pix_fit_win(1) = find(ncal.WL >= 217,1,'first');
end

if ~isempty(ncal.max_fit_WL)
    d.pix_fit_win(2)  = find(ncal.WL <= ncal.max_fit_WL,1,'last');
else
    d.pix_fit_win(2)  = find(ncal.WL <= 240 ,1,'last');
end

% SAMPLE PIXEL RANGE INDEX IDENTIFIED FROM SAMPLE LINE
saved_pixels = d.spectra_pix_range(1): d.spectra_pix_range(2); % index

% SIZE OF SAMPLE INTENSITY SPECTRA MATRIX(samples X wavelengths)
[rows,cols]  = size(d.UV_INTEN);

% SUBSET CALIBRATION DATA OVER PIXEL RANGE IDENTIFIED FROM SAMPLE LINE
% THEN TRANSPOSE TO ROWS & REPEAT ROWS SO SAME SIZE AS UV_INTEN
REF  = ones(rows,1) * ncal.Ref(saved_pixels)';  % Reference intensities
WL   = ones(rows,1) * ncal.WL(saved_pixels)';   % Wavelengths
ESW  = ones(rows,1) * ncal.ESW(saved_pixels)';  % Extinction coefs SW

ENO3 = ncal.ENO3(saved_pixels)'; % Extinction coefs NO3 (array)

if isfield(ncal,'EHS') % NOT USED YET BUT CARRY ALONG
    EHS  = ncal.EHS(saved_pixels)';  % Extinction coefs HS  (array)
end

% SUBTRACT DARK CURRENT INTENSITIES
UV_INTEN_SW = d.UV_INTEN - DC * ones(1,cols); 
% Get UV reference intensity for pixel range from cal file

% CHECK FOR UV_INTEN_SW <= 0 (DC >= UV_INTEN_SW) & SET TO NaN
tINTEN_SW = UV_INTEN_SW <= 0;
if sum(tINTEN_SW(:)) > 0
    UV_INTEN_SW(tINTEN_SW) = NaN;
    disp(['UV Intensities <= DC found for ',d.file_name, ...
        '. Setting these low intensities to NaN']);
    %clear tINTEN_SW
end
    
% ABSORBANCE SPECTRA FOR ALL SAMPLES IN PROFILE
ABS_SW = -log10(UV_INTEN_SW ./ REF);

% ************************************************************************
% CALCULATE TEMPERATURE CORRECTED ABSORBANCES
% ************************************************************************
A = 1.1500276; % Coefficients from Sakamoto et al., 2009
B = 0.02840;
C = -0.3101349;
D = 0.001222;

% Build matrices for calculations all size of spectra - avoid some loops 
% calculate corrected absorbances for all samples all at once
ctd_temp  = d.T * ones(1,cols);           % measured sw temperature, C
ctd_sal   = d.S * ones(1,cols);           % measured sw salinity
cal_temp  = ones(rows,cols) * ncal.CalTemp;  % temp isus calibrated at, C


%disp('Using cal temp from cal file not from message file')

% Opt_wl_UV = ones(rows,1) * (isus_cal(:,1))';  % matrix of wavelengths
% E_swa_tcal = ones(rows,1) * (isus_cal(:,2))'; % sw extinction coefs, cal t

% Refer to Computing DAC Nitrate v0.2
% Sakamoto temperature correction "f function" at ctd temperature
F_sw = (A + B .* ctd_temp) .* exp((C + D .* ctd_temp) .* ...     % Eqn. 3
       (WL - WL_offset));

% Sakamoto temperature correction "f function" at calibration temperature
% WL_offset
F_cal = (A + B .* cal_temp) .* exp((C + D .* cal_temp) .* ...    % Eqn. 3
        (WL - WL_offset));

% Calculate in situ sea water extinction coefficients (mostly due to Br)
ESW_in_situ = ESW .* F_sw ./ F_cal;      % Eqn. 2

% Add pressure corrrection if flag is on
if Pcor_flag == 1
    pres_term   = (1- d.P ./ 1000 .* ncal.pres_coef) * ones(1,cols);
    ESW_in_situ = ESW_in_situ .* pres_term;
end

% SW (bromide) absorbances
ABS_Br_tcor      = ESW_in_situ .* ctd_sal; % in situ bromide absorbances

% Calculate predicted absorbance of nitrate + baseline
% The 2nd term here is removing the temperature corrected absorbance
% of seawater from the spectra
ABS_cor = ABS_SW - ABS_Br_tcor  ;    % Eqn. 4

clear ctd_temp cal_temp ctd_sal REF ESW pres_term

% ************************************************************************
% CALCULATE THE NITRATE CONCENTRATION & THE SLOPE AND INTERCEPT OF THE
% BASELINE ABSORBANCE
% This calulation occurs over the "fit window" which is usually a subset
% of the "spectra window" which is a subset of the calibration spectra
% window. The calibration window is the full wavelength range of the
% spectraphotometer
% ************************************************************************

% SUBSET TO FIT WINDOW
t_fit     = saved_pixels >= d.pix_fit_win(1) & ...
            saved_pixels <= d.pix_fit_win(2); % fit wind pixel indicies
      
Fit_WL    = WL(1,t_fit)';
Fit_ENO3  = ENO3(t_fit)';
Ones      = ones(size(Fit_ENO3));

M     = [Fit_ENO3 Ones/100 Fit_WL/1000]; % /1000 numerical stability??
M_INV = pinv(M);
% /100 & /1000 may be in here for numerical stability. Need to /100 and
% /1000 at the end to get back to actual baseline and intercept
colsM = size(M,2);
% Columns may need to be added for more solutions
% # samples x (# fit parameters + QC values)
NO3 = ones(rows, colsM +3)* NaN;

for i = 1:rows
    samp_ABS_cor = (ABS_cor(i,t_fit))'; % NO3 + BL
    
    % NO3 = [ baseline intercept, and slope]
    NO3(i,1:3) = M_INV * samp_ABS_cor;
    NO3(i,2) = NO3(i,2)/100; % baseline intercept
    NO3(i,3) = NO3(i,3)/1000; % baseline slope
    
    ABS_BL  = WL(1,:) .* NO3(i,3) + NO3(i,2); % Baseline absorbance
    ABS_NO3 = ABS_cor(i,:) - ABS_BL; % Calculated nitrate absorbance
    ABS_NO3_EXP = ENO3 .* NO3(i,1);
    FIT_DIF = ABS_cor(i,:) - ABS_BL -ABS_NO3_EXP;
    RMS_ERROR = sqrt((sum(FIT_DIF(t_fit).^2)./sum(t_fit)));
    ind_240  = find(Fit_WL <= 240,1,'last');
    ABS_240  = [Fit_WL(ind_240), samp_ABS_cor(ind_240)]; % WL ~240 & ABS
    
    NO3(i,colsM+1:colsM+3) = [RMS_ERROR ABS_240]; % RMS WL~240 ABS~240
    % ********************************************************************
    % TEST -- PLOT SPECTRA
    % ********************************************************************
    %     if RMS_ERROR > 0.003 || ABS_240(2) > 0.8 % ONLY PLOT IF NOT GOOD
    if fig_flag == 1
        clf
        figure(1)
        set(gcf, 'Units', 'normalized');
        set(gcf,'Position',[0.17 0.23 0.52 0.55])
        a1 = axes('Position', [0.13 0.12 0.775 0.59]);
        plot(WL(1,:), ABS_SW(i,:), 'k-','LineWidth', 2) % sw spectrum
        hold on
        x_lim = [min(WL(1,:)) max(WL(1,:))];
        xlim(x_lim)
        y_lim = ylim;
        set(gca, 'Linewidth',3)
        xlabel('Wavelength, nm', 'FontSize', 14)
        ylabel('Absorbance','FontSize', 14)
        set(gca, 'FontSize', 14)
        
        plot(WL(1,:), ABS_BL, 'g-','LineWidth', 2) % BL
        plot(WL(1,:), ABS_Br_tcor(i,:)+ABS_BL, 'b-','LineWidth', 2) % Br spectrum
        plot(WL(1,:), ABS_cor(i,:), 'ro','LineWidth', 1) % NO3+BL observed
        plot(WL(1,:), ABS_NO3_EXP +ABS_BL, 'r-','LineWidth', 2) % NO3+BL exp
        
        y_lim = ylim;
        ylim(y_lim)
        
        fit_patch = fill([Fit_WL(1),Fit_WL(1), Fit_WL(end),Fit_WL(end)], ...
            [y_lim(1), y_lim(2), y_lim(2), y_lim(1)],'y');
        set(fit_patch,'FaceAlpha',0.20, 'LineStyle','--','EdgeColor', ...
            [160 160 160]/256, 'LineWidth',2)
        
        %     s1 = {['FIT WINDOW (',num2str(Fit_WL(1),'%4.1f'),' - ', ...
        %           num2str(Fit_WL(end),'%4.1f'),')'],['NO3 = ', num2str(NO3(i,1),...
        %           '%2.1f'), ' µM'], ['RMS error = ', num2str(NO3(i,4),'%2.3e')]};
        
        s1 = {['Float:', d.float, '   Cast:',d.cast, '   Depth:', ...
            num2str(d.P(i),'%0.1f')],['FIT WINDOW (',num2str(Fit_WL(1), ...
            '%4.1f'),' - ', num2str(Fit_WL(end),'%4.1f'),')'], ...
            ['NO3 = ', num2str(NO3(i,1),'%2.1f'), ' µM'], ...
            ['RMS error = ', num2str(NO3(i,4),'%2.3e')]};
        
%         text (Fit_WL(1) + (Fit_WL(end) - Fit_WL(1))/2 , ...
%             y_lim(1) + (y_lim(2) - y_lim(1))* 0.9,...
%             s1,'FontSize', 16, 'HorizontalAlignment','center', ...
%             'VerticalAlignment', 'Top')
        text (WL(1,1) + (WL(1,end) - WL(1,1))* 0.4, ...
            y_lim(1) + (y_lim(2) - y_lim(1))* 0.95,...
            s1,'FontSize', 16, 'HorizontalAlignment','center', ...
            'VerticalAlignment', 'Top')
    
        
        set(gca, 'box', 'on')
        
        
        set(get(get(fit_patch,'Annotation'),'LegendInformation'),...
            'IconDisplayStyle','off'); % Don't show patch
        legend('Sample','Baseline', 'Br+BL', 'NO3+BL obs', 'NO3+BL fit', ...
            'LOcation','NorthEast')
        legend('boxoff')
        
        hold off
        
        % PLOT RESIDUALS
        a2 = axes('Position', [0.13 0.75 0.775 0.20]);
        FIT_DIF_E3 = FIT_DIF*1000;
        max_y = max(abs(FIT_DIF_E3));
        a2_ylim = [-max_y - 0.1*max_y, max_y + 0.1*max_y];
        ylim(a2_ylim)
        xlim(x_lim)
        ylabel({'Obs - Fit','x 10^{3}'},'FontSize', 14)
        set(a2,'Xticklabel','','box', 'on', 'FontSize', 14, 'Linewidth',3);
        hold on
        
        dif_patch = fill([Fit_WL(1),Fit_WL(1), Fit_WL(end),Fit_WL(end)], ...
            [a2_ylim(1), a2_ylim(2), a2_ylim(2), a2_ylim(1)],'y');
        set(dif_patch,'FaceAlpha',0.20, 'LineStyle','--','EdgeColor', ...
            [160 160 160]/256, 'LineWidth',2)
        
        plot(x_lim,[0 0],'--','color','k' ,'LineWidth', 2)
        plot(WL(1,:), FIT_DIF_E3, 'ko','LineWidth', 1, ...
            'MarkerFaceColor','r') % Resid * 1000
        
        hold off
        pause
    end
end

% ************************************************************************
% Create an output matrix and plot the data
% [SDN, DarkCur,Pres,Temp,Sal,NO3,BL_B,BL_M,RMS ERROR,WL~240,ABS~240]
NO3 = [d.SDN, d.DC, d.P, d.T, d.S, NO3];







