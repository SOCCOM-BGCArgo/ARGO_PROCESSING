function outfilename = write_PREDEPLOYMENT_CALIB_forBGCmeta(wmo,float_num,save_dir)
% Function to generate txt file for AOML that lists the following fields for input into meta files.  Each field will have a listing corresponding to PARAMETERS in the Bfile
% PREDEPLOYMENT_CALIBRATION_EQUATION
% PREDEPLOYMENT_CALIBRATION_COEFFICIENT
% PREDEPLOYMENT_CALIBRATION_COMMENT
% PARAMETER_ACCURACY
% PARAMETER_RESOLUTION

% TANYA MAURER
% MBARI
% DEC 18, 2020
% WRITTEN TO ASSIST AOML WITH POPULATION OF SENSOR ACCURACY, RESOLUTION AND
% CALIBRATION INFORMATION IN THE META.NC ARGO FILES.
%


% TESTING:
% float_num = 'ua12540';wmo = 5905105;
%%%%%

caldir = 'Z:\ARGO_PROCESSING\DATA\CAL\';
%isusdir = 'W:\floats\';
isusdir = 'S:\UW\';
calfile = ['cal',float_num,'.mat'];
floatnumer = char(regexp(float_num,'\d*','Match'));
load([caldir,calfile])

FLOAT_INFO.FLOAT_Id = float_num;
FLOAT_INFO.WMO_Id = cal.info.WMO_ID;

%%
%DOXY-------------------------------------------------------------------
%--------------------------------------------------------------------------

if cal.info.O2_flag == 1 % has O2?
    RESOLUTION.DOXY = 0.1; %umol/kg
    if strcmp(cal.O.type,'4330')==1
        %TPHASE:
        PREDEPLOYMENT_CALIBRATION_COMMENT.TPHASE_DOXY = 'Phase measurement with blue excitation light; see TD269 Operating manual oxygen optode 4330, 4835, 4831.';
        PREDEPLOYMENT_CALIBRATION_EQUATION.TPHASE_DOXY = 'NA';
        PREDEPLOYMENT_CALIBRATION_COEFF.TPHASE_DOXY = 'NA';
        ACCURACY.TPHASE_DOXY = 'NA';
        RESOLUTION.TPHASE_DOXY = 'NA';
        %TEMP_DOXY:
        PREDEPLOYMENT_CALIBRATION_COMMENT.TEMP_DOXY = 'Optode temperature returned.';
        PREDEPLOYMENT_CALIBRATION_EQUATION.TEMP_DOXY = 'NA';
        PREDEPLOYMENT_CALIBRATION_COEFF.TEMP_DOXY = 'NA';
        ACCURACY.TEMP_DOXY = 'NA';
        RESOLUTION.TEMP_DOXY = 'NA';
        if isfield(cal.O,'SVUFoilCoef') == 1 % FROM ARGO DOXY PROCESSING DOC: This is Case 202-204-305
            ACCURACY.DOXY = 10;
            PREDEPLOYMENT_CALIBRATION_COMMENT.DOXY = 'See Processing Argo oxygen data at the DAC level cookbook, http://doi.org/10.13155/39795';
            if isfield(cal.O,'ConcCoef') == 1
                PREDEPLOYMENT_CALIBRATION_EQUATION.DOXY = 'Phase_Pcorr=TPHASE_DOXY+Pcoef1*PRES/1000; CalPhase=PhaseCoef0+PhaseCoef1*Phase_Pcorr+PhaseCoef2*Phase_Pcorr^2+PhaseCoef3*Phase_Pcorr^3; MOLAR_DOXY=ConcCoef0+ConcCoef1*MOLAR_DOXY; MOLAR_DOXY=[((c3+c4*TEMP)/(c5+c6*CalPhase))-1]/Ksv; Ksv=c0+c1*TEMP+c2*TEMP^2; O2=MOLAR_DOXY*Scorr*Pcorr; Scorr=exp[PSAL*(B0+B1*Ts+B2*Ts^2+B3*Ts^3)+C0*PSAL^2]; Pcorr=1+((Pcoef2*TEMP+Pcoef3)*PRES)/1000; Ts=ln[(298.15-TEMP)/(273.15+TEMP)]; DOXY=O2/rho, where rho is the potential density [kg/L] calculated from CTD data';
                % Currently MBARI processing does not incorporate the "A" term
                % in the Scorr.  We need to add this and reprocess (it results
                % in a very small difference in [O2]), but until we do, specify the
                % equation as we apply it.
                % Same goes for 4330 with poly coeff cal below, and sbe63.
                %             PREDEPLOYMENT_CALIBRATION_EQUATION.DOXY = 'Phase_Pcorr=TPHASE_DOXY+Pcoef1*PRES/1000; CalPhase=PhaseCoef0+PhaseCoef1*Phase_Pcorr+PhaseCoef2*Phase_Pcorr^2+PhaseCoef3*Phase_Pcorr^3; MOLAR_DOXY=ConcCoef0+ConcCoef1*MOLAR_DOXY; MOLAR_DOXY=[((c3+c4*TEMP_DOXY)/(c5+c6*CalPhase))-1]/Ksv; Ksv=c0+c1*TEMP_DOXY+c2*TEMP_DOXY^2; O2=MOLAR_DOXY*Scorr*Pcorr; Scorr=A*exp[PSAL*(B0+B1*Ts+B2*Ts^2+B3*Ts^3)+C0*PSAL^2]; A=[(1013.25-pH2O(TEMP,Spreset))/(1013.25-pH2O(TEMP,PSAL))]; pH2O(TEMP,S)=1013.25*exp[D0+D1*(100/(TEMP+273.15))+D2*ln((TEMP+273.15)/100)+D3*S]; Pcorr=1+((Pcoef2*TEMP+Pcoef3)*PRES)/1000; Ts=ln[(298.15-TEMP)/(273.15+TEMP)]; DOXY=O2/rho, where rho is the potential density [kg/L] calculated from CTD data';
                PREDEPLOYMENT_CALIBRATION_COEFF.DOXY = ['Pcoef1=0.1, Pcoef2=0.00022, Pcoef3=0.0419; B0=-0.00624523, B1=-0.00737614, B2=-0.010341, B3=-0.00817083; C0=-4.88682e-07; PhaseCoef0=',num2str(cal.O.PCoef(1)),', PhaseCoef1=',num2str(cal.O.PCoef(2)),', PhaseCoef2=',num2str(cal.O.PCoef(3)),', PhaseCoef3=',num2str(cal.O.PCoef(4)),'; c0=',num2str(cal.O.SVUFoilCoef(1)),', c1=',num2str(cal.O.SVUFoilCoef(2)),', c2=',num2str(cal.O.SVUFoilCoef(3)),', c3=',num2str(cal.O.SVUFoilCoef(4)),', c4=',num2str(cal.O.SVUFoilCoef(5)),', c5=',num2str(cal.O.SVUFoilCoef(6)),', c6=',num2str(cal.O.SVUFoilCoef(7)),'; ConcCoef0=',num2str(cal.O.ConcCoef(1)),', ConcCoef1=',num2str(cal.O.ConcCoef(2))];
            else
                PREDEPLOYMENT_CALIBRATION_EQUATION.DOXY = 'Phase_Pcorr=TPHASE_DOXY+Pcoef1*PRES/1000; CalPhase=PhaseCoef0+PhaseCoef1*Phase_Pcorr+PhaseCoef2*Phase_Pcorr^2+PhaseCoef3*Phase_Pcorr^3; MOLAR_DOXY=[((c3+c4*TEMP)/(c5+c6*CalPhase))-1]/Ksv; Ksv=c0+c1*TEMP+c2*TEMP^2; O2=MOLAR_DOXY*Scorr*Pcorr; Scorr=exp[PSAL*(B0+B1*Ts+B2*Ts^2+B3*Ts^3)+C0*PSAL^2]; Pcorr=1+((Pcoef2*TEMP+Pcoef3)*PRES)/1000; Ts=ln[(298.15-TEMP)/(273.15+TEMP)]; DOXY=O2/rho, where rho is the potential density [kg/L] calculated from CTD data';
                % Currently MBARI processing does not incorporate the "A" term
                % in the Scorr.  We need to add this and reprocess (it results
                % in a very small difference in [O2]), but until we do, specify the
                % equation as we apply it.
                % Same goes for 4330 with poly coeff cal below, and sbe63.
                %             PREDEPLOYMENT_CALIBRATION_EQUATION.DOXY = 'Phase_Pcorr=TPHASE_DOXY+Pcoef1*PRES/1000; CalPhase=PhaseCoef0+PhaseCoef1*Phase_Pcorr+PhaseCoef2*Phase_Pcorr^2+PhaseCoef3*Phase_Pcorr^3; MOLAR_DOXY=[((c3+c4*TEMP_DOXY)/(c5+c6*CalPhase))-1]/Ksv; Ksv=c0+c1*TEMP_DOXY+c2*TEMP_DOXY^2; O2=MOLAR_DOXY*Scorr*Pcorr; Scorr=A*exp[PSAL*(B0+B1*Ts+B2*Ts^2+B3*Ts^3)+C0*PSAL^2]; A=[(1013.25-pH2O(TEMP,Spreset))/(1013.25-pH2O(TEMP,PSAL))]; pH2O(TEMP,S)=1013.25*exp[D0+D1*(100/(TEMP+273.15))+D2*ln((TEMP+273.15)/100)+D3*S]; Pcorr=1+((Pcoef2*TEMP+Pcoef3)*PRES)/1000; Ts=ln[(298.15-TEMP)/(273.15+TEMP)]; DOXY=O2/rho, where rho is the potential density [kg/L] calculated from CTD data';
                PREDEPLOYMENT_CALIBRATION_COEFF.DOXY = ['Pcoef1=0.1, Pcoef2=0.00022, Pcoef3=0.0419; B0=-0.00624523, B1=-0.00737614, B2=-0.010341, B3=-0.00817083; C0=-4.88682e-07; PhaseCoef0=',num2str(cal.O.PCoef(1)),', PhaseCoef1=',num2str(cal.O.PCoef(2)),', PhaseCoef2=',num2str(cal.O.PCoef(3)),', PhaseCoef3=',num2str(cal.O.PCoef(4)),'; c0=',num2str(cal.O.SVUFoilCoef(1)),', c1=',num2str(cal.O.SVUFoilCoef(2)),', c2=',num2str(cal.O.SVUFoilCoef(3)),', c3=',num2str(cal.O.SVUFoilCoef(4)),', c4=',num2str(cal.O.SVUFoilCoef(5)),', c5=',num2str(cal.O.SVUFoilCoef(6)),', c6=',num2str(cal.O.SVUFoilCoef(7))];
            end
        else % FROM ARGO DOXY PROCESSING DOC: These are Case 202_204_303 (with ConcCoef), or Case 202_204_302 (no ConcCoef)
            ACCURACY.DOXY = 15;
            PREDEPLOYMENT_CALIBRATION_COMMENT.DOXY = 'See Processing Argo oxygen data at the DAC level cookbook, http://doi.org/10.13155/39795';
            if isfield(cal.O,'ConcCoef') == 1
                PREDEPLOYMENT_CALIBRATION_EQUATION.DOXY = 'Phase_Pcorr=TPHASE_DOXY+Pcoef1*PRES/1000; CalPhase=PhaseCoef0+PhaseCoef1*Phase_Pcorr+PhaseCoef2*Phase_Pcorr^2+PhaseCoef3*Phase_Pcorr^3; deltaP=c0*TEMP^m0*CalPhase^n0+c1*TEMP^m1*CalPhase^n1+..+c27*TEMP^m27*CalPhase^n27; AirSat=deltaP*100/[(1013.25-exp[52.57-6690.9/(TEMP+273.15)-4.681*ln(TEMP+273.15)])*0.20946]; MOLAR_DOXY=Cstar*44.614*AirSat/100;ln(Cstar)=A0+A1*Ts+A2*Ts^2+A3*Ts^3+A4*Ts^4+A5*Ts^5; Ts=ln[(298.15-TEMP)/(273.15+TEMP)]; MOLAR_DOXY=ConcCoef0+ConcCoef1*MOLAR_DOXY;O2=MOLAR_DOXY*Scorr*Pcorr; Scorr=exp[PSAL*(B0+B1*Ts+B2*Ts^2+B3*Ts^3)+C0*PSAL^2]; Pcorr=1+((Pcoef2*TEMP+Pcoef3)*PRES)/1000; DOXY=O2/rho,where rho is the potential density [kg/L] calculated from CTD data';
                %PREDEPLOYMENT_CALIBRATION_EQUATION.DOXY = 'Phase_Pcorr=TPHASE_DOXY+Pcoef1*PRES/1000; CalPhase=PhaseCoef0+PhaseCoef1*Phase_Pcorr+PhaseCoef2*Phase_Pcorr^2+PhaseCoef3*Phase_Pcorr^3; deltaP=c0*TEMP^m0*CalPhase^n0+c1*TEMP^m1*CalPhase^n1+..+c27*TEMP^m27*CalPhase^n27; AirSat=deltaP*100/[(1013.25-exp[52.57-6690.9/(TEMP+273.15)-4.681*ln(TEMP+273.15)])*0.20946]; MOLAR_DOXY=Cstar*44.614*AirSat/100;ln(Cstar)=A0+A1*Ts+A2*Ts^2+A3*Ts^3+A4*Ts^4+A5*Ts^5; Ts=ln[(298.15-TEMP)/(273.15+TEMP)]; MOLAR_DOXY=ConcCoef0+ConcCoef1*MOLAR_DOXY;O2=MOLAR_DOXY*Scorr*Pcorr;Scorr=A*exp[PSAL*(B0+B1*Ts+B2*Ts^2+B3*Ts^3)+C0*PSAL^2]; A=[(1013.25-pH2O(TEMP,Spreset))/(1013.25-pH2O(TEMP,PSAL))];pH2O(TEMP,S)=1013.25*exp[D0+D1*(100/(TEMP+273.15))+D2*ln((TEMP+273.15)/100)+D3*S]; Pcorr=1+((Pcoef2*TEMP+Pcoef3)*PRES)/1000; DOXY=O2/rho,where rho is the potential density [kg/L] calculated from CTD data';
            else
                PREDEPLOYMENT_CALIBRATION_EQUATION.DOXY = 'Phase_Pcorr=TPHASE_DOXY+Pcoef1*PRES/1000; CalPhase=PhaseCoef0+PhaseCoef1*Phase_Pcorr+PhaseCoef2*Phase_Pcorr^2+PhaseCoef3*Phase_Pcorr^3; deltaP=c0*TEMP^m0*CalPhase^n0+c1*TEMP^m1*CalPhase^n1+..+c27*TEMP^m27*CalPhase^n27; AirSat=deltaP*100/[(1013.25-exp[52.57-6690.9/(TEMP+273.15)-4.681*ln(TEMP+273.15)])*0.20946];MOLAR_DOXY=Cstar*44.614*AirSat/100;ln(Cstar)=A0+A1*Ts+A2*Ts^2+A3*Ts^3+A4*Ts^4+A5*Ts^5;O2=MOLAR_DOXY*Scorr*Pcorr;Scorr=exp[PSAL*(B0+B1*Ts+B2*Ts^2+B3*Ts^3)+C0*PSAL^2]; Ts=ln[(298.15-TEMP)/(273.15+TEMP)];Pcorr=1+((Pcoef2*TEMP+Pcoef3)*PRES)/1000; DOXY=O2/rho, where rho is the potential density [kg/L] calculated from CTD data.';
            end
            % In MBARI cal file:
            % PolyDegT ==> "m0...m27"
            % PolyDegO ==> "n0...n27"
            % FCoef ==> "c0...c27"
            % A == PA_A in calc_o24argo
            % B == PB_A
            % Note that coeffs in cal file go to n=1:21, the rest are =0 and
            % not stored.
            Ccoeff_str =[];
            Ncoeff_str = [];
            Mcoeff_str = [];
            for i = 1:28
                if i > length(cal.O.FCoef) %we only store up to 21 entries because typically the rest are zeros...place checks on this, and add the zeros so as not to confuse users, since 0:27 is listed in equation.
                    cftmp = 0;
                else
                    cftmp = cal.O.FCoef(i);
                end
                Cstr = ['c',num2str(i-1),'=',sprintf('%7.6e',cftmp),', '];
                Ccoeff_str = [Ccoeff_str Cstr];
                %----------------------------------------------------------
                if i > length(cal.O.PolyDegO) %do this check separately on each coef variable (lengths should be the same, but...just in case)
                    nftmp = 0;
                else
                    nftmp = cal.O.PolyDegO(i);
                end
                Nstr = ['n',num2str(i-1),'=',sprintf('%7.6e',nftmp),', '];
                Ncoeff_str = [Ncoeff_str Nstr];
                %----------------------------------------------------------
                if i > length(cal.O.PolyDegT)
                    mftmp = 0;
                else
                    mftmp = cal.O.PolyDegT(i);
                end
                Mstr = ['m',num2str(i-1),'=',sprintf('%7.6e',mftmp),', '];
                Mcoeff_str = [Mcoeff_str Mstr];
            end
            if isfield(cal.O,'ConcCoef')
                PREDEPLOYMENT_CALIBRATION_COEFF.DOXY = ['Pcoef1=0.1, Pcoef2=0.00022, Pcoef3=0.0419; B0=-0.00624097, B1=-0.00693498, B2=-0.00690358, B3=-0.00429155; C0=-3.11680e-07; PhaseCoef0=',num2str(cal.O.PCoef(1)),', PhaseCoef1=',num2str(cal.O.PCoef(2)),', PhaseCoef2=',num2str(cal.O.PCoef(3)),', PhaseCoef3=',num2str(cal.O.PCoef(4)),'; ',Ccoeff_str,Mcoeff_str,Ncoeff_str,'; ConcCoef0=',num2str(cal.O.ConcCoef(1)),', ConcCoef1=',num2str(cal.O.ConcCoef(2)),'; A0=2.00856, A1=3.22400, A2=3.99063, A3=4.80299, A4=0.978188, A5=1.71069'];
            else
                PREDEPLOYMENT_CALIBRATION_COEFF.DOXY = ['Pcoef1=0.1, Pcoef2=0.00022, Pcoef3=0.0419; B0=-0.00624097, B1=-0.00693498, B2=-0.00690358, B3=-0.00429155; C0=-3.11680e-07; PhaseCoef0=',num2str(cal.O.PCoef(1)),', PhaseCoef1=',num2str(cal.O.PCoef(2)),', PhaseCoef2=',num2str(cal.O.PCoef(3)),', PhaseCoef3=',num2str(cal.O.PCoef(4)),'; ',Ccoeff_str,Mcoeff_str,Ncoeff_str,'; A0=2.00856, A1=3.22400, A2=3.99063, A3=4.80299, A4=0.978188, A5=1.71069'];
            end
        end
    elseif strcmp(cal.O.type,'3830')==1
        % FROM DOXY PROCESSING DOC: This is Case 201_202_202
        ACCURACY.DOXY = 25;
        PREDEPLOYMENT_CALIBRATION_COMMENT.DOXY = 'See Processing Argo oxygen data at the DAC level cookbook, http://doi.org/10.13155/39795';
        PREDEPLOYMENT_CALIBRATION_EQUATION.DOXY = ['UNCAL_Phase=BPHASE_DOXY; Phase_Pcorr=UNCAL_Phase+Pcoef1*PRES/1000; DPHASE_DOXY=PhaseCoef0+PhaseCoef1*Phase_Pcorr+PhaseCoef2*Pcorr_Phase^2+PhaseCoef3*Pcorr_Phase^3; MOLAR_DOXY=c0+c1*DPHASE_DOXY+c2*DPHASE_DOXY^2+c3*DPHASE_DOXY^3+c4*DPHASE_DOXY^4; ci=ci_0+ci_1*TEMP+ci_2*TEMP^2+ci_3*TEMP^3, i=0..4;O2=MOLAR_DOXY*Scorr*Pcorr; Scorr=exp[PSAL*(B0+B1*Ts+B2*Ts^2+B3*Ts^3)+C0*PSAL^2]; Pcorr=1+((Pcoef2*TEMP+Pcoef3)*PRES)/1000; Ts=ln[(298.15-TEMP)/(273.15+TEMP)]; DOXY=O2/rho, where rho is the potential density [kg/L]calculated from CTD data'];
        %       PREDEPLOYMENT_CALIBRATION_EQUATION.DOXY = ['UNCAL_Phase=BPHASE_DOXY; Phase_Pcorr=UNCAL_Phase+Pcoef1*PRES/1000; DPHASE_DOXY=PhaseCoef0+PhaseCoef1*Phase_Pcorr+PhaseCoef2*Pcorr_Phase^2+PhaseCoef3*Pcorr_Phase^3; MOLAR_DOXY=c0+c1*DPHASE_DOXY+c2*DPHASE_DOXY^2+c3*DPHASE_DOXY^3+c4*DPHASE_DOXY^4; ci=ci_0+ci_1*TEMP+ci_2*TEMP^2+ci_3*TEMP^3, i=0..4;O2=MOLAR_DOXY*Scorr*Pcorr; Scorr=A*exp[PSAL*(B0+B1*Ts+B2*Ts^2+B3*Ts^3)+C0*PSAL^2]; A=[(1013.25-pH2O(TEMP,Spreset))/(1013.25-pH2O(TEMP,PSAL))]; pH2O(TEMP,S)=1013.25*exp[D0+D1*(100/(TEMP+273.15))+D2*ln((TEMP+273.15)/100)+D3*S]; Pcorr=1+((Pcoef2*TEMP+Pcoef3)*PRES)/1000; Ts=ln[(298.15-TEMP)/(273.15+TEMP)]; DOXY=O2/rho, where rho is the potential density [kg/L]calculated from CTD data'];
        mycoefstr = [];
        for ia = 1:5 %pC*
            xtmp = ia-1;
            eval(['TMP = cal.O.pC',num2str(xtmp),';'])
            TMP=flipud(TMP); %because of the way we store and access the coeffs for use in polyval
            for ib = 1:4
                TMPb = TMP(ib);
                tmpstr = ['c',num2str(xtmp),'_',num2str(ib-1),' = ',num2str(TMPb,'%9.8e'),', '];
                mycoefstr=[mycoefstr tmpstr];
            end
        end
        PREDEPLOYMENT_CALIBRATION_COEFF.DOXY = ['Pcoef1=0.1, Pcoef2=0.00022, Pcoef3=0.0419; B0=-0.00624097, B1=-0.00693498, B2=-0.00690358, B3=-0.00429155; C0=-3.11680e-07;; PhaseCoef0=',num2str(cal.O.pP(2),'%8.7f'),' PhaseCoef1=',num2str(cal.O.pP(1),'%8.7f'),' PhaseCoef2=0, PhaseCoef3=0; ',mycoefstr];
        %TEMP_DOXY:
        PREDEPLOYMENT_CALIBRATION_COMMENT.TEMP_DOXY = 'Optode temperature returned.';
        PREDEPLOYMENT_CALIBRATION_EQUATION.TEMP_DOXY = 'NA';
        PREDEPLOYMENT_CALIBRATION_COEFF.TEMP_DOXY = 'NA';
        ACCURACY.TEMP_DOXY = 'NA';
        RESOLUTION.TEMP_DOXY = 'NA';
        %BPHASE_DOXY:
        PREDEPLOYMENT_CALIBRATION_COMMENT.BPHASE_DOXY = 'Phase measurement with blue excitation light; see TD218 operating manual oxygen optode 3830, 3835, 3930, 3975, 4130, 4175.';
        PREDEPLOYMENT_CALIBRATION_EQUATION.BPHASE_DOXY = 'NA';
        PREDEPLOYMENT_CALIBRATION_COEFF.BPHASE_DOXY = 'NA';
        ACCURACY.BPHASE_DOXY = 'NA';
        RESOLUTION.BPHASE_DOXY = 'NA';
    elseif ~isempty(strfind(cal.O.type,'SBE63')) %use strfind because sometimes SBE63 entry has trailing chars?
        if strcmp(cal.info.float_type,'APEX')==1
            disp(['APEX FLOAT WITH SBE63!!! --->  ',cal.info.UW_ID])
            pause
        end
        % FROM DOXY PROCESSING DOC: This is Case 103_208_307
        ACCURACY.DOXY = 10;
        PREDEPLOYMENT_CALIBRATION_COMMENT.DOXY = 'See Processing Argo oxygen data at the DAC level cookbook, http://doi.org/10.13155/39795';
        PREDEPLOYMENT_CALIBRATION_EQUATION.DOXY = ['V=(PHASE_DELAY_DOXY+Pcoef1*PRES/1000)/39.457071; Ksv=C0+C1*TEMP_DOXY+C2*TEMP_DOXY^2; MLPL_DOXY=[(A0+A1*TEMP_DOXY+A2*V^2)/(B0+B1*V)-1]/Ksv; O2=MLPL_DOXY*Scorr*Pcorr; Scorr=exp[PSAL*(SolB0+SolB1*TS+SolB2*TS^2+SolB3*TS^3)+SolC0*PSAL^2]; Pcorr=1+((Pcoef2*TEMP+Pcoef3 )*PRES)/1000; TS=ln[(298.15–TEMP)/(273.15+TEMP)]; DOXY[umol/kg]=44.6596*O2/rho, where rho is the potential density [kg/L] calculated from CTD data'];
        %       PREDEPLOYMENT_CALIBRATION_EQUATION.DOXY = ['V=(PHASE_DELAY_DOXY+Pcoef1*PRES/1000)/39.457071; Ksv=C0+C1*TEMP_DOXY+C2*TEMP_DOXY^2; MLPL_DOXY=[(A0+A1*TEMP_DOXY+A2*V^2)/(B0+B1*V)-1]/Ksv; O2=MLPL_DOXY*Scorr*Pcorr; Scorr=A*exp[PSAL*(SolB0+SolB1*TS+SolB2*TS^2+SolB3*TS^3)+SolC0*PSAL^2]; A=[(1013.25-pH2O(TEMP,Spreset))/(1013.25-pH2O(TEMP,PSAL))]; pH2O(TEMP,S)=1013.25*exp[D0+D1*(100/(TEMP+273.15))+D2*ln((TEMP+273.15)/100)+D3*S]; Pcorr=1+((Pcoef2*TEMP+Pcoef3 )*PRES)/1000; TS=ln[(298.15–TEMP)/(273.15+TEMP)]; DOXY[umol/kg]=44.6596*O2/rho, where rho is the potential density [kg/L] calculated from CTD data'];
        PREDEPLOYMENT_CALIBRATION_COEFF.DOXY = ['Pcoef1=0.115, Pcoef2=0.00022, Pcoef3=0.0419; A0=',num2str(cal.O.PhaseCoef(1),'%9.8f'),', A1=',num2str(cal.O.PhaseCoef(2),'%9.8f'),', A2=',num2str(cal.O.PhaseCoef(3),'%9.8f'),'; B0=',num2str(cal.O.PhaseCoef(4),'%9.8f'),', B1=',num2str(cal.O.PhaseCoef(5),'%9.8f'),'; C0=',num2str(cal.O.PhaseCoef(6),'%9.8f'),', C1=',num2str(cal.O.PhaseCoef(7),'%9.8f'),', C2=',num2str(cal.O.PhaseCoef(8),'%9.8f'),'; SolB0=-0.00624523, SolB1=-0.00737614, SolB2=-0.010341, SolB3=-0.00817083; SolC0=-4.88682e-07;'];
        %PHASE-DELAY:
        PREDEPLOYMENT_CALIBRATION_COMMENT.PHASE_DELAY_DOXY = 'Sensor output phase-delay.';
        PREDEPLOYMENT_CALIBRATION_EQUATION.PHASE_DELAY_DOXY = 'NA';
        PREDEPLOYMENT_CALIBRATION_COEFF.PHASE_DELAY_DOXY = 'NA';
        ACCURACY.PHASE_DELAY_DOXY = 'NA';
        RESOLUTION.PHASE_DELAY_DOXY = 'NA';
        %TEMP_DOXY:
        PREDEPLOYMENT_CALIBRATION_COMMENT.TEMP_DOXY = 'Optode temperature.';
        PREDEPLOYMENT_CALIBRATION_EQUATION.TEMP_DOXY = 'TEMP_DOXY=1/(TA0+TA1*L+TA2*L^2+TA3*L^3)-273.15; L=ln(100000*TEMP_VOLTAGE_DOXY/(3.3-TEMP_VOLTAGE_DOXY)';
        PREDEPLOYMENT_CALIBRATION_COEFF.TEMP_DOXY = ['TA0=',num2str(cal.O.TempCoef(1)),', TA1=',num2str(cal.O.TempCoef(2)),', TA2=',num2str(cal.O.TempCoef(3)),', TA3=',num2str(cal.O.TempCoef(4))];
        ACCURACY.TEMP_DOXY = 'NA';
        RESOLUTION.TEMP_DOXY = 'NA';
        %TEMP_VOLTAGE_DOXY:
        PREDEPLOYMENT_CALIBRATION_COMMENT.TEMP_VOLTAGE_DOXY = 'Thermistor voltage, in volts.';
        PREDEPLOYMENT_CALIBRATION_EQUATION.TEMP_VOLTAGE_DOXY = 'NA';
        PREDEPLOYMENT_CALIBRATION_COEFF.TEMP_VOLTAGE_DOXY = 'NA';
        ACCURACY.TEMP_VOLTAGE_DOXY = 'NA';
        RESOLUTION.TEMP_VOLTAGE_DOXY = 'NA';
    end
end

%%
%NITRATE-------------------------------------------------------------------
%--------------------------------------------------------------------------

% currently no difference in accuracy specification between suna and isus.
if cal.info.isus_flag == 1
    ACCURACY.NITRATE = 3; %umol/kg
    RESOLUTION.NITRATE = 0.01; %umol/kg
    if strcmp(float_num,'6091SOOCN')==1 || strcmp(float_num,'7567SOOCN')==1 %This float has Nitrate sensor onboard that died right away (no isus files).  Need to populate fields for consistency with PARAMS listed in its corresponding Bfiles.
        PREDEPLOYMENT_CALIBRATION_COMMENT.NITRATE = 'The nitrate sensor on this float failed upon deployment.';
        PREDEPLOYMENT_CALIBRATION_EQUATION.NITRATE = 'NA';
        PREDEPLOYMENT_CALIBRATION_COEFF.NITRATE = 'NA';
        PREDEPLOYMENT_CALIBRATION_COMMENT.UV_INTENSITY_NITRATE = 'NA';
        PREDEPLOYMENT_CALIBRATION_EQUATION.UV_INTENSITY_NITRATE = 'NA';
        PREDEPLOYMENT_CALIBRATION_COEFF.UV_INTENSITY_NITRATE = 'NA';
        ACCURACY.UV_INTENSITY_NITRATE = 'NA';
        RESOLUTION.UV_INTENSITY_NITRATE = 'NA';
        PREDEPLOYMENT_CALIBRATION_COMMENT.UV_INTENSITY_DARK_NITRATE = 'NA';
        PREDEPLOYMENT_CALIBRATION_EQUATION.UV_INTENSITY_DARK_NITRATE = 'NA';
        PREDEPLOYMENT_CALIBRATION_COEFF.UV_INTENSITY_DARK_NITRATE = 'NA';
        ACCURACY.UV_INTENSITY_DARK_NITRATE = 'NA';
        RESOLUTION.UV_INTENSITY_DARK_NITRATE = 'NA';
    else
        
        % Only way to get pixel_start and pixel_end returned from float is from
        % an isus file returned from float!?  This info isn't in cal file cuz
        % cal file contains all 256 pixels... so, need to parse one isus file.
        for iis = 1:5 %try up to first five isus files; should have a full file in there somewhere!
            if strcmp(float_num,'9761SOOCN')==1
                isusfile = 'W:\floats\alternate\f9761\9761.001.isus'; %only 1 cycle for this float and it's in the alternate...
            else
                isusfile = [isusdir,'f',floatnumer,'\',floatnumer,'.00',num2str(iis),'.isus']
            end
            [SUCCESS,MESSAGE,MESSAGEID] = copyfile(isusfile,'C:\temp\');
            if ~SUCCESS
                isusfile = [isusdir,'n',floatnumer,'\',floatnumer,'.00',num2str(iis),'.isus'];
                [SUCCESS,MESSAGE,MESSAGEID] = copyfile(isusfile,'C:\temp\');
            end
            fn = ['C:\temp\',floatnumer,'.00',num2str(iis),'.isus'];
            spec = parse_NO3msg(fn);
            if spec.EOT == 1
                break
            end
            delete(fn)
        end
        if cal.N.pixel_base == 0
            spec.spectra_pix_range = spec.spectra_pix_range +1;
            spec.pix_fit_win = spec.pix_fit_win + 1;
            disp('Pixel registration offset by +1')
        end
        PIXEL_START = spec.spectra_pix_range(1);
        PIXEL_END = spec.spectra_pix_range(2);
        if isnan(sum(spec.pix_fit_win)) && isfield(spec,'WL_fit_win')
            pfit_low = find(cal.N.WL >= spec.WL_fit_win(1),1,'first');
            pfit_hi  = find(cal.N.WL <= spec.WL_fit_win(2),1,'last');
            spec.pix_fit_win = [pfit_low pfit_hi];
            clear pfit_low pfit_hi
        end
        if ~isempty(cal.N.min_fit_WL)
            spec.pix_fit_win(1) = find(cal.N.WL >= cal.N.min_fit_WL,1,'first');
        else
            spec.pix_fit_win(1) = find(cal.N.WL >= 217,1,'first');
        end
        
        if ~isempty(cal.N.max_fit_WL)
            spec.pix_fit_win(2)  = find(cal.N.WL <= cal.N.max_fit_WL,1,'last');
        else
            spec.pix_fit_win(2)  = find(cal.N.WL <= 240 ,1,'last');
        end
        
        PIXEL_FIT_START = spec.pix_fit_win(1);
        PIXEL_FIT_END = spec.pix_fit_win(2);
        UV_INTENSITY_REF_NITRATE_NTRANS = cal.N.Ref(PIXEL_START:PIXEL_END);
        %A-D are coeffs for the temperature correction.
        A=1.1500276;
        B=0.028400;
        C=-0.3101349;
        D=0.001222;
        OPTICAL_WAVELENGTH_OFFSET = cal.N.WL_offset;
        OPTICAL_WAVELENGTH_UV_NTRANS = cal.N.WL(PIXEL_START:PIXEL_END);
        TEMP_CAL_NITRATE = cal.N.CalTemp;
        E_SWA_NITRATE_NTRANS = cal.N.ESW(PIXEL_START:PIXEL_END);
        E_NITRATE_NTRANS = cal.N.ENO3(PIXEL_START:PIXEL_END);
        %%% MAY NEED TO ADD CALIBRATION DATE??  IS THIS A META FILE FIELD?
        PREDEPLOYMENT_CALIBRATION_COMMENT.NITRATE = 'Nitrate concentration in umol/kg; see Processing Bio-Argo nitrate concentration at the DAC Level, Version 1.1, March 3rd 2018.';
        PREDEPLOYMENT_CALIBRATION_EQUATION.NITRATE = 'The sensor returns UV_INTENSITY_DARK_NITRATE and UV_INTENSITY_NITRATE(Ntrans), a subset of continuous pixels of UV_INTENSITY_NITRATE(N), N = 1 to 256. The Ntrans indices span the interval [PIXEL_START, PIXEL_END] subset of the original array (1 to 256). Thus Ntrans(i) refers to pixel N = (PIXEL_START+i-1). PIXEL_START and PIXEL_END are defined from calibration data so that the [PIXEL_START, PIXEL_END] interval is the smallest interval of pixels that correspond to the [217 nm, 250 nm] interval of wavelengths. Only a subset of the [PIXEL_START, PIXEL_END] interval is processed to compute nitrate concentration. This subset is defined as the [PIXEL_FIT_START, PIXEL_FIT_END] interval which is the smallest interval of pixels that correspond to the [217 nm, 240 nm] interval of wavelengths (thus PIXEL_FIT_START = PIXEL_START). In the following equations the data are computed for each pixel R = PIXEL_FIT_START to PIXEL_FIT_END; ABSORBANCE_SW(R)=-log10[(UV_INTENSITY_NITRATE(R)-UV_INTENSITY_DARK_NITRATE)/UV_INTENSITY_REF_NITRATE(R)]; F(R,T)=(A+B*T)*exp[(C+D*T)*(OPTICAL_WAVELENGTH_UV(R)-OPTICAL_WAVELENGTH_OFFSET)]; E_SWA_INSITU(R)=E_SWA_NITRATE(R)*F(R,TEMP)/F(R,TEMP_CAL_NITRATE); ABSORBANCE_COR_NITRATE(R)=ABSORBANCE_SW(R)-(E_SWA_INSITU(R)*PSAL)*[1-(0.026*PRES/1000)]; Perform a multilinear regression to get MOLAR_NITRATE with estimated ABSORBANCE_COR_NITRATE(R) with ABSORBANCE_COR_NITRATE(R)=BASELINE_INTERCEPT+BASELINE_SLOPE*OPTICAL_WAVELENGTH_UV(R)+MOLAR_NITRATE*E_NITRATE(R); NITRATE=MOLAR_NITRATE/rho, where rho is the potential density [kg/L] calculated from CTD data.';
        PREDEPLOYMENT_CALIBRATION_COEFF.NITRATE = ['PIXEL_START=',num2str(PIXEL_START),...
            'PIXEL_END=',num2str(PIXEL_END),...
            ', PIXEL_FIT_START=',num2str(PIXEL_FIT_START),...
            ', PIXEL_FIT_END=',num2str(PIXEL_FIT_END),...
            '; UV_INTENSITY_REF_NITRATE(Ntrans)=',sprintf('%8.3f,',UV_INTENSITY_REF_NITRATE_NTRANS'),...
            '; A=',num2str(A),', B=',num2str(B),', C=',num2str(C),', D=',num2str(D),...
            '; OPTICAL_WAVELENGTH_OFFSET=',num2str(OPTICAL_WAVELENGTH_OFFSET),...
            '; OPTICAL_WAVELENGTH_UV(Ntrans)=',sprintf('%6.3f,',OPTICAL_WAVELENGTH_UV_NTRANS'),...
            '; TEMP_CAL_NITRATE=',num2str(TEMP_CAL_NITRATE),...
            '; E_SWA_NITRATE(Ntrans)=',sprintf('%8.7f,',E_SWA_NITRATE_NTRANS'),...
            '; E_NITRATE(Ntrans)=',sprintf('%8.7f,',E_NITRATE_NTRANS')];
        %UV_INTENSITY_NITRATE:
        PREDEPLOYMENT_CALIBRATION_COMMENT.UV_INTENSITY_NITRATE = 'Intensity of ultra violet flux from nitrate sensor.';
        PREDEPLOYMENT_CALIBRATION_EQUATION.UV_INTENSITY_NITRATE = 'NA';
        PREDEPLOYMENT_CALIBRATION_COEFF.UV_INTENSITY_NITRATE = 'NA';
        ACCURACY.UV_INTENSITY_NITRATE = 'NA';
        ACCURACY.UV_INTENSITY_DARK_NITRATE = 'NA';
        RESOLUTION.UV_INTENSITY_NITRATE = 'NA';
        RESOLUTION.UV_INTENSITY_DARK_NITRATE = 'NA';
        
        %UV_INTENSITY_DARK_NITRATE:
        PREDEPLOYMENT_CALIBRATION_COMMENT.UV_INTENSITY_DARK_NITRATE = 'Intensity of ultra violet flux dark measurement from nitrate sensor.';
        PREDEPLOYMENT_CALIBRATION_EQUATION.UV_INTENSITY_DARK_NITRATE = 'NA';
        PREDEPLOYMENT_CALIBRATION_COEFF.UV_INTENSITY_DARK_NITRATE = 'NA';
        ACCURACY.UV_INTENSITY_DARK_NITRATE = 'NA';
        ACCURACY.UV_INTENSITY_DARK_NITRATE = 'NA';
    end
end
%%
%PH------------------------------------------------------------------------
%--------------------------------------------------------------------------

% currently no difference in specification between sbe and mbari pH stems
if cal.info.pH_flag == 1
    ACCURACY.PH_FREE = 0.05;
    RESOLUTION.PH_FREE = 0.0004;
    ACCURACY.PH_TOTAL = 0.05;
    RESOLUTION.PH_TOTAL = 0.0004;
    ACCURACY.VRS_PH = 'NA';
    RESOLUTION.VRS_PH = 'NA';
    ACCURACY.IB_PH = 'NA';
    RESOLUTION.IB_PH = 'NA';
    
    PREDEPLOYMENT_CALIBRATION_COMMENT.VRS_PH = 'Voltage difference between reference and source from pH sensor (in volt).';
    PREDEPLOYMENT_CALIBRATION_COMMENT.PH_FREE = 'See Processing BGC-Argo pH data at the DAC level, https://doi.org/10.13155/57195';
    PREDEPLOYMENT_CALIBRATION_COMMENT.PH_TOTAL = 'See Processing BGC-Argo pH data at the DAC level, https://doi.org/10.13155/57195';
    
    PREDEPLOYMENT_CALIBRATION_EQUATION.VRS_PH = 'NA';
    PREDEPLOYMENT_CALIBRATION_EQUATION.PH_FREE = 'k0T=k0+k2*TEMP; pcorr=f1*PRES+f2*PRES^2+f3*PRES^3+f4*PRES^4+f5*PRES^5+f6*PRES^6; k0TP=k0T+pcorr; Tk=273.15+TEMP; Cltotal=(0.99889/35.453*PSAL/1.80655)/(1-0.001005*PSAL); ADH=3.4286e-6*TEMP^2+6.7524e-4*TEMP+0.49172143; IonS=19.924*PSAL/(1000-1.005*PSAL); log10gammaHCl=[-ADH*sqrt(IonS)/(1+1.394*sqrt(IonS))]+[(0.08885-0.000111*TEMP)*IonS]; deltaVHCl=17.85+0.1044*TEMP-0.001316*TEMP^2; log10gammaHCLtP=log10gammaHCl+[deltaVHCl*(PRES/10)/(R*Tk*ln(10))/2/10]; PH_IN_SITU_FREE=[(VRS_PH-k0TP)/(R*Tk/F*ln(10))]+[ln(Cltotal)/ln(10)]+2*log10gammaHCLtP-log10(1-0.001005*PSAL)';
    PREDEPLOYMENT_CALIBRATION_EQUATION.PH_TOTAL = 'lnKhso4fac = (-deltaVHSO4 + 0.5*KappaHSO4*(PRES/10))*(PRES/10)/(R*10*Tk); Khso4 = e^(-4276.1/Tk + 141.328 - 23.093*ln(Tk) + (-13856/Tk + 324.57 - 47.986*ln(Tk))*IonS^ 0.5 + (35474/Tk - 771.54 + 114.723*ln(Tk))*IonS - 2698/Tk*IonS^1.5 + 1776/Tk*IonS^ 2 + ln(1 - 0.001005*PSAL)); Khso4TPS = Khso4*e^(lnKhso4fac); Stotal = (0.14/96.062)*(PSAL/1.80655); PH_IN_SITU_TOTAL = PH_IN_SITU_FREE - log10(1 + Stotal/Khso4TPS);';
    
    PREDEPLOYMENT_CALIBRATION_COEFF.VRS_PH = 'NA';
    % NOTE THAT FOR SOME OLDER MBARI FLOATS, F(P) HAVE FEWER THAN 6 COEFFS, IE
    % 9031SOOCN, 7672HAWAII.  SO, FIRST DO A CHECK ON THE 6TH COEFF.
    LFP = length(cal.pH.pcoefs);
    pHFcoefs = [];
    for ifp = 1:6 %want to fill 6 coeffs
        if ifp<=LFP
            pHFcoefs(ifp) = cal.pH.pcoefs(ifp);
        else
            pHFcoefs(ifp) = 0;
        end
    end
    PREDEPLOYMENT_CALIBRATION_COEFF.PH_FREE = ['R=8.31446; F=96485; k0=',num2str(cal.pH.k0,'%7.6f'),', k2=',num2str(cal.pH.k2,'%9.8f'),', f1=',num2str(pHFcoefs(1),'%6.5e'),', f2=',num2str(pHFcoefs(2),'%6.5e'),', f3=',num2str(pHFcoefs(3),'%6.5e'),', f4=',num2str(pHFcoefs(4),'%6.5e'),', f5=',num2str(pHFcoefs(5),'%6.5e'),', f6=',num2str(pHFcoefs(6),'%6.5e')];
    PREDEPLOYMENT_CALIBRATION_COEFF.PH_TOTAL = 'NA';  %The way I've specified PH_CALIBRATION_EQUATION_TOTAL is that it's just a conversion from PH_CALIBARITON_EQUATION_FREE.  The latter depends on the sensor calib coeffs, so already listed...
    
    %DO THE IB_PH SPECIFICATION SEPARATELY, AS THIS IS ONLY POPULATED ON OUR
    %APEX FLOATS (DURAFET CALIBRATED AT MBARI).  DO A SIMPLE FLOAT-TYPE LOOKUP
    %IN THE CALFILE INFO (MAY NEED TO GET MORE SPECIFIC AS PH SENSORS ON NAVIS
    %START TO RETURN MORE DIAGNOSTIC DATA)
    if strcmp(cal.info.float_type,'APEX') == 1 %APEX float, IB_PH is populated
        PREDEPLOYMENT_CALIBRATION_COMMENT.IB_PH = 'Base current of the ISFET pH chip.';
        PREDEPLOYMENT_CALIBRATION_EQUATION.IB_PH = 'NA';
        PREDEPLOYMENT_CALIBRATION_COEFF.IB_PH = 'NA';
    end
end

%%
%CHLA----------------------------------------------------------------------
%--------------------------------------------------------------------------
if cal.info.chl_flag == 1
    ACCURACY.CHL = 0.08; %mg/m3
    RESOLUTION.CHL = 0.025; %mg/m3
    PREDEPLOYMENT_CALIBRATION_COMMENT.CHLA = 'See BGC-Argo processing chlorophyll-A concentration at the DAC level, http://dx.doi.org/10.13155/39468.';
    PREDEPLOYMENT_CALIBRATION_EQUATION.CHLA = 'CHLA=(FLUORESCENCE_CHLA-DARK_CHLA)*SCALE_CHLA';
    PREDEPLOYMENT_CALIBRATION_COEFF.CHLA = ['SCALE_CHLA=',num2str(cal.CHL.ChlScale,'%6.5f'),'; DARK_CHLA=',num2str(cal.CHL.ChlDC)];
    PREDEPLOYMENT_CALIBRATION_COMMENT.FLUORESCENCE_CHLA = 'Uncalibrated chlorophyll-a fluorescence measurement.';
    PREDEPLOYMENT_CALIBRATION_EQUATION.FLUORESCENCE_CHLA = 'NA';
    PREDEPLOYMENT_CALIBRATION_COEFF.FLUORESCENCE_CHLA = 'NA';
    ACCURACY.FLUORESCENCE_CHLA = 'NA';
    RESOLUTION.FLUORESCENCE_CHLA = 'NA';
end

%%
%BBP700------------------------------------------------------------------------
%--------------------------------------------------------------------------
if isfield(cal,'BB')
    fn = 'C:\Users\tmaurer\Documents\AOML\BBP_processing_update\MBARI_Bbp_fixed_list_WMOonly.txt';
    MBBPfix = load(fn);
    xx = find(MBBPfix(:,1) == cal.info.WMO_ID);
    ACCURACY.BBP700 = 'NA'; %no specification currently
    RESOLUTION.BBP700 = 'NA';
    if strcmp(cal.BB.type,'FLBBAP2') == 1
        khi = 1.097;
        thedeg = 124;
    elseif ~isempty(strfind(cal.BB.type,'MCOMS'))
        khi = 1.142;
        thedeg = 150;
    end
    if ~isempty(xx)
        PREDEPLOYMENT_CALIBRATION_COMMENT.BBP700 = 'The particle backscattering in m-1, estimated from the "BETA_BACKSCATTERING700" counts. Sullivan et al., 2012, Zhang et al., 2009, BETASW700 is the contribution by the pure seawater at 700nm, the calculation can be found at http://doi.org/10.17882/42916. Reprocessed on May 10, 2018 from updated coefficients provided by Andrew Bernard, SeaBird (http://doi.org/10.17882/54520).';
    else
        PREDEPLOYMENT_CALIBRATION_COMMENT.BBP700 = 'The particle backscattering in m-1, estimated from the "BETA_BACKSCATTERING700" counts. Sullivan et al., 2012, Zhang et al., 2009, BETASW700 is the contribution by the pure seawater at 700nm, the calculation can be found at http://doi.org/10.17882/42916.';
    end
    PREDEPLOYMENT_CALIBRATION_EQUATION.BBP700 = 'BBP700=2*pi*khi*((BETA_BACKSCATTERING700-DARK_BACKSCATTERING700)*SCALE_BACKSCATTERING700-BETASW700)';
    PREDEPLOYMENT_CALIBRATION_COEFF.BBP700 = ['DARK_BACKSCATTERING700=',num2str(cal.BB.BetabDC),'; SCALE_BACKSCATTERING700=',num2str(cal.BB.BetabScale,'%4.3e'),'; khi=',num2str(khi,'%4.3f'),'; BETASW700 (contribution of pure sea water) is calculated at ',num2str(thedeg),' angularDeg'];
    PREDEPLOYMENT_CALIBRATION_COMMENT.BETA_BACKSCATTERING700 = 'Uncalibrated backscattering measurement.';
    PREDEPLOYMENT_CALIBRATION_EQUATION.BETA_BACKSCATTERING700 = 'NA';
    PREDEPLOYMENT_CALIBRATION_COEFF.BETA_BACKSCATTERING700 = 'NA';
    ACCURACY.BETA_BACKSCATTERING700 = 'NA';
    RESOLUTION.BETA_BACKSCATTERING700 = 'NA';
    if strcmp(float_num,'0565SOOCN')==1
        ACCURACY.BETA_BACKSCATTERING532 = 'NA';
        RESOLUTION.BETA_BACKSCATTERING532 = 'NA';
        PREDEPLOYMENT_CALIBRATION_COMMENT.BETA_BACKSCATTERING532 = 'Uncalibrated backscattering measurement.';
        PREDEPLOYMENT_CALIBRATION_EQUATION.BETA_BACKSCATTERING532 = 'NA';
        PREDEPLOYMENT_CALIBRATION_COEFF.BETA_BACKSCATTERING532 = 'NA';
        ACCURACY.BBP532 = 'NA';
        RESOLUTION.BBP532 = 'NA';
        PREDEPLOYMENT_CALIBRATION_COMMENT.BBP532 = 'The particle backscattering in m-1, estimated from the "BETA_BACKSCATTERING532" counts.  See BGC-Argo processing particle backscattering at the DAC level, http://dx.doi.org/10.13155/39459';
        PREDEPLOYMENT_CALIBRATION_EQUATION.BBP532 = 'BBP532=2*pi*khi*((BETA_BACKSCATTERING532-DARK_BACKSCATTERING532)*SCALE_BACKSCATTERING532-BETASW532)';
        PREDEPLOYMENT_CALIBRATION_COEFF.BBP532 = ['DARK_BACKSCATTERING532=',num2str(cal.BB.BetabDC),'; SCALE_BACKSCATTERING532=',num2str(cal.BB.BetabScale,'%4.3e'),'; khi=',num2str(khi,'%4.3f'),'; BETASW532 (contribution of pure sea water) is calculated at 150 angularDeg'];
    end
end


%%
%CDOM----------------------------------------------------------------------
%--------------------------------------------------------------------------
if isfield(cal,'CDOM')
    ACCURACY.CDOM = 'NA'; %no specification currently
    RESOLUTION.CDOM = 'NA';
    PREDEPLOYMENT_CALIBRATION_COMMENT.CDOM = 'NA';
    PREDEPLOYMENT_CALIBRATION_EQUATION.CDOM = 'CDOM=(FLUORESCENCE_CDOM-DARK_CDOM)*SCALE_CDOM';
    PREDEPLOYMENT_CALIBRATION_COEFF.CDOM = ['SCALE_CDOM=',num2str(cal.CDOM.CDOMScale,'%6.5f'),'; DARK_CDOM=',num2str(cal.CDOM.CDOMDC)];
    PREDEPLOYMENT_CALIBRATION_COMMENT.FLUORESCENCE_CDOM = 'Uncalibrated fluorescence from coloured dissolved organic matter sensor.';
    PREDEPLOYMENT_CALIBRATION_EQUATION.FLUORESCENCE_CDOM = 'NA';
    PREDEPLOYMENT_CALIBRATION_COEFF.FLUORESCENCE_CDOM = 'NA';
    ACCURACY.FLUORESCENCE_CDOM = 'NA';
    RESOLUTION.FLUORESCENCE_CDOM = 'NA';
end

% keep save_dir ACCURACY RESOLUTION PREDEPLOYMENT_CALIBRATION_COMMENT PREDEPLOYMENT_CALIBRATION_EQUATION PREDEPLOYMENT_CALIBRATION_COEFF FLOAT_INFO
outfilename = [save_dir,wmo,'_',float_num,'_meta_fields_MBARI.mat'];
save(outfilename,'ACCURACY','RESOLUTION','PREDEPLOYMENT_CALIBRATION_COMMENT','PREDEPLOYMENT_CALIBRATION_EQUATION','PREDEPLOYMENT_CALIBRATION_COEFF','FLOAT_INFO');

%%-------------------------------------------------------------------------
% % % % MAY NEED TO MODIFY PROCESSING AND ACCOMPANYING META FILES SPECIFICATIONS
% % % % TO INCLUDE THE "A" FACTOR IN THE DOXY SCORR AT SOME POINT....need to check with Henry about this...
% % % D0 = 24.4543
% % % D1 = -67.4509
% % % D2 = -4.8489
% % % D3 = -5.44e-4
% % % % T = 5;
% % % % TK = T+273.15;
% % % % S = 35;
% % % EXP1 = (D0+D1.*(100./(T+273.15))+D2.*log((T+273.15)./100)+D3.*S);
% % % EXP2 = (D0+D1.*(100./(T+273.15))+D2.*log((T+273.15)./100));
% % % pH20_1 = 1013.25.*exp(EXP1);
% % % pH20_2 = 1013.25.*exp(EXP2);
% % % A = pH20_2./pH20_1;
% % %
% % % pH20_3 = exp(52.57 -(6690.9./TK) - 4.681 .* log(TK));
%%-------------------------------------------------------------------------