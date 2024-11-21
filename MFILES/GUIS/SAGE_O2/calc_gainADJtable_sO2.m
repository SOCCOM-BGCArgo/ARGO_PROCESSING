function [tableDATA,brkRESIDS] = calc_gainADJtable_sO2(DATA,cyST,cyEND,is_floatQC)
% ************************************************************************
% calc_gainADJtable_sO2.m
% ************************************************************************
%
% Function to calculate slopes and intercepts for breakpoint segments
% identified by user within "oxygen gain adjustment' table.
%
%
% USE AS: [DATA.tableDATA,DATA.brkRESIDS] = calc_gainADJtable_sO2(DATA,cyST,cyEND,1)
%
% INPUTS:
%    DATA: input DATA structure for SAGE_O2 gui
%    cyST, cEND: cycle start and cycle end variables.  Either 1x1 vectors or
%       1xn vectors containing start and end cycles for each breakpoint
%       segment.
%    is_floatQC: logical (1 or 0), identifies whether the inputs are user
%       defined (0) or defined by the floatQC doc (1).  If =1, then slopes and
%        intercepts within tableDATA are used to calculate residuals, as
%        opposed to recalculation of slopes and intercepts.  This is needed for
%        appropriately calculating AIC for floatQC that uses an average gain
%       over the whole time period instead of the automatic slope&int
%       calculations in the new method.**
%
% OUTPUTS: tableDATAP: updated tableDATA array for gui display
%          brkRESIDS: vector of residuals (after regressions are applied to
%          each segment) for use in calculating AIC value.
%
%
% AUTHOR: Tanya Maurer
%         MBARI
%         tmaurer@mbari.org
%
% DATE: 03/14/2017
% UPDATES:
% NOTES: ** is_floatQC may be removed in future updates.
% ************************************************************************
%
% ************************************************************************

RES=[];
tableD=nan(size(cyST,1),3);
for i = 1:size(cyST,1) %do calculation for each adjustment row
    regrX = DATA.refdata(DATA.refdata(:,2)>=cyST(i) & DATA.refdata(:,2)<=cyEND(i),1);
    if DATA.iswoaref == 0
        if ~isempty(DATA.rawGAINS{3})
            regrY = DATA.rawGAINS{3}(DATA.refdata(:,2)>=cyST(i) & DATA.refdata(:,2)<=cyEND(i)); %if new aircal method exists, use it
        else
            regrY = DATA.rawGAINS{1}(DATA.refdata(:,2)>=cyST(i) & DATA.refdata(:,2)<=cyEND(i));  %otherwise use older
        end
    elseif DATA.iswoaref == 1
        regrY = DATA.rawGAINS_WOA{1}(DATA.refdata(:,2)>=cyST(i) & DATA.refdata(:,2)<=cyEND(i));
    end
    
    regrX = regrX-nanmin(regrX); %subtract minimum x-value, so x variable is 'time since cycle at specified breakpoint'
    N = isnan(regrY)| isnan(regrX);
    rX = regrX(~N);
    rY = regrY(~N);
    
    if is_floatQC == 1 %&& DATA.tableDATA(1,3) == 0 %for calculating preloaded mean gains
        tbld = DATA.tableDATA; %will already exist
        mCy = tbld(i,3)./365;
        bCy = tbld(i,2);
        newY = rX.*mCy+bCy;
        res = newY-rY;
    else
        if i == 1
            [mCy,bCy,~,~,~]=lsqfity(rX,rY);
            newY = rX.*mCy+bCy;
            newoff = newY(end);
            res = newY-rY;
            ISSIG = test0slope(rX,rY,0.05);
            if  is_floatQC ~= 1
                if ISSIG.test == 0
                    msgbox('SLOPE IS INSIGNIFICANT.')
                else
                    msgbox('SLOPE IS SIGNIFICANT.')
                end
            end
        else
            if isempty(rX) && isempty(rY) %small fix to prevent code barf if missing air data at last cycle (but NCEP present) 8/26/21
                bCy = newoff;
                mCy = 0;
            else
                RY = rY-newoff;
                mCy = (sum(rX.*RY))/(sum(rX.*rX));
                newY = rX.*mCy+newoff;
                bCy = newoff;
                res = newY(2:end)-rY(2:end);
                newoff = newY(end);
            end
        end
    end
    %     if is_floatQC == 1
    %         mpYR = mCy;
    %     else
    mpYR = mCy.*365;
    %     end
    %     mpYR
    tableD(i,:) = [cyST(i) bCy mpYR];
    RES=[RES;res];
end
tableDATA = tableD;
tableDATA(isnan(tableDATA))=0;%replace nan with zero for consistency...
brkRESIDS=RES;
