function [gui,DATA] = apply_QCadj_sO2Argo(gui,inputs,DATA)
% ************************************************************************
% apply_QCadj_sO2Argo.m
% ************************************************************************
%
% Function to apply QC adjustments, as listed in the Oxygen Gain adjustment
% table.
%
%
% USE AS: [gui,DATA] = apply_QCadj_sO2Argo(gui,inputs,DATA)
%
% INPUTS:
%    gui, inputs, and DATA area all structures of handles and inputs to the
%    SAGE_O2 gui.
%
% OUTPUTS: updated structures for gui and DATA.
%
%
% AUTHOR: Tanya Maurer
%         MBARI
%         tmaurer@mbari.org
%
% DATE: 03/14/2017
% UPDATES:
% NOTES: 
% ************************************************************************
%
% ************************************************************************

        if inputs.rorq ==1 
            DATA.RAWorQCprof = DATA.O2data{1};
            DATA.RAWorQCair = DATA.O2air;
            DATA.RAWorQCwoa = DATA.SURF_SAT;
            gui.whichAX = gui.t1ax;
        else
            gui.whichAX = gui.t2ax;
            QCadjs = DATA.tableDATA; %QC adjustments, from Sage_O2 adjustment table. 
            CYC=[];
            G=[];
            for i = 1:size(QCadjs,1);
                cyST = QCadjs(i,1);
                if i == size(QCadjs,1) %last segment
                    cyEND = max(DATA.refdata(:,2)); %end cycle is last cycle
                else
                    cyEND = QCadjs(i+1,1)-1; %else, end cycle is cycle before next breakpoint
                end
                m = QCadjs(i,3);
                b = QCadjs(i,2);
                xx = DATA.refdata(DATA.refdata(:,2)>=cyST & DATA.refdata(:,2)<=cyEND,1); %time
                X = xx-nanmin(xx);
                c = DATA.refdata(DATA.refdata(:,2)>=cyST & DATA.refdata(:,2)<=cyEND,2); %cycle
                g = m.*X./365+b;
                G = [G;g];
                CYC = [CYC;c];
            end
            %sometimes cycles are missing aircal readings, interpolate gain
            %values for these (same for profile data)
            Gint = interp1(CYC,G,DATA.SURF_SAT(:,1));
            DATA.RAWorQCwoa(:,2) = DATA.SURF_SAT(:,2).*Gint;
            if ~isempty(DATA.O2air{1}) && get(gui.rb3(2),'Value')~=1;
                DATA.RAWorQCair = DATA.O2air;
                for  s = 1:DATA.howmanyairs;
                    airGAIN = repmat(G,1,5);
                    DATA.RAWorQCair{1}{s}(:,7:11) = DATA.RAWorQCair{1}{s}(:,7:11).*airGAIN;
                end
            end
            %multiply profile data by gain
            DATA.RAWorQCprof = DATA.O2data{1};
            profgains = nan(size(DATA.RAWorQCprof,1),1);
            profCYC = unique(DATA.RAWorQCprof(:,2));
            GintP = interp1(CYC,G,profCYC);
            % ie for float 7558 there are more profiles beyond the last
            % profile with air cal.  So, interp function results in NaNs
            % for the last 3 profiles.  If more profs beyond last aircal,
            % propagate the last gain value to the end profiles.
            ipmax = nanmax(CYC);
            GintP(profCYC>ipmax) = G(CYC==ipmax);
            
            for j = 1:length(profCYC);
                cyclenum = profCYC(j);
                r = find(DATA.RAWorQCprof(:,2)==cyclenum);
                profgains(r) = repmat(GintP(j),length(r),1);
            end    
            profGAIN = repmat(profgains,1,5);
            DATA.RAWorQCprof(:,7:end) = DATA.RAWorQCprof(:,7:end).*profGAIN; 
        end
end %end function