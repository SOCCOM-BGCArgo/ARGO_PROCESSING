% find gridded SOCAT pCO2 at each float profile position.  No interpolation,
% just match the closest date and grid point.

% calculates pCO2 from CO2 sys without WIlliams  corr and with or  withoout pK1 adj  (-0.01)
% KJ
% 
clear
close all

savefile = 1;  %  create mat file if 1

%SOCAT gridded netcdf file
pCO2nc='SOCATv2023_tracks_gridded_monthly.nc';  %  x/360, y/180, t/636 months from 1/1/1970,

% kj=ncinfo(pCO2nc);
% return

    pCO2ave=ncread(pCO2nc,'fco2_ave_weighted');    % read pCO2 _weighted or _unwtd
    pCO2std=ncread(pCO2nc,'fco2_std_weighted');    % read max pCO2 _weighted or _unwtd
    
    pCO2ls=pCO2ave ;   % pco2 is avg + 1 SD
    X=ncread(pCO2nc,'xlon'); % read lon -180 to 180
    Y=ncread(pCO2nc,'ylat');  % read lat
    mnth=-ncread(pCO2nc,'tmnth');   % mnth - as days - negative values
                                    % and data resolution is monthly
    base=datenum(1970,1,1);
    newdte=datenum(1970,1,-mnth);   % use datestr(double(newdte)) to convert to string date

    
    %     Dte=-ncread(pCO2nc,'date');   % date monthly as 6 element array with year month day hh mm ss, starts 01/15/1982
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  read floats
dir = 'C:\Users\johnson\Documents\MATLAB\ARGO_PROCESSING\DATA\FLOATVIZ\QC\';
floats=ls([dir '*txt']);   
nfloats = size(floats,1);
pressurelimit = 10;  % only data above pressurelimit

% will save output in nx12 "out" array with WMO, lon, lat, date, cycle, temp,
% salt, press, pHfloat, tafloat, pCO2 float, pCO2 landschutzer
out=zeros(1,12);




for k= 1: nfloats  % loop over all floats
    if strfind(floats(k,:),'_')>0   % skip blanks
        continue
    end
    if nfloats(1)>1
        fl = [dir floats(k,:)];
     else
        fl = [dir floats];
     end    
    
    
    %%%%   read floatviz qc file
    disp(floats(k,:));
    Str=floats(k,:);
    Str(Str < '0' | Str > '9') = [];
	WMO = sscanf(Str, '%d');
 
    if isempty(WMO); continue; end    %  if its not a WMO number skip

    data=get_FloatViz_data(fl);    % need Josh's get_Floatviz m file

    % find columns needed
    Dteind = find(ismember(data.hdr,'SDN'));
    Lonind = find(ismember(data.hdr,'Lon [°E]'));
    Latind = find(ismember(data.hdr,'Lat [°N]'));
    Tind = find(ismember(data.hdr,'Temperature[°C]'));  
    Saltind = find(ismember(data.hdr,'Salinity[pss]')); 
    Pressind = find(ismember(data.hdr,'Pressure[dbar]')); 
    pHind = find(ismember(data.hdr,'pHinsitu[Total]')); 
    TAind = find(ismember(data.hdr,'TALK_LIAR[µmol/kg]')); 
    pCO2ind = find(ismember(data.hdr,'pCO2_LIAR[µatm]')); 
    
    if isempty(pCO2ind); continue; end   % no pCO2 data, skip
    
    cycle = data.data(:,2);
    PressAll = data.data(:,Pressind);    % pressure
    shallow = PressAll < pressurelimit;  % get obs at depths above pressurelimit
    Press=data.data(shallow,Pressind);
    cycle = data.data(shallow,2);
    cycleuniq = unique(cycle);           % cycle numbers that are available
  
    % salt
    Salt=data.data(shallow,Saltind);
    Salt(Salt == -1e10) = NaN;
    Saltflag=data.data(shallow,Saltind+1);
    Salt(Saltflag ==8) = NaN;
    
    Temp=data.data(shallow,Tind);
    Temp(Temp == -1e10) = NaN;
    Tempflag=data.data(shallow,Tind+1);
    Temp(Tempflag == 8) = NaN;
    
    pHflt=data.data(shallow,pHind);
    pHflt(pHflt == -1e10) = NaN;
    pHflag=data.data(shallow,pHind+1);  % quality flag, good data = 0
    pHflt(pHflag ==1) = NaN; % was NO3flag == 8    
    pHflt(pHflag ==4) = NaN;  
    pHflt(pHflag ==8) = NaN; 
    
    TAflt=data.data(shallow,TAind);
    TAflt(TAflt == -1e10) = NaN;
    TAflag=data.data(shallow,TAind+1);  % quality flag, good data = 0
    TAflt(TAflag ==1) = NaN; % was NO3flag == 8    
    TAflt(TAflag ==4) = NaN;  
    TAflt(TAflag ==8) = NaN;     
    
    pCO2flt=data.data(shallow,pCO2ind);
    pCO2flt(pCO2flt == -1e10) = NaN;
    pCO2flag=data.data(shallow,pCO2ind+1);  % quality flag, good data = 0
    pCO2flt(pCO2flag ==1) = NaN; % was NO3flag == 8    
    pCO2flt(pCO2flag ==4) = NaN;  
    pCO2flt(pCO2flag ==8) = NaN; 
    
%%% calculate pCO2 in situ from TA and pH without Williams corr
%      Csys=CO2SYSv3(TAflt,pHflt,1,3,Salt,Temp,Temp,Press,Press,0,0,0, 0,1,10,1,2,1);
      Csys=CO2SYSv3K1(TAflt,pHflt,1,3,Salt,Temp,Temp,Press,Press,0,0,0, 0,1,10,1,2,1);
      pCO2CS=Csys(:,21);
      pCO2CS(pCO2CS==-999)=NaN;
    
     flat = data.data(shallow,Latind);
     flat (flat<-1000) = NaN;  % replace missing value -1e10
     flon = data.data(shallow,Lonind);
     flon (flon<-1000) = NaN;  
     fdate = data.data(shallow,Dteind);

     mx=size(cycleuniq,1);
     if mx == 1; continue; end   % only 1 cycle
     
     bx=cycle(1);                % first cycle available - not always 1
        fSalt=nan(mx,1);         % intialize output data
        fTemp=nan(mx,1);
        fPress=zeros(mx,1);
        fpHflt=zeros(mx,1);
        fTAflt=zeros(mx,1);
        fpCO2flt=zeros(mx,1);
        fpCO2CS=fpCO2flt;
        fLat=zeros(mx,1);
        fLon=zeros(mx,1);
        fDate=nan(mx,1);
        fWMO=zeros(mx,1);
        fcycle=cycleuniq;
        pCO2lsf=zeros(mx,1);
     
     for i=bx:mx                      % this is a loop over all cycles for a given float
        fSalt(i)=nanmean(Salt(cycle==i));
        fTemp(i)=nanmean(Temp(cycle==i));
        fPress(i)=nanmean(Press(cycle==i));
        fpHflt(i)=nanmean(pHflt(cycle==i));
        fTAflt(i)=nanmean(TAflt(cycle==i));
        fpCO2flt(i)=nanmean(pCO2flt(cycle==i));
        fpCO2CS(i)=nanmean(pCO2CS(cycle==i));
        fLat(i)=nanmean(flat(cycle==i));
        fLon(i)=nanmean(flon(cycle==i));
        if fLon(i)>180; fLon(i)=fLon(i)-360; end
        fDate(i)=nanmean(fdate(cycle==i));
        fWMO(i)=WMO;
        fcycle=1:max(cycle)'; %nanmean(cycle(cycle==i));      
     
 %%%%%%   NOW FIND THE CORRESPONDING pCO2 in the SOCAT data
        refdate =datenum(2022,12,15);  % this is max date in SOCAT 2023
        if fDate(i)>refdate   % Dec 15 2022  last date in SOCAT  if data in 2023  shift back 1 or 2 years
%             rDte = fDate(i)-ceil(double((fDate(i)-refdate))/365.25)*365;  %  move back to date in 2020
            continue  % exit for loop
        else
             rDte = fDate(i);
            
        end    

 
 %%%% find index value in SOCAT array of closest lat, lon, and date       
        [minlat(i) Blatind(i)] = min(abs(fLat(i)-Y)); 
        [minlon(i) Blonind(i)] = min(abs(fLon(i)-X)); 
        [mindte(i) Bdteind(i)] = min(abs(rDte- newdte));
%    use the indexes to find SOCAT lat, lon, date
        blat(i,1)= Y(Blatind(i));
        blon(i,1) = X(Blonind(i));
        bdte(i,1) = newdte(Bdteind(i));
%    get the SOCAT pCO2 if there is a valid float pCO2       
        if ~isnan(fpCO2flt(i))
            pCO2lsf(i) = pCO2ls(Blonind(i),Blatind(i),Bdteind(i));
        else
            pCO2lsf(i) = NaN;
        end

     end    % end of loop over float cycles
     
 
 % plot of all float and SOCAT pCO2 versus time    
    fig1=figure(1);
    set(fig1,'position',[100 300 500 400]) 
    hold on
    plot(fDate,fpCO2flt,'b+') 
    plot(fDate,pCO2lsf,'r+')
    ylim([150 600])
     ylabel('pCO2, blue=Float, red=SOCAT')
     xlabel('Year')
    datetick('x','yyyy')
    
 
% plot of all SOCAT pCO2 versus float pCO2, will get a 1:1 and model II regression later    
    fig2=figure(2);
    set(fig2,'position',[600 300 500 400])
    hold on
    plot(fpCO2flt,pCO2lsf,'b+')
     xlim([150 600])
     ylim([150 600])
     xtickformat('%.0f')
     xlabel('pCO2 Float')
     ylabel('pCO2 SOCAT')      
     fpCO2flt(fpCO2flt>1600) = NaN;
     fpCO2flt(pCO2lsf>1600) = NaN;
 
 
%  plot for each float of float and SOCAT pCO2 vs time 
      if sum(~isnan(fpCO2flt))>0
        fig3=figure(3);   
        set(fig3,'position',[1000 300 500 400])         
        plot(fDate,fpCO2flt,'b+')
        hold on
        plot(fDate,pCO2lsf,'r+')
        ylim([150 600])
        xlim([min(fdate) max(fdate)])
        ylabel('pCO2, blue=Float, red=SOCAT')
        xlabel('Year')
        datetick('x','mm/yy') 
        info=append(string(fWMO(1))," ",string(nanmean(fLat)), " ", string(nanmean(fLon)));
        title(info)
        hold off
        gname=append('figs\',string(fWMO(1)),'.png');
        saveas(fig3,gname)
      end
 
        
 % create an array with all data 
 % will save output in nx12 "out" array with WMO, lon, lat, date, cycle, temp,
 % salt, press, pHfloat, tafloat, pCO2 float, pCO2 landschutzer
     mx = size(fWMO(~isnan(fpCO2flt)),1);
     
    if size(out,1)==1    % if first float add data
        out(1:length(fWMO(~isnan(fpCO2flt))),1)=fWMO(~isnan(fpCO2flt));
        out(1:length(fWMO(~isnan(fpCO2flt))),2) = fLon(~isnan(fpCO2flt));
        out(1:length(fWMO(~isnan(fpCO2flt))),3) = fLat(~isnan(fpCO2flt));
        out(1:length(fWMO(~isnan(fpCO2flt))),4) = fDate(~isnan(fpCO2flt));
        out(1:length(fWMO(~isnan(fpCO2flt))),5) = fcycle(~isnan(fpCO2flt));
        out(1:length(fWMO(~isnan(fpCO2flt))),6) = fTemp(~isnan(fpCO2flt));
        out(1:length(fWMO(~isnan(fpCO2flt))),7) = fSalt(~isnan(fpCO2flt));
        out(1:length(fWMO(~isnan(fpCO2flt))),8) = fPress(~isnan(fpCO2flt));
        out(1:length(fWMO(~isnan(fpCO2flt))),9) = fpHflt(~isnan(fpCO2flt));
        out(1:length(fWMO(~isnan(fpCO2flt))),10) = fTAflt(~isnan(fpCO2flt));
        out(1:length(fWMO(~isnan(fpCO2flt))),11) = fpCO2flt(~isnan(fpCO2flt));
        out(1:length(fWMO(~isnan(fpCO2flt))),12) = pCO2lsf(~isnan(fpCO2flt));
        out(1:length(fWMO(~isnan(fpCO2flt))),13) = fpCO2CS(~isnan(fpCO2flt));
    else  % if later float append data
        clear new;
        new=zeros(mx,12);
        new(:,1)=fWMO(~isnan(fpCO2flt));
        new(:,2) = fLon(~isnan(fpCO2flt));
        new(:,3) = fLat(~isnan(fpCO2flt));
        new(:,4) = fDate(~isnan(fpCO2flt));
        new(:,5) = fcycle(~isnan(fpCO2flt));
        new(:,6) = fTemp(~isnan(fpCO2flt));
        new(:,7) = fSalt(~isnan(fpCO2flt));
        new(:,8) = fPress(~isnan(fpCO2flt));
        new(:,9) = fpHflt(~isnan(fpCO2flt));
        new(:,10) = fTAflt(~isnan(fpCO2flt));
        new(:,11) = fpCO2flt(~isnan(fpCO2flt));
        new(:,12) = pCO2lsf(~isnan(fpCO2flt));
        new(:,13) = fpCO2CS(~isnan(fpCO2flt));
        out = cat(1,out, new);
    end
%  


% 
     
end

% add 1:1 and model II regression to figure 2
figure(2)
plot([150 600],[150 600],'k--')
[mf,bf,rf,smf,sbf]=lsqfitgm(out(~isnan(out(:,12)),11),out(~isnan(out(:,12)),12));
plot([150 600],[150*mf+bf 600*mf+bf],'r-')


%save out array
if savefile==1
    save('float_SOCAT_pCO2_Co2sysK1','out')
end

return

% subset plot to No. Atl.
