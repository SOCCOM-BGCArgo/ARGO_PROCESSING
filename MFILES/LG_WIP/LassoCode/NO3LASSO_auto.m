% Calculates NO3 by Lasso regression and add the Lasso value and raw NO3 to QC ODV file 
% 
% Uses calc_FLOAT_NO3_Mod or calc_FLOAT_NO3Mod250 functions to calculate
% raw NO3 from each .isus file for a float and then uses spec and cal data
% to run a Lasso regression.  Results for both are added the the ODV QC
% file.
% 
% The ABS_cor,ENO3,WL, and fit_WL parameters for
% use in the LASSO method come from the FLOAT_NO3_Mod function.
% 
%
% parse_NO3msg, parseNO3cal, 
%
%

clear
close all
 warning('off','all')

 
 load('\\atlas\Chem\ARGO_PROCESSING\DATA\CAL\MBARI_float_list.mat'); 

% select float and institution
inst='ua';  % 'wn'  'un' 'ua' 'ss' whoi or uw or sio      
float='21810';          %  19881 17NAtl 20804 WNP  8474 TPOS SO  18608 EqPSo 165 cycles 90/90 4 0 15 1 1 0 17862 EqPNo 165cycles
                   %19018 30S 1 year
                   %bad baseline SoPac 20704  19751  and 20532 20592 no of Eq also odd
                   %ETNP 20043  20137  20175  20162 20520
                   % 6960 60, 400, 10, DC 0 tder 0
                   %8482  60 500 mod2 10 sumlim 15
                   %19129 NoATL, 19531 NoATL  19970 - amazon
                   % so of eq example 21142 

fltdesc = cat(2,float,'\');

%ADJUSTABLE Parameters for LASSO fit
% uses a maximum WL of 250 nm for NO3 fit
% nparms=3;   %  now set automatically based on min O2 on profile, 
%nparms = 4 not enabled.

alpha = 1;  % LASSO parameter, 1=LASSO regression, almost 0 (breaks on 0
              % is Ridge Regression.  Between 0 and 1 is a mix.  About 0.25
              % seems good
              %20175 = 0.25, 20137 = 0.025, 20520 = 0.001
%  lambda = [.01 0.005 0.001 7.5e-4 5e-4 2.5e-4 1e-4 7.5e-5 5e-5];
%   CV = 3;  % cross validate with cv folds. Set to 1 for no cross validation
%nlambda = 10;
 lambda = 0.001;     % 0.0005 is a good value. 0.0025 will fit only NO2 ususally. if lambda = 0 then fit lambda
 CV = 1;  % cross validate with cv folds. Set to 1 for no cross validation 
 reltol = 0.0001;   % Lasso relative tolerance 

% Adjustable parameters 
maxWLto250 = 0;    % 1 to set max WL to 250 and use calc_FLOAT_NO3_Mod250, default is value in cal setup
maxWLto230 = 1;

ttype='LNE';    % LNE or Saka temp correction algorithm in nitrate calculation
 pause_on = 0 ;    % 1 to pause after each profile
 writefile = 0 ;   %  set to 1 to append NO2 to ODV QC file  
                   % doesn't yet append the LASSO values

 savefloatvizdir='C:\data\isus\argo\argo_float_files\no2\';  % directory to save floatvizqc files with NO2 added
if inst == 'wn'
    datadir= '\\seaecho\floats\WHOI\wn';
    caldir = '\\atlas\Chem\ARGO_PROCESSING\DATA\CAL\NO3_CAL\';
    qcdir = '\\atlas\Chem\ARGO_PROCESSING\DATA\FLOATVIZ\QC\';
    cal_file = join(cat(2, caldir, 'wn', float, '.cal'));
elseif inst =='un'
    datadir= '\\seaecho\floats\uw\n';
    caldir = '\\atlas\Chem\ARGO_PROCESSING\DATA\CAL\NO3_CAL\';
    qcdir = '\\atlas\Chem\ARGO_PROCESSING\DATA\FLOATVIZ\QC\';
    cal_file = join(cat(2, caldir, 'un', float, '.cal'));
elseif inst == 'ua'
    datadir= '\\seaecho\floats\merged\UW\f';
    caldir = '\\atlas\Chem\ARGO_PROCESSING\DATA\CAL\NO3_CAL\';
    qcdir = '\\atlas\Chem\ARGO_PROCESSING\DATA\FLOATVIZ\QC\';
    cal_file = join(cat(2, caldir, 'ua', float, '.cal')); 
elseif inst == 'ss'
    datadir= '\\seaecho\floats\SIO\ss',;
    caldir = '\\atlas\Chem\ARGO_PROCESSING\DATA\CAL\NO3_CAL\';
    qcdir = '\\atlas\Chem\ARGO_PROCESSING\DATA\FLOATVIZ\QC\';
    cal_file = join(cat(2, caldir, 'ss', float, '.cal'));      
end

 ncal = parseNO3cal(cal_file);    % read ISUS calibration file
%  clear ncal; load('ncal19531');   % ajusted Io for 19531

summary=zeros(237,10);
f=1; fltnum=float;

clear out 

%get all the isus files

if inst == 'ss'
    idir=[datadir float '\NO3\'];    % directory for .isus files
    ifiles=ls([idir '*.no3']);
else
    idir=[datadir float '\'];    % directory for .isus files
    ifiles=ls([idir '*.isus']);    
end

%%%%%%%%%%%%
%%%%%%%%%%%%  read the floatviz qc file, 

 
    f=d.list(:,2);
    floc=strfind(f,float);
    fidx=find(not(cellfun('isempty',floc)));
    WMO=d.list(fidx,3);
    
    qc_file = join([qcdir WMO 'qc.txt']);
    qc_file=strrep(qc_file,' ','');
    fid = fopen(string(qc_file));
    qcdata = get_FloatViz_data(string(qc_file));
 %   fclose('all')   % close all files
    
    
%     Oxyind = find(ismember(qcdata.hdr,'Oxygen[µmol/kg]'));
%     NO3ind = find(ismember(qcdata.hdr,'Nitrate[µmol/kg]'));

%     Oxy = qcdata.data(:,Oxyind);
%     Oxy(Oxy==-1.e10)=NaN;
%     NO3 = qcdata.data(:,NO3ind);
%     NO3(NO3==-1.e10)=NaN;
%     OxyNO3=Oxy;
%     OxyNO3(isnan(NO3))=NaN;   % OxyNO3 will hold oxygen values at depths with NO3 obs.
%     
    Cycleind = find(ismember(qcdata.hdr,'Station'));
    Pressind = find(ismember(qcdata.hdr,'Pressure[dbar]'));
    Cycle = qcdata.data(:,Cycleind);
    Press = qcdata.data(:,Pressind);
    
  
% structure of NO2 output data
allno2 = struct('cycle', [], 'press', [], 'slope', [],'inceptg', [],'no3raw',[]...\
    ,'no3lasso', []);

if size(ifiles,1)>0    % if there are isus files then go to work
 
 
 % loop over all isus files.  for testing use i= 1:4 
 
 lfit = cell(size(ifiles,1),70,2);        % lfit will hold the detailed LASSO fit info. for each ISUS scan 
                                          % indices are over station,
                                          % depth, 1 is conc info, 2 is fit
                                          % params
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    LOOP over all profiles 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    i is the profile index,
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    j will be the observation
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    index in each profile.
                                          
    for  i= 100:103 % size(ifiles,1)   %   loop over .isus files
i
       if inst=='ss'
           verbosity = 0;
           spec=parse_NO3sio([idir ifiles(i,:)],  verbosity);
       else
           spec=parse_NO3msg([idir ifiles(i,:)]);
       end
       
       if isempty(spec.P)
           continue
       end
       
       if isnan(spec.P)
           continue
       end

      Pcor_flag = 1;  % pressure correction on = 1 for NO3 calc
 
      %%%%%%%%% this is how to shift the fit window
%       ncal.min_fit_WL = 217;
%       ncal.max_fit_WL = 235;
      
      
     if maxWLto250==1    
        NO3all = calc_FLOAT_NO3_Mod250(spec, ncal,Pcor_flag);      % compute nitrate using max WL of 250 nm, good for noisy baseline
     elseif maxWLto230 == 1
         NO3all = calc_FLOAT_NO3_Mod230(spec, ncal,Pcor_flag);
     else
         NO3all = calc_FLOAT_NO3_Mod(spec, ncal,Pcor_flag); 
     end
     NO3new=cell2mat(NO3all(1));
     ABS_cor=cell2mat(NO3all(2));
     ENO3=cell2mat(NO3all(3));
     WL=cell2mat(NO3all(4));
     t_fit=cell2mat(NO3all(5));
     ESW=cell2mat(NO3all(6));
     ABS=cell2mat(NO3all(7));
     

if size(NO3new,1)>30                   %  just check for good data
          clear cycle xtemp xpress xtemps yvar NO2 NO3Lasso NO2Lasso S2O3Lasso SO3Lasso

         xpress=NO3new(:,3) ;   % pressure
         
         
 
 %do the LASSO fit here 


 %%  i is profile index, j is observation in profile index
 %% this will use the predicted baseline slope and will fit the intercept as a LASSO regression constant
 bl=NaN(size(ABS_cor,1),size(WL,2));
 
 
 for j=     1:size(ABS_cor,1)   %  loop over each observation in profile todo  LASSO

    
    nwl=sum(t_fit);  % number of wavelengths used
    Y=ABS_cor(j,t_fit); %%%% 8/10/24 moved ' in one )

     FW=WL(1,t_fit);
        X= NaN(2,size(FW,2));
        X(1,:)=FW;
        X(2,:)=ENO3(t_fit);
   
        X=X'; 
 %  [CONC,FitInfo]=lassoglm(X,Y,'normal','Alpha',alpha,'Lambda',lambda,'CV',CV); 
 
 
     if CV ==1
        [CONC,FitInfo]=lassoglm(X,Y,'normal','Alpha',alpha,'Lambda',lambda,'RelTol',reltol); 
     else
         [CONC,FitInfo]=lassoglm(X,Y,'normal','Alpha',alpha,'Lambda',lambda,'CV',CV,'RelTol',reltol); 
     end     
   

 
 %     
 %   [CONC,FitInfo]=lassoglm(X,Y,'normal','Alpha',alpha,'NumLambda',nlambda);

    lfit{i,j,1}=CONC;        %lfit 1 holds the LASSO concentration estimates
    lfit{i,j,2}=FitInfo;     %lfit 2 holds the LASSO fit parameters

    
    if CV==1
        index1se = 1;
    else
        index1se = FitInfo.Index1SE;   % Lasso best fit if lambda is an array
    end
 
        NO3Lasso(j)=CONC(2,index1se);
        BLslope(j) =CONC(1,index1se);
 end
 
  
        
        
        % save the data
         cycle(1:size(NO3new,1))=str2num([spec.cast]);  
       if i == 1
            allno2.cycle=cycle';  
            allno2.press=xpress;
            allno2.slope=BLslope(j);
            allno2.inceptg=NaN;
  
            allno2.no3raw=NO3new(:,6);
            allno2.no3lasso=NO3Lasso';

 
        else     
            allno2.cycle=[allno2.cycle ;cycle'];
            allno2.press=[allno2.press;xpress];
             allno2.slope=[allno2.slope BLslope(j)];
            allno2.inceptg=[allno2.inceptg NaN];

            allno2.no3raw=[allno2.no3raw; NO3new(:,6)];
            allno2.no3lasso=[allno2.no3lasso; NO3Lasso'];
 
       end 

    end
    end



     fig4=figure(4);
    set(fig4,'position',[600 200 500 400]) ;
     plot(allno2.press, allno2.no3raw, '+')
    hold on
    plot(allno2.press, allno2.no3lasso, 'or')
    xlabel('Pressure');
    ylabel('Nitrate \mumol/L');
    title(float);
    legend('NO3 Original','NO3 LASSO')
    hold off
    
        fig3=figure(3);
    set(fig3,'position',[400 200 500 400]) ;
    plot(allno2.no3raw,allno2.no3lasso, '+')

    xlabel('NO3 Original \mumol/L');
    ylabel('NO3 Lasso \mumol/L');
    title(float);
  
 

end  % end of if 

if writefile == 1   % write the odv file with NO2 appended

    % need to get WMO number

%     if fid ==-1; continue; end
    qcout = join([savefloatvizdir WMO 'lasno3.txt']);
    qcout=strrep(qcout,' ','');
    fidout=fopen(string(qcout),'wt');
    tline = fgetl(fid);
    i=1;
    printflag=0;
   
    npress=[allno2.press];
    ncycle=[allno2.cycle];

    while ischar(tline)
     %   disp(tline)


        if size(tline,2)>6 & tline(1:6)=='Cruise'    % data headers
            printflag=1;
            tline = [tline char(9) 'NO3raw[µmol/kg]' char(9) 'QF'];
  
            tline =[tline char(9) 'NO3Lasso [µmol/kg]' char(9) 'QF'];

        else if printflag ==1  %   these are data lines, add NO2

            sline = split(tline, char(9));
 
            tf2 = ncycle==str2double(sline(2));
            [tf1 tind] = nanmin(abs(npress(tf2==1)-str2double(sline(10)))) ;
            
            no3out=-1e10;
            no3outqc=1;
            no3temp = allno2.no3raw(tf2==1);
            no3tempqc = ones(size(no3temp,1));
   
            if abs(tf1)<2 
                no3out=no3temp(tind);
                no3outqc=no3tempqc(tind);
            end
            clear no2temp no2tempqc;
            tline = [tline char(9) num2str(no3out,'%7.2f') char(9) num2str(no3outqc,'%2.0f')];            

            no3out=-1e10;
            no3outqc=1;
            no3temp = allno2.no3lasso(tf2==1);
                if abs(tf1)<2 
                    no3out=no3temp(tind);
                end               
                clear no3temp
                tline = [tline char(9) num2str(no3out,'%7.2f') char(9) num2str(no3outqc,'%2.0f')];
  
            end
        end    

        fprintf(fidout,'%s\n',tline);
        tline = fgetl(fid);
        i=i+1;    
    %      if i==50; break; end

    end
%     fclose(fid); fclose(fidout);
end

%end
