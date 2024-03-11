function SBE = parse_cal_pdf(pdf_file)
%
% This function tries to find various calibration coefficients for
% profiling floats stored in pdf files
%
% A file exchange function by Derek Wood is used to extract info:
% https://www.mathworks.com/matlabcentral/fileexchange/
% 63615-read-text-from-a-pdf-document
% Before pdf reader can work add a couple lines to start%up.m to
% switch to the dir where the mfile is found and execute the following line:
% javaaddpath('iText-4.2.0-com.itextpdf.jar')
% it only needs to be executed once upon the start of a matlab session
%
% Several possible pdf calibration file types: 
% APEX:  FLBB, pH(if SBE stem on APEX)
% NAVIS: MCOMS, SUNA, pH, SBE63, SBE83(soon?)
% SOLO:  ECO TRIPLET, SBE63, pH

% CHANGES:
% 05/10/23 - JP - Added compatibility for SBE higher order f(P) sheet
%     format as well as k2 f(P) for pH calibration files as well as capture
%     further MCOMS variablity in cal files & exclude new oil cal in multi
%     page pdf
% 05/23/23 - JP - Amended to extract FLBBFL CHL435
% 01/22/24 TM, more edits to different cal sheet flavors!!!

% ************************************************************************
% TESTING
%pdf_file = 'C:\temp\1114-12157-10374 Float Calsheet.pdf'; % pH
%pdf_file = 'C:\temp\SBE B5 CalSheet.pdf'; % pH

%pdf_file = 'C:\temp\SBE 10399 CalSheet.pdf'
%pdf_file = 'C:\temp\11018 Stem Calsheet.pdf'

%pdf_file = 'C:\temp\MCOMS-0286 FDOM Char Sheet.pdf'; % FDOM
%pdf_file = 'C:\temp\MCOMS_034_CDOM.pdf'; %CDOM

%pdf_file = 'C:\temp\FLBBAP2-5534_(CHL)_CharSheet.pdf';
%pdf_file = 'C:\temp\MCOMS-0286 CHL Char Sheet.pdf';    % CHL
%pdf_file = 'C:\temp\MCOMSB-129_(CHL)CharSheet.pdf';

%pdf_file = 'C:\temp\FLBBAP2-5534_(700nm)_CharSheet.pdf'; % BBP
%pdf_file = 'C:\temp\FLBBAP2-6333_CharSheets.pdf'; %CHL & BBP in one pdf
%pdf_file = 'C:\temp\MCOMS-0286 700nm Char Sheet.pdf'; %BBP
%pdf_file = 'C:\temp\10952 Stem pH Calibration Sheet.pdf';
%pdf_file = 'C:\temp\10954 Stem pH Calibration Sheet.pdf';

%pdf_file = 'C:\temp\FLBBAP2-6121_CharSheets.pdf'

%pdf_file = 'C:\temp\11045 Stem Calsheet.pdf';
%pdf_file = 'C:\temp\11022 Stem Calsheet.pdf';

% SOLO CAL PDF PARSER TESTING
%pdf_file = 'C:\temp\FLBBCDRT2K-6610_CharSheets.pdf'; % SOLO w/ ECO triplet
%pdf_file = 'C:\temp\14319-10908pHCalsheet.pdf'; % SOLO w/ pH
%pdf_file = 'C:\temp\83-9013OxygenCalibration.pdf'; % SOLO w/ O2 phase coef
%pdf_file = 'C:\temp\83-9013TemperatureCalibration.pdf'; % SOLO w/ O2 T coef
%pdf_file = 'C:\temp\SBE_T54_CalSheet.pdf'; % SOLO w/ pcomp pH
%pdf_file = 'C:\temp\1114-12157-10374 Float Calsheet.pdf'; %ph for un1114
%pdf_file = 'C:\temp\15144-11042pHCalsheet.pdf';
%pdf_file = 'C:\temp\FLBBCDRT2K-6612_CharSheets.pdf'; %ss0002 ECO
%pdf_file = 'C:\temp\FLBBCDRT2K-6610_CharSheets.pdf'; %ss0001 ECO
%pdf_file = 'C:\temp\FLBBCDRT2K-7206_CharSheets.pdf'; %ss0003 ECO
%pdf_file = 'C:\temp\FLBBCDRT2K-6612_CharSheets.pdf'; %ss0002 ECO
%pdf_file = 'C:\temp\FLBBAP2-6925_CharSheets.pdf';
%pdf_file = 'C:\temp\18374-11964 Float Calsheet.pdf'; %un1506 pH
%pdf_file = 'C:\temp\1116-12196-10391 Float Calsheet.pdf'; % un116 pH

% differnt parsing between chem & local for Chl Scale??
%pdf_file = 'C:\temp\MCOMSC-0401_CharSheet.pdf' % un1510
%pdf_file = '\\seaecho.shore.mbari.org\floats\UW\n1510\MCOMSC-0401_CharSheet.pdf';
%pdf_file = 'C:\temp\MCOMSC-129__CHL_CharSheet.pdf';

%pdf_file = 'C:\temp\FLBBFLAP2-7343 700nm CharSheet.pdf';
% pdf_file = 'C:\temp\FLBBFLAP2-7714_CHL CharSheet.pdf' ;
%pdf_file = 'C:\temp\FLBBFLAP2-7343 CHL2 CharSheet.pdf'; %FLBBFL w\ CHL435 ua21290
%pdf_file = 'C:\temp\FLBBFLAP2-7714_CHL2 CharSheet.pdf'; %FLBBFL w\ CHL435 ua21910
%pdf_file = 'C:\temp\17682-11764 Float Calsheet.pdf'; %un1473 9th order pH
%pdf_file = 'C:\temp\17680-11218 Float Calsheet.pdf';%un1472 6th order pH
%pdf_file = 'C:\temp\16512-11568 Float Calsheet.pdf'; %un1363 k2f(P)
%pdf_file = 'C:\temp\18118-10783 Float Calsheet.pdf'; % un1502 12th order f(P) & k2f(p)
%pdf_file = '\\seaecho.shore.mbari.org\floats\SIO\ss4002\cals\SBE83-0041-Temperature.pdf';
%pdf_file = '\\seaecho.shore.mbari.org\floats\SIO\ss4002\cals\SBE83-0041-Oxygen.pdf';
% ************************************************************************

tmp = pdfRead(pdf_file); % get text from pdf file
for tmp_ct = 1:size(tmp,2) % step through "pdf pages" if greater than 1
    
    % THE NEXT FEW LINES ARE IMPORTANT! The pdf reader pulls in some odd
    % hidden characters. This reg exp line only grabs common visible
    % characters & leaves the hidden ones behind via indexing
    pdf_txt = tmp{1,tmp_ct};
    inds    = regexp(pdf_txt,'[\w.-() /=]'); % this removes hidden chars in str
    pdf_txt = pdf_txt(inds);


    % Now also including oil cal for FDOM - ignore
    page_skip_filt = 'Oil\s+Calibration';
    if ~isempty(regexp(pdf_txt, page_skip_filt,'once'))
        continue
    end

    
    % *************************************************************************
    % pH
    % *************************************************************************
    % k in k2 k0 can change case ie "11045 Stem Calsheet.pdf"
    
    if size(regexpi(pdf_txt,'f\(P)|k0|k2'),2) == 3
        SBE.pH.SN = '';
        SBE.pH.k0 = [];
        SBE.pH.k2 = [];
        SBE.pH.fP = [];
        
        % "SBE" not always in SN string block so do "zero or more",ie "11045 Stem Calsheet.pdf"
        %SN = regexp(pdf_txt,'(?<=Serial\s+Number\:*\s*)\d+','once','match');
        SN = regexp(pdf_txt,'(?<=Serial\s+Number[\:\sSBE]*)\d+','once','match');

        if ~isempty(SN)
            SBE.pH.SN = SN;
        end

        % GET pH COEFFICIENTS
        poly_order1 = str2double(regexp(pdf_txt,'\d{1}(?=th\s*order)','once','match'));
        poly_order2 = unique(regexp(pdf_txt,'f\d+|K0|K2', 'match'),'stable');
         if size(poly_order2,2) ~= poly_order1+3
             fprintf(['WARNING: stated f(P) poly order (%0.0f) is not constent ',...
                 'with meta data coefficient count (%0.0f)\n'], poly_order1, ...
                 size(poly_order2,2)-3);
             fprintf('Reseting poly order using meta data info to: %0.0f\n', ...
                size(poly_order2,2)-3);
             poly_order1 = size(poly_order2,2)-3;
         end


        coef_str   = regexp(pdf_txt,'(?<=f0\s*)\-*\d{1}.+\d{1}(?=\s*Sea)', ...
            'match','once');
        coefs      = str2double(regexp(coef_str,' ','split'));
        
        % JP fix for non standard SBE pdf file 09/01/21
        tnan       = isnan(coefs);
        if sum(tnan,2) > 0
            disp('WARNING: non standard SBE pH Calibration file!')
            disp('Removing NaNs in coefs array & trying to procede.....')
            coefs = coefs(~tnan);
        end
        
        if size(coefs,2) == poly_order1 + 3 % f(p) f0 + k0 + k2 ALL GOOD
            [S,idx] = sort(poly_order2); % k0, k2, f0, f1, ....
            coefs   = coefs(idx);

%             tf0     = strcmp(S, 'f0'); % loose f0
%             S       = S(~tf0);
%             coefs   = coefs(~tf0);

            SBE.pH.k0 = coefs(1);
            SBE.pH.k2 = coefs(2);
            SBE.pH.fP = coefs(3:end); %f1, f2.... fn (low to high poly order)
            %SBE.pH.fP = flip(SBE.pH.fP);

%             SBE.pH.fP = coefs(1:poly_order1);
%             SBE.pH.k0 = coefs(poly_order1+2);
%             SBE.pH.k2 = coefs(poly_order1+3);
        else
            disp('pH coeficient count is ambiguous')
        end
        
        % In case higher order reported with coefs == 0
        t0        = SBE.pH.fP == 0;
        SBE.pH.fP = SBE.pH.fP(~t0);

    % deal with yet another pdf format higher f(p) & k2 f(p) format
    % JP 05/09/23
    elseif size(regexpi(pdf_txt,'f\(P)|k0|k2'),2) >= 3 % i.e. 1363
        SBE.pH.SN = '';
        SBE.pH.k0 = [];
        SBE.pH.k2 = [];
        SBE.pH.fP = [];

        SN = regexp(pdf_txt,'(?<=FET\s+Serial\s+Number\:*\s*)\d+', ...
            'once','match');
        if ~isempty(SN)
            SBE.pH.SN = SN;
        end

%         poly_order2 = unique(regexp(pdf_txt,'(f\d+|K0|K2f\d)(?=\s*\=)',
%         'match'),'stable'); %ug, difficult to automate this.
        poly_order2 = regexp(pdf_txt,'(f\d+|K0|K2f\d)(?=\s*\=)', 'match');
        coefs       = regexp(pdf_txt,'(?<=\=\s*)[\deE\.-+]+','match');
%         keyboard
        if size(poly_order2,2) ~= size(coefs,2) 
            fprintf('WARNING: pH meta & coef array sizes are not equal!\n')
        else
            coefs   = str2double(coefs);
            [S,idx] = sort(poly_order2); % k0, k2, f0, f1, ....
            coefs   = coefs(idx);

%             tf0     = strcmp(S, 'f0'); % loose f0
%             S       = S(~tf0);
%             coefs   = coefs(~tf0);

            SBE.pH.k0 = coefs(1);
            tK2       = strncmp(S,'K2',2);
            SBE.pH.k2 = unique(coefs(tK2),'stable'); %this is where unique should be placed...?  If K2(P), the first coef will be repeated (as an option for K2(P) constant), ie see cal for 4003.
            tFP       = strncmp(S,'f',1);
            SBE.pH.fP = coefs(tFP);

            % if polyorder > 10 need to do some sorting for proper order
            % (un1502) so do for all

            coef_ct_str  = regexp(S(tFP),'(?<=f)\d+','match');
            coef_ct      = cellfun(@str2double, coef_ct_str);
            [jp,ia]       = sort(coef_ct);
            SBE.pH.fP = SBE.pH.fP(ia);

            % remove any coef == 0 (i.e. un1503)
            t0 = SBE.pH.fP == 0;
            SBE.pH.fP = SBE.pH.fP(~t0);
        end


%         k0_filt = '(?<=K0\s+\=\s+)[\d\.-]+';
%         k2_filt = '(?<=\s+K2f\d+\s+\=\s+)[\d\.E-]+';
%         fp_filt = '(?<=\s+f\d+\s+\=\s+)[\d\.E-+]+';
% 
%         k0 = str2double((regexpi(pdf_txt, k0_filt,'match')))';
%         k2 = str2double((regexpi(pdf_txt, k2_filt,'match')))';
%         k2(k2==0) =[];
%         fP = str2double((regexpi(pdf_txt, fp_filt,'match')))';
%         fP(fP==0) =[];
% 
%         SBE.pH.fP = fP';
%         SBE.pH.k0 = k0;
%         SBE.pH.k2 = k2';
    end


    % *************************************************************************
    % CHLOROPHYLL
    % SBS has lots of differnt pdf formats & even if the look the same
    % there can bee hiiden differences when the string is extracted from te
    % pdf SN 7714 & SN 7343 - 2 new sensors with chl435
    % *************************************************************************
    if size(regexpi(pdf_txt,'chlorophyll'),2) > 0 % chl calibration

        chl_text      = pdf_txt;

        chl_name = 'Chl';
        if ~isempty(regexp(chl_text,'435\s+nm','once')) % FLBBFL 435 CHL
            chl_name = 'Chl435';
        end

        SBE.(chl_name).Model = '';
        SBE.(chl_name).SN    = '';
        SBE.(chl_name).Scale = [];
        SBE.(chl_name).DC    = [];

        % Get serial number
        SN = regexp(chl_text ,'(?<=S/N.+-)\d+','once','match');
        if isempty(SN) %un1116
             SN = regexp(pdf_txt,'(?<=MCOMSC*.+-)\d+','once','match');
        end

        if ~isempty(SN)
            SBE.(chl_name).SN = SN;
        end

        % Try & get model
        Model = regexp(chl_text,'(?<=com)\w+(?=[\w\s]+Characterization Sheet)', ...
            'match','once');

        if isempty(Model) % Try a different format% ie ss0002 cal sheet
                    Model = regexp(chl_text,'(?<=Digital)\w+(?=[\w\s]+Characterization Sheet)', ...
            'match','once');
        end

        
        % Try again with another format diff (ie FLBBAP2-6925_CharSheets.pdf)
        if isempty(Model) 
            Model = regexp(chl_text,'FLBB\w+', 'match','once');
        end

        if isempty(Model) % try another flavor % un1510
            Model = regexp(pdf_txt,'MCOMSC*(?=\-\d+)','match','once');
        end

        
        if ~isempty(Model)
            SBE.(chl_name).Model = Model;
            if strcmp(SBE.(chl_name).Model,'E')  & ... % SN 7714
                    ~isempty(regexp(chl_text,'FLBBFL','once'))
                SBE.(chl_name).Model = regexp(chl_text,'FLBB\w+', 'match','once');
            elseif strcmp(SBE.(chl_name).Model,'E') % ss0003
                SBE.(chl_name).Model = 'ECO';
            end  
        end

        % Get dark counts
        DC = regexp(chl_text ,'(?<=Date\s*)\d+(?=\s*counts)','once','match');
        if ~isempty(DC)
            SBE.(chl_name).DC = str2double(DC);
        end
        
        % Get Scale factor
        if regexp(chl_text,'ECO','once')
            SF = regexp(chl_text ,'(?<=counts\s*)[\dE\.\-]+(?=\s*.g/l/count)', ...
                'once','match');
            if ~isempty(SF)
                SBE.(chl_name).Scale = str2double(SF);
            end
        elseif regexp(chl_text,'FLBB','once')
            SF = regexp(chl_text ,'(?<=counts\s*)[\dE\.\-]+(?=\s*.g/l/count)', ...
                'once','match');
            if ~isempty(SF)
                SBE.(chl_name).Scale = str2double(SF);
            end
        elseif regexp(chl_text ,'MCOMS','once') % 2 flavors for MCOMS (at least!)
            SF = regexp(chl_text ,'(?<=equation)[\dE\.\-]+','once','match');
            if isempty(SF)
                SF = regexp(chl_text,'[\dE\.\-]+(?=\s*.g\s+CHL\s+l-1)','once','match');
            end
            if ~isempty(SF)
                SBE.(chl_name).Scale = str2double(SF);
            end
        end
        clear SN Model DC SF
    end
    clear Dc SF SN chl_txt
    
    % *************************************************************************
    % BBP
    % *************************************************************************
    if size(regexpi(pdf_txt,'Wavelength 700'),2) > 0 % bbp700 calibration 
        SBE.bbp700.Model = '';
        SBE.bbp700.SN    = '';
        SBE.bbp700.Scale = [];
        SBE.bbp700.DC    = [];
        
        % Get serial number
        SN = regexp(pdf_txt,'(?<=FLBB.+-|MCOMSC*.+-)\d+','once','match');
        if ~isempty(SN)
            SBE.bbp700.SN = SN;
        end
        
        % Get model
        Model = regexp(pdf_txt,'(?<=com)\w+(?=[\w\s]+Characterization Sheet)', ...
            'match','once');
        if isempty(Model) % try another flavor
            Model = regexp(pdf_txt,'MCOMSC*(?=\-\d+)','match','once');
        end
        % Try again with another format diff (ie FLBBAP2-6925_CharSheets.pdf)
        if isempty(Model)
            Model = regexp(pdf_txt,'FLBB\w+', 'match','once');
        end
        if ~isempty(Model)
            SBE.bbp700.Model = Model;
            if strcmp(SBE.bbp700.Model,'E') || strcmp(SBE.bbp700.Model,'FLBBCDRT2K')% ss0003
                SBE.bbp700.Model = 'ECO';
            end
        end
       
        % Get dark counts 1st occurance of pattern is DC
        DC = regexp(pdf_txt,'\d+(?=\s*counts)','once','match');
        if ~isempty(DC)
            SBE.bbp700.DC = str2double(DC);
        end

        % Get Scale factor
        SF = regexp(pdf_txt ,'[\dE\.\-]+(?=\s*\(m-1sr-1\)|\s*m-1 sr-1)', ...
            'once','match');
        if ~isempty(SF)
            SBE.bbp700.Scale = str2double(SF);
        end
        clear SN Model DC SF
    end
    clear bbp_txt
 
    % *************************************************************************
    % CDOM or FDOM
    % *************************************************************************
    
    if size(regexpi(pdf_txt,'FDOM|CDOM'),2) > 0 % FDOM calibration
        SBE.FDOM.Model = '';
        SBE.FDOM.SN    = '';
        SBE.FDOM.Scale = [];
        SBE.FDOM.DC    = [];
        
        % Get serial number
        %SN = regexp(pdf_txt,'(?<=MCOMS.+-)\d+','once','match');
        SN = regexp(pdf_txt,'(?<=FLBB.+-|MCOMSC*.+-)\d+','once','match');
        if ~isempty(SN)
            SBE.FDOM.SN = SN;
        end
        
        % Get model
        Model = regexp(pdf_txt,'(?<=Output)\w+(?=[\w\s]+Characterization Sheet)', ...
            'match','once');
        if isempty(Model) % try another flavor
            Model = regexp(pdf_txt,'MCOMSC*(?=\-\d+)','match','once');
        end

        if ~isempty(Model)
            SBE.FDOM.Model = Model;
            if strcmp(SBE.FDOM.Model,'E') % ss0003
                SBE.FDOM.Model = 'ECO';
            end
        end
        
        % Get dark counts 1st occurance of pattern is DC
        DC = regexp(pdf_txt,'\d+(?=\s*counts)','once','match');
        if ~isempty(DC)
            SBE.FDOM.DC = str2double(DC);
        end
        
        % Get Scale factor
        SF = regexp(pdf_txt,'[\dE\.\-]+(?=\s*ppb)', ...
            'once','match');
        if ~isempty(SF)
            SBE.FDOM.Scale = str2double(SF);
        end

        clear SN Model DC SF
    end


    clear fdom_txt DC SF

    
    
    % *************************************************************************
    % OXYGEN: SBE63 or 83
    % FOR SIO THIS COMES AS 2 CAL FILES, ONE FOR [O2] & ANOTHER FOR O2
    % TEMPERATURE
    % *************************************************************************  
    % O2 PHASE CONC. CAL SHEET
    O2sz = size(regexpi(pdf_txt,'OXYGEN\s+CALIBRATION'),2);
%     s83sz = size(regexpi(pdf_txt,'83\s+CALIBRATION'),2);

    if O2sz > 0 % O2 PHASE TO [O2]
        SBE.O2.Model = '';
        SBE.O2.SN = '';
        SBE.O2.A  = [];
        SBE.O2.B  = [];
        SBE.O2.C  = [];
        SBE.O2.E  = []; % Collect but not currently not used
        
        % GET MODEL & SN
        SN = regexpi(pdf_txt,'(?<=SERIAL NUMBER.+)\d+','once','match');
        SBE.O2.SN = SN;
        SBE.O2.Model = regexp(pdf_txt,'SBE \d+','match','once');

        % GET A PHASE COEFS %TM change to regexpi !!  (for SBS83)
        A        = regexpi(pdf_txt,'(?<=A[012]\s+\=\s+)[\d\.-+e]+','match');
        SBE.O2.A = str2double(A);
        B = regexpi(pdf_txt,'(?<=B[01]\s+\=\s+)[\d\.-+e]+','match');
        SBE.O2.B = str2double(B);
        C = regexpi(pdf_txt,'(?<=C[012]\s+\=\s+)[\d\.-+e]+','match');
        SBE.O2.C = str2double(C);
        E = regexpi(pdf_txt,'(?<=E\s+\=\s+)[\d\.-+e]+','match');
        SBE.O2.E = str2double(E);
        
        clear SN A B C E
    end
    
    % O2 FOIL TEMPERATURE CAL SHEET
    SzOXYGEN = size(regexpi(pdf_txt,'OXYGEN\s+TEMPERATURE\s+CALIBRATION'),2) ;
    Sz83 = size(regexpi(pdf_txt,'83\s+TEMPERATURE\s+CALIBRATION'),2);
    if Sz83 > 0 || SzOXYGEN > 0
        SBE.O2.Model = '';
        SBE.O2.SN = '';
        SBE.O2.TA  = [];
        
        % GET MODEL & SN
        SN = regexpi(pdf_txt,'(?<=SENSOR SERIAL NUMBER.+)\d+','once','match');
        SBE.O2.SN = SN;
        SBE.O2.Model = regexp(pdf_txt,'SBE \d+','match','once');
        
        % GET A FOIL TEMPERATURE COEFS
        TA        = regexp(pdf_txt,'(?<=TA[0123]\s+\=\s+)[\d\.-+e]+','match');
        if isempty(TA)
            TA        = regexp(pdf_txt,'(?<=a[0123]\s+\=\s+)[\d\.-+e]+','match');  %sometimes they are listed more simply as 'a0 = ...', for sbs83.
        end
        SBE.O2.TA = str2double(TA);
        
        clear TA SN
    end
%     clear pdf_txt

end

% DEAL WITH MORE SBE INCONSITENt FORMATING!!!!
if isfield(SBE,'Chl') && isfield(SBE,'FDOM') && ...
        isempty(SBE.FDOM.Model)
    SBE.FDOM.Model = SBE.Chl.Model;
end

if isfield(SBE,'Chl') && isfield(SBE,'bbp700') && ...
        isempty(SBE.bbp700.Model)
    SBE.bbp700.Model = SBE.Chl.Model;
end







