function SBE = parse_cal_pdf(pdf_file)
%
% This function tries to find various calibration coefficients for
% profiling floats stored in pdf files
%
% A file exchange function by Derek Wood is used to extract info:
% https://www.mathworks.com/matlabcentral/fileexchange/
% 63615-read-text-from-a-pdf-document
% Before pdf reader can work add a couple lines to startup.m to
% switch to the dir where the mfile is found and execute the following line:
% javaaddpath('iText-4.2.0-com.itextpdf.jar')
% it only needs to be executed once upon the start of a matlab session
%
% Several possible pdf calibration file types: 
% APEX:  FLBB, pH(if SBE stem on APEX)
% NAVIS: MCOMS, SUNA, pH, SBE63, SBE83(soon?)
% SOLO:  ECO TRIPLET, SBE63, pH

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

% ************************************************************************

tmp = pdfRead(pdf_file); % get text from pdf file

for tmp_ct = 1:size(tmp,2) % step through "pdf pages" if greater than 1
    
    % THE NEXT FEW LINES ARE IMPORTANT! The pdf reader pulls in some odd
    % hidden characters. This reg exp line only grabs common visible
    % characters & leaves the hidden ones behind via indexing
    pdf_txt = tmp{1,tmp_ct};
    inds    = regexp(pdf_txt,'[\w.-() /=]'); % this removes hidden chars in str
    pdf_txt = pdf_txt(inds);
    
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
        SN = regexp(pdf_txt,'(?<=Serial Number (SBE)*\s*)\w+(?=\s*Calib)', ...
            'once','match');
        if isempty(SN) % if SN is empty try pcomp cal sheet flavor
            SN = regexp(pdf_txt,'(?<=Serial Number SBE\s*)\w+(?=\s*p-comp)', ...
            'once','match');
        end
        
        if ~isempty(SN)
            SBE.pH.SN = SN;
        end
        
        % GET pH COEFFICIENTS
        poly_order = str2double(regexp(pdf_txt,'\d{1}(?=th\s*order)','once','match'));
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
        
        if size(coefs,2) == poly_order + 3 % f(p) f0 + k0 + k2 ALL GOOD
            SBE.pH.fP = coefs(1:poly_order);
            SBE.pH.k0 = coefs(poly_order+2);
            SBE.pH.k2 = coefs(poly_order+3);
        else
            disp('pH coeficient count is ambiguous')
        end
        
        % In case higher order reported with coefs == 0
        t0        = SBE.pH.fP == 0;
        SBE.pH.fP = SBE.pH.fP(~t0);
    end
    


    % *************************************************************************
    % CHLOROPHYLL
    % *************************************************************************
    if size(regexpi(pdf_txt,'chlorophyll'),2) > 0 % chl calibration
        chl_text      = pdf_txt;
        SBE.Chl.Model = '';
        SBE.Chl.SN    = '';
        SBE.Chl.Scale = [];
        SBE.Chl.DC    = [];
        
        % Get serial number
        SN = regexp(chl_text ,'(?<=S/N.+-)\d+','once','match');
        if ~isempty(SN)
            SBE.Chl.SN = SN;
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
        
        if ~isempty(Model)
            SBE.Chl.Model = Model;
            if strcmp(SBE.Chl.Model,'E') % ss0003
                SBE.Chl.Model = 'ECO';
            end  
        end

        % Get dark counts
        DC = regexp(chl_text ,'(?<=Date\s*)\d+(?=\s*counts)','once','match');
        if ~isempty(DC)
            SBE.Chl.DC = str2double(DC);
        end
        
        % Get Scale factor
        if regexp(chl_text,'ECO','once')
            SF = regexp(chl_text ,'(?<=counts\s*)[\dE\.\-]+(?=\s*.g/l/count)', ...
                'once','match');
            if ~isempty(SF)
                SBE.Chl.Scale = str2double(SF);
            end
        elseif regexp(chl_text,'FLBB','once')
            SF = regexp(chl_text ,'(?<=counts\s*)[\dE\.\-]+(?=\s*.g/l/count)', ...
                'once','match');
            if ~isempty(SF)
                SBE.Chl.Scale = str2double(SF);
            end
        elseif regexp(chl_text ,'MCOMS','once') % 2 flavors for MCOMS (at least!)
            SF = regexp(chl_text ,'(?<=equation)[\dE\.\-]+','once','match');
            if isempty(SF)
                SF = regexp(chl_text,'[\dE\.\-]+(?=\s*.g CHL l-1)','once','match');
            end
            if ~isempty(SF)
                SBE.Chl.Scale = str2double(SF);
            end
        end
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
        SN = regexp(pdf_txt,'(?<=FLBB.+-|MCOMS.+-)\d+','once','match');
        if ~isempty(SN)
            SBE.bbp700.SN = SN;
        end
        
        % Get model
        Model = regexp(pdf_txt,'(?<=com)\w+(?=[\w\s]+Characterization Sheet)', ...
            'match','once');
        if ~isempty(Model)
            SBE.bbp700.Model = Model;
            if strcmp(SBE.bbp700.Model,'E') % ss0003
                SBE.bbp700.Model = 'ECO';
            end
        end
        
        % Try again with another format diff (ie FLBBAP2-6925_CharSheets.pdf)
        if isempty(Model)
            Model = regexp(chl_text,'FLBB\w+', 'match','once');
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
        SN = regexp(pdf_txt,'(?<=FLBB.+-|MCOMS.+-)\d+','once','match');
        if ~isempty(SN)
            SBE.FDOM.SN = SN;
        end
        
        % Get model
        Model = regexp(pdf_txt,'(?<=Output)\w+(?=[\w\s]+Characterization Sheet)', ...
            'match','once');
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
    end
    clear fdom_txt DC SF
    
    % *************************************************************************
    % OXYGEN: SBE63 or 83
    % FOR SIO THIS COMES AS 2 CAL FILES, ONE FOR [O2] & ANOTHER FOR O2
    % TEMPERATURE
    % *************************************************************************  
    % O2 PHASE CONC. CAL SHEET
    if size(regexpi(pdf_txt,'OXYGEN\s+CALIBRATION'),2) > 0 % O2 PHASE TO [O2]
        SBE.O2.Model = '';
        SBE.O2.SN = '';
        SBE.O2.A  = [];
        SBE.O2.B  = [];
        SBE.O2.C  = [];
        SBE.O2.E  = []; % Collect but not currently not used
        
        % GET MODEL & SN
        SN = regexp(pdf_txt,'(?<=SENSOR SERIAL NUMBER.+)\d+','once','match');
        SBE.O2.SN = SN;
        SBE.O2.Model = regexp(pdf_txt,'SBE \d+','match','once');

        % GET A PHASE COEFS
        A        = regexp(pdf_txt,'(?<=A[012]\s+\=\s+)[\d\.-+e]+','match');
        SBE.O2.A = str2double(A);
        B = regexp(pdf_txt,'(?<=B[01]\s+\=\s+)[\d\.-+e]+','match');
        SBE.O2.B = str2double(B);
        C = regexp(pdf_txt,'(?<=C[012]\s+\=\s+)[\d\.-+e]+','match');
        SBE.O2.C = str2double(C);
        E = regexp(pdf_txt,'(?<=E\s+\=\s+)[\d\.-+e]+','match');
        SBE.O2.E = str2double(E);
        
        clear SN A B C E
    end
    
    % O2 FOIL TEMPERATURE CAL SHEET
    if size(regexpi(pdf_txt,'OXYGEN\s+TEMPERATURE\s+CALIBRATION'),2) > 0
        SBE.O2.Model = '';
        SBE.O2.SN = '';
        SBE.O2.TA  = [];
        
        % GET MODEL & SN
        SN = regexp(pdf_txt,'(?<=SENSOR SERIAL NUMBER.+)\d+','once','match');
        SBE.O2.SN = SN;
        SBE.O2.Model = regexp(pdf_txt,'SBE \d+','match','once');
        
        % GET A FOIL TEMPERATURE COEFS
        TA        = regexp(pdf_txt,'(?<=TA[0123]\s+\=\s+)[\d\.-+e]+','match');
        SBE.O2.TA = str2double(TA);
        
        clear TA SN
    end
    clear pdf_txt

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







