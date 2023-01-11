function [Y,Xf,Af] = ESPER_phosphate_4_Other_1(X,~,~)
%ESPER_PHOSPHATE_4_OTHER_1 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:27.
% 
% [Y] = ESPER_phosphate_4_Other_1(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timesteps
%   Each X{1,ts} = 7xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 1xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-0.999999984769129;-0.999998793562952;-77.961;0;29.715;-2.18999082656054;-2.02];
x1_step1.gain = [1.00003298043486;1.00000060321889;0.0138606862425759;0.000266418009857466;0.193068829037552;0.0582073296585987;0.00819772922900357];
x1_step1.ymin = -1;

% Layer 1
b1 = [2.4426922438517237879;2.1082568389349196636;2.5916916302870744993;3.6236987212824498528;-2.332838662811775432;3.5332013389581393703;-0.4629639448043169625;-2.710537109223527974;-2.3989084151175799065;-4.8832044594830152917;-0.87135386904725153911;-1.5440065741451496351;-3.9644881646128591335;-0.1518704348293893891;2.4448941697068446999;3.8669944659206434423;2.3482270350749980814;1.8513685398358823431;-0.37088272673814554059;-1.6785286967540762326;3.1323720296479966585;0.15002031399838924863;0.93209565993655496552;-1.2408755963589606086;-3.4705929101435173401;2.5127363743168777965;-0.8156147064481106268;0.9874910385443954608;2.8111286466264502337;2.4548225319606244987;1.3645924889991793982;3.8978497333041643991;-3.1484489521929681644;4.2245734750162382554;-4.0594007063660049184;-2.6172783716202512316;2.5984432196472107002;-3.5370510887157959878;0.21386288106254100483;2.7559665383270082195];
IW1_1 = [-0.49055243886782984308 0.07292040030077505397 0.80515616391792632633 -0.33679115484047378315 -5.9819763853523948427 -0.078662922506969504277 1.1082589897084862418;-0.26283193607579169448 0.40775987287814535742 -0.65831306826901814233 0.94185618817845162543 0.58062962162092346396 -0.29691713657894774903 0.77361352102510116957;-0.35393247525439514023 0.08863869072768525148 -0.3881295569702696846 0.53750461241122338141 -2.4672465918841464116 0.48473542275641201549 1.2549060914110261589;-0.66381472549153885954 0.089105650788799603546 3.6090381711695171951 -0.74218326314642579788 1.9046916875046233653 -2.2653495930276004167 2.040780067395529862;1.8291227565486305462 1.5377001057175516685 0.65246451738960753453 -1.2420018997548167849 -0.0082329792338318268097 0.26393699807745796271 -0.28978471037907277408;-0.038136721652269972571 -0.27141972261759605978 1.5657414969803169491 0.014039437174632428407 1.4473361686084988964 3.390047159186115433 0.65880151006747067477;0.030103831822748172709 -0.093900254091692669789 0.84701388589014858788 -0.20064653315651317933 -0.67537398698430184396 -0.49985755812714560697 0.57676895597309807417;0.0097063138714640721794 -0.55913889241875480707 -3.2344087148803035703 -0.26080742136596851699 -1.660308356137395025 0.4968488873427596153 -2.1648568740773335151;0.42560855028743566075 0.41999238380591263908 0.94671719174375001948 -0.99928145947890434364 -0.43723878047481085085 -0.19112303430633453827 -1.2992631987722251807;0.28482656074996126616 -0.065146889945119329512 -2.4224304869189010603 0.10835076886806020713 -9.369619842896213413 -3.0545964957804438988 1.3119804041098723779;0.54348647788046455798 0.92261806449268857211 1.0901623768852559682 0.1820903269383880585 -1.9371492407488553056 1.3291777438557521585 -1.2663713801216027743;0.45861561892941121732 0.75945322984441132608 0.86969997983525881846 0.19260215001029049042 2.1293157234916133547 -0.13177664796908664102 0.026733750544085124745;-0.019234275397680834774 -0.27180452016156253503 0.18966890244921288944 -1.2041499054541602298 -1.2793597182564060422 1.1092503558559954069 -2.290672097366081239;0.074858160892550026921 -0.13861473037540508746 -0.241797627087279976 -0.24167088638782313614 -2.4416516141407615947 1.413924688577881783 -0.57722038396031516427;-0.3840412124562837537 -0.052182943963062208836 1.6285744002088291538 -0.56337929086786797672 -1.7924415495679881793 1.5526466141942698584 -1.5758711880702140906;0.32462275928563483207 0.17502224226438772847 -0.10490974227047032696 -0.90390654517469481988 -0.41271643693766291161 5.5803182767123882257 1.5207572542151970296;-0.48323187205727063942 -0.68275056748376217186 -1.3145216886130710598 -0.16320563399612061928 1.0696433719338398483 2.3335541763143168303 -0.73121935992405495686;0.43264177555143312759 0.69109251076559985449 0.85072119912175092615 0.53526383699604329625 -1.1357358349214601301 1.2396403185272768521 1.1397413937605180845;0.32623286356016101273 0.19814476791419002866 -1.2765210911299009755 0.86209013010322743664 6.1306500215142065713 0.22274868749754844766 0.74006230352398960637;0.94440498122287375171 -0.41169509573678902958 -1.1562089563889528243 -0.40852728563624535196 -1.3435084273367203789 1.6980980878782780952 0.30317173755409015934;0.097950033948826037422 0.0028390819795553906037 1.5021573623774111539 0.071239728054910481103 1.5323682879314690641 3.3428685647980604934 0.36488590977843976537;0.33028932696602508745 -0.32642696086344169526 0.77241344775235887266 -1.2839686467163109462 1.0839427513783645995 0.42956532505613076456 0.0028827112081526613807;0.16513865139873334931 -0.15254710663646170454 5.4457207485912428879 0.9954571566290916218 -0.25016243068366600388 1.0487376858076735697 0.75335564130007337535;0.042880319491126514531 0.42946181881433226346 -0.77684210449991120129 -0.048223745811586116172 1.3737656976588628055 -1.5295374929706280209 -1.647338952183171612;-1.423425218099530376 -0.24597011812525776397 -4.9597303379375912158 0.79268580811346367465 -1.1514037767826521552 1.1709565380782029909 -2.599051027362516475;1.2803792047469941817 2.1803540644322079167 1.0065260409702827538 -0.48092465661161037938 1.2219321988457618744 -0.86420707314169642466 1.8328802385539517683;-0.75019527587298051596 -0.40986549621455764525 0.5627127109961646001 0.39846169513315343735 -1.915553586532906083 -1.2596872162765364589 -0.060071131951512542202;0.14023101199858681309 -0.094095134335701435813 5.743483666018095235 1.0664056878714791576 -0.12675958304591028458 1.1457176287492962796 0.77993743734286213876;0.54660928788750706708 -0.50807799183810886667 -3.8788380428100248665 0.040250251085836276366 -2.5809466064262966967 -0.032933151177715155877 0.66649791253596146756;0.079171215490840815554 0.60037783060048099948 -1.8614720085882536793 -0.040819475848400683937 3.4437824607160703039 -1.278991733813900078 2.8989385121438759896;-0.0590885725941714679 -1.0510763733919730178 -0.83478223940287443749 0.67812546845111731919 4.4435748379558326349 -0.15724413230272468733 0.67493872197835902416;1.6145072772225579349 0.35644025523187278237 5.05312598405885538 -0.69121495959855094871 -0.84769889476157433972 -0.81005933270391627943 3.1274466822277280365;-0.36009473905638395008 -0.53174638098876947367 -1.1123992583934028122 -0.32774966147901962721 -1.8338866833177394078 -3.1376198416930303203 -0.3412376066796118379;1.3311796684551768788 -0.58491410714111491398 1.3960574027780030537 -0.90943097833006503183 -0.58175440392118304977 -1.4671857395082918973 3.2064052360472818926;0.24588250524505456118 -2.2054651432797678368 2.3331403609110381581 1.4000551164374754887 -0.44189221310320470071 -0.7002959016155789751 -4.2910953390217363435;-1.0909539941649237793 -1.4164696426375209537 2.9607095167217698517 -1.3227592814352171491 -2.3806604887400206749 -2.0703834835436096817 -1.0402778938098873329;0.3382169604202734936 0.087029216970178235879 -0.94576930237804579527 0.02005407146810493646 3.2464809088378454582 0.78542837105756113658 -0.43683781964275680831;0.066390833937270021803 -0.26208436281668484602 1.5025250870112882851 -0.3427940349200785386 -4.9119115511966882082 -1.0031591520898355263 0.90127620003093478207;0.33449785996555275736 -0.31560456758173394132 0.70800362599551103227 -1.0165537570013924995 -0.57627855808732031395 0.52449178229774429116 -0.013825276380795383999;0.99322224221334320937 2.0563407005477443157 1.0209629576224930059 0.12268692561320858336 0.12796096591395839104 -0.4754775237817017719 1.7314767275175642069];

% Layer 2
b2 = -1.8930414659343632877;
LW2_1 = [-1.6335692491234827806 -3.6430307605064449206 4.2358931153491976573 -0.54790310404378916864 -0.20605896553141886396 -1.045073642194823238 -3.7206559705840636099 0.25416845992140452148 1.0004507001349509654 -0.38518565133880400264 -0.55262823243670111673 1.2874880193333053668 -2.4767778224579073409 -2.0416139711155345537 -0.84300982439244820377 0.28421916840114774816 -0.4683103673330591965 0.7607822264723738126 -1.2505092332254819087 -0.96051062797344355371 1.2697772967947291711 -2.0234988986642408193 -1.7724591325802450115 -0.91197873067477264186 -0.48678537221169188021 -0.37423294250414818052 0.44919092107930125302 1.6633627667770380931 0.21853916074904872202 0.29238997336565336616 0.25099513903223596989 -0.38780624000291508136 0.56585964355628137756 -0.48528331664373841603 -0.19422589314206722966 0.080139961475654539202 -2.0816553015283627559 -1.3986862499099170787 3.1024377083319532034 0.52180064848995200144];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.568038194888224;
y1_step1.xoffset = -0.04979;

% ===== SIMULATION ========

% Format Input Arguments
isCellX = iscell(X);
if ~isCellX
  X = {X};
end

% Dimensions
TS = size(X,2); % timesteps
if ~isempty(X)
  Q = size(X{1},2); % samples/series
else
  Q = 0;
end

% Allocate Outputs
Y = cell(1,TS);

% Time loop
for ts=1:TS

    % Input 1
    Xp1 = mapminmax_apply(X{1,ts},x1_step1);
    
    % Layer 1
    a1 = tansig_apply(repmat(b1,1,Q) + IW1_1*Xp1);
    
    % Layer 2
    a2 = repmat(b2,1,Q) + LW2_1*a1;
    
    % Output 1
    Y{1,ts} = mapminmax_reverse(a2,y1_step1);
end

% Final Delay States
Xf = cell(1,0);
Af = cell(2,0);

% Format Output Arguments
if ~isCellX
  Y = cell2mat(Y);
end
end

% ===== MODULE FUNCTIONS ========

% Map Minimum and Maximum Input Processing Function
function y = mapminmax_apply(x,settings)
  y = bsxfun(@minus,x,settings.xoffset);
  y = bsxfun(@times,y,settings.gain);
  y = bsxfun(@plus,y,settings.ymin);
end

% Sigmoid Symmetric Transfer Function
function a = tansig_apply(n,~)
  a = 2 ./ (1 + exp(-2*n)) - 1;
end

% Map Minimum and Maximum Output Reverse-Processing Function
function x = mapminmax_reverse(y,settings)
  x = bsxfun(@minus,y,settings.ymin);
  x = bsxfun(@rdivide,x,settings.gain);
  x = bsxfun(@plus,x,settings.xoffset);
end
