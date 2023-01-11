function [Y,Xf,Af] = ESPER_phosphate_14_Other_1(X,~,~)
%ESPER_PHOSPHATE_14_OTHER_1 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:31.
% 
% [Y] = ESPER_phosphate_14_Other_1(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timesteps
%   Each X{1,ts} = 6xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 1xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-0.999999984769129;-0.999998793562952;-77.961;0;29.715;-0.9];
x1_step1.gain = [1.00003298043486;1.00000060321889;0.0138606862425759;0.000266418009857466;0.193068829037552;0.0415843642790311];
x1_step1.ymin = -1;

% Layer 1
b1 = [-1.7015569442873224038;2.518195762054127318;-2.0707999284862022016;1.0602605384710754866;-1.6082744260514980628;6.1827492723982082623;-3.2378249173741524558;-5.5591300042196429843;1.8857814745923342503;12.904971696986963892;-12.934429275422292704;1.1916657864541828626;1.9307448929545407701;2.0029929367719141453;2.3010146618632885129;0.56205228168367626918;-0.49360831569223051485;0.015894550062008974689;1.4141576340825217084;-0.63690085629668624634;0.8520335957716455022;-1.1106322169349354123;1.0259974430716174432;-2.7963130246019232139;-0.37319240815254606236;1.9723189361512265538;-0.69124907292622783661;-2.3189219321368752524;1.3198725690507009389;1.9111323312425723131;0.0057717718081127349083;-1.4488844151368127644;1.7819545015973725643;0.99386187876988318912;-1.3638973843238666461;-49.870569650176008736;-0.78214165016774483963;-1.3565536153085262505;11.720462059467887883;-1.5178380290401392028];
IW1_1 = [-0.57994645157632773724 0.87502415754585249896 -1.1441505898459651824 0.94277806360167926147 0.74854089461197803601 -0.84641809164234815821;-0.12128997582628070484 0.012195988114732683044 -0.96019090086357550273 0.39564474357696594042 -0.29977900305738214293 1.2804041653630189312;0.10549085725872135721 -0.020262115242623508771 0.93872879832703381986 -0.27042856938771892006 0.31917973099749541444 -0.87153439992873815712;-0.24214972220547772985 -0.31282595853827255272 -0.66985965307471595498 0.35040728311891577995 -0.33988783699274610628 1.1600108370894204324;0.0092354580091084347016 0.19297680455710719549 -0.36936655308984622792 -0.039884064592396302473 -0.30209087339785584136 -1.4404039460732560585;-3.1063580034937934293 -1.3912555477558168349 3.0957313480909176562 0.60128731459313822238 -3.3714147576104651627 -7.1365871240131619757;0.32709471417334096133 0.47270838629511596762 -5.9833914850555389719 -2.1813852046093682446 11.565488115710282457 -0.65407656974914762937;14.73968834902002456 14.346134471768877461 2.5653331470866911346 2.2167525767612428211 -2.9175678001330940248 8.4740751205234143839;0.059334095511230443287 -0.39568514791468123448 3.759156435382278616 -1.595338293452140066 0.22951531266612834603 -1.1693221828054125577;-15.378394232068041347 16.784385386561030629 -2.5318740826632883945 -2.4591926937269334985 3.746537791347478219 -0.11741314377118181045;15.362081825511630484 -16.752907227350760877 2.4918535492564606315 2.470736647637643113 -3.947059736248927031 0.16928239088664381873;-0.01961100166874846229 0.4466069848202559589 -2.2349530669289445761 1.5698006036604785773 -4.0521786798975965027 -0.46477291554286898068;0.014917542566642259805 -0.39281577492530761475 -4.3349767687338705002 0.49998635650769374106 1.6252826623138092454 1.3796060471636859557;0.73340418729660261654 -1.0465100854800395425 -0.63650063151018010377 1.3476263493683349193 7.0596455005701814756 0.7912864961831926891;-0.10355440009796045009 0.035048757676561077823 0.5099112624137192018 3.4940469469694233062 3.9724257896540509094 1.9077861183787283572;-0.1187807748298435051 -0.066813932982689283602 1.4044566203394168724 -0.22977040235137946578 2.3691488048385784815 0.55170684360211097452;-0.16648064407120338748 0.022523152635189568754 1.6523664861598803544 -0.17671774639941936758 -3.7565015715059253409 0.013763209149391367539;-0.17517064436521784243 -0.0066102086864434309726 -2.555248657364539433 1.0030642141747463381 -2.5658894660710136471 -1.2210054466057593814;-0.9543639035547982763 1.1435010430960261818 -0.82350432668320394658 1.386765584360493353 -5.4966340301206555452 -1.0058296655722238633;-0.14812067735437256522 0.043266220933384558212 1.533345847791013794 -0.13933830231284849432 -3.4883448901443294687 -0.010336453110154889337;-0.29812133360226095169 0.46517763315752758135 0.47775822668301132223 0.024253410840996504627 1.8710591381730561356 -0.10925319167063056447;0.032098728993464567538 -0.54043023884755969632 2.0971872086071163821 -1.5790672214776722004 4.1578835246684979765 0.51833394394779952474;-0.13664797791578525699 -0.25929011681548241874 -0.59123037053287286202 0.33613087518600415216 -0.36478955386089684243 1.0900653592963513994;-0.14715671673058841606 0.3854740486249030873 3.1193835780432195293 -0.25121344431035769196 -5.9512404069557094388 -0.39349062775849397156;0.096994733835203267303 0.90728980131330305703 -3.9522307400468186245 1.9047945038366467241 -1.4353909742507764413 2.3704902842299930832;0.25832226300354271809 -2.5666859294310753192 16.007761204810730504 -0.075035872877424275584 0.051737688231804676575 -2.1328860828741209588;-0.04676751359725794116 -0.8975889748031128379 0.90447630371977227881 -0.20705827856395700448 0.87741240179568147362 1.0678250821011754379;0.18124962158366603182 1.0142021656349862191 1.5024875165438233982 -2.3899092628817357742 -3.0029835293158537723 0.45247378839935697981;-0.0054847949297368529115 0.38457580681818881452 -2.4424172143301681182 1.6255931626752813113 -4.1594397272174150615 -0.44132380104570251289;0.018019456954739460769 -0.37678231172940279681 -4.3745601275376460038 0.50237522645771004459 1.4600294524248675909 1.3568357228613436849;-0.16422806697844238855 0.02462765854062761578 -2.6604366711992031291 1.0687488382754084526 -2.7766326734826503042 -1.2646373505973453888;-0.068608688034413126133 -0.026210439170342183374 1.400195098081598255 -0.033991989736319574988 0.55876074567166633678 0.79028087027819382548;0.094329600088723902784 -0.68686684025265387632 -1.060507638510391093 1.5892997403974953308 -1.6659951258676990449 -0.39874040601103000192;0.14249516680704552019 -0.066630373157192598899 -1.6823638639836981934 0.081456436419293273832 3.7752212157947067084 0.044219742651856776272;-0.34581055232057938964 0.61892547225217575679 -0.90295355140293898089 0.68818764259981912712 -0.77383872322025681889 -0.82792947206719491593;-9.7276220206185506356 -26.685709454713862243 3.8281358760998900337 -22.947810037969642138 9.6149036379205359282 -0.09904291362128309073;-0.24744336106210096382 0.49732418133802980265 0.64947565828032383095 -0.06842038967462682475 0.9513101030108486178 -0.99576528779658468604;-0.032776132837728282066 -0.039305468579660281236 1.1956574104019546123 0.10101378121832460533 0.15821610962345938889 0.92767700176996525929;1.4291717152235396959 9.0499843845441905188 -2.9973672811037550545 3.2111581677100446974 -4.879895237881050285 0.3114763064527851788;-0.059435937095075754011 0.34768785793850687238 -3.9402284526611781601 1.6178679781852634001 -2.1322441459706835687 -0.85148073686985414898];

% Layer 2
b2 = 4.073225417096367984;
LW2_1 = [1.2749747482645163998 14.429513542729846876 20.187780573725689237 -8.8527972239222894046 2.7151913849667206335 -0.023957864023534478776 -0.032298958413436551496 0.0065438184446674069195 -0.29328650751978485589 -2.8894347511738893708 -2.8996469778001352502 -10.513259181467049075 3.9794045734798548786 0.10079782836547772351 0.062778562482542182854 -1.6310222021266511661 -16.885736990456479134 5.2513959666091576395 -0.15901439158121535855 28.776101993599649376 0.90210400681722768823 -4.1242653842607408521 14.624213297629692576 0.45309106279673005968 0.09118638004915249895 0.011322127261495281564 0.34229386656264759736 -0.14180854307303178419 5.9482904971782151193 -4.077332036719870878 -5.2009683249307263608 -7.1047433343101324965 -0.50768180231591131246 11.567083049748955403 -2.7164591319920150525 -0.20643624748150277148 2.1134824489811188108 6.3812826309182906215 -0.3157270590291603729 -0.38506592612874240844];

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
