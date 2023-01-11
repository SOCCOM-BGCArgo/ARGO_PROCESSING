function [Y,Xf,Af] = ESPER_talk_8_Atl_1(X,~,~)
%ESPER_TALK_8_ATL_1 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 12-Jan-2021 17:03:47.
% 
% [Y] = ESPER_talk_8_Atl_1(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999999390765;-40;0;11.29;-1.99149636107278];
x1_step1.gain = [1.0000000017546;1.00031699124092;0.0153846153846154;0.000310077519379845;0.0705602618858232;0.0604107502911417];
x1_step1.ymin = -1;

% Layer 1
b1 = [-4.4314806654442664;-0.43217886267366551;3.2313314617241269;-2.4923282636774182;-4.6256153710713734;-2.3168783811308868;-6.5713866155071425;0.25200120709109058;6.1502005360040224;-4.3359888817007306;1.2545194289074879;-1.6822213072810439;-0.57833302485583449;-4.5403757549624215;1.4695457653393045;-0.93184758801794432;0.17328110818452083;-4.4733428790426313;1.0299650999428132;19.235459337122144;7.5633355328071588;1.1083162697486459;4.725389187064736;-1.8128283113553969;5.9866319164567612;-1.6799044410659869;-4.3373359244891221;-0.77864092381526373;2.8656819014605883;-3.436772144515849;0.2891566225896966;2.826621923901226;2.6027238945127524;-1.4583633214902572;2.9831382525907206;-2.4309836672413714;4.113111726334953;0.17232392888135242;5.0939772935911476;-2.7687616318653796];
IW1_1 = [1.6868083537679428 -3.7832613474707699 2.6486593992825567 1.5703296068372061 -2.0331358886580002 -1.9162816510258898;-0.34791783589542941 0.5519510474480962 -0.28878670894652259 -0.2982248255463788 -0.2138667089717548 3.2607257985075337;-3.5450063663437255 0.25505730327950465 -0.94570857110338657 0.94154973923473539 1.0004855609935579 -3.7514017764279792;1.2772060643061913 -0.083855961672005694 0.46781519517883646 -0.11213186675641329 -0.12649041031126945 -1.740711348453416;2.8953495570826724 3.1687585316471756 -0.66706806848021483 -0.50230919530827789 2.5652291700038754 0.26679615225710696;1.9663119340389792 0.40826439193928549 1.7253080181712463 -0.5494056305013697 4.9902880710980879 -1.7988966110660642;0.017584435485410797 -1.5145276775800294 0.93839245809225436 -0.3841877040611073 6.2529152864052042 -1.7514321945150468;-5.8928640401756782 -1.891367961089784 -0.37645058271078313 -1.0598671092137797 0.31235166800470548 -1.5708343763747195;-1.9096319518359999 3.6420254811497998 2.5841843289448634 -0.62330425508075449 1.2338396167473087 5.6033271972417378;1.7596636615296526 -1.6276178078987271 -1.4265161754479219 0.77159888084372807 3.3321744046092161 -0.93349737326727045;-0.82392231612552558 -1.6783666480434327 -1.1141099080946582 0.53378930616210407 0.48159548103973732 0.21474552331815089;0.79043400684499932 1.6471539909904371 1.1598461667133189 -0.45409981532143134 -0.38332912952250126 -0.4813294138933793;-2.4675364616215392 0.42782130734956908 -2.3337652769068447 0.96035435127328295 0.0093169815869527887 1.8179652964354767;-0.68340554063616732 -2.4524514081960342 -2.2365661235661651 3.9790004576077123 4.7956328453192745 -3.8409628199488934;4.4407511838045695 1.0486340999204746 4.3841795482406738 -5.7991319352431709 -6.5352211143324688 -2.0234356839921306;3.6902625225189021 0.23262628993679174 -3.05389019996603 1.6712695712498353 0.93078115422397867 1.4135147200182552;1.2426105014319022 -0.97499272199924891 -0.93823823581708599 0.41640464436474189 -2.4327683163932594 -0.26983339464778927;-2.2364849620366871 -4.3159105855357556 2.4239187572682015 -3.0970079813502949 1.6123121359974908 -2.3550427240698535;0.42260191836952571 2.725350494130399 5.9019035841149527 -0.96116267802711275 -4.0665156174850097 4.1783525870223945;5.4728199572032326 5.4810041559989893 -28.870107245335582 -4.156944059346813 2.9755287072016596 -9.8786387912207232;2.821288098280156 1.4313903475568772 -1.1847563262571716 -2.6250509849629693 -1.0422917043526987 2.287958241383707;1.7576235911014646 1.5276018066409938 1.6331439828861458 -2.1500586791283269 -5.6348758832335024 -6.3912953855983252;-3.5297859750165719 -0.40864029054849477 -2.7004280013677677 4.213475556733262 1.0462510057714087 2.2832302797601569;0.54203413408722212 -2.2927964059311865 -2.5820282055862767 -0.18046882862518637 -1.5702809647494507 -3.0815115259842338;-1.6479811357801883 -1.7392341085919421 1.6780844471799263 2.9029754314071901 -9.1258829548146654 2.1382074636089268;1.3980017459200507 -0.11781833258998795 0.85256996081826542 -0.20330359798419409 -1.630731962298646 -1.5912446006361445;-0.38740464760055193 -2.1438465306557495 1.5385300591571476 0.013766630826709949 1.8828585748322788 -1.2079933474291067;-2.0367819659017434 2.7725622272433661 -3.6390416706260917 2.5833455882985046 -0.22982436139954276 0.96466752851323301;-0.045586132281231977 -0.86503427863500726 -1.385438208634024 -0.015338271858497857 -3.5022964005986466 -1.0278170084281826;-0.88156176504998252 -0.083728477975923463 -0.65945562212631248 1.9237776548230243 4.1557935378097266 -0.2546513921769723;-2.3785985800634126 1.9296016438582981 0.34334219459011417 1.1464569448616391 2.8213217488362647 3.9623046636671897;0.99920770959332517 -2.2489194748495072 -0.88377255398821408 0.82047954003253654 -0.038878158112288418 1.6966446004184372;-8.8700649779457166 -2.1982923836416139 0.078407172247324158 0.14904925107274233 1.6971166374456426 0.77225930804773735;2.0865180699608059 0.63702287826475779 0.81750795070839022 0.54511550509190088 2.998245108539924 -1.9208187879214538;3.5017668314678057 4.9216560934020412 1.620559826832322 2.2622557683793336 -3.5367561659195457 4.6390061107606009;0.5521979642846998 2.7821401185062511 -3.7264337817604161 -0.54550254837283019 -0.32520571461724623 0.25585407090180351;-0.021470367197873427 0.65695422953002147 -0.58923758157696082 -0.067546874137377574 -5.5073881925046626 0.26639622399139246;1.2012181245569149 2.6664346368051448 0.67859166420601347 -0.60252709259739223 -1.0691898986596242 -0.017503029786332829;7.5832164210570525 0.7293778632367055 -0.26925524493783903 -0.096559288106522673 -4.5026448644216828 0.84317783640962685;-2.1329322474141641 -1.1697467837723003 -1.6775097190168713 -0.088098061671652531 2.3515256865401852 1.4740810302440577];

% Layer 2
b2 = 0.30594124423349228;
LW2_1 = [0.024868330981178983 0.030002029858486438 0.012603940639056052 0.38992985182671847 0.045852888206829256 0.13292666413371326 0.086727233926532518 0.013373986098861278 0.011208498764016996 0.091464332480353405 -0.55020643484339726 -0.52385896517523733 0.10993705559707881 -0.0079245132845549762 -0.0060452207005803173 0.021109273077987205 -0.22965877558859821 -0.012672081688374782 -0.013050721005481016 -0.0078452332786131986 0.06611455327594562 -0.010198593268428134 0.012337804537035119 -0.016326680224359322 -0.022539829905326186 -0.37929799888877175 -0.05906463682315552 -0.024209050365954034 -0.11574247378155138 0.042657034198150111 -0.0066239510195582156 0.057242818409092303 -0.005912833210191662 0.079862689118352384 -0.021368624085835338 -0.0095511507945620199 0.15588406677287303 -0.11675483512292588 0.0080178406034849704 0.077847212468007299];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.00111018595614765;
y1_step1.xoffset = 1025.5;

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
