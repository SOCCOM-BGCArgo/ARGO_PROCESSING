function [Y,Xf,Af] = ESPER_TA_3_Other_1(X,~,~)
%ESPER_TA_3_OTHER_1 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:00.
% 
% [Y] = ESPER_TA_3_Other_1(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timesteps
%   Each X{1,ts} = 8xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 1xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-0.999999984769129;-0.999540494077341;-75.751;0;29.715;-2.05132667284751;-133.803853032615;-0.1386];
x1_step1.gain = [1.00000001035699;1.00022983623563;0.0140762793578401;0.000302617642608564;0.251952632905014;0.0584431850247228;0.00432440760873659;0.00826143816767911];
x1_step1.ymin = -1;

% Layer 1
b1 = [-8.2941962720667827824;8.5982444334425363053;0.22474549093968959612;-0.72556303477749417574;1.4596248572022625378;0.75621345875367085299;1.21621338197635942;-1.1952038616945772631;-0.37682382440934664158;0.51486312944621348553;-2.3525453856960378829;0.750424564893550472;1.3251168254506970889;-1.0914064019401070738;6.4922411485311908308;2.1103777252980799872;-3.2284179940416715127;11.049201778187873302;-1.4704437453653451673;1.2574394308838972023;1.9484108048827650794;-1.0244701771431794501;6.3146108978074009244;0.13133500103440678974;-0.68461206866564872087;-0.50195242627123382029;-0.1613793472539784557;-2.5057962657804693229;-3.0139414052297004432;-1.2448919884276874459;0.5280486340914523069;0.99190264681757145393;0.97442140014933664727;1.3761495700180588209;0.7325544735718551026;0.22126388828899770655;-5.6504586630839552797;1.6387590183760485374;-7.2717525143933920972;-0.38089808916581258114];
IW1_1 = [-2.112619273395891728 2.3711229518176342168 -3.8309835038847177557 -1.298455729220997279 0.46872796569729019422 -5.3709592214438934832 5.8119179513694962935 -2.7970104680821523679;-1.254065021982371686 -0.080633214276320849812 3.0514812063865668357 -1.6255264924546213301 -14.327163546115448511 4.3226852731350113856 -2.1865891827767494071 -1.4640047203619341865;-0.066950100061980488686 0.41082761199782980022 -0.25153024137852242337 0.42248163681805039893 0.024582634043050399092 1.1766102264064550731 0.58847207225499276184 0.89934458074787815285;-0.081878902994115307368 0.095888263620824767952 -0.074918730214396753175 -0.075354564472684654186 -0.4868976424913107337 -0.44338501786875594846 -1.2283453925524479011 -0.085169378808806775472;0.51724135703557605481 1.9123999139690608384 -0.076302945510270170359 4.6670589428219688699 2.081914066845293565 4.5709203462084682101 1.1381794702536720543 -4.0440972289124568206;-0.24249563521970632363 0.24401155708179858572 0.19142090210576034881 0.47188863365829403218 -0.26143689016824639415 1.560402318848565173 1.5305167495646225806 1.2532217227821071948;-0.40869403003670790042 0.31071450263194777541 0.33321522125161140426 0.098240555106295354704 0.60947162373347052444 -0.36790652692476943919 0.41398425223874596446 0.61255874914411378995;0.38655122437132882762 -0.23653211804582729982 -0.27510295452632166358 0.026720352872609761524 -0.24826822233005035878 0.39008785697174985785 -0.42419662904814398052 -0.65710979810992298233;-0.01127707768458135508 0.10044091031972181283 2.1860936624575084686 -1.3327007292390984361 -1.4636798042268885922 -2.2967339403760238525 0.013512955345042756625 -2.1474385590347004893;-0.65902332118120632298 0.12725421947925782962 -0.42911950932797721103 0.20369750939756312014 0.9469434188897538629 -0.87721564427557485732 -0.76163714799181925219 -0.61702139086958573966;0.14089654819863059454 0.36803148465128576738 -0.69677390924664961425 -0.45927708064404065835 0.24056306368195387591 -3.1369550702866466629 -1.7457062849466236099 0.12633013920592894852;0.45761867675816059187 -0.043128650616357143421 1.8977233960080459596 0.48857527306246062437 -2.2244976117456807252 0.076834784989102969388 -0.82253120334533591507 0.30910525131543986621;-0.062037167094886165597 0.35050128081800635549 -2.0209464495301450171 2.8195383858295359758 3.6976568878886295622 1.4590130514717034593 0.06170001757362335304 -0.28754253656287237506;-0.0483665563730857756 -0.12761928482522336137 -1.9006718367795067959 0.95070036987317807498 4.015421247109278724 0.81156913688082699121 -1.5152566143830974443 2.8970882609584975853;-2.4248676658126306904 2.9796664782784505476 -3.8918388864922155435 8.08561958535464953 4.138459565640245863 -0.65353959478192669152 3.2766253307824975138 -3.803637174294820511;0.45024257601895317782 0.63759300315145639448 -0.34073277956666458977 0.9459521289367665009 -0.58958618097295922844 0.9850763683904580903 0.46102325567149821106 0.83070741502646849863;-0.34846653751184147874 0.090269723043269134433 0.2873040046496127542 0.88871771548549272079 0.33824917647037361368 -1.6115802217215360326 1.652682995460425186 0.8623231295164593968;3.3287445397376105305 3.6331652518953267972 -6.8945660159140533807 0.28891274741038669438 -0.32612174741590904681 0.76067637948745214782 -0.20547166814555295478 0.86096260076673103434;1.6971548730703547303 5.1647981920043193327 -0.61443854108864870422 -3.5124044561823937194 -5.1908040902305270237 -9.6898236549177561017 -3.6639379631370521118 4.6102668839077836083;-0.10502510538473316948 0.35192637232635837208 -0.28346494343898381851 1.131261731636639345 3.1711766151960341098 0.1779219306084507124 -0.40859049437908029345 -0.47669987826371706952;0.25231018658250914388 -0.33252403903262800489 1.339693799987133227 0.24955656690339908521 -6.371347015834120775 1.147842600249233902 -1.4669848972477255167 -1.1754587481867788679;-0.85179327530710058358 -1.404728203834464173 -0.64773857987050664775 1.3095610808490891763 0.29382426454486143541 -1.2153768512256712953 0.94459142502233262029 -1.2577922949715598655;-2.0721181009655884608 0.20957924027961163937 2.5215651932768694543 -1.1914442505867439692 6.5158693924818082976 4.1742226603896641635 -8.586800337684781681 3.4066570830743421183;0.92040902967849025362 0.15799013238187167851 -0.30197962277564527334 -0.94906288839350194664 4.4971850980613767845 -0.32582392301978396754 0.00074847204909572750964 -2.4730240041169189169;-0.72659886427360664829 1.4748312575973514527 0.58176500899267069578 1.4081424972892420122 1.621010395113719138 2.4556964234232654398 2.2028092187396404178 1.8969654308114400898;-0.098208892653502494197 1.1982022538068088924 1.2347484228940308615 0.78181835549914124872 1.7166677520052004979 -1.1563900836205491807 0.076898722515181552883 0.11053390882138108831;1.2846411504254089397 -0.8794421419247800964 -0.81719872789717706851 1.8677618103278554873 2.1421444833167964106 -0.16587603886019308685 0.21598165226941878903 -1.5968700174875400943;-0.028637864330167885074 2.6226757330963978987 0.95580653585150987528 -0.071489594708845166648 -0.35761235509225963725 -0.29804594208625040386 1.0945167692170780782 0.45366870941697967634;-1.2335219455917600939 0.21732607007527801612 1.8969810956810353453 2.2200234700479080097 6.1112822560606412381 -4.2113711134597160424 -0.71108274359185952296 -1.1906955004606076809;0.93565744340990064654 -0.22408473549710403172 -2.1821431665899604369 -3.8912586872114434122 -0.5408317319077010632 -1.6569848984539135017 -0.7940378493931727899 0.84946498969146644331;0.097726764063874663191 -0.13226883027908947765 0.012645022777855924323 0.12139286431682196066 0.5173291392640106956 0.43395751818899791497 0.80757912346787574798 -0.099945737883478630303;0.40750019976252005449 -0.590484139105845518 0.92580654312453769528 0.40057206014797752669 1.9890797086761600188 0.99042263135766117621 0.55963121779849100967 0.51327745825305270966;0.41760519753833391654 -0.5686584257730241676 0.79974879857087788348 0.36463746232358984489 1.7710636232071284013 1.0566482564349695128 0.61892404106282639109 0.51398067387574741804;0.41549439775276764131 -0.41161168864144132007 -0.12327154397152299825 0.094026480936431747204 0.23227501365049285642 1.1811223447302756639 0.44896843314947965098 -0.26474115060171016234;-0.011348847511795698284 -0.16170364109873583414 0.27807164817842416893 0.3730596754358816991 0.50972314164032772599 1.6759854912907155011 0.24596235774396352225 -0.78653604849105651642;0.29947052424930908909 -0.57621270398530177737 0.37864602856587076252 -0.30397130637560582089 -0.88224996936318056395 -0.82361043536340450277 -0.25987096443096424636 -0.93331734075352812408;0.54280231216533958438 -3.8408252311536346113 3.6786234626675069492 -0.31734193925641979162 4.6058146255049905449 -1.1298092442081533182 -1.9054729410512825094 -0.32618686100358279534;-0.53176459143926779749 -1.120028201564707615 -4.314625148872152316 0.53384773390924555869 -5.2908572975282961437 4.8888405850607439262 2.7732607781244174028 -1.0079657284138652251;-5.2505993195728759559 -2.4217884459215310855 -1.583471118763578378 0.66698631710191502719 6.9974019433545944224 -1.9645637058090530402 -2.2860232489943923539 -3.32061288037124136;0.097984523864235617552 -0.46352375392401651277 -2.7105575601838318711 0.95071432029390245955 -0.20609382303193424235 -0.6064278731981043169 -2.6780647093647012014 -0.58579727874879394811];

% Layer 2
b2 = 1.1342788502252303395;
LW2_1 = [-0.011660013506319035631 -0.031002927932010270612 -1.6394378530045761178 -1.8941512582726593994 0.0096850327132229577926 0.40402544555571717755 2.2261987441404880173 2.7236464374836497626 -0.23054922158710663971 -0.58683022543413609906 -0.13348343644107762973 -0.1722911315588725456 0.062988660820144704022 -0.244951054896030751 -0.030408717091692627488 -0.25086059073749217285 0.62126308579643907404 0.05490084125934878434 0.0075283578636073326257 0.32664049212018181345 -0.056767861839405658886 0.038053519204732544345 -0.015897024748817496775 0.23020438644280699636 0.059157218296389944834 0.12034881334498327699 0.040199252194925767789 0.076210908213974401515 0.052513647881553784935 -0.028876517822022610554 -3.5573441796399589698 -1.6414120011986819136 1.9118265514900361079 0.93233915652113508177 -0.33717891756017942528 -1.046276040647508232 0.08331009519882667802 0.034200228060098841754 -0.013610255869589266836 -0.06792941087186626925];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.00492610837438424;
y1_step1.xoffset = 2075.6;

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
