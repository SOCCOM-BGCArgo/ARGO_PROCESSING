function [Y,Xf,Af] = ESPER_oxygen_8_Other_1(X,~,~)
%ESPER_OXYGEN_8_OTHER_1 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:46.
% 
% [Y] = ESPER_oxygen_8_Other_1(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999998793562952;-77.961;0;29.715;-2.18999082656054];
x1_step1.gain = [1.00003298043486;1.00000060321889;0.0138606862425759;0.000266418009857466;0.193068829037552;0.0582073296585987];
x1_step1.ymin = -1;

% Layer 1
b1 = [-16.814191841458644205;-2.3717483287991871777;-0.28300872782377056724;-1.679565929534555746;-2.5180210637954969499;-0.64393070585993550381;0.92560643554195642935;-0.94690114989674056822;-0.45672748264140694552;5.5323175490613980543;3.7019414612161258837;-0.29902907991230631701;7.9840811306652437196;1.8165983424983118688;1.4441259177714609319;4.695751757987296493;0.066217615788517619757;-8.6208020454034404167;-5.8618335929261897732;-4.1769838855991698878;3.3082629135502914508;-3.7527365425407590571;16.008392761415350947;0.38662178274872704975;1.1054379869603490061;-2.3654812223719248188;1.5416512670771076454;4.5837093620612439082;-0.83992136717288934378;1.7421134634336716829;2.331064536862852421;-2.064361198713645873;4.9849625845532132828;9.3378067292805262412;0.94880858377394639014;6.5493822891248409945;-0.66830839183127854763;0.2034138746296249578;-0.7753450324496485857;-4.7905932948577971686];
IW1_1 = [-0.43353042306541056394 9.3053214914700355109 3.2796458370280969774 -6.6188900920758939606 7.6261168220469794932 -3.3293762582155275176;-0.37178879020677563227 0.54071747972392714132 -9.3225011470922947154 -3.3538102735858532633 -17.379944601234711854 4.2547721426142235401;0.37311268433385724741 -1.0176504668878576521 2.7602683788695578926 -1.2854151807255884687 -0.5074129822348261154 1.4848053763967203889;0.14878011070715760344 0.052951352002886045711 -1.2241996453064196704 0.68995506936139372467 7.3046826642575570787 -4.0368499636898897975;1.4223851806327170166 0.94037156578782921112 -2.5466047876852111465 -0.11305232308951983222 0.17725985080635009439 0.97201413607754827417;-0.37356100005629822025 0.49698412592308294222 3.1795487807222437837 1.1816184660975108667 -5.1766876248007411832 2.8447020494176373795;-0.78792354585988111371 0.010271766137513189954 1.0115959931028353225 0.36912981832629199186 -2.7114044621021426984 0.31023836809644905754;0.10098332431663377406 0.48079463358506657755 2.0909729592926691311 -5.9495164865795731046 11.654750865221648937 3.8392938940844825169;1.4552902764189754592 -0.12136555879132704394 -1.6607607691251125015 -0.92292899652944437872 9.1212922657575763452 -5.1826857398940378019;9.1006440612292109904 -5.7099585846114466747 13.65392293703049198 0.48585260623925674572 -7.089719163098758159 -2.4130500200362696717;-0.64065256979794749093 -0.7883634015309283205 2.112410399428918506 0.68876830425490298548 -0.55212049339396473968 3.1272957462315327248;-3.0085313492160281612 -3.1472531805532817373 6.0128324852285892632 -2.0097791661648498618 -5.0831663027311000391 -1.0154999482734285987;0.26992373575456335422 -0.43575080247482445017 13.328108655952783934 6.7837465266995273794 -2.5475142164048274473 0.1428115515810788605;0.10632696958688946098 1.0348698718055617807 -0.6338099395814852377 0.34024017816499252387 -1.9065415994617107565 1.6787552476458122985;-1.0632253760269536613 1.4385633186742063927 -4.5981086488854474581 1.9248427576367701519 -3.549516449509827698 2.2273589109627645044;-0.73674473505953386177 1.2456796149389888217 -4.2399565651445589864 2.1515826526540640096 5.3790274173331011909 5.2240461254580452533;0.174867483058987222 -0.12206428023595798649 -0.85136648617059029753 0.11948084894341573503 -9.581645876665858097 3.1566371513109934277;-0.26156146900044618953 -0.38622576329653246496 4.1221383465970182769 -8.2663090396565195306 25.815873283284844319 -2.9447257871303738064;0.70593984845508850334 -0.53306087977683513923 4.5628985300684856696 -2.4533377474427573439 2.7192814344190274767 -2.1317745556821026121;-0.40123168898999683751 -1.730157357431060694 0.46134212256445300548 -1.0717033999130356481 0.97599875865013008269 0.017481898389994855042;-2.1876104148946664552 -0.37928932096301420929 3.6431850306789383964 0.84299393326669735949 -11.31247791749696141 4.6749449755512149096;-0.8195700500142685252 2.2207366982588081505 12.031547624441028432 -3.8019930865295239464 0.89379470788458059705 0.011248620076070910828;0.37500442817985624178 -8.8688196683641038476 -3.1266605412466592639 6.174535548472259272 -7.2565043147293364356 3.170389936968987854;-1.3448823602485513451 0.063968031432581701212 1.6363891064672311249 0.78661193040412380473 -8.3593368467332691552 4.9320320507291013357;-0.010039340560989672146 0.94701668570811237124 -0.48186848857125624823 0.10556239383317921054 -1.5782140194070171102 1.090288455898881681;1.4822597509167059204 1.0484087739208127754 -2.2926645670383622999 -0.047490120698540132238 -0.021089096131137102796 0.94182303664486488604;-1.5299097236425198698 0.77285851812593764354 2.8601396076806553381 0.43020582328933570926 -6.6657621409446177196 0.76903418851379679033;1.514768709603030894 0.93337658226598962941 15.884959282877817444 -6.6204404442058724101 -15.754999638203379675 -1.2652100158432368904;-0.39910310959237388362 0.80302728776114551401 1.6205069072571605115 0.045523982565292608893 1.4004860143783179893 -1.9733920666620006212;0.54583566269719896447 -0.78872704246151481566 -6.6867968186924571228 0.32697517337041759333 9.1655447590047049999 -4.0922536553240203006;0.15728690983196025321 1.0111806584545997012 -0.66502545842610050197 0.46983425338027001716 -2.0377882400726674561 2.1763236550893023669;-0.2696301382253674106 1.0938727977745605902 -1.0318350129827118522 0.18675712698102037157 1.1508172292036267237 -1.2560593640184971598;-1.4233636047712077133 -0.45578470992480413226 -0.30791779258436891009 -0.78147092348754443325 -3.5305244140080613846 7.3316133944280883483;-6.600733649111267809 7.9606472122448082018 -0.020222708363471258997 0.94009939048666590899 -1.6548533840780599302 -0.23245652374574488075;-5.5244407776938491139 -0.59277414917423543184 2.5438763220180007352 1.5261000126140762401 -10.343658665613046921 2.4744695828910674429;2.8319533028623817827 -0.19505914210211006132 0.27852986789142741841 1.271738904150056193 -9.5425705149467745514 5.9355105452759691076;-0.3720321719423876905 0.47966995853963184659 3.1300150395465027309 1.0564661696001074187 -4.8042237590135004055 2.7231532146615040624;-0.1449767684943717716 0.11426064711585545708 0.75035540466200689735 0.20910368770847265196 8.8469191515893790267 -2.8111767731456223096;1.9031568842888162418 0.19464266773736352278 -1.2296516938824935128 -0.32874431405516252314 0.46991125927142957108 0.48160159171889987961;4.8027040669721259292 0.53645708103644296205 5.6035595167068574796 -1.0085345543985118066 -2.6272880349690130863 -2.7222050629763567819];

% Layer 2
b2 = -2.0252578131885576873;
LW2_1 = [-4.4278429058622013059 -0.081033352197574390341 -0.39886600017258150297 -0.56757873776810330302 -4.1542495059237616317 -4.9367617891975275057 -1.5301281924235148324 -0.10588535485715003248 3.1805401838826066196 0.031400122617037597161 0.39546161058019929646 0.089617985344769057754 -0.083078865617277980471 12.02427759176017652 0.19380078812994075421 -0.13242857582711162201 -3.6469049614210975463 0.082689215642827676445 -0.57281224712473743921 -4.1922460515251929181 -0.18217314144397792064 -0.054263164015346997127 -5.0911966689609862158 3.3984664498631920182 -6.2542702652736865687 4.4236688559953663002 0.35816498026116211362 0.030833269575422646985 0.19926261956031690592 0.2698339841157376684 -7.0874944265201289539 0.85657676796530946017 -0.17643077325885506146 -0.15332365750342935584 0.09701671302467151925 -0.11899658139261481971 5.6208701555286779694 -4.2398743263897342715 -0.72487840091821142963 0.17009283334456373749];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.00454442172233583;
y1_step1.xoffset = 0.3;

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
