function [Y,Xf,Af] = ESPER_TA_16_Other_1(X,~,~)
%ESPER_TA_16_OTHER_1 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:05.
% 
% [Y] = ESPER_TA_16_Other_1(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timesteps
%   Each X{1,ts} = 5xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 1xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-0.999999984769129;-0.999540494077341;-75.751;0;29.715];
x1_step1.gain = [1.00000001035699;1.00022983623563;0.0140762793578401;0.000302617642608564;0.251952632905014];
x1_step1.ymin = -1;

% Layer 1
b1 = [1.151899937676278185;-0.46896510169269411961;2.1641615568847467621;1.0912440017713371798;-0.50210898482510379903;6.5129725815126935373;-15.133501133452510246;-0.17087289257144194776;-2.6918396092473217962;-1.0084344964041731529;-2.2479228690243853173;18.211952946364124983;0.21059371527247003342;-0.61560468131123136093;5.4247567434864940594;3.0176069839648036464;-7.2904037493024862471;9.5574624929988960531;-1.1878750782543368025;-17.02049275863872424;-0.058478942353391592801;-0.85053508733025229294;-14.076961370804619733;4.1129930989639351679;-1.4351307625949105695;2.9984979786158096893;0.51477506616240809034;-0.44730014970480957981;1.3132549212052182508;-0.52496611998011633027;1.6603408631292666531;-0.31179064109734944132;-3.0442861106700163099;7.9674400669057723334;-1.4213523922113031617;4.0314448177767481241;-1.0781315895320318088;4.0971694479617122298;19.623472674013630268;4.7767747716050141804];
IW1_1 = [-3.9832018134160787071 0.21246535375580657967 0.70140201834917337287 1.2748637447647621634 0.72002497518716768177;1.4370279999274422345 -0.13502069419375289927 -0.093571288206000038934 -0.32269511729833671998 0.44388364810148289719;-0.1814540976599762212 -0.42237437026493679371 2.7409991729085203538 0.7962091970855934564 0.97568645030173584143;-2.0240517654573157991 0.73778003722816987153 0.02431091662127884967 -0.62793222543633187005 1.2154523096150482431;1.3690704395717556707 -0.29772679692333514634 -0.088194469784423273895 -0.27652791390861419796 0.57364924460958199859;-0.077989697925543707835 0.70227429795798812506 6.2526775153266749285 7.5897329114368892178 -0.38182836028196653899;-0.96755167887382897973 0.24098470309418992819 -1.7979420159193841844 -16.528045953357533904 7.1382753075291613953;-5.3859396057706563354 -3.2327760243241034921 -1.0544131787129285005 0.30119114391913548445 -3.4952599244441526771;-0.14124788447794411894 -0.069895605493910883954 0.48211793414490083931 -0.23462785954934395338 12.616030138472531519;1.8990681748840578802 -0.68831421867223152589 0.071940118827573146243 0.62065940192631718908 -1.3192479691122946939;0.18250893187835010023 0.39828529357128944488 -2.8356379586971498519 -0.77233274490304060222 -0.92576971092333426583;0.346943123135273368 -4.5290167937847565227 12.090853974046853025 23.835378751222194182 -13.734004878677342987;-0.18688015784170561306 -0.25975752919845962463 0.84313980525489962137 -1.3204933690549787162 0.43339352376066064831;5.0828882591888460141 -1.5021226909507499059 -8.1191237522501538137 1.5844233422722560078 -3.8294074351212081808;-3.8749709606217352942 1.6237576140780707323 5.2177413306530633363 -0.68547999289503136477 1.9782605959912686888;0.12700944535482394659 0.065079976637258174499 -0.45365733680356512547 0.54725356611432385279 -12.693295120885105476;-1.2308761377614347854 -0.013155279613713601397 -2.3126096624555239067 -6.7954535452086490821 14.874109890053222927;-4.9499216336503879532 -8.4972992345763831423 -6.441793180641291805 -0.18156140985237295693 -1.6710612738895440632;-0.80812570374382752814 2.1341161748513224339 -0.42909284568473998744 -0.97750419979887859245 2.029727383970516108;-1.7338870645672588289 -0.24855356602322209447 -4.5055009421536480474 -19.501814600358876817 13.943447928734549635;-0.14996803972405733929 0.14152992197510350092 -2.4643525191935964358 -5.4721767170690096549 -4.2006444932493511857;-0.4552737817982064894 1.7188441336365194445 -0.25756293185988521355 -0.96497961990782787911 1.7578127822340585684;-11.513291623580672152 16.20695249357155987 -11.631206103043913558 -0.53930843540957185134 -13.562124520478869627;0.95806749847560190858 -0.039446910031875864755 1.0103895331795382795 4.2960795951317907182 -8.4017504660321584709;3.53080158176039971 0.86672221382517766752 0.11058394652885225606 -0.21056953080982485083 3.7902362895458319159;0.022643435460580060153 -0.27612595710439929997 0.50332435615105086679 2.1373424716383198074 2.2295975508719290126;-1.2891465680263398674 -0.22685885815994183412 -0.16261769563322583942 -0.11048775260930257525 -1.6123191420495039949;9.3391685378356275038 8.6734860224691541219 -10.332831981215653272 5.2857457789796118419 -18.541048679428556056;-1.6610410813923490725 9.6192354269160276203 9.1692213935107886869 -5.8568969721236072701 -8.6924432676322371805;0.086620703984762181804 1.0953401830025044639 1.4968599554328005929 -0.53271512479678584384 -0.32554606917036676972;3.3666536478588007242 -4.9391244726059069237 -3.3029629186425997212 -2.8917391463092481096 0.46778018448651215166;-2.451615677770620394 -0.70239233055981609155 -0.27058979241943831306 0.0087333951580123320851 -2.0231753737200115495;-1.3312629596480947125 -3.1432668716659279973 -3.2298545818963133591 0.80436193724042492903 0.15160790290139994352;-1.7509770080178816976 6.7581520360381865942 8.9757329618739998267 -0.46250293905087763724 3.3842550119030301303;2.1575117885544528029 0.55198222014794684487 0.23106081518725823942 0.10232187870912479166 2.3677868955005836327;4.4314739228385819914 -3.0533828979684347615 0.55494832145460148976 0.42012952994294877618 -2.2842646087729869819;-2.3700076222051253971 -1.116016210480349713 -3.5829785050270772651 1.3035399968113354774 -0.48403853686470071915;-0.67219203659151949282 0.19340116772098439601 -4.8047896354405956743 3.4323154888990035083 1.7732225080448222876;-7.4036850388370885057 -16.084759612813058993 5.6407238870332054859 -2.7864240597160745416 -0.48231428412571819564;-0.46682259224470723513 0.33452493175496689215 -5.3607986197020043306 2.8564097743638057203 -0.30051615582246998803];

% Layer 2
b2 = -0.40376954569246420457;
LW2_1 = [0.15859089304067969239 5.3006079739937677076 4.5747736321674841165 3.1269717189832388904 -4.9064017796896788326 0.051058545518501162463 -0.049410310052850799079 -0.07598945044271573912 1.9913650411030989762 3.0616981473747446607 4.751911895770515315 0.024478382058322541159 0.73878934665203810361 -0.046310511871946291595 -0.15747470176143965603 1.9548431855826267078 -0.12083290016056831562 0.085400133257734642256 -0.48145121811457697314 -0.044080799297562686667 -0.066179948478742783857 0.59645941340725561197 0.015656465795238234257 -0.15933565446117506359 -0.62815797823948660294 0.32665244942461674205 -3.6895274222488918525 0.022126643429985999589 0.027706877859218295018 -0.37906144709991373976 -0.017786738277593895641 1.1488216169323179905 0.097602766754019923412 -0.075002099244751666007 -1.6666286093679627545 -0.070709858480851511819 -0.079560393449057131843 0.29367340245836798251 0.039143328021204827627 -0.29343602304896204691];

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
