function [Y,Xf,Af] = ESPER_talk_11_Other_1(X,~,~)
%ESPER_TALK_11_OTHER_1 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 12-Jan-2021 17:03:48.
% 
% [Y] = ESPER_talk_11_Other_1(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999540494077341;-75.751;0;29.715;-133.803853032615;-0.1386];
x1_step1.gain = [1.00000001035699;1.00022983623563;0.0140762793578401;0.000302617642608564;0.251952632905014;0.00432440760873659;0.00826143816767911];
x1_step1.ymin = -1;

% Layer 1
b1 = [-1.5510528735695142;1.5038658486551828;7.2699703675851612;2.3293497367545544;-0.7373536134792148;-5.1177949464505641;4.0121299255835172;1.3280231450132267;0.88991881235605741;-1.7342344895442241;-0.57966039211334652;-2.2052233869494726;-0.27832026332101401;-0.43688591083479894;1.0650997767300505;0.39107114737876225;-0.84878516596535958;-2.8929626031916587;0.53383432878645898;-0.29042320326233834;2.0196944900275136;0.99033243153546524;-0.70991782537102144;-0.52246243190013919;1.0739552963983001;-0.58358310771691824;0.18700309744077151;-2.7269399709683162;-1.4884289641794581;3.2790443681597239;-2.4714426678313814;-1.3479140550396229;-2.0663207191700326;-2.6826327536307715;0.68453878135409729;-1.6683312554871081;1.5354224134362735;2.1190955744648221;-3.1135898446560857;5.4783844489176179];
IW1_1 = [-0.70736582600962483 0.15359669187009034 0.41903673027769883 0.034040916581251815 2.0097819573262243 1.4412888909416131 -0.28282371230322173;-0.84057868291136395 2.7369403503596565 1.8654053069510488 -0.95426908096311935 0.18267853726733527 -2.4669659197778739 1.2995115988219046;-2.5360227529418276 -5.4777814635837139 -0.42527987074920265 -0.98177690179263799 1.3918125967371757 -1.3082868412447619 -3.1938423658167996;-1.6715350548176817 0.60641259280621007 1.4580392320429152 -2.7260436294250394 -1.8260326010324865 5.9039361982491467 0.099216918193145584;0.12347952223942117 0.10044296395523028 -1.2030693030465698 -0.35536354920611996 0.96775998461146862 -0.16771397503950514 -0.77142164370537347;-2.9190921940818475 0.64234745243615254 -3.5374362854311325 1.4169760893720453 -0.93333170670032728 -0.42719935505185247 -1.7979557246791147;-0.50565763028902089 -1.0162767882783412 -7.8313500722031026 0.40952755711872008 3.3206423234508908 -4.256452567423441 0.41082240290555672;0.61118345668357021 -0.15943468491004528 0.40200463232438066 0.70304033780820385 2.2799587175650422 -0.23274576344035713 0.49802377074658427;-0.20509542944936304 -0.27715868296140805 1.7894188483210927 0.39341173802084251 3.1505518038284368 0.66569046984972913 -1.9987094670530088;-1.7491166251389807 1.1733341164839364 -0.92077607557344632 -2.0150554390569591 3.0957579685280554 2.1423192454592308 0.34922803781227157;1.0863452684524881 0.70364988773736292 0.60576925110874724 0.54722104675223759 -0.53515031302890126 -0.18349283382643128 -0.23996732001183183;-0.43647972657664258 -2.6635660757023802 2.9025727612826464 -0.15851172038350581 -4.3688471884113058 0.79683194044708983 2.0759806300189054;-0.2491813909719108 -0.097534641177032352 0.57082792454924036 0.30856896688617441 1.4918687711139433 2.3134307389970341 0.89058879926348422;0.36421418630313213 -0.49852774107638231 0.0064782860270860172 -0.28566858426746311 4.335949021029144 0.23385931559067913 0.0071496404515493262;-0.85023265597930153 0.5648196251350025 -0.35297435325371923 0.10497976987160526 -1.334507553430436 1.2776118068623852 0.43378500271817799;0.24345392844278846 0.093389537423772831 -0.54269717931631845 -0.18794635618465219 -1.4593759099532035 -2.265444768641538 -0.84767102672835126;-1.1331232571122405 0.62024429039005768 -0.92763926448760459 1.3836263222435288 0.46145084563038513 0.080564533407124578 -0.87846310003395456;-4.6846941847303301 2.8538728008159437 1.2450911900601214 -0.3752159642873083 -0.82643527380306792 -0.49433981227364038 0.90787563078017297;-1.2497735246514479 -0.73488555109267739 -0.72672796153417463 -0.63266738770188435 0.86514626060784849 0.15932233068522461 0.34321326361456833;-1.3138005600175291 0.85138782765515097 -0.096114755158822854 -1.7272985233616873 1.7094626294588184 0.51394257932847121 1.1380707817136102;-0.16727627976597473 -0.51300497093566055 -4.9299500617999845 -0.55258254544552698 -1.020660326929369 -0.62856608867221353 0.82853851673831802;1.1428016390802935 -0.63928369475330504 0.99858083601425418 -1.480658673353032 -0.6850471532689214 -0.26603737516404086 0.95624302728076971;-0.074982475492211889 -0.047487948374007473 -0.13873715918189106 0.7177632847226415 3.5528845997827383 0.34133733070505945 0.6301537739538603;0.18599477496122174 0.21867271044501368 -1.4817081320184573 -0.32660802993847632 0.15409546568398874 -0.17678290382472056 -0.9048049773015967;0.6172539033740464 -0.23492943669561936 0.37737314735295135 0.70816603184653815 2.6270789608305094 -0.070794000772477883 0.26313884742148264;-0.76890907694854671 -0.28033099486215113 -0.83108793829446614 0.25633317671984451 1.1617894809851588 0.32558026861254308 0.13682085759595716;-0.13426523465684229 0.13911144773478135 -0.44589704066734387 0.78566668650426097 -0.73017286914595347 -1.9813192324943312 0.28047770844466419;-0.13610557653110258 0.82202626387373046 4.1762947222804296 0.46577731378257725 1.646226611588046 0.22821986857112483 0.3759630808675769;-0.38115851839416431 -0.027817185096486851 1.5267141660244727 -1.3454332479025626 -1.1764308593640918 -0.17342043989121753 0.35853981085270736;2.4372421652178184 -0.043957696806626131 -1.4506493373193112 2.0489074804122889 1.7790701484416318 1.126579539162275 -2.3501450398013159;-0.15394486536944427 0.81834688312790738 3.8675926697967422 0.42792862420739625 1.6719290907471647 0.13035998403977175 0.38153242359698436;0.24734425683328187 0.24796919122859015 -1.7643854643042542 -0.21802617211819819 -2.0448768146513796 -1.0145740113756379 0.83220069124761253;0.19402644184170412 0.47939200736852483 5.0268887751362357 0.54937807093560664 0.92518144221957244 0.5903964490291429 -0.98612955953940684;-3.4595781567772534 0.69624101348840439 2.4371894435899826 -1.8196449042180312 6.915597800585032 -5.2213809877961506 -1.3777753350029152;0.86376883070610022 0.25249003454491065 -1.548072065523225 1.1029155239697244 -3.0972117227012927 -1.8078184324614619 -0.50170933276411334;-0.33211526998820601 -0.078938064250287307 1.5700873235836257 -1.4617617001155521 -0.56796396511219571 -0.22333786912987058 0.5867491463197424;0.49645351419245259 -0.14825291905102481 -0.41089555641396003 0.31816965357367905 -0.42840924041062184 -0.32488520985138991 0.8699685977926187;-6.2053813683937893 -6.715174136250889 0.11966421381348456 0.20413368257944262 -3.675716153570816 1.0042459476245182 -0.61332212422582577;0.66842279356182543 0.49286438651289383 2.4231227032594527 -1.5067807321579021 -0.86842712638280328 -0.35895869115118267 -1.0635167125150256;-1.0760163082010368 1.7562258290160171 -1.7007625701318896 2.566669637631009 -0.40812678776023203 1.9474392528936986 0.83880248493832077];

% Layer 2
b2 = 0.85447483280724046;
LW2_1 = [-0.28974639240200767 -0.030610864752021644 -0.11585090794064787 -0.038919081693363097 -1.3140510638031342 0.034475448233443881 -0.02674632298430895 1.6824493457311771 0.43527418090309028 -0.073590328866084964 1.6020664158259459 -0.028169238461242144 -2.2796309033352098 0.20788509278289632 -0.15860624628877854 -2.4966967606776991 1.242543546751032 0.035630028497590581 1.2880379178767478 0.13393628184335246 0.91085314209500634 1.1589932491796848 0.25498674229904184 1.0136976472299835 -1.4744509671936288 0.27183757157592597 -0.30526379384126001 -0.95014182394505242 -1.5616872592443094 0.086173468823149357 1.0222365646301306 0.66878424711751328 0.8166840444384994 -0.013063965225228867 0.08973108789661427 1.6175194156876047 -0.50022777180779943 0.012504976583633417 -0.17636031986042364 -0.13739107762076913];

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
