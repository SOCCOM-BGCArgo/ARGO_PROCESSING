function [Y,Xf,Af] = ESPER_tco2_1_Other_3(X,~,~)
%ESPER_TCO2_1_OTHER_3 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 12-Jan-2021 17:03:52.
% 
% [Y] = ESPER_tco2_1_Other_3(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timesteps
%   Each X{1,ts} = 9xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 1xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-1;-0.999998793562952;-77.961;0;29.715;-2.08184650291172;-0.9;-133.803853032615;-0.1386];
x1_step1.gain = [1.00000000274156;1.00000063368067;0.0138606862425759;0.0002963841138115;0.248478071810163;0.0588359962254402;0.0419287211740042;0.00432440760873659;0.00826143816767911];
x1_step1.ymin = -1;

% Layer 1
b1 = [-0.83024371453096679;0.063482079712840631;-0.97714670059844344;0.53404852215356091;4.3736450455660174;1.9785035225104746;0.26153755225480846;0.46563120988079976;0.22708727263439979;-0.84281563698464113;3.1836308843668242;-0.16630048796909691;-0.3129628133876346;1.4894559188770422;0.49419612844344291;0.22485751798221434;-1.2730946184580658;0.11257463509019691;1.0316479323477628;3.688531875696941;5.4406682809529414;0.098721107899187796;-5.261601514076764;0.72624589400495532;3.4513267371260334];
IW1_1 = [-0.23580811256750406 -0.27163719007161219 -0.90325817499639971 -0.054710897771987684 -2.082519744469324 0.73528215722957935 0.42528069960522219 0.94650243248534127 -1.3785951518937958;-0.054733970404800041 -0.031738377742341853 -0.10794795638479152 -0.020884327020505567 0.53446950313136077 0.47590659421464354 0.40618080807106 -0.49573281919162532 0.16890747454715782;0.056467331403761797 0.28166743270820294 0.12551035823108533 0.34296168704476671 1.0973007827416681 -0.10044567879399917 0.44552212869327246 0.42963094275603358 -0.60113409329746803;-0.17271507256155855 0.15886137672816975 0.40156682182935183 0.29859252923921537 1.0339887674726767 1.3808210788109352 1.7676022433511407 -0.41538047462866168 -0.083056017570505417;0.095185579425212094 -0.0098278280052114114 0.058436108529883117 -0.42129693183715988 -2.7155125500996786 3.3413689711013586 0.37729470871765924 -1.506344857899707 0.098320586711721766;0.0042653293138463031 0.21541493934877262 1.9134714013086325 0.35068513521084604 -0.45760210804047058 1.0570645486189483 1.2953181002718881 -1.8120931755983389 0.93153788009226279;0.14243339295596755 -0.11661472456722924 -0.099217206642344083 -0.38184214889483986 -0.57273725660685471 -0.12354577746773944 0.01203077061830217 -0.65782331463305188 0.39751752029319021;-0.057959495316824414 0.12077969848238153 -0.57221361667293746 -0.19996345246249289 -0.57625568953729245 -0.55455949758629997 0.20086018525779373 -0.65779113456678284 0.4305073337873741;0.094225894794233739 -0.10624237175186942 0.68867596771494954 0.10461522586756625 1.0440948190820416 -0.49940430629448407 -0.54666620542124833 1.4885822816449861 -0.66774324511502048;0.44144542103266771 -0.26838657989545039 -0.12001343060609755 1.3295470237652431 0.69805047065643944 -0.53740300570944255 -3.5590001150010289 -1.960559292029562 3.2428407421640881;0.026342558398669839 -0.034825402045310844 2.015150836668528 1.5654017539619451 -0.45989033759130521 -0.98254369048599732 0.55458494852208107 1.9376360896563414 -1.7233114744435731;0.6094968466745504 -0.41300613638619627 1.3026332087476513 -0.18634842504478644 2.5121134842667474 -0.34748353662535636 0.65833094820973015 -1.879217262093033 1.4175424647768549;-0.0058391791137250119 -0.02588917028599335 0.26262523389377129 0.17288023728952792 0.96834924209837525 0.51653502834503862 0.16027186229037113 -0.044360983384120201 -0.22104233468324339;0.14203832692814003 0.43217694713243854 1.2985397498710598 -0.12819183015726421 -0.44198003584836715 -1.4455514562016565 -0.68365024743784064 -1.5859156253770845 1.6180877958892357;-0.37947905241652324 0.36630685675758895 -0.88992411564846241 -0.34797883344539093 0.28897125973123583 -0.77373228543371919 1.2423848948199059 -1.0786540247647916 -0.83605513702432122;-0.0073758738127014806 0.16079943269388486 1.6435315430479185 0.22166932144384152 0.73155359885691928 -0.53065380318600908 0.51883316697114434 -0.55857921448541692 0.27085427977332049;0.022867746454465467 -0.2805569727066054 0.38175356490798695 -0.094418849379514988 -0.63558174967379655 0.42831330062973172 -0.38030877171136374 -0.16226740032540732 0.11266852732735555;0.20894025127424762 0.091480696420893237 -0.46822353080495693 0.37221165096208497 2.1654455502068224 0.38834990982782847 -0.19130557478984256 2.1928653951350463 -1.2281762886875052;-0.45354586726208496 0.033073799383529744 -0.41822573599042473 -0.088550205293739914 -0.33903434256382242 -0.49593387718484971 1.0085038150195966 -0.51049166978132454 -1.065498005065642;-0.58229332048056981 -0.17554418264050689 -0.18005988968240852 1.2115867153672859 -0.15123979664146942 -0.49112599913629745 0.088350445217343365 0.8100606291018817 1.4240291326344283;0.071848395774593884 0.068424152178090308 0.57159901356998133 0.41147828128010994 -6.6902754305487013 4.4795154116231437 -0.10952150087384108 1.9024280696897142 0.31093419707666714;-0.16772455286305024 -0.078060376875873305 0.014248277917674072 0.20335832111839749 -0.7332758644386147 -1.3565957166151175 -1.0708869905340834 0.91358536143179303 0.08367558159336376;1.4893090565276546 2.3606967086302255 1.5860234462909661 -0.55084900186793662 0.74999653322345083 -1.5056675125166303 -2.3462259365941316 -0.31824902587464032 1.4220524812918629;0.13453099417174177 -0.15414012230601562 0.04694064067473111 -0.28507646981590873 -0.27809025699608148 0.42489473735805089 0.30069547840372818 -1.1604845418091791 0.49490026392788139;-0.78386762945854627 0.9140777415458502 -1.552062210016828 1.7573297210939198 0.73018814035578994 -1.6778279355732897 -0.2734069961845807 0.75191215131375422 0.95749351298667928];

% Layer 2
b2 = [6.3973860926449753;-1.7142857008623222;8.2662786496210749;4.9558126675334977;-4.0287191205884261;0.35914269745053418;-1.9482111491513707;2.97304041863184;-1.0817372037466364;-2.0149870496219835;-1.1343435922571574;-5.7684540445099977;-0.25171130087625909;2.3925644421809138;-4.5044679755098773];
LW2_1 = [-0.057139091597295459 1.7271112871864522 5.609242704684525 -3.1234248389788393 -0.028085800107221934 -1.8839793990399405 -3.7782811390212032 3.3826976892871183 3.1995764923797756 -0.43530223489042819 -1.9362796798976456 3.1481070010660681 -10.670904934592048 0.58492353930611174 -1.0298271566249084 1.4755449436384132 2.7944159422424262 -0.69021636120102225 -4.4203199199533083 -1.3776575801771118 -0.56114970333211311 -5.4163346378529083 -1.7882753546083214 1.6413933708150636 -0.60589201328912645;0.26436275435710899 -1.910878855127117 -4.44894329337908 2.423587696004085 -3.6323980974757091 2.5454210765345979 -2.2993141953760916 -0.38027123737901891 -0.64492764944579517 -0.33769576081322811 2.4301961714269322 -0.34660931794547661 0.70157258359225494 -2.5598056437116057 3.6801548995597599 0.051980219241001303 -1.573807115846231 -0.64853115204016198 -2.3462543387337407 2.1638333856563432 -6.1600630531271161 -3.3078478339706638 -2.983173062744032 -0.21883337299530731 -2.3330976700274739;7.6105743024809023 3.701682131828389 4.1731979337901661 -5.2190202395724601 -4.2845959439909951 -4.866024277626221 5.8153164683696694 -0.95133899993724802 7.2029878806136303 -0.4995533077281587 -0.5122770642845863 2.6964119228906309 2.766351142414297 7.775027004025028 6.7879586894029682 3.8640703893257884 1.3109334393914278 1.8418167225052802 -0.99547599807564724 -2.2979908938670865 -1.3967848914231673 -1.0272241556706061 1.2449659106034077 2.4236841814843419 -0.22533284889322372;9.0852466704788295 -2.9339928862215787 6.8371477951784598 2.2004562243042192 1.1899234190900998 -1.7063920215112345 -0.94434926454225609 -5.5785898383772352 3.7849878276314626 -7.5649989857139897 -4.9817768435369922 14.13627888916813 -7.7294153210674894 -0.20981654138627559 5.5471701702807286 0.55114933070385874 -8.5481029283408532 -0.68373508504585567 2.1104029806814806 -19.431489067307787 -11.467337584197669 -3.9845200240825669 3.3172689029636553 -10.736112401594568 21.230118376862322;-0.50676302614870261 1.7401958622806433 0.54315581087821385 0.23393139199168844 -0.010782047948752878 0.60988857733451773 0.097033032256270788 -2.41285747930716 0.45550932979967135 0.79636370849148508 2.0057600893036458 0.027661960834820285 -2.9382644441811214 0.027603537436187833 -0.15047231829002683 -1.2453966117311874 -3.5496957665375177 0.04405876734285067 -0.75952422974940037 0.73876929675545222 -0.22833460287991691 0.62357926962141408 0.11638022134491503 0.61837487651397094 -0.16607184753766052;2.6067955996585024 0.8252109700011262 4.0116487076178489 -0.30296944928792885 2.6336942770739893 -0.61926820940317362 4.05588143164176 0.85950758596226251 -1.7269239716211751 3.6783000574284452 2.2676909338175437 1.3526189410133753 7.2550757349676145 2.2012730098046753 2.1669257794504868 1.8080853784072732 0.16837190268760679 -2.1683467210106073 -2.4852479592852266 5.5722792133466568 -1.7641244423227098 -1.6727590115208051 -1.4052726863732552 1.1374440578283636 1.3213067557206726;-0.70823732666480976 -2.4255844505362845 -0.55863352342199324 1.4121705479895568 -1.2357170195014753 -1.2561088130305809 0.81023657919274572 -0.051390782331124776 -1.4560012276157652 0.079432358842361342 -0.91041535061302303 -0.55476190797981373 4.2909672844398035 0.053108587362091189 3.5152545406366014 0.89487778618507252 -2.4835518216698294 -0.21259527448377349 -0.072329220125640845 4.6963166960018974 0.93344009361343616 0.28387693463827446 2.7423025522846203 4.3390545317870766 -3.4393820167445077;-0.16820391051912342 7.50164044405438 0.012825034748425097 -0.89667402508354765 3.0972714785181643 3.1242847411695305 -2.9809600926462565 1.0748658801891278 4.0051230882760178 -0.033124733570328968 1.2370656949123928 -0.44727690048453683 0.50537159282874178 -0.45956477716250171 3.4565907200679691 -2.6010929314607267 -0.89118546553804612 -1.6767151611343658 -3.0422253328444104 -0.12768279267602231 -0.55601804523892628 1.9623046211550459 9.133387337339764 0.65123861040924502 -0.29370438868289167;-1.7711392539250819 -3.6258160726486568 -5.8437689156426158 5.2179747772382408 -0.76681033521952713 -3.4015577538140969 2.5865678431505321 -0.5769505471823887 -1.713542483274767 4.3362639329176398 0.2912584245788048 1.0968399516642535 0.11875109487857914 -2.3275515561957083 -0.61518823935920286 2.3760241288323596 2.317470679286652 1.4165882110187513 3.4375307414030201 4.0142762234328728 1.6761761972357596 -0.33358836227967276 1.904324066154093 -3.5946220347590958 0.1826373626534496;-4.4990349226974358 4.0516561403972089 1.3077906764468088 3.737579493243457 2.0325302398156633 -0.24562510522542633 -8.6145896035310638 -1.635727731763976 4.7703621672433165 7.5111230452568565 2.2283779488461164 -5.8126733857753452 3.5196556825532284 0.81827455106402391 -1.2525309269985669 -2.3823996620808088 1.9374611182829895 4.4021179805103099 -4.3223567389628768 -1.8012631358201145 0.64408652419886669 5.8808662574220678 -5.4818317743800336 24.793765225299495 -3.006429611549823;2.0476599977083794 -0.187430242457922 2.7805406364499712 0.22017930066563107 -0.59554339777612597 -0.9905485831279055 0.9180807287843491 -1.3732064958138142 2.0462092380901926 2.4441063223232424 1.6482764363617914 -0.39453169338514438 2.3814241299366246 3.2780322006765856 -1.1218430519871283 -1.5328703035466384 3.2791226746384194 -1.9258199695953628 1.4083943580999396 -3.5241864886772083 -0.40989811052127711 0.43545699116901182 -6.1589923607625039 1.5079493941608908 3.163996603622802;-0.97318141533510061 -0.0049358014368507857 0.16258526199108311 0.33201651940550447 0.64293210462828665 1.1244124738474326 2.8019078035483913 0.49624893883973409 0.29955844119949204 -0.80545371932107135 2.2808874784380126 0.18112877016675955 -0.32241253561083749 -1.2945610954838926 2.1379223604151862 -0.68958243257225893 4.3956242974874353 0.59510383314010173 -3.4901836684944225 6.1237592472797981 1.2158973092982699 1.9206864378864692 1.1818812202496172 -2.8190543830248029 -0.57378070652222768;-0.86335339662764454 -1.190204487018796 -1.6274880742319178 -1.1378940151389716 1.7407867746673298 0.48370736261958619 -8.4562546897519546 -0.86443038747064338 0.33384462518332875 -0.035206327481327876 0.95824974177082645 0.4397773865062875 -6.5815561354900405 -0.76215398325879602 3.0699935422887519 -0.083155241674080907 3.4749981249868513 0.75184860284444732 -5.0745726338879544 0.30258300848824371 -0.0037887602306562549 -2.1454246654515305 0.17848666109826092 3.4356348485691326 -0.67596397520107898;0.81502463897002309 4.6267942301221492 0.74852477323176769 -1.3403447977758638 -2.032538503325374 -2.846011602917033 -2.9692186199130362 -1.2560526686914444 1.2002009659395569 -2.6907108398839199 -5.1258996459410202 -0.19176562721256851 -6.1399590612104973 0.64811582897596531 -1.2367850907380162 3.9933881352076424 -2.0768674614326077 -0.50799992135946359 2.7632151067479946 -0.77118223048983303 1.1773940087550787 -2.2845045114097267 0.65826950517173077 1.9027621427251356 -1.5303756066479368;-14.226755025499861 -12.973481027504608 -13.383777043800468 17.980631323226394 -18.411493452449747 0.21501394986321384 12.40558136066568 -9.7967937757691939 18.538449527028551 7.5740423077302461 -4.5819186093633295 -13.030248669382804 7.628800314482195 -2.0193860934946288 4.676283180685239 -6.8268878501617785 2.3127167872809187 -3.9323380329028983 -3.9570959974874249 -21.919883103927194 4.8052875141172215 20.107468966492828 -0.13510554386154688 15.44282974633381 2.9699285537820064];

% Layer 3
b3 = 0.010001322678575229;
LW3_2 = [-0.014491920892337678 -0.013849490522395943 -0.017851788861504961 0.0034741659986962086 0.54956272938189632 0.016089908106878919 -0.078100914526883947 0.052325656741536899 0.016343314596707285 -0.010932713382549177 -0.056295521835997495 0.16177478152145591 -0.3950612491152336 0.17873716434980735 0.0072963232970834771];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.00314743155299998;
y1_step1.xoffset = 1774.26765577007;

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
    a2 = tansig_apply(repmat(b2,1,Q) + LW2_1*a1);
    
    % Layer 3
    a3 = repmat(b3,1,Q) + LW3_2*a2;
    
    % Output 1
    Y{1,ts} = mapminmax_reverse(a3,y1_step1);
end

% Final Delay States
Xf = cell(1,0);
Af = cell(3,0);

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
