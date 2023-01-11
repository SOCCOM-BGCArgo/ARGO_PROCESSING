function [Y,Xf,Af] = ESPER_pH_10_Atl_3(X,~,~)
%ESPER_PH_10_ATL_3 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:23.
% 
% [Y] = ESPER_pH_10_Atl_3(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999025224415;-0.999999994516886;-40;0;11.29;-0.2178;-0.12];
x1_step1.gain = [1.00000048738803;1.0003169936794;0.0153846153846154;0.000310077519379845;0.0705602618858232;0.0479365703301392;0.0145655815308426];
x1_step1.ymin = -1;

% Layer 1
b1 = [2.20636427040742511;2.1712014180433873989;1.8017292149646664257;1.9557713853122045133;-1.5980321290498573639;0.74378204220184673812;0.47474161801737913002;-0.27061379570602506206;-1.2128088696346202013;1.4289984284528636227;-0.83709998523398110315;0.10375181539850651669;-0.015182240735084855457;0.52670436708449597241;-0.48965504879114724046;-1.3120500854566896631;0.76965866253268666597;-1.4987627385255157808;-1.0685409214699903568;0.43659793658850759579;-1.1589639688593214029;1.5881775292540636357;1.3181573457518998005;1.6992187509042946392;1.8920525428440104765];
IW1_1 = [-0.85597847234786461623 -0.9583607872327634114 0.85293803365497167235 0.0014593679465920330988 0.38477124001942575315 -0.077338924354531796146 -1.491370311392345549;-0.4109031633924613014 -0.28369672191371081427 1.3701826563024885886 1.4991145799096938784 1.5994686569013922828 -1.496398099045456842 -0.49635160401176120493;-0.41523635371523220927 1.4851219969903299845 1.9791968188030995002 -0.52696825468871710818 -1.6473863270305084683 0.37950715869431272953 0.41126170049915206084;-0.66618144524480427116 -1.1079197528469379996 -0.98036072831421006146 -0.37432208906037744978 0.86852992008876550578 -0.8209255855619094211 -0.47086026356430527029;0.20762923397650187574 0.76192049138972262057 0.82892544817011815095 0.72281069867016067843 -2.0221384460696132379 0.2777572556095467271 -1.0849156615241948298;-0.25976447573827388116 -0.42869580372303789373 -0.088373710448587550204 -1.3112000826218359606 -2.5622437067627403806 -1.7316524901116265855 1.4681098006383916932;-0.62969438354699003302 1.2655109122140248967 -0.94496567663593844522 0.82714811184142167377 -0.59857132494793141575 -1.5850065597622366997 0.44587481354782643184;0.32788014050317249026 1.2639487476171382685 -3.1471510645345261992 -0.22023333705996567256 0.48321666312412125954 0.0032065551170566098868 -0.77683710924025661804;-0.24794001418950226134 1.0675443967810467338 -0.96178689981052278846 -0.5800058425897137715 2.3635318190112584169 0.87705172284212451483 0.37256840960301867982;-1.3419211177228778098 -0.079496329025653456046 2.1849425722986337561 0.59989917765646849634 -2.965580254891621248 0.51487683628731617524 -0.40504861024514871026;0.94822069888144666638 0.073393041963282376616 -1.5944361712229713213 -0.17231603684329266368 0.84242676516573222756 0.91043265708690257476 -1.7019398972021766436;-0.48626607482895894963 0.0099587810385651058531 1.6562412070782503992 1.0984293190272951257 0.6941952036104969137 0.64986112160326514608 0.1307741168119137265;-1.0058928661472865684 0.79423442932229926239 -2.1726311762871763911 -1.1287095775334516912 1.3327938086737953771 -0.62921860662859718794 1.7083476070073473352;0.33939173871030225982 0.18831959233848205582 1.4243767913459124408 -0.64182818157086674038 -1.4745269412103461093 0.41317770741718712335 1.3554668938751190943;0.45412216417682593761 1.2311785964099530855 1.6514756629117450792 -0.62770862632174484474 1.1968354223524582114 -0.40887811495169829046 0.71044794675557598751;-0.98990419741641300089 -0.59388447110058084633 -1.3405114637906885289 -0.63431766156824964664 -1.7703608078675927562 -1.7066264605256604714 0.35078397906264169315;1.9271237908716256637 -0.93478775451173035105 -1.9560616404765127641 1.9736818558859952066 2.1319579761851072952 0.3378002191866249504 -0.251318207691838158;-0.39768800494010719815 0.34436749526994669957 0.34452444076070776191 -0.82541786033691944002 0.20002929733119551026 1.2089492679091708993 -1.2690457310664424728;-1.9885968530250330844 1.623728267188286889 -1.3642978557051597122 0.7936884639343984027 -1.3058817209370581214 -1.2737647884312510715 -0.3889392440695939035;1.4083525659254179097 0.2113476836327926045 -1.3866847525194787583 -0.15439160045191818726 -0.0012556041671167479656 -0.67734418825335962477 -0.70440650278781735683;-0.58863486422783939389 -0.81452892549932387478 -0.3117384929764375201 1.0040728689486875247 -0.38651805506536257573 -0.97324479390126183986 0.31171531240013833353;0.3977345357069808629 -0.40783984859242328191 1.725316331039588702 0.88282723642648397444 -2.0847215303706665424 -0.012990202305954377457 -1.2503626136610230812;0.97110445875972872809 -0.0080343462100835559558 0.84748373327721915338 -0.95480376641002440152 -1.2549650575204798209 -0.60315420186581980566 0.63336512534750566417;1.4072796037986718964 -0.78823277285902948464 2.3180564013560589487 0.87325982503738275309 1.9775575213918099582 -0.17393065853479330873 0.23285434312369482357;1.8681584976765517325 -1.9565751949775862517 -0.30711584471944652908 0.50669462810720378165 0.58385920664208823538 0.41277308535238010911 -1.290056498302525112];

% Layer 2
b2 = [1.5474455289861932528;1.4396198165376723477;-1.2534380873871417972;-0.70657734070322142372;0.48192824679531676146;-0.92902759641828702719;0.48658943890352329698;0.19952325194626283067;-0.37589556012063291002;0.31396081305927420191;1.0068728480266935232;1.0994778700513330172;1.0489574203302094446;-1.6831881932060215323;1.4332279412484234538];
LW2_1 = [0.0056155658646638659942 0.13285177140242987459 0.66878180757904504006 -0.1260462444850388164 0.3455752782297726311 0.12433519985748099634 -0.6120373974572537179 0.37031878183821725914 0.78646440088431346371 -0.42140779102803771483 -0.86928418078074864983 1.2830760731782411277 -0.012656132444049719363 -0.68970298701280285947 0.40443982060487482233 0.63850188684160014141 -1.0046399596912525265 0.38868646404377382719 0.31103127940159142062 0.0652982465586814248 0.35001019314859749265 0.46080534395472211306 -0.32838382873944477192 -0.41805702634769825732 0.19410308828098421507;0.29618960910276848075 0.915713998054422329 -1.5227246860137906292 -0.19663710440726306361 0.25244461360887265888 0.38741457987406929808 0.24357368507274276825 -0.43152572942853811844 0.043156193373476571529 0.79683188595889076211 0.38911070199363601008 -0.34549669297014840019 -0.50520887119739221749 -0.089272026037807275523 1.1523773053252419629 0.62316662152794499718 0.0064605754119087133863 0.059710866943858624534 0.16470269634013062765 -0.089394016614659690956 0.04159086865010431372 -1.0706788906517483717 -0.12229808044099943487 0.057448982815554557047 0.4473781646037886528;0.094624226802370495215 0.16466657699524475666 0.33912996763307190484 -0.11513984208328643943 0.30973633673238443187 -0.16768031680513909021 -0.86263765038635975557 0.21250642076784384549 0.40191656145170362446 -0.74041911160158835337 -1.133068756358008411 0.53225775309765055443 -0.26693422176846937077 -0.66489080265638922373 -0.64972499246078951884 0.62661654767965380142 -0.087923240097334653953 0.57673323158713152292 0.90148005879794546047 0.13131117260162789018 0.038751839641236691947 0.67554005136714623259 0.44544738892243457284 0.20389553989117301103 -0.083460269188366603688;0.52916824811528795358 0.054202795152401998802 0.67270884147368958939 0.69499199517602683329 0.071020340159481595843 -0.16247687939113025779 0.69506084335588058654 0.50868126745924069354 -0.59108976553977066182 -0.93094191288625194325 0.32531396967904396833 0.058401154882614006447 -0.16189222352637805402 0.72460196400964793639 -0.56560777397322159477 -0.11257270109021290383 -0.10737222504345055041 -0.12559346315951844253 -0.96070270363282328496 -0.57628821713493394441 -0.27910080864574970949 -1.1830311015397656504 -0.14593680314683613042 -0.089415825808062232571 -0.38428638357884975507;-0.1734857810071256945 -0.092799733970548581308 0.35583738878206905287 -0.32906375173425095326 0.069452046893590649312 0.44351959251688760633 -1.0572049448016986872 -0.048673796393931288706 -0.1323951818138252301 0.97538076796062000984 0.017231999508420119821 -0.62093647320017875124 -0.15116388277328413059 -0.083552327628752767752 -0.12350689604124569232 0.26964655225289402285 0.16276230954332837619 0.015057744910006392799 -0.66400959059142217189 0.6344748393465451386 0.30940350580406539294 0.027002945728432602851 -0.11861465440894715384 -0.25534159093088903969 -0.71235367438659813466;0.079580238652571341773 -0.23932213916478084492 -0.10390078900816472618 -0.90242406729269930565 -0.14280961931297009304 -0.042801701173622176688 0.079430874564369993718 -1.1640638569348993325 0.75490393677809630724 0.44583754719944196365 -0.62417377360778947271 -0.25757132621051526122 0.17534647419472801366 1.2641989987383308591 0.0372692755537804149 -1.3345051288654636856 -0.67847878462669475308 -0.34780939403685545708 0.72522728215005782815 0.71192585584790601949 0.43465128689168219811 -0.34564748468922756874 0.96229597876383199395 0.14398419530376282971 0.8190722027065800992;-0.22986892371846306959 0.94260971657582504069 -0.30259912754065976248 -0.023172160261583012869 -0.12119848279409138614 0.6253681086419848878 -1.0545295781645334099 -0.23068825272947973692 0.09891326678018322327 0.46074033628915300742 -0.31468565862189912252 0.24967361222800124598 -1.0404880669717677666 -0.24729731290036238711 -0.14840060970633256465 -0.47149198117321899915 -1.0070992895748072105 0.91185906709448205643 -0.63890656010299518019 0.82638867537444660183 0.44764810279594585696 -0.69811240083804604328 -0.70930759301690060603 0.72205321092260144056 -0.33934392283541481294;0.18827377108977505493 -0.19045712706334849695 0.068313448049685462915 0.15145699888097879926 -1.6141823940441233365 0.74059986234588359721 -0.2770589685196426788 -3.6478221369704846508e-05 -0.69650208055535700602 0.15079312745248479866 0.5542719022012806418 -0.11532764667826669869 -0.048711538245278643289 -0.15919513597523571313 -0.14629132768136257625 0.31154767733339955305 0.43519682604381942914 0.81537932252865075178 0.24957644347953555908 -0.73383096581460294239 0.45508232417736271813 0.11431697683104322105 -0.37351866906691544301 0.56022803413141641915 0.54953665164821374756;-0.27493610844242638569 -0.20046534289983372612 -0.4026090644381906114 -0.48679957904783460743 0.30510114490150974653 0.50135775638724788816 -0.77210969850076793808 -0.90538054727847983028 -0.012061791174002975749 0.18629416876796420222 -0.30114083126224289177 0.17575456583866333471 0.1565178962583370359 -0.30022847764181881391 0.15562755630329444534 0.22488166686166091712 0.84486339010993549703 0.04531535271432313855 -0.71489810996685987377 -0.18264935088896053217 -0.41470579509686900099 -0.39583560976820597599 -0.51267727971537624398 -0.38538873183527361244 0.33839139645091947939;-0.22924234862715708871 -0.26619718469965797381 1.2008487871869280017 0.20711307445354534829 0.20553137008281766041 0.61851888451701941829 -0.50757288948095702352 0.43768112281030063793 -0.53458459256206702914 0.26639063230798620818 -0.39224218154060086494 0.19625551024604209749 1.3001446581755813359 1.6956836735507923031 -1.3429900862812786944 -0.001765007454566838458 0.23411626677519084549 -0.59585089866974683126 0.044527570511899282335 -0.10856639757074373898 0.081797911393153038828 1.104546747772035209 0.72074743197560098285 -0.81321404682553255494 -0.17008098224946094357;0.69827280290595561407 0.76561295723672606606 0.11430756412124132293 0.78677705222785132566 -0.23354912212570314023 -0.76403642348412004193 0.63600519507405983699 0.0021246095371196317171 0.71860096657079275406 0.74508676749561741914 -0.23102172085022879 0.14064555433509362548 -1.1938745688973726455 -0.5858525842013186935 1.1971092703622436826 1.080798665562448857 -0.058900314852318211289 -0.57294284630338498943 0.20004697792776215959 -0.24379996298004708022 -0.65661423514860473016 -0.96475392448288277425 0.75446852429927635519 0.20433567917628295274 0.26961996434736024142;0.55351418450371525637 -0.2932165678518794949 -0.51869287582773337952 0.69457298951539403475 0.088021902561484455241 -0.66726110563114815477 -0.54561767672827787834 0.80425224811756423904 0.48388047706908965395 0.12158220975202266301 -1.1454487352393984967 0.98226556213474824908 -0.1268519248054884585 0.13611946892126516873 -0.91060744234619916959 -0.88641531626915215369 0.71781939253980819871 -1.4094951784890723534 0.15204725827633211455 -0.18078883182249724793 0.046853911097909470285 0.013592425960924076414 0.23618459034413671893 -0.17486378745451647276 0.86032973544054547332;0.69546176141036086804 -0.15289745596946122119 0.39435130863355605246 -0.56758861324422227312 0.22848322136134527338 -0.046660492879509171682 0.57009436775547972509 0.37495744027717498037 -0.028319711466545282458 -0.7824302551853776766 -0.59049435314011289311 0.40460622759169695417 -0.55270926252834828851 -0.66777067215079222162 -0.41453440254044138236 0.32991032784032831016 -0.30677726422571394149 -0.08764508416371721844 -0.72410338273943541942 -0.48700859920890599009 -0.40503909565127038661 -1.0572725416730812675 0.068404220737138463071 -0.47002731012655551979 0.11838608136621441014;-0.11496696785289131637 -0.63000847060482723094 0.86505271496621849447 -0.31910647117778978554 0.24485890381640515456 0.02359252047250188361 -0.99653072290185717996 0.32308840981024694594 0.5458928184912359427 -0.45108304353941469644 0.59053045336618792405 -0.34888401993646561694 -0.49953023017098729941 0.20132021193921786906 -0.59256000458654467522 -0.40032467332035809315 -0.025636823030734133544 0.10150622351275929989 0.50491101476188460229 -0.19273210040478827576 0.39806008141487225771 -0.83971240517140255388 0.19302942866746669592 0.6289879479132577389 -0.63854066234468198626;0.27423367591330571891 -0.35523058763208187338 0.46644289019236417904 -0.051482044807781686402 -0.01684711103598952181 -0.016352399524499780514 0.51977361342927730359 1.0109238767662103164 -0.5269486159350351695 0.38441378294367717361 -0.44848167844173281171 0.06973473277475174259 -0.46578913937507182563 -0.56531007913491171379 -0.25505621800298083901 0.20435874680382454249 -1.040899877898820991 0.35189361879237152086 -0.13028294939583812129 0.52291683878047345857 0.049587308384107083026 0.14054795270144712638 0.28763449978291155684 0.51256469697669326813 -0.046331243373864017254];

% Layer 3
b3 = -0.26183340845494212923;
LW3_2 = [-0.78866060402832471876 -1.6444407138809349345 -1.7802616456819115864 -1.7498681106159901244 1.192719784513214254 0.92420991829592591404 -0.81800081599312324254 1.1050053847238472393 -1.055280712674263377 -0.72112346052556386411 -1.3161468075885607298 -1.2635091307847874376 1.7561236570963838499 -1.9440685466364693301 -1.2879852716746125196];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 2.03830178109291;
y1_step1.xoffset = 7.47181919175367;

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