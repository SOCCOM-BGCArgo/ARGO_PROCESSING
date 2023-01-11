function [Y,Xf,Af] = ESPER_silicate_11_Atl_1(X,~,~)
%ESPER_SILICATE_11_ATL_1 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:41.
% 
% [Y] = ESPER_silicate_11_Atl_1(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999994516886;-40;0;9.1619;-177.232103709389;-0.28];
x1_step1.gain = [1.0000000017546;1.00003626070655;0.0153846153846154;0.000309885342423303;0.0656325902039552;0.00425400504498441;0.0470499670650231];
x1_step1.ymin = -1;

% Layer 1
b1 = [0.10190086942604288656;-2.4677416516608188246;-1.2304021976252135229;-0.45813362521563411045;-5.1041212093887065038;-1.129895406530097679;-2.3384798709561258612;2.2082118398683174831;-0.68255400877531013926;2.2773354545739015187;2.4593818605363626517;-0.26884025774413239374;-0.82683634244293402471;1.3678891334127349921;-2.3115974112208617086;-0.55466779795191989422;0.059898716498264781283;0.23529069703608651576;-1.0682629104144563126;1.6655815838272236551;-0.22911419460051091201;-1.4232561364086318001;1.3088806879507601799;-0.37016429445723425529;2.667399610106186536;-0.61228923480738939311;-1.1419858066423231868;-0.19458096399515012465;1.4780818406268143583;-0.26185569405523401176;1.6029164184087301148;4.2585165874323633872;-0.73210167319117358531;0.95229382327277212195;-2.1313575482968185604;3.1773630281339220183;-3.7319335983408317148;2.6507991011909406431;-3.5063562431614014869;1.7456170985004675344];
IW1_1 = [-0.44862608957871036308 -0.75658863886528782139 0.70565783431871298426 -1.8584859080951732224 3.1600016672176449717 -1.8640201429975478131 -0.496881662884584685;-0.62874035661155303067 -0.26419759974767931254 -0.62341263236144284221 0.84079465947879394871 2.3952562809159831936 0.59869970641900516206 0.5424466345015727109;-0.31867674649291627409 0.061321382988745683218 -0.87205991239957592054 -0.45312943239829023279 -0.46626861134209030002 1.0778920621762335674 1.8702562829371029984;2.3980377064452098601 0.8212399160549005428 -0.52085266702326604715 0.67659346544104048693 0.6399717314948437652 -1.3166987458910195397 2.0123004334212621913;3.1647479940223375294 -1.708742624197793436 0.10397193336032790889 1.6029146990750533508 2.1522396127540357469 -1.5773574611438849757 0.062055435111252575664;0.94988977265657181803 1.8724538891222222148 -2.4403901887290890826 1.0851287454907254659 -1.1779948285348882386 -0.1610792097170182513 1.7027659972311093295;3.6071955882057702425 -1.3394014015825197461 3.8947851708746323851 -1.5957898413667475079 0.86723017384162659482 2.8183935265472670473 -1.3522909036293853102;-1.0694864836406086095 -2.4855954697595104008 -0.51031135111097292167 1.0008384092329440929 -2.8192847617381175951 -0.32049164770797627932 -1.269169291000153077;0.062356943682319455891 0.18483138513435282135 1.1788453580087756567 -0.70197218316363596546 -1.6994905572511957992 1.7267304770480051967 -0.93579839227357319853;-0.21329334463779564546 -0.079115660631522324864 -0.64151446588959248896 0.89021049936873242192 -3.5184322155266087861 3.0192123710854703589 -0.14628432040760180888;-3.7540945726014687267 1.3312648760481391452 -4.0107983878507909026 1.6709269349734039434 -1.0139387796974270728 -2.9277217898839471566 1.2971972610224786759;0.10535876377677254323 -0.1141554192925603034 -0.25695448565004341823 1.2071458182138625759 0.45936253617760181012 0.72735096880484739756 -1.160718134707788618;0.6466569238935750974 0.29329065970858630985 0.40865548133754991955 -1.1881879457500708064 -0.11436698364202699252 -0.80404890587212374253 -0.91394346781301671889;0.16181963552555359009 -0.19666171644799135265 -1.0713133960132512623 1.3499353428647080655 -0.329572077770584837 1.0180616486521008035 1.5242905441372165409;0.50288807252839973216 -0.22564012097322749817 -1.2445932413909983616 0.9914692160908185059 2.1482443546813634327 0.74101059272473790429 1.875345891436719592;1.5668710170774553525 -1.4993067671491553394 0.35282612618366454571 0.77517368119052743047 -3.1908283608387475461 1.3734027666072157192 1.0279218114810844131;-1.6755773134092049048 3.6704480852597360929 2.3193127060855642796 1.0080768617942421628 4.1517105791029162987 1.4176568223041576822 -1.4218988466948343508;-0.13284158985091903959 -0.24412487564242224547 -0.18194353874389124215 1.533158651278406559 -0.73105585571427633695 3.7572769037641622347 0.094720370998550240671;0.61527526494138973234 0.31261092306937654506 0.35204219675666431622 -0.88652029030747647731 0.69204450555544450019 -0.92487771284160069651 -0.81282522976711979457;0.41278689372010207714 -0.28514134440116950531 -1.21967291793956778 -0.11485950736973125907 -0.76805717621924218808 -4.9338480934417106027 0.19483177590247230393;-0.75298328966931804551 -0.29655211153196492369 -0.74523377958699754675 1.1832446911185694471 -1.1851684909213049401 0.58024631576543483913 0.51246547542635678685;0.52666695239666916351 0.5146163537423602552 1.4744359920233727212 1.8705375448224776935 3.5109842242418478619 1.9240388494416325038 -0.96605186240048346313;-0.17472622711024801689 0.55379432167968090095 -0.22366522623838114292 -1.6288154532942868968 -0.95699503653047779395 0.19711029676983612968 1.7133991748691297374;-1.4950118316225573523 1.4541641475165192343 0.96404883487258474783 -3.16848329827496622 5.1468895807765999706 1.4871988149994672668 -1.9948467966314640254;0.57302208188666359234 -0.030099433295627080198 0.65530435197763592381 -1.250132562929065605 -4.0537522957268263823 -0.78369814695013761963 -0.46409356627256870276;-1.8752264137492182794 1.4055761495895795488 -2.6668786990631141443 0.758484918771975547 -2.1598732665520627982 -1.594995638620962719 -1.0069123548941678603;-1.9610502055450214076 -1.5113117156204918246 0.52807114208687744306 1.1043689396234368516 -1.0638308641309701663 2.2724963485736280333 -1.2817927537804409965;-1.1622508280326362406 -0.55631889616059837511 -0.081113479051071801784 0.52297962806748266029 -0.89546559880512832663 1.6306152569452574408 -0.32071299305967243543;-0.063093943110368463922 1.0753602237324748625 -0.57861445537228717306 -0.97660110378276232534 -1.5304164547191245749 2.0679121222286171999 -0.060609612743719905315;-0.8904481602118884398 -3.4058019506176862023 0.79271132580776326737 0.089155289065579568231 1.3433220552439919526 -1.4583792907562103114 -0.093068930029570534623;-0.24898127372845824112 -0.071845163055925545526 -0.41488609688576116863 0.15503325861228625659 -2.3047868236781909879 1.543720030664891274 0.75732915511417508991;0.30261770099778262377 2.9368021827693548964 -0.66394809925525977157 1.8696210867972078429 -1.5773183244267579006 0.84298957306611987939 0.57399322272945985102;-0.17091066060335913246 0.29216015899920488863 1.2231416865997630783 -0.81217582047911773824 -1.9063021780055344756 1.7251667164622332962 -1.1326139062317008843;1.6270956401673115632 0.26839750040821463983 2.3650750650217498006 1.3881248770371130696 0.22205267766120151407 -0.54510280290663837288 0.19718221778937047861;-1.7308271150587215548 -0.52255334242245699006 -2.965868374214515768 -2.0045640601840291062 0.28640400245427921622 0.68428642480749335419 -0.12884492346008671881;0.63097145916101005092 -0.26991538162000905388 0.39765567484101121964 -1.2351360810843698967 -4.4904767233222093736 -0.74191486763128899895 -0.17323304429982885888;-0.35314267646518748833 1.5465041138354949091 1.0173690114926365879 0.16523647456691156754 1.0238149594171108703 -0.096425241114472037829 0.082545822185203607257;-0.017592510094373903229 2.1286027052667217951 0.012430792602618027931 1.3135745037834789528 -0.93201298640585394928 -0.31491201584144945524 0.52866108434624148327;0.44381788697707047486 0.33701064000326308445 1.8442638854819290994 -0.92629483103263166033 1.3570138920356655809 2.0926998346771505943 -0.83624691300778575265;-1.8597790524117878697 1.578764209964695242 0.10349806096380984033 0.26982235092028394874 -1.6941515059118072983 0.58499417386874963398 1.147664323370244821];

% Layer 2
b2 = -0.18855910846884180176;
LW2_1 = [-0.6665849533336012378 -5.830468043981691828 0.45711125137630959081 0.10156210059373041443 -0.3786626913289735108 -0.14058984204277308327 -1.1916949442792368341 0.088612422390053610655 -1.2085213632351923252 -0.82312646225312680048 -1.1181141901967339702 0.70906571288085573634 -2.2867533960662269621 -0.43713934846446839133 -0.52556766125621490193 0.34988986542047451067 0.10447471466175031751 0.37151398742147473397 3.4259961317105163126 0.13638423439163527284 4.1432111962898527935 -0.15152817465795809238 0.42660921698084908193 -0.25594050085230135583 -2.6417965128540528141 -0.3193442554496106478 -0.28651294542952471245 0.97034817404146367359 -0.28670765524542274338 0.067320999895140579139 2.1226734608871331567 -0.11702952390797131954 0.9337358523548756617 0.38544449851379708472 0.2386433897101466195 1.9516756598325766259 0.61287985051078586185 0.18921476040723947398 0.32032659936348573115 -0.05503882168605552655];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.0145232735458572;
y1_step1.xoffset = -0.52;

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
