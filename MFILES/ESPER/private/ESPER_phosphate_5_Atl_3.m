function [Y,Xf,Af] = ESPER_phosphate_5_Atl_3(X,~,~)
%ESPER_PHOSPHATE_5_ATL_3 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:28.
% 
% [Y] = ESPER_phosphate_5_Atl_3(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999994516886;-40;0;9.1619;-1.99149636107278;-0.28;-177.232103709389];
x1_step1.gain = [1.0000000017546;1.00003626070655;0.0153846153846154;0.000309885342423303;0.0656325902039552;0.0604107502911417;0.0470499670650231;0.00425400504498441];
x1_step1.ymin = -1;

% Layer 1
b1 = [-3.430735538199228607;0.94657882172530960485;0.33823174661237115357;-2.8957386993549198984;-1.6207575575524382216;2.6816807322410101122;-0.90166785572309215979;0.89050338006690221704;1.5299093752028323223;-0.1959411254978609962;-1.0658895232112304274;2.211384811844814724;0.63113887878194052661;-1.1017134827963008803;-1.0944411432557259456;0.036670955675984530375;-1.0359575986590050523;1.2629288405392697925;0.11605726432714501273;-1.6375430785509610399;-0.91654935396289305416;1.4965424668477844694;1.8807277949207608092;1.4420560238672870401;-3.3627142810682566498];
IW1_1 = [0.4090536802428537122 -0.6706765486455512626 2.5367419162250199882 -0.04941221860453252207 -1.0788545225130852234 -2.3047537028742626219 1.1309652160709242619 -0.73176182484552565466;-0.47920090916332297448 -1.6799701498473622596 -0.98958716058340623345 -0.051140361437722721416 0.35110278665895994221 0.10755607389080928948 -1.4622720268500009055 -0.51199179725562604659;-0.74323734467087976263 0.85703244806217004914 0.35428012158589955449 -0.19436949925783642001 0.86606342097713695516 -0.70253447539800673027 -0.78067702799992955143 -1.1146137092825942894;0.29181748106998911529 0.12893292805947506596 1.2663828127842295945 -0.066612924733930989341 3.5599574488487695234 -0.44521892687954661705 0.82704935826114567288 -0.40794874825555926812;1.3069280817956401286 1.1590528576412681971 -0.91410919090629394468 0.49062286376908942342 -1.435081576645427015 -1.2874013918148230484 1.0055570233299486738 0.33488163984579061294;0.15184364606770148209 -0.96031989627608216331 -1.3267522251580068371 0.085448649113860949966 -1.322434442992059811 1.3612019235350609758 -0.65917227758194008036 -0.82008562524482342404;0.33356111209817546071 0.17406718166505000256 -2.0224623302430222438 -0.76776509343654231188 -0.12146202999051775029 -0.48892200315832679003 0.25618351190886790025 0.34490624450225609854;0.3699181845875650354 -0.19616715373453968985 0.78914903579656148214 -0.8349461053293509849 -0.78224824315382057716 0.93171949731926906324 -1.6622356408637197767 -0.31689333430569693872;-0.06599935984556537516 -0.36822975843898325676 0.12817890630423239839 -0.05868132911668789975 -3.722538400308716966 -0.56935546427073291387 1.0232371629903873611 -0.23245444094537964785;0.87491006314714359071 -0.88076070843843690117 1.1234967323692173213 0.61541285857593119868 -1.0592042929178415722 1.3479483628464385969 0.43970169380129575654 -0.51782035137694426563;0.55960136405900651191 -0.074123158589239290794 -0.23368102193778850229 0.48442001401404466421 1.9177610663682134184 -0.20483599498471821709 0.82875984103678668191 0.061363304058722020762;0.10623031899041923687 0.35217207047122750518 -1.7116321423870499707 -0.36437068614420625723 -2.6033104360566099444 1.9544818381428492593 -1.7883735366355737462 0.29980228966948258007;0.55787493504226726149 -0.75398178060255682809 1.0399605539537888355 -0.81230856177623722569 -1.4869630048548594381 0.40238358643558308492 -0.98007402930948561703 -0.9611467265072004551;2.2525663597505101166 0.52755414962169477899 1.3460981515811991649 -0.36120739690500236474 2.1862789647508216184 -0.21271608678940878701 -1.3685147228450029733 0.73291549271221512196;-0.26315065061835757332 -0.49647432334841418156 -2.0096348036495266243 0.23326898284756827917 0.46269359182907388028 -0.82965403030436701837 -0.67650242568814045185 1.1987917364873730275;0.68281869794590888123 -0.032866902929626920504 1.7019762795418358348 1.4450411964334295423 1.0742981283069619192 -2.2080483193626312399 2.1638529399253898511 -0.40630620203930040146;-0.033449937303047301562 0.68725478775907655749 -1.6113318149649733524 -0.22283556725502134954 0.39072501269793347145 -0.8202943351269607053 -1.4145086734848288934 0.2932061851726158741;0.94984352432748619943 0.67011415485948522264 -0.38304995890586240881 1.0766098003932760285 0.9543773216057350739 -0.16728319283710027521 0.60500880012464741675 -0.50440048627370159817;-0.65065921119781988669 -0.79623574767210170133 0.55816833300587886946 -0.12150592799328696625 1.3004421209144860772 2.0419653597694091829 -0.41618582580650492764 0.9254915791920605983;-0.76740416448697490459 0.078062896164512063368 0.47977775504172676113 0.22178607573928721397 0.10363193560013274352 0.79295717712972368929 -2.0471815872262228986 -0.36986488327169125823;0.4044962367197293962 0.08192626189333253317 -0.53363101189081685316 1.1410223133378991633 -0.78662909660014945779 1.3589885371316277052 -0.55364274289813786289 -0.45782896299826730813;1.7392631264734836893 1.0313395559216755881 1.1998923138393065013 -0.19279626478287059665 1.3766576692990162289 0.16840187802271502915 -0.96189699552565133267 -0.88983542225728051722;-0.12074272837602301189 0.23523742110784978676 0.37150281897274006759 -0.18168231879497298564 -2.5677728287249403927 -0.97652522812417330567 -0.22932330125000288668 -0.16377785781593512393;1.2563762743927253673 -0.44898294964479712466 -0.6800146247349394546 -0.9435910447995953243 1.5831113499138060874 1.4061392089807869255 -0.70917179248676431857 -0.6404279977931851553;-1.8481315202354293703 -0.26871072502203879884 -0.044369559972042380003 -1.0912827054344012012 0.19111825608792204489 -0.27528022708564436893 -0.095835920225238735437 0.36496277004835764313];

% Layer 2
b2 = [-1.0840234219334881072;0.93964943115348853464;0.94544695415127466553;-1.0422602689024396128;-0.64834615242753512732;0.052704828775006228037;-0.78118716011569488966;0.21468511020639552411;0.74501801548721080515;0.35501798974844617218;0.67045532564182475088;-0.95027510482086596788;-1.2600678415792310183;1.3137712267015164702;-1.3975576134319926958];
LW2_1 = [-0.43212954869434178073 0.066746813520764836092 0.063915205138010564023 0.54681878636380065739 0.20248339034329698727 -0.35435483391764127648 -0.019843284026366786454 0.57117083796933043605 0.61604874630193862117 0.1395181762420334981 -0.39751801364591310417 0.62473620935833673595 -0.56294467571164630293 0.5492461220274748035 0.70204051222534014087 -0.019538177437310174506 -0.69993277181504676143 0.090510572041472892035 0.33030184058219252963 -0.12605492670522966403 -0.47408021183694443179 -0.55676967297009194446 0.94503758916567204285 0.60484317399318321407 -0.31798266171026762228;0.52621322250707236456 -0.88235519516382876493 0.14673834943647950935 -0.50058835808163515413 0.77125708743829524128 1.1225549297094390422 -0.37366832302661273646 -0.35035193117874979807 0.72753734456731933289 0.35282818248749642231 0.2638674025159301606 0.14339748378796068895 -0.49168314743864105498 0.18648214523511058571 -0.52055848358450151903 -0.92988377258516929924 0.67698993711350385905 1.1278844131388445771 -0.41727681190562287972 -0.66572402234081573535 0.3078474291713484412 -0.075346045970725042396 -0.63582331548422232892 0.3203777394878516982 0.24011591123635300793;-0.94758619171406799353 0.27690291653427095975 -0.2651466350892436874 0.38292455784344830994 0.64973026128409472335 -0.28945860872483586412 -0.71991475437295271433 0.069719301044491921449 0.26865221857887444656 0.022621220346498424614 0.18722607568771199027 0.14544939317091490349 -0.16919870713398188489 -0.098574638730434857581 0.54275090413029170033 -0.037300415531223561627 -0.34397300091506388675 -0.53681159111070719803 -0.10876732349785336873 -0.5379884470683325981 0.40941927487049328827 -0.26776265568294982389 0.94497267546231877855 0.37511590844869369121 -0.29514676583971855006;0.35499812922813883675 -0.18733229144368143682 -1.1843755555667898616 -0.1471866151082991403 0.041560830639075452275 -0.16557885998623186885 -0.23079455177100649532 0.12470756469246691012 -0.2761224377149188447 -0.58608146527719973129 -0.21918592732643032983 -0.28378891882938805935 -0.4845383582310232029 -0.27052228500219316354 0.87029111023747807163 -0.36140796616413234377 0.55841522920862673995 0.17904101791259247056 -0.42828460267479595558 -0.13968725877196838669 0.66980088200625431671 -0.1714399944551272259 -0.20212744024035844692 -0.65969619716477234661 -0.37884631111693944927;0.31136825282706037665 0.57247414237712479501 -0.89578092289824917671 0.43620873818472993833 -0.45855897679277946022 -0.34791416391283191967 -0.098420631534030106335 0.028128757180961250217 0.75334918046447452156 0.4918939430085146336 0.443223638382303331 0.13232393677695461487 0.21606874366486270023 -0.65291073080293993147 -1.1521644362551730634 -0.53058379052385384878 0.64494788260649205114 -0.91683525748066341254 -0.6885912303556259717 0.59505741967701786255 -0.73560747974299456331 0.97438939522343914756 0.17109392885333438161 -0.20540970713828654581 -0.88425368597644216084;0.089451016684015643987 -0.1695511420986504092 0.40596047051112804649 0.6421641967184070765 0.15773344853707890745 0.20532874632364131462 0.53511078358822461443 0.43702456668131600681 -0.29110475147312309119 -0.27826480037763678954 0.18586998364497431258 -0.31738333272900559612 -0.031066007851231272618 -0.039434298094301177551 0.32896454307086553381 0.41029280868556244899 -0.25482867641367284994 -0.68847987616867167837 0.51221536834290404272 0.51453889267179131473 -0.2726976528703884961 0.18938477827360197803 0.18939862216835251507 -0.0040556641216533094452 -0.13721145597432604646;0.41475990388501227102 -0.081752627589223803017 -1.1022185242322473098 0.33280129457951301841 0.10196905289000653128 0.15979684681996450535 -0.36746486442015890983 0.28217118394805473791 0.054799331054003586061 -0.44716819034611215988 -0.76994950973006448525 0.11415761876400265518 -0.69948305474129746084 0.058942210293125220366 0.1758082329815526812 -0.5131004038249360466 0.28305549820556102913 -0.39124871794007709536 -1.0313181513921738652 -0.75498438090356001506 0.54494914034509278622 -0.52318281112269926503 -0.084131403147504682516 0.010017127097039807415 -1.094337574140940017;-0.16204078924307779852 0.086482172756444680717 0.14180108109835920516 -0.15750926912441551297 -0.081591436281147497467 -0.1424481991358903521 -0.37490993563784991149 0.30207402043999559416 -0.78457745464948891456 -0.26652705072546767351 0.061463066149930231652 -0.35243009357593924058 -0.17569601310881932665 -0.61745056185910263125 0.96610454386913158231 -0.046054990776277650311 -0.68353566248864494614 0.65936082242085047866 0.098334708320733280051 -0.41454567505716649611 0.71974421546222056012 -0.087603281561222906548 0.18487793702691041053 0.32526925775678705577 -0.082849822759644115022;1.041575684073517527 0.44267266856359532845 0.26317290318257197113 0.42589242653822634299 -0.35640319316396151805 -0.93573509719795566753 1.5175570181341664355 -0.81106710321265651498 0.72017102195460791503 0.34382592459916716532 -0.27994781053259210069 0.73568694290918024148 -0.57021857820746224554 1.1578844151013354402 -1.2587048502826301277 -0.33729123983573722212 0.40571341471978417115 -1.1984432527267683266 -0.33642931688724292405 0.57592073496744367045 -1.1979703022477101726 -0.52949291978350843113 -0.21867588651740949013 1.0098369965857580954 0.89075854387149766556;1.4385232315812432891 0.6855963237031790225 1.0103732221955032955 0.99990772951120565626 -0.28877958330010528032 -1.3133487942804000781 1.3529981085562852172 -0.63428146691852815842 -0.23015776286705863507 0.94045479962656419737 -0.46974253692639145186 0.55566365618773172486 -0.50592325419224026817 1.4424602475661776335 0.95985176802782634642 -0.27104945891870146912 0.48484656296516881469 -0.77507108782636424671 0.12523582388143267297 0.37244945068575302516 -0.023366534877242411961 0.018544058007567160856 0.86835722887854127627 1.2044172235351813249 1.3675265826080418297;-0.31269564974165187099 0.10314700018981434737 0.36584157390087468364 -0.18471001441901407913 0.11480444818796829709 0.02084415984856285986 -0.068326609627916590206 -0.16500444785587944918 1.5386111390845134128 1.1548614436328956678 -0.044830318874282973352 1.3655636568711591039 -0.52275441433233515998 0.24019227129171394353 -0.28828703566611546538 0.46516044515208848642 -0.6312317480693966365 -0.75107459514163210645 -0.12430685554099268897 -1.1328563413284638006 -0.21406446121976274721 -0.2619724031271934428 0.34684784846175176209 -0.31690954808785404273 -0.0028778344397017499579;0.027529898074452933965 0.087558848586120383017 -0.10103893202501748605 -0.094351297544893381053 -0.31527514655444854874 0.34897399412472157598 0.44597957288849038315 -0.30389560815082522494 -0.67084283930837018683 0.14342779795483723282 -0.8742917850107989608 0.26024184592458593501 0.22511451316779410026 0.077310256011213121519 0.10268792257184515548 -0.032968243843586778064 -0.20478440408111237003 0.18588614132519870603 -0.26413732741995143982 -0.17122288586444367464 -0.16900708596252528304 -0.1533336125239080916 -0.37566365329976786347 0.59959392198723560075 0.59638757048462820443;-0.21537743819983387605 0.27514007960817271625 0.44817628769447465098 0.15528552796182981677 -0.71490542745439622507 -1.0107550255069799316 -0.4220060283892365649 0.30348888806720325517 1.2101377691114514334 0.64205177172229443983 -0.0097559056641218673545 0.79835727621235164086 -0.35584145058092753189 0.47534232272248200246 0.38410871332752027918 0.30445507403293481374 -0.62215876167317862855 -0.022933932663316110995 0.2659984150984456841 -0.23514331684087413921 -0.20835502353710422252 -0.46095865776603378583 0.74825600412220860669 0.54286440405042257762 -0.31098909622152942633;1.1463713141960083508 -1.5513149540136006443 -0.63763332458193811014 -0.94635459130720223087 -0.63239170519906262591 1.3638016893018229947 0.60283210773369144686 -0.34733009848666646091 -0.94767081566462363362 0.090917494085232891665 -0.79326392051968186436 -0.46454270949629472565 -0.54136320492789802294 -0.43023754259978697068 -0.26352062738211651238 -0.4422890194926075158 0.72961436528005729762 -0.81766412467932159114 -0.93117019994322325349 -1.0412543598610231044 0.64853778326911637464 0.5577152245026439914 1.1840106738888640514 -0.080892989584930394265 1.6793005111011267694;0.28645607161089731152 0.081689218685691983302 -0.5362168699954550366 -1.0510317833598596948 0.073296605691805655214 -0.93394924990316685598 0.77563750914021911687 0.63398088682559416984 -1.4333390245573320865 0.15233453077191874692 0.85645149354691008092 -0.39770367089911906611 0.0002803455751670836521 0.39215044269138649158 0.6102223219447812097 -0.17646101902614311219 -0.62821181621128407979 -1.8997400628608511752 -1.0346426534865724722 0.35617255600335107069 -0.34438135234917999217 -0.12497818032374094499 0.67158833727148525838 -0.1535154091985841962 -0.25126423146266835262];

% Layer 3
b3 = 0.54321046665678640597;
LW3_2 = [1.0903248724323741925 -0.37396330709241887291 -1.0906016335365067427 1.1559736591748626111 0.40000626492901542042 -0.57601787099641454937 -1.1099227942646487932 1.8019116224387570213 0.77714678750131460649 -0.55779439595995405909 0.36601431035744380615 -1.0755855709470905079 -0.93278223638153578978 0.37401376103781636129 1.1533556387866037429];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.604229607250755;
y1_step1.xoffset = -0.03;

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
