function [Y,Xf,Af] = ESPER_talk_11_Other_3(X,~,~)
%ESPER_TALK_11_OTHER_3 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 12-Jan-2021 17:03:49.
% 
% [Y] = ESPER_talk_11_Other_3(X,~,~) takes these arguments:
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
b1 = [0.03479285192079451;0.094546511506201705;-0.65199845618563934;0.28958034402310889;0.35931047073574529;-1.1198883599978271;0.47201969557305756;0.036418411209768214;-0.51289923308981544;-0.36818133333563136;0.26812348218647186;0.28871227188183518;-0.75920336059109761;-1.7988416526289079;0.01331184079080203;-0.58161025531300714;0.8445646869022353;-0.65240193446919137;-0.3647862853556319;0.29151989478627449;-2.1377361877390153;0.48042017473297888;0.0091077871248508337;-0.70904668394585713;-0.082922040651872919];
IW1_1 = [-0.15450318060645596 0.18163540541221607 -0.44028466230417912 -0.31647768805264731 -0.13349055395732898 -0.18569866888644812 -0.27533678695668939;0.019619637323558047 -0.092349277977469585 0.20083606598748893 0.11898383831869178 -0.97045645415393522 0.72125127329786842 0.23771443117887556;0.83053185476333957 -2.0337222654647142 1.7025002050959301 0.68926665848201396 0.66940911945264969 2.3081571306310575 0.27258067129597435;-0.69044427066596925 -0.1126826072834374 -1.9559117774993664 -0.35529394722977903 -1.1858354250477117 -0.40313246405685943 0.25352276410421781;-0.10369188476545925 0.040936714174877764 -0.26697681391110756 0.095107840077190836 0.70376509877838334 -0.31743936883161711 0.22242199757380704;-0.6872358264937819 -0.65609599347817738 -1.5212445376635086 0.057204746536978109 0.58117778590231439 0.37564383587177702 -0.77343530505385572;-0.71353868196026604 -0.14319331061740928 -1.6398444568212367 -0.16890068434502523 -1.4470522791441409 -0.40835106266171933 0.03769545292211731;-0.043083061562932672 -0.14329516791822336 -0.34073157616707345 0.56039466436582674 0.93093020134874471 0.47351342596903678 0.41615584532426864;0.20391961878387044 0.2464128526686982 -0.49081888486574793 -0.20919807535783186 0.59617281737476469 -0.0075978909655980958 -0.89409365062941537;-0.30126706466721803 0.096391052061748361 -2.7891340818743533 0.36904130143188146 -0.57203440852929666 0.045536559023839934 -0.42456437371146949;-0.38502621803996545 -0.54879167230823866 0.31183102535941987 -1.152143019844335 -1.0740917730442623 0.27603567175284255 0.92180647317954834;0.14286929038164203 0.02346411035388412 -1.1420898444876046 -0.096938507362063275 -0.28420180577973986 -0.79313081268206687 -0.27296194736022422;0.12984818896519948 -0.30198982395606833 -0.041069514351821268 0.10464323777819995 -0.49245025786388724 -0.34733218027192364 -0.12712526332523896;-0.21637082739788999 0.23065103089763006 -2.4751625206235204 -0.17436468892253995 0.13993750922074022 -0.64390390442557932 -2.0852692525028891;0.15482330471895811 0.086577274319057096 0.027661170956973673 0.48387913095020263 1.0600435476926131 0.44916318425175894 0.36321903972865377;-0.12532543627223153 0.26668420040770224 1.0961110085688819 0.036958666747662371 0.83188892344564513 1.2853372676784989 0.14975144124946416;-0.2307477088281113 -0.074494147923434909 0.019583132012624408 0.078301835255643687 -0.22662831521890039 -0.58012024262346429 0.94297940959969784;-0.41339531626872389 0.13405540151591075 -2.6864117847400215 0.34545455392801333 -0.12787401700142695 0.046505879302704944 -0.75215615013929482;0.13265025956825596 0.035110688552683816 -0.4477109082142679 0.18889465645158851 -1.4493938594135249 0.78529734660541228 -1.0996137898248914;0.071668554521581718 0.025365940839381589 0.011772527716212105 -0.10176087914612708 2.3977455123491342 0.14873351071568672 -0.99989293499551124;-1.476429541542289 0.063946718403238442 -0.56107737349690756 1.221923198495344 -3.5419745237791824 1.7762539085680564 -0.87034655770858882;-0.049046080991389913 0.037175744197983064 1.163315460597385 -0.13623195218231732 0.6517862413395934 0.86281526895781702 1.0243796212976173;-0.23280307853887913 0.31257646571247955 -0.53404979661281959 -0.55614962256816147 -0.14883091032127727 0.29115583977411363 -0.36338124730231414;-0.18495412185709487 -0.070327806120204603 -1.0673839216299061 1.2624444184325101 1.1851134726521031 -0.39121327595726552 -0.40475991549560353;-0.089554285443099693 -0.13453276842322684 0.071702999345365601 -0.022062675857744383 2.3059546101872863 0.31897277461217161 -0.44794429814234377];

% Layer 2
b2 = [5.825304307223079;-10.926268281132828;-2.4185560578828103;-0.10032564737193377;5.389847159604348;3.2653756900540136;6.6566088299415398;-1.1590699388896288;-6.531250684154605;2.3704025834328779;-0.82444437440879115;-4.0174083201523274;-0.2850859048752285;-3.7946179017295227;-2.8419947543003587];
LW2_1 = [0.5663188798047375 -1.6124957423442778 -0.23546736449191205 -2.0209533781119027 -4.7664612211434179 -0.68341828406898397 1.3009323284337082 4.5596084851136176 3.5518057810760881 -0.95318136398524145 0.25274674230708227 -2.7079186754603186 0.36653277407926232 -1.1382366906641819 -3.6732185226584062 -1.2408588810160468 0.16917926415202508 1.4902148640074073 -1.4123243446994924 -5.4396887896652366 0.37760998770974485 -0.81407142969977275 0.69965140278652116 0.29907319020486139 0.50081178576114604;8.0273861514921609 8.0109846344828917 -0.39673788161729867 3.5660060204679871 4.3113880881003404 -3.2948504861194507 1.5655098554269014 1.4372608531500326 -5.6091490882962969 -7.9342352953896782 -3.2657820016838985 -4.9845902965063127 -11.139366760423483 -0.23586062348793246 -6.5486542731374051 -0.72666669939351491 -8.7996181308797468 9.6523981440316735 -4.5253918175810268 2.6031641766737592 -2.9069390826694175 -1.1225513703581584 -8.8026859551411114 -0.91157358486764806 -2.9941120045168796;11.47853932075575 0.23448472316995139 -0.27740296187213109 -3.9290128416400756 -2.5017948911467456 1.959756504688869 1.270426132636963 1.6901392467210188 3.7408136778237293 0.39952393774598277 -1.1549475349932998 -13.323928162756406 -8.19724273509493 0.3920722071759824 -2.1847493004894787 -12.063816758187089 7.151063667532112 1.8859660121065513 3.4042581250005512 -9.9371556413512359 -2.2800606141486179 4.472273401014994 -3.8464938226738044 -0.75304903998414496 3.1136525788953633;7.5350936228554373 3.3152761607091641 0.68500779680329493 -4.9741188911304342 -0.44066695273761686 -0.37766971267596322 4.8754137580855943 -2.7487833168942362 -0.82913544575265996 -3.5404692358537599 -0.033869940104663682 1.3765277536387748 2.7454583780601358 -1.96615806826217 5.5021219524553269 -2.9144078822376418 -0.059289920089519194 1.5971515467477115 1.3329307650492757 -2.8996320469296997 -0.31785946701027945 -1.4060232264292416 0.11061432480656096 -0.69042877836313232 4.3637041686835207;-1.6655388814375176 -4.2532881115061247 10.240615319070461 2.2080816054450509 -5.3838988581849829 -7.0676209325710113 2.2769177782552639 1.0645365309837365 6.9852917290663861 1.9442538478853943 -9.3452462401917575 0.97575471344662934 0.39119820598739913 -5.4369315236430813 0.9897762139334475 -5.3059782421553736 -8.1547033878005717 -4.6415866282467571 2.3413882791749567 3.6989771701677268 12.008809096453364 -4.8708911472003624 -9.4841591324507952 -6.806016392597992 7.9898133347358291;7.6497157537250891 5.5309006037774084 -0.05037733625651615 -3.6513459536449053 -2.3410647185492475 0.71041940401448767 3.6006057998014676 2.9984476492643473 2.8963785455126514 -0.060265506710991401 -1.2756362614057351 -2.7324466978568673 -2.8761607787934915 -0.10288478740442622 0.086727414206689474 -2.1308377592309351 -1.9024744984771007 0.28292727907524395 -5.3290538002848118 2.4879211260506024 1.4307852695649481 -1.403210026469544 -10.832646891140204 -1.5775224396555276 -1.0783237611429652;15.041025177018456 -9.8125892697158701 -0.68210974630928134 5.9274034620282654 -0.30766194893193638 -0.42085981046963283 -10.142780399595118 9.2467050509720981 -13.328260134869533 -5.0853690298746921 1.2852188573165855 -1.3985215396630324 4.4192522912572327 -1.121448468386685 -2.1330364943762548 0.56331669370088588 -7.1417204541625283 6.50045872583775 4.1194323071865702 4.6427569599755278 3.2486346435140736 -3.0664849419741835 -4.4532702707021805 -3.2539191124498892 -3.9641056696033936;1.1111875494251358 3.5855437047576282 -0.13618980888262622 -5.9634299725874911 0.98214593338203504 1.5219254447532671 10.19282360377276 1.4540238076296581 -0.96498499239674251 4.6292891668127156 -0.339030413982852 8.7340500839359958 1.340288135986105 -0.31401638317196395 -0.94272220724857503 8.2915892383279761 -5.7304799840650418 -8.7831351391371086 -3.8109859222924491 -3.7393488942760182 -4.521641092980337 -4.4017485030810004 -0.49370359738088543 -0.67750500883692588 -1.1326134246573978;-2.46487662018043 1.2571817715059457 0.017773600310566122 -7.4162856667093831 -3.2860697214920602 2.9840905304075505 8.7954186857075705 -2.2334343343782312 -0.33798268751303528 -2.1311533154675089 -1.5494407225136573 1.1416457090398127 -4.314362881603012 -1.31977572897422 6.4114942118607292 -2.8115179581701084 2.4141088844129581 1.3874195084758898 2.2762976212189212 -0.6545453166891364 -1.3489667620968779 1.2336037769296972 2.3495980139869443 -2.1677714473561061 1.4337214348150327;-2.3401723749435441 -4.2512662538493213 -1.3192725600876809 6.4124994326938225 9.519640050707288 -1.3121241461633883 -7.9728382853523758 8.560974895785403 -12.688969365026955 -1.0981855607079309 2.292131097436104 1.2780386747597092 -8.50323281531284 0.83344310303105495 -8.4646170470748547 -5.061041531562438 -18.704458067271105 0.22997625252432402 -0.59168366297652042 4.5702898235388334 7.4558258952142635 2.1089119575001578 -1.0406222301764547 -0.08433953892415412 -5.2856550441118832;0.37688028154488756 -1.4543798282632103 -0.088806327070689958 1.3796227144724691 -1.1102935648942218 -0.21057362417283854 -1.214942200817914 -0.33604423513054738 -0.83994769138060077 0.019433519763882381 -0.058889634904145803 0.15447589157745947 -1.5278784136229258 0.49718181563226926 -0.25405650543046188 -0.18918051171927072 0.90771496626364967 -0.123283295748675 1.3714564150439834 -0.58111423012815033 -0.16061480233669609 0.80104985144887508 -1.2452266864401789 0.071777415122132107 0.71824617062578722;4.345622996275142 -2.7003511989175788 -0.92618588127062629 -3.7737910165070758 5.9538657134042321 -1.0616164345124117 2.180057771795362 -0.46344779969413308 -2.1108207017116918 2.4514267549351763 0.60162870557309722 -9.2748433050475949 -2.0528516692623531 -1.0804746767551898 2.5050183366661347 -5.5269244187113609 0.80226701359569619 -1.1091295142132362 5.2570752848807798 -3.7649450336579058 -2.1048537060544517 -0.6524620883817891 0.40341523095161635 -0.44046637928957111 2.7605219403497481;1.4250042678590349 -2.3095215621214877 -0.025058979463892548 -4.4325591708058907 -5.2782137325441321 0.16292327910086241 2.8179652803913284 1.9649432119005519 2.53351885067958 2.475117335434077 1.1733571008933641 -2.0738164969123556 -1.6883322036581134 -0.54342303436444195 1.6918402114867481 -0.71969589542251056 3.6152959000338338 -1.5863753770434648 0.49175487329808887 -0.74513413511750237 -0.2404997761879093 -1.2922560837081198 2.1279047570714775 1.5457805066092882 0.58210693637696098;5.7600318433142865 12.282198001064758 2.9307805772701863 -8.2723273062392391 8.5489902962182089 -2.2038642860993787 7.7438191773982243 2.4586430367877368 -1.188135712684713 -3.7928242312822933 1.3602802480219702 5.0509579164732079 -6.5069340647157468 -0.26422114677781977 -5.9019554455741021 0.24728693912771127 -1.7740717790983629 3.4568782907876048 -1.5514831802403632 -3.0860552834988324 6.869116904725753 0.67093028696615598 -1.5983414471215012 -2.8547721930581247 4.7857085977784806;-3.6376163331030593 -4.8094810499195582 -0.31764873928899467 -3.6557291596261448 -3.7306739080969238 -0.4900519941696711 2.3773622094424844 3.063627149659526 2.5929964009620323 0.46011263907771044 1.4130609089292765 -1.4165660023289786 -4.0105255444411192 0.10831243200351914 -0.75648014496518023 -1.6382160730711941 5.4287601729599784 -0.24940565302811002 3.0224774622014516 -2.8427184016400298 -1.1383881333896753 0.24711174842628866 2.9101556914846638 1.5760522694094268 2.0367052221213235];

% Layer 3
b3 = 0.2929624359278496;
LW3_2 = [-0.31651551816997059 0.037890057106051513 -0.062157812468828164 -0.12811834657259927 -0.0060096619956404706 -0.045739835720013784 -0.044424613044727196 -0.06373854331722266 -0.075261563510004945 0.033901614605488788 -1.1871481732736324 0.12156609396964456 -0.43145191423980256 0.024548348567208122 0.32824067324344902];

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