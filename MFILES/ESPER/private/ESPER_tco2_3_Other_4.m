function [Y,Xf,Af] = ESPER_tco2_3_Other_4(X,~,~)
%ESPER_TCO2_3_OTHER_4 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 12-Jan-2021 17:03:53.
% 
% [Y] = ESPER_tco2_3_Other_4(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-1;-0.999998793562952;-77.961;0;29.715;-2.08184650291172;-133.803853032615;-0.1386];
x1_step1.gain = [1.00000000274156;1.00000063368067;0.0138606862425759;0.0002963841138115;0.248478071810163;0.0588359962254402;0.00432440760873659;0.00826143816767911];
x1_step1.ymin = -1;

% Layer 1
b1 = [1.6475588843989881;1.2128911201452686;0.55133702148819641;-0.06449967008958199;-0.042236531714878678;0.54324930883095079;0.1889161822379069;3.7492686659136152;0.54998221558659821;0.58465282775578653;0.20237538556569989;-1.3581716422446306;-0.31331967650162562;0.085424184176045595;-1.2696969554612554;3.2067160728565334;-0.61242441654201663;0.25964057611426472;0.058543063516089913;-0.14302598868062227;-0.54373047267436658;-4.3637064611947221;0.78115322486587924;1.0073460789648356;0.2652099502336685;-2.3284976617632096;1.7807349987648249;-6.7340954269485351;-2.339238736461434;-0.53462176479247192];
IW1_1 = [-0.35394097951130843 0.5618058227240541 0.10325678955486045 0.12535684875232173 -0.37538600501059238 1.9974613716997847 2.229433025773119 1.6436003587338508;-0.022321991245551637 -0.27038284555221831 -1.186905936405753 0.26371201400174443 0.53313469308868577 0.34483239143262839 -0.44320582699948907 0.074546201511203766;-0.25388618819034497 -0.17969824568311354 0.50069890974250564 0.71867876118186713 -1.2261482158780117 0.49603208222239925 -1.4002809060184169 0.69490717370796407;-0.17746693895569812 0.5535202644170476 -0.16792134932791089 0.13316004092187167 -0.22216952024894521 0.045879984123649097 1.1549684034558105 -0.52270915403190144;-0.37169120296843844 -0.17950833535285685 -0.98076517143341846 0.46002255386075219 0.29527787403087447 -0.46081972547282507 0.33793577222529064 -0.19023726064754679;0.29970168830140592 -0.78782508865094414 -0.45240030408711862 -0.010343650486514394 -0.57538744428208632 1.3477548585658183 1.2102962357938398 -0.017472512463626154;-0.18698258593404551 -0.33105221802981344 0.53771343664170579 -0.14530626548206019 0.75242379214769328 0.92437099213249074 0.57243619764706666 0.4649212887726607;-0.58956270238223807 -0.1608536441296152 2.2812585259294615 -0.29095821700930347 -2.3612524770278749 -0.68045995692275696 -2.3522548238598846 0.81576070083375918;-0.37191031361990012 -0.42152709502795871 0.48016202693130694 -0.09116649959400791 0.99100041502909886 1.3133657510638088 0.28013997067611351 0.75419474410261467;-0.38241007822604856 0.22400104487588018 -4.5166557538663437 -0.68248815606098956 1.229504956972828 -1.6433652439429596 -0.36896901578807739 0.069720959514155773;-0.14130052013771766 0.2224288055787294 -0.96084072287632882 0.055560685731424503 0.15702781887956349 0.34308158292732471 -0.84416235963937736 -0.055741269774505553;0.044513626723747757 0.52483924231157553 -1.2599139178317 0.8994506079551351 0.11512016773583386 0.34828019683911943 1.3854484486100733 -1.6962786680774609;0.51896532039829191 0.37597199041489388 0.2506953046273464 -0.174290118343988 0.79912371067885568 0.60366869405344536 -1.052383225439143 0.3355702959556599;0.17347076313136628 0.14831949392690422 -0.70409535004417523 0.072559756973896564 -0.53968576710182115 -0.80292903980435348 0.16389233851198465 0.57423269437987223;-0.856807432042753 0.22472941851924444 -0.50877346445698801 0.091486499401663868 2.9168567148367166 -1.8399753414900593 -0.33047210746746131 -0.83365842622069519;-0.26332766928486107 -0.62406493396742491 -0.84699408275213839 0.68609040966528989 -5.3812964250992312 1.9847008600422074 0.26661451801335401 0.12823365102215864;-0.28063975973044691 0.31611051684592478 0.10720505615602215 0.44173086206989964 -1.4904356665492051 -0.17573434750434683 1.6354849096361423 0.16480959405889101;0.0086632588234119446 -0.15917272774946223 -0.13689517304182214 -0.22897391130366931 0.10749408983262766 -0.37044905521213622 -0.081562496706997903 0.87344580179234377;-0.11290637268920772 0.55862863773419236 -0.55837309204478836 0.1550913021189245 -0.7994396626791026 0.042386196476912136 -0.30511053565975843 -0.159805883872532;-0.14848863785957159 -0.13551558870281968 0.44295006645645474 0.3362381732668614 -0.60161757797848958 0.15586688068452612 0.8282005811921076 0.24576459043842377;0.26141991312659624 -0.92133107806305303 0.49978918724391119 0.04119504574126362 0.099750710160978764 -0.34681422632183384 -0.080064754919210473 -0.055193136515916029;0.056613440543034764 -0.21952036760950736 -1.1339199533682121 -0.5994402353606908 0.49920795859142297 2.0017361599048904 -0.29543834683031422 -1.6028986562796714;-0.344088691901549 0.010616456021716986 -2.7665285663401233 0.76996873284075196 0.55719099203039546 2.2999577727427645 1.298971795848914 1.326808398336327;1.3106453513804091 1.0979924335865789 -1.6178995067289266 -1.6059847328253571 -1.8752686995932257 -0.78718876793115988 -2.50713225727768 0.50810807745495334;-0.32457348983959311 -0.029721745344359854 -2.4999073076398379 -0.70827603546506734 0.40002079289800707 -0.97759005005130573 -0.66375956191936669 -0.12710089262570282;-0.056581895861515859 -0.072296718452854389 -0.99332204423944703 -1.1046430125268181 -3.84608977927883 -1.1924438712576613 0.47944041186326986 0.53116173941657352;-0.10055274331096611 1.0926248990603791 -1.4259080803342206 0.31922859546363902 -0.53944150619488429 0.76446573395055806 0.43893226325207046 0.4813365154767496;-0.73447483491840682 2.703146425883161 0.71532575666462217 -0.17406468004214842 0.94703370716039581 -2.8772213926586239 0.44344576170773509 -1.2239785838454886;-0.35523049754876551 -1.1599465472848804 -0.0039129665409612111 -1.9793235556449562 -1.3480262269317429 -0.040271093860804991 -3.2515563395385185 3.4974721316888693;-0.2300274472812795 0.19545843085664782 3.5294502833214629 0.049184504865885145 0.8404525053380425 3.2554095084852741 0.36061678632774691 1.332684159565497];

% Layer 2
b2 = [0.28547372735864668;-0.80804978760497337;-0.78207509242297335;-2.4022947743937806;1.2198318954090612;-1.4632837175871969;-1.8050568584019107;3.20848428733825;1.5913396284593959;-3.6124526584882459];
LW2_1 = [-0.14657090589793331 0.33025399639961267 0.079835109869120541 -0.21438570144780345 -0.33113481594153987 -0.11597286036200793 0.65082650988654611 0.20734062577303752 -0.17628043625866521 0.023140940140167678 -0.91607478416658206 0.48329534502291993 -0.077807688090650584 -1.1083363604829939 0.072019495956351751 0.039901458240230739 0.77341760858580599 1.6658453574844438 1.0754697696345048 -1.8450339428467706 0.9072513669881902 0.2708513292171929 0.029819271437285955 -0.19862819391970948 -0.24839903950516351 -0.15622062279076818 0.77033237189907644 0.40345426141476798 -0.48999160495414762 -0.091810497667889721;-0.1578295918753844 5.5577272639402082 -0.61858211381494455 -1.3280804687551502 -4.2554058319232713 0.4144473320968633 -1.131634281044043 -1.3762930061477565 0.0097523948743139521 1.6824642779198289 -2.8036320591475827 2.3670768502136257 -2.7506084693581809 1.9777714409173592 1.1612808997556747 0.40055057129541805 1.4083375467599695 0.12759084133840251 -1.1784514413381557 -1.4677225705437409 -2.9224762205802453 -1.7396804615256813 0.26187952388010244 0.62837777718653853 -2.9014436238592891 -1.522189781624445 -1.3838723875266905 -0.23054276192864165 1.5433403630880487 1.2082841376955944;0.49177540273269932 -2.9376048641589021 -0.58718054223328342 -3.0449078302431971 1.7778869407899918 2.5587989788647878 -0.81674484836542838 2.4359523345184053 2.3642685660365808 0.22914581813958843 -2.6403657626595116 0.1396436582105389 -0.89478012249210137 0.45535127154848543 0.20881367964392025 -0.81757380197486007 -0.68624635069684681 -1.3715651603078773 7.4471098846231243 -0.0025160740347284979 1.1398941865120047 1.8697916227997884 -0.46060667312601039 0.023311283869387821 -0.91348742751072365 -0.070957574089570791 1.3857983514117744 -2.6572116480688854 1.3593994923368717 -0.74613694053970858;0.65659453299197279 -1.5404852938222202 0.84634364332439571 1.4658827174853795 1.605434119495631 -0.064802991106937213 3.0365869732364317 -1.0706782838990454 -2.2929599882334233 0.26547893332621758 3.4894621266329455 -0.69363904656306907 0.39695436400964507 2.0697785308814955 -0.84346861053317368 -0.22024987724729775 -0.2590968214178927 -1.4964123160341094 -1.8943127464524561 0.10326415231992998 0.9574643980035662 -2.942717003451131 -0.47591007251823469 -0.069669327912947498 0.65979665195212456 -1.0884286763392763 0.22375285356132255 -1.0789423949409518 0.34620368777702076 0.46962211365806439;0.20139988602490094 1.4539443257337232 -2.0976502951067797 -2.5883299153282864 -0.38963818339355805 0.47751033485313826 -1.1882728587284976 3.2478649650201672 -2.1700686759782499 0.14822627966624938 -2.3786715500027209 1.7504724198215973 -3.1059327173927715 -3.376264581480803 -0.41153370421898677 -0.84890693861020572 1.4064213753193493 4.0758879499149927 3.5031778889175023 1.9571234039738141 -2.6456841580885353 -1.1343535748299329 1.9463282400514588 -0.083863116679679092 0.7125322007273116 -1.7293572980456995 -0.19438108393704884 2.4713619789818826 0.74830726889617605 2.6734447459299981;0.480755842504347 -2.1730406979092129 -0.34016656465427003 -3.1918235426421382 0.46045188183314045 0.39546406981145343 1.0962992174958857 0.76843002879938982 -0.045616481064233438 -0.26600407823815103 1.6462053139465476 -0.77286270418730785 -2.4340176725912372 0.024942471267345834 0.7048657493859849 -0.76587592085004208 0.1948264081000258 -3.2279022049210928 -0.96090356860335324 0.92963787112039842 -1.7822813248206613 -2.0858793121863299 0.479272627957609 -0.1530721891162424 0.0027156374884944969 0.25681005104661758 -0.7381146371765458 -0.033453371438083503 -0.46171508499981062 0.14430881072473264;1.1395002908013563 0.32031392472474507 0.26242288944167513 0.51439145872354397 -0.71329062548424738 -0.87868451306736028 1.0878538589302136 -1.5780293360823758 -1.0475876057739872 -0.05690713964114541 1.3786428254620531 0.054095379575175219 -0.30358869519438031 0.25776773497639927 -0.25474140557571251 0.33005506665251061 -0.32416522038335005 0.40839605441621191 0.85776471510599084 0.10374299114296359 3.8410681647437341 0.59329993098033829 -0.64950599416431265 -0.44434416288350509 1.4680260263005269 0.43709578276114464 2.3182466848070571 -1.4156384306197156 -2.1249505331823899 0.70310024971233898;-0.310933599980312 0.072995644600655024 0.46711873269635662 -0.0061063325235463806 -0.09888888006123954 -0.26889591017767134 0.78014908487266066 0.40906893297275404 0.12460759382585396 0.20756615578758408 -1.35321499293189 0.76401521844061937 0.2098633744320727 -2.0106409269527763 -0.051389815173228465 0.096207807124078032 1.5748945264411944 2.6642140938597989 1.2926872669372915 -4.3205751304871427 1.6015042056651949 3.3565218392629039 0.067507802852600371 -0.37212689734659243 -0.54679414239400281 -0.43305225840878081 1.255301834803503 0.84854351357210633 -0.83707707164074774 -0.060201331109790882;-1.3658054287410273 -1.7446041489605135 1.2580383218388238 0.93545157478774399 2.5603085879032674 1.539507203646266 0.37740081275532544 -1.382177424755537 2.0499433864644727 -1.34735620886406 0.31017696737517259 0.65010797513407159 -0.49696212272896967 2.3883899324203663 -2.2873840736450037 1.0554255172396672 -0.0055023222654801817 0.45500876564323595 0.92329464202969391 -3.6634898712696486 -0.40175280443861128 -2.3890735921088972 -0.25994119102908664 1.304470370218985 1.5916825867484896 -1.8663025411061633 0.074751563685003894 3.2890600959762448 -2.4505681301792293 0.96904933893028733;1.5326725389825113 -0.073683963708348135 2.4669499754005342 0.56425773268230361 0.64062260549889316 3.64321737679338 1.9850371513386305 0.1568908607125668 -0.48237978712671142 1.2190332193212485 0.54453339067563267 -1.4836908876380004 1.0738169841179581 1.4126708207215968 1.3602835543950109 -0.17233713182659144 0.1272618637413433 -2.9335919113976945 -2.4296813638543968 -0.054303643733255885 -0.22353359739123704 0.1890744036072726 -0.16420503521295771 -0.56627442354483271 -2.6288642289972159 0.098881840532217644 1.7294366820643949 0.57106271458123681 0.6696794026144951 1.1701193118842272];

% Layer 3
b3 = -0.56352424226463271;
LW3_2 = [2.4031131679932001 0.10975522411197802 -0.081419392144948186 0.25038977832084602 0.069303431242541513 0.06636784267480586 -0.13687563357416527 -1.0734696055071837 -0.039088918396055071 0.032517212787419668];

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
