function [Y,Xf,Af] = ESPER_nitrate_10_Atl_4(X,~,~)
%ESPER_NITRATE_10_ATL_4 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:35.
% 
% [Y] = ESPER_nitrate_10_Atl_4(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999994516886;-40;0;9.1619;-0.03;-0.52];
x1_step1.gain = [1.0000000017546;1.00003626070655;0.0153846153846154;0.000309885342423303;0.0656325902039552;0.604229607250755;0.0145232735458572];
x1_step1.ymin = -1;

% Layer 1
b1 = [6.5307088047972614575;1.9966686601451253402;1.1550694963150371919;6.8888546772115208938;-2.9352894487537888324;-2.7196159171058535975;-2.3334837795399159255;0.85288305682837572963;4.0993294400683399914;-3.3435905353995440947;0.046713246509972994325;-3.8132055492432539445;-2.5521282925657549079;3.7406293990872314659;-6.3842774030645550454;0.08885845488764676503;-1.1416672055839274247;-0.04980261507089494416;-2.8983252685603084764;-0.59780032313642195163;-0.50742899828619070224;-3.1552496816574278071;1.8683769931524676355;-2.8868419920849102134;-2.7536569776340407145;-4.9093308890558926549;-5.0098924216167262102;0.63731029658121440828;-7.2277055029772565575;0.67590667208582522196];
IW1_1 = [-1.9784170099670399523 -0.31290713167431644104 5.0578289609885871414 1.0419161539024748286 -1.4913595817992080228 0.072311302027858070929 -0.92139395115437849881;-0.36550418911935300681 2.1915980564668346986 2.2348419034522426507 -0.072421680382117936681 -0.031550150449865137225 0.18804886661772285406 0.34450066821471730938;1.6834226369682665503 2.2331967134599830516 1.0568072338767962037 -0.033365893094663846841 0.29112233884260463324 -0.91662588098223696242 -0.76457626281802548718;-4.9214701233116544543 1.9208782834260766936 0.70937375554484749518 -0.27694629123891145728 -0.65329730300248589625 -0.9488061203184381398 0.66095628863916366935;0.87756162339733534417 -0.78778746500726370972 0.64403806288732756702 -1.1941143308690900238 2.7776348040597111932 2.9298122881363228842 -1.4649350500295763577;-0.44428676898165436748 -0.55967854946762385371 0.36011073175860996987 -2.6088736284896274498 1.2099142238182605258 0.80694668797328861221 1.7065573801371580753;0.51313963391232653155 -1.1013820079891891712 -0.90086496075669886263 1.0520207482842931501 2.4398721478338707414 -0.10306056646821229461 -0.48911591169277962576;-1.0275313451350231997 2.0275868199210016485 -0.70898790181073312322 -0.17214488926346291353 2.6015529816604781033 4.1471015690844659218 -1.7457654427053175272;-0.4855219422276709218 0.45385231521545837952 -0.428994482146920042 1.7039949948614532271 -0.034665424973854608304 3.5271478325254070896 0.060303631339129970534;-0.099814372262353470711 -0.28461797912419800616 0.07024325861449308428 -2.5071298161182893338 2.7422572197856696086 0.46615857847765440614 1.5785686251208383979;1.4350495875721656081 1.2691186581243325548 -0.6850665457143784165 0.70120407505931892089 0.16274996357327814467 -0.56555211320634302119 1.7620171698450708409;0.55060008310347674509 -0.40393077486178696756 1.0563287725711565379 -1.9409263560213525501 -2.3829412729586132436 0.18693200758337666367 -3.9138432094477555445;-0.072835546499722653158 0.15783652476730217207 0.081372477812750518278 -1.0918750451294441106 -1.5361413634942753248 -0.84091388962782187555 -2.3698198081174082219;-0.18496726393367457497 0.21074497575103184888 0.12741402669676063497 1.818285751058553279 -2.9781224813003133178 0.58565145683079367078 -1.4381491570646292022;-0.46879254412925364104 0.92425638863982961801 2.2930419313069609899 0.38863700758171043637 10.076466218208917525 1.2567997093024563249 -0.99793124130408927819;-0.90632740307312864569 0.69763154296079243899 0.22721367274913359391 -0.16992630358101137866 1.0742436135490174287 -0.64998113942065438575 0.10874964718723803547;-1.254138255018246495 -0.9785334708552053673 0.3582098066859605856 -2.4316092884090219073 -0.41538733613250700127 1.8451821738905009873 0.24709485058505525901;-0.42474459917035967127 0.15834763698760528916 0.33309464013987000186 -0.048007486617406947005 -0.020092462370836221275 -0.3609776782481048385 0.84594901509715658694;1.9388075715236841035 1.21062878958228759 -2.0689841602593315173 -1.1165575967711189875 -0.021133030251119573617 -0.0097954810925255948045 0.79318177975609416297;1.3070986211786306264 -0.52963108977490247486 -1.0665152499650618623 0.60534434709227580296 0.83198530265925196847 2.6796032516861920669 -3.2256360036387548007;0.85751271056401656701 -0.57450073120186484488 0.15405404021270210912 0.3550183666504679314 -0.33240146219474880906 0.16403916919881161851 0.18987221922386840478;0.23319568388559730465 0.68281891724859267168 -0.581878150933390903 -0.77559618918349226835 4.2595519061797650906 1.4960372353099642417 -0.6514062299942906531;1.923386967913762291 -5.4850341927014074628 3.1160552723336065739 -0.72128911974210441205 1.5864870738054126509 0.49580397109834012381 -0.45468174303378477985;0.10435357394651488538 0.49563922376963459593 0.69149319510506623843 -0.80356050462774186727 0.58004478414550741938 -2.151555448126872605 -0.42853236408765332799;-0.070681589880859402308 -2.3829402737118057232 -1.1199083876428050921 0.18869200793196982957 1.242867819882715219 1.6471548379469720391 -2.1123353703536338344;-1.1950805176258199047 0.42650951016036686703 0.066614532761260641935 -5.404308848032658652 -1.2960865440529467651 1.3101151438805522798 1.4096568959978390545;-3.1779872928326904358 0.074633945877971757943 -4.2675741730143972319 1.5716593405649319504 5.3531626822421110035 0.3629774607532498254 -1.2813696808644090552;0.20069601976582418512 -1.1271540714873513611 -0.77244290177379959506 0.63303118350853748719 2.1214513658751719127 -1.3707078538886803276 1.0759074705419806595;-0.26990344594339865658 -0.30553402878008117938 2.0597670220825268217 -1.0368429591347272378 -1.8767780806338174582 0.81229476490645535769 -8.8030078751786504654;1.5680226647494306391 -1.4959342931564039159 1.4102978718224299381 1.1109270315313608535 1.5715907339270323018 0.72683520706478743634 2.569132988754189828];

% Layer 2
b2 = [-1.3326759972319865977;3.5413750080904877748;3.5879686436577427067;-0.71439918441939331828;1.3578625085428135399;-5.2932334474791957746;-1.1645605183695877383;1.3103413376935735268;-0.013589108287158686503;-4.2318107635883865925];
LW2_1 = [0.1336450609117040278 0.13873353134080534943 0.443471651599784078 0.052761859854757321708 0.19548692959347718778 1.7679771984841825549 0.54368942857816948688 0.41909972162157660414 -0.4606956356691808141 -1.597088323847460023 -0.24834660900444929799 -0.47946188688903990105 1.0129273039555402036 0.80804738730908354238 0.1474711993896982154 0.69704179918329778332 -0.2105384506935722988 1.9063130969910784263 0.34731880809924542408 0.074157606329619596419 0.19733079195926475635 0.63969473937618392956 0.45828087076069662364 -0.25727759373097364159 0.54891630271841362898 -0.26219781941325392971 0.23793184629786268136 0.54484794147219606231 0.39564298237165884675 0.23554179986787929235;1.1034489428798768262 -0.00018798330895337250768 -0.42686448437081414564 0.073641459616458349036 -1.405245746119446526 0.40240262370917589951 -2.4169393798399814166 -2.0689344986836228912 0.65347896390072479633 -2.2667826656899991988 0.2206319004115425142 1.6803191972702697932 2.3757873615227622643 -0.7309714559751839813 0.52507057211096108684 3.2025466487967113594 1.1046760311198604576 -2.6101703000313305481 1.4676194681764676009 -0.1755680772970682757 3.9181221045541954062 0.96498150017160599923 0.6347564006390857827 -5.1499711868614443944 1.5104577903159106533 0.46251242214826032706 0.87896150088436197922 -2.3244044379870851813 0.41135102421304142739 1.5388189790246118438;-0.39683733010454869117 -0.48098633555127956463 -0.53029938929206621445 -3.092736320791670046 -0.7207383303471527336 -0.091087352807346810835 -2.0432054529310073043 -0.15777436760851926145 0.96712430387149406208 -0.67981554785251607775 -0.57078348617943031051 -0.30936139861879391466 -0.50123024299635754453 -5.3094986204269911667 -0.59769885637810238332 -4.2657588019175864957 -0.46174966939338668581 -4.4200239862070542785 -1.3923174766916588396 -0.23197925705743421831 -0.63462390378666722057 -1.2559733870975697823 1.4688198918500421364 0.58130959562741890689 -0.21360332850680774364 -0.026115798826778606478 0.25747786189367610765 -0.77911957904673789788 0.74012816431744976597 -1.3818352077657314325;-0.45983979173785860528 0.052873855050872907435 0.52304398354699632723 -1.1283595059911046832 0.44289219503875115747 -1.9097029108182999035 0.54798957145247284561 -0.20339598059311844724 -0.060312126738326372888 2.348888090041058252 0.69747781994472624056 -0.3531223059995469149 0.1347761449283534918 4.212477853044875431 -0.029791622032993857888 1.3407894748208035729 0.063134990165187476752 4.3157362125270228859 -0.57579307232335485978 0.7784071581919352889 1.0954241085004219336 1.2021021198061638824 -0.19119720558536681954 1.1885780863449519895 0.82640803087362257884 -0.92516979352483941224 0.36617239205403689928 1.4578709307401822581 -0.26277143639992023294 0.057988533272378912919;-1.0659037503954942316 1.5441886766392651786 0.30361037852673627979 2.5699731747254732639 -0.17292459728931236507 -2.5750716310871091252 -0.87637751924916129642 -1.2435852755335039532 1.840540718854137836 2.9331872984489826273 0.48637649139111438501 -0.97939606825959046876 -3.4278699900011111446 2.6602558467639387807 0.072174676141079294189 -4.0602695076183046297 -0.8389003727435100366 2.4236979105886087282 -1.6136168395093786199 0.92543044374859595447 -0.97325046833782680178 2.9472805892445430587 -4.3822933755926145594 4.8553475832478509844 0.60722849475508600126 3.560675253355543024 -1.3727998726028767251 4.210630498573845415 0.90013255971300798475 -2.8658481534945829949;-0.29764443334715784406 0.71193960006991785505 0.30211413977686152332 0.22080451275971363367 1.7796260908989234295 -6.4778354753972031688 1.573009046997682292 0.21139165639677282016 -2.8972709457926559828 6.8166341729942212524 0.081323485229299039001 -0.88147155139964272763 -2.5717506827352716314 7.6616945483195380007 0.27640889168205556858 -2.9855815148324031938 0.75386221468940495249 3.6575010947280368612 0.27073486236027466401 1.8143224559035708854 -4.1382593581668993821 -0.389550105682682557 -0.68722779790902532859 1.802499224465343719 -0.39467944536853494242 -2.5764001324537919579 0.4801525394069183772 -0.74726956433754188058 -0.54615545162675482871 -0.34282362621987366724;0.0077140664400638441975 -0.34127332581286717872 0.12996851493297864599 0.24516031007034355671 0.27575391585372843739 0.63014434727637391731 0.55817857000739845574 0.55092859438006291928 -1.0874936740643614819 0.35587144208039650506 0.055296092599327620332 0.065014992858305434287 -0.26481778792669430356 2.4704729722713896223 0.28683691898198487857 -0.36304657530088946249 -0.015844998230958844826 2.2338523360033981646 0.25259759573472828631 0.085572593772448884986 -0.46280945843777165916 0.27893522594759734279 -0.10843651155307179601 -0.67346251613622998278 0.050892534109185379176 -1.042890652190436418 0.14725975009715266961 -0.2790348284852397498 0.34268322068301243499 0.25720017571743020168;0.71669872411577850713 1.1439097695349267347 0.54517381430241951268 0.037769533091836240968 0.43447682878632259973 1.8586829911021975903 0.55470391378538896898 0.25289522427650951153 0.21972784715211904349 -2.0773434881086374482 -0.50443578281040368783 -0.0653365661042945578 -0.42170007880763393793 -2.6069243123160217301 0.16959548098421595164 -3.3683891785066335522 0.082417716217715056293 3.0743450548527548882 0.97909209787290418792 1.5079634470953060532 -5.2180773369318673716 -1.6244433113459133544 -1.9963418026004422501 -0.032114027239258789714 -0.84157936338249084152 -3.5617621966024355373 0.010794294162544113003 -2.1421634290855178229 -0.79784416505149613474 0.20260678514418309359;-0.93436718085353687169 -2.0301681859641527161 -0.65933575435328939385 -3.1871716626959778473 -0.33897573935424551994 -2.0162512826054581616 -2.1934860061521499475 -2.7474859305959316025 1.4695751478479743835 2.4328087693233269029 -0.50218120216625683483 -1.2299723759079317364 -1.3503014624655809151 -3.4159686029196190127 1.8280459877964805315 -2.2348593418217301831 -0.91341965024398286133 -2.8097655961038738859 -2.0636638527325135151 1.1265434789754242839 -0.70827781446990267789 -0.52892514997009554278 0.80906078996578190132 3.2521167533918289294 -0.76354496905335844659 -2.176893343113906365 -0.91408220604983303481 2.0644839549133346956 0.7294443774830751126 -1.4882184022880966534;5.457326353338571856 -0.67455828408421325548 -1.8980613433300044246 3.6288285468589376315 0.79174053000828326088 0.52744900040540443342 1.5016126626741068684 -0.095384753124138105251 -0.44869451760635004067 -0.78690072256184562249 3.2649171104124006249 -0.38572518767581071009 -0.42357786958515031284 -1.9033728804437590743 -0.036161363144052841023 0.51503796203235929774 2.6200439918986373478 -6.3125142132621716584 3.4318287099891398029 1.455448164137235878 -3.6449653862306639063 -2.5190991469923615931 -3.0300392408146143985 1.229625002522556132 -2.6365864944589532115 -2.1524543780806366655 0.11134505188989148039 -3.6278309735990954543 1.4142824819595873898 1.3006395306080447138];

% Layer 3
b3 = -0.85147670475774039911;
LW3_2 = [-0.64807947069425708708 0.096635175301701833894 0.19641839578095046659 0.53432283219560516851 0.057894615938094827479 -0.19225951605410948608 0.76388606609519138502 0.2938514978830722435 0.12841280819810371727 -0.099571851767030081226];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.0470499670650231;
y1_step1.xoffset = -0.28;

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
