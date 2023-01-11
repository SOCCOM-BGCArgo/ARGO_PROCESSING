function [Y,Xf,Af] = ESPER_tco2_7_Other_2(X,~,~)
%ESPER_TCO2_7_OTHER_2 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 12-Jan-2021 17:03:55.
% 
% [Y] = ESPER_tco2_7_Other_2(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-1;-0.999998793562952;-77.961;0;29.715;-2.08184650291172;-133.803853032615];
x1_step1.gain = [1.00000000274156;1.00000063368067;0.0138606862425759;0.0002963841138115;0.248478071810163;0.0588359962254402;0.00432440760873659];
x1_step1.ymin = -1;

% Layer 1
b1 = [1.9445585912447463;0.73325737578994066;0.79244793557291282;-1.0422218661111915;-0.38947496675527121;-2.1360841337380099;-3.4736258719104618;-0.79803376452084551;-0.90537976233567174;-0.26142166302858588;-3.2767598273160612;-0.21444472654665744;-0.00099540122087697397;-0.50947192610737591;0.72335152937061231;1.1533130680863317;-1.8917584205925364;-0.88351719445830479;4.6286507985248813;-0.96925723103902395];
IW1_1 = [-0.27821364356915734 -0.2396377361648023 0.91271610675576542 0.0034119882651295441 0.60330213903022178 1.6088400885791998 -0.22429943282015646;-0.092802775425881737 -0.35015254546772073 -0.22951466333714818 0.10552712910929055 -0.56213322628979701 -0.96393234340891532 -0.59292941687325895;-0.10240278812799855 -0.033809441443148049 -0.23965650038415226 0.02093276002837308 1.4970424510788283 0.25946837878230322 0.45389019547464882;0.11068583444349564 0.14806651246407551 0.5288699665114176 -0.12692398214039372 0.30445377218987207 0.83393146502564419 1.3656587439329289;0.13207624575178933 -0.13334703668687614 -1.5328991191038785 -0.22765434880159319 0.27630156573154963 -0.55456743506881467 -0.73041645079771123;-0.44906140401315492 0.3723976980499572 1.5436038622794057 -0.059109203590559801 0.3580841030008825 -0.33628074986844297 0.018254650307639293;0.39667941791368755 -0.050312036316387337 -2.2940285420936997 -0.55043927345663468 0.17039776991143699 2.6348955028708203 -2.4406886714272598;-0.1189366064739534 -0.6372408054002372 -0.70022533797653996 0.097340977699924836 0.35629886849486375 -1.0943788460852746 -0.10930700196133546;0.13260252080534679 0.0088404842663757817 -0.39991189250415854 0.29432897015256426 0.6716534310864466 -0.97928429237217263 1.0262303382212654;0.16173216240183758 -0.17066523476253223 -0.59219405094303024 0.25289254002709155 0.036229972480133169 -0.55065104307173751 -1.0183043758120225;-0.3647247997396818 0.28108006801787644 1.8373928565188136 -1.0862097665887143 0.42395479403529346 1.3150969164657635 -0.23195427711108116;0.077074702123736433 0.34788799681576071 0.59641665619236339 -0.32359973451219459 1.5698999320595248 0.11288376703114227 -0.12535813253277722;-0.96709944967826911 0.36255494457369219 -0.89586391474483029 0.18789668282202945 -1.3094712746406854 0.60376896401622915 -0.99762687584322352;0.051579419405797494 -0.16658722058645936 -0.63105560145719941 0.32013482627568218 -0.17566492997541042 -0.54652830680979947 -0.13434206945624325;-0.13697582893402485 0.099588653525914783 0.18408580759253912 -0.14100885464310212 3.5322600625425782 0.3791288311613592 -0.65105664265857099;0.00042389620629296057 0.020691145015814458 0.88412019839487821 -0.15942653795983319 2.2515422910174743 1.6109484854716831 0.21425721763601985;-0.5751313693743425 0.54357323901470156 1.866843496560348 -0.038589565938614082 -0.30514850975478225 0.026074612705399319 -0.50610420291734048;-0.30720699236232496 0.024805684196673006 1.6597480787423098 0.15867659014645635 -0.66323973925478608 -0.19119585730001581 -1.0785281861055733;-0.73473783618670074 -1.4182524974330548 3.8396727383566591 -0.65269180655503456 -1.227050943855581 -1.79914775943455 2.7834809631896289;0.0029104900156310589 -0.10930202897919754 0.19881679757002738 -0.095616044189328481 0.11469085493119971 -0.12160091892962258 0.94587559283791889];

% Layer 2
b2 = [-3.3770592896708536;5.040178669082974;5.5361563735105577;3.8586013596431545;-2.2261997876220923;4.529002320259421;2.5842212956086197;0.97081254796038674;4.4422134613637532;-4.0200011496434946;-3.432810466609912;1.6146583790193534;0.20896588128929275;10.251900361331122;-4.1579156855331796;-0.66316673864482811;1.5844483643064291;-3.2384234513700232;-0.84011686864281199;1.0306184747015854];
LW2_1 = [1.4359550814747424 -0.89094443502941034 0.53129770392713738 -1.229254928706919 -3.7310213703812733 -3.3683272439289285 -0.28040198503362168 5.7392226965764008 -5.3667911840002072 -1.0043947424946555 -0.74332056738564534 7.8414233320634574 -1.3009300638064385 -1.242603456977714 -1.4654614562093082 -6.0499371601253538 1.9925563995701767 0.48759194247769494 -2.5930864372226377 2.457933232618744;-1.2879210862887307 -0.1147314758480999 -0.63566342320483693 -2.9854523679310554 -0.36469196269868587 2.0029887750099151 -1.0982051057923901 -2.0651773103634694 0.56626963396607055 1.3434354325182705 5.3256350300776996 -3.6328243076399898 -0.76881980370923864 -3.29175582083082 1.5620110031405681 0.12387998094988888 -4.1233846628241064 1.5186869687966722 -2.6008224025310152 2.6683942737193651;-4.0110370018077424 -4.5766121953322161 4.7442061645490092 -1.3845906310425828 -0.034716265878011854 1.5264536735750445 0.21813107040403384 -1.2508416640775037 -5.523914835012449 3.6743955521029927 0.18419244841281557 0.046876212775496168 0.2521979958057739 5.5609172853831739 -0.23089569983930769 2.3659938592733902 -2.3951840825245805 3.9967923787329553 -0.33208666324474445 6.7106589447532254;-3.2262415586487392 -3.3673658515870719 -4.8119959529045762 1.299657900739225 0.11181483135054307 0.26541243129737341 0.70267282961517641 4.7910127648864549 -1.1077560299753724 -1.9159757400642738 -2.0360375536723079 4.8323235371393336 1.0804081697464314 1.8554200236674585 -2.7444323583944468 3.4742490386092495 -4.2112470790005512 3.5935420708755088 1.1044908428844045 -1.9875132433550009;-2.1751337791940646 8.607837312690064 -0.23266440828617674 4.6777864421751199 2.4166189938466007 2.3509520357681568 0.073953365601158713 0.082161870214579866 -0.079599718840812406 -0.96142462208464208 0.93492362986433142 -2.7680999548934517 -0.25471092098801706 -3.1792438834072207 0.54717422601266419 0.75701168836112009 3.4358721427597834 -2.4968969249361765 -2.4589049321944114 -3.2626007028458619;-0.11324478435122995 0.45624311448645288 5.7301363771417293 3.4515977954329804 -1.1487099922239739 7.878781305413467 -0.84708235412805533 0.30728995040823653 -1.9260501979215521 2.0002215577219435 -3.085457179557288 -3.7198965278640217 -1.4632749507547849 -1.6136585957827037 -5.7297655969002834 0.89442557134634482 3.2957065416324056 -4.0717337771019846 0.99872536385580524 -1.8035583392199255;-2.812735364775552 -7.0675892989563058 5.2432355919629421 2.7511693408856437 -3.1713205992912794 4.5433306996499354 -2.5333965113601415 1.8957143578208675 -4.9134695457178248 3.1304907852508896 -5.0581199563032957 -6.4227609654792488 5.8678339604248633 3.2290532878200482 -3.7803082748025223 10.370933201737856 9.0760244143463904 -6.7387119311295205 4.5142418979197734 1.8744236787264918;0.6939625164981722 0.13494168445786495 0.84097737157542196 0.52237518533765981 0.61405631862630061 -1.11938344237131 0.36018373972763967 0.055147262788941095 -0.051660694718357394 -0.30113424849037318 0.27845254477557368 0.89587741964388201 -0.13147810737165752 0.27286419395929667 -0.24635136446306471 -0.64645299106830612 -0.34135253007452571 1.0851192906753211 -0.21660266444465454 0.26480509004861019;-4.764881229126849 -6.0092402827458837 1.3340577615947722 -5.19744583152763 -1.4611314181418324 -4.2070618171802332 -2.0898288914126111 -1.6062851061931973 -1.7780827521362501 3.1074860258737513 1.0218111541694193 -4.1538784416125436 0.84898849434805956 0.4256321545807662 2.3325722683465999 3.1111623209284835 4.2450870013598045 -1.0219461684477422 -1.2276369734594037 9.1159904416133433;1.0694009923322809 0.69200119126507786 3.6840384760918701 1.7135538492086151 -3.3033642664315503 2.3363376287232507 -1.2777657459650007 -0.053243068311833217 1.0305844073350221 1.9280361451922992 2.1547639979298241 -2.864543516963074 2.0089898788545097 -4.6161312187078574 -2.8162391498849848 1.776055301708426 -2.4471933825631544 -1.790108776686834 -1.2961362364694142 -3.6135511847477257;1.8558227139223038 0.92654468548131741 -1.4644101083293668 -1.4041427446481327 0.86631505079754523 0.66621769999640013 -1.6267267676282349 -0.65627515388862279 1.5205636813957693 -3.6320466524557564 1.204518725307687 -2.0566539743836887 -1.7131423493355207 -1.6102547549920507 1.0533174938867356 -2.5497672249460255 -2.1655168178463242 0.77895222895500726 -3.2134485651765257 -4.265162983284565;-1.1739390043708702 -3.1776508113113549 -0.50167068881361176 -3.5453171365774687 -1.6661886821124074 2.2352446480921873 -0.34874976081971248 1.2753352721639359 -0.78569109066400433 0.15976017848326574 -2.03984967499321 2.4038586795671404 0.39791436756421639 2.459768961864869 -1.2267934398292442 -0.14763529890371971 -1.0126294368350821 -0.98432260535919791 0.24785992392479958 2.8445467320320366;-1.0031258120525557 -0.3222505323610706 1.7088382536681068 1.2088834487677236 -0.87526937917724501 0.86951684606222113 -0.21388210369711252 0.2427976989819115 -1.6245794442062815 2.2368080527713432 -2.4317239696953332 1.2155235816474201 1.3031164756885953 1.2591052924958559 -0.60600985985815292 1.0264499104081681 1.0888032615917476 -1.2087752344641649 1.1828019415907387 3.0116194053711758;-2.5713743654366241 -0.019623166532195122 2.1075206050383453 -2.4480441870772709 -1.5640560938383499 2.9945046219170033 -0.95166407888679128 -1.0606811024209506 -4.5329221229864247 -1.2588545418194586 1.6554870428757416 0.33321526391598943 0.072783149710545378 9.692561744022413 -1.5352500785667815 1.5507433610302293 -4.6656665069071321 3.9869266824334235 0.27032110737693871 9.0829136102596308;0.74329260292036337 6.7325425998962869 5.6471369065386012 3.2905102157690851 3.0058350344826446 -4.2281870794057301 -1.4639405898403512 -1.6302146821738304 -2.023140796266651 -1.1905522546449645 1.3148697058916803 2.79523472260867 0.90304236238688029 2.1053752394009875 0.0015196441564594848 -0.26227289086246225 1.8095518612355437 2.2176683140317208 -2.3744765149096265 0.057089537374390063;-0.37147000758616344 5.531327547511391 -3.0144780120246986 0.74660331890109188 0.025866330923045907 1.4221008868390441 0.76872240693891492 -1.6603242041519375 -2.3227573216062787 -3.1490009189899681 1.9200787555144978 1.5528245018493434 -0.77226091269078523 6.7250922195687659 -1.4338599835269323 1.5012717360265613 -2.3364476858248002 -0.88289600350423292 1.7882518457245982 2.0336353963726652;-0.99232784112923689 -2.6629596221176941 -3.4944928878455856 -0.36719670057158971 3.6939566522672633 2.8448944312288473 -1.7795227394160733 -1.1957319391948764 0.75326727349088507 -2.8593389368351838 3.3863766338910102 -2.6092554255351716 -0.30098410664974545 2.4573572910273387 3.7930718307607014 2.232276499884422 -5.0278061878476992 3.2447133494800946 1.2817817117645098 -0.25142219034526164;1.1939029217941959 0.80328398857228378 0.087924610353188826 -1.1682513148983531 0.75496787566566337 0.072630682638718455 -0.89217821016358867 -0.30644727305948177 -1.0650223582281784 -2.115403594845751 -0.97938988735205812 2.1442883146294451 -1.0507427062842292 1.5751958554972321 -0.61083792069265797 -1.1569546817955201 -0.77911726754561261 0.209723832819078 -1.3796056413080939 2.0262123815944135;-1.0574899362962031 -0.11892960748858974 -1.4809541981605303 -0.58939664063584218 -0.8708972133822005 1.4038401657219903 -0.41748920451146548 -0.085249050725757697 0.083954116454900418 0.24575018989758868 -0.63186286682958104 -1.256868354493335 0.23849705667946075 0.011919935676851293 0.5048400089688676 1.1304080801473679 0.79974317189770794 -1.6473974954337713 0.38892737296042851 -0.28643193976251419;0.90401033595659674 -3.7033163721487945 1.3878055479366469 -2.3503039135064565 -1.5579478043588522 4.8566760143305139 -1.0582816796579959 -0.47082049338354776 -0.24466610323680435 1.6058323281041222 -4.3027072924207816 1.8060097565153166 0.36848263043875429 2.9023406890256651 -1.1187038136069487 1.0322650802157134 -1.6589162685763119 0.92067319478171694 2.4956120134287252 0.98182516925864471];

% Layer 3
b3 = -5.0620811758127013;
LW3_2 = [0.092229664064935643 0.17411440735071196 -0.18695470506814699 -0.12702776674082189 0.92207216227277955 0.046749102957397398 -0.007922755965326593 8.8979173717277309 -0.085855860947825866 -0.066037264676550647 0.27042443580200048 0.35828899845124246 0.39365806300075584 0.17177625279777969 -0.04160702251987531 0.16042345646957873 0.10437285026592036 -0.33805676236917859 4.4096524081372301 2.4908950646017645];

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
