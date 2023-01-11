function [Y,Xf,Af] = ESPER_TA_3_Atl_4(X,~,~)
%ESPER_TA_3_ATL_4 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:00.
% 
% [Y] = ESPER_TA_3_Atl_4(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999999390765;-40;0;11.29;-1.99149636107278;-143.241523403244;-0.12];
x1_step1.gain = [1.0000000017546;1.00031699124092;0.0153846153846154;0.000310077519379845;0.0705602618858232;0.0604107502911417;0.00458552991639557;0.0145655815308426];
x1_step1.ymin = -1;

% Layer 1
b1 = [-0.63958144270253414199;3.0226463518903039507;0.75904293989353988081;1.3154081753970983115;1.6746513736060910649;4.1659054667671124861;-1.8154202262558036818;0.78328293906540125136;0.11857898497111968306;0.83369127274963183982;1.0562029013158806823;-1.5881056432198905615;-2.3046046550093279848;1.5040422807831381746;-1.2662271038964414416;-2.496671260756820665;-0.17420414763877672959;0.49677551916480516825;-1.9265899951508504984;-3.2566741127820582058;-2.5727342152045911128;0.71722115318386625571;-0.18444845081804137066;2.9630828627302867595;1.0016092631832487214;-3.0479698435611930485;2.9547490243795895282;-0.59822097207890490012;2.7895962536963301304;-0.67147148876992102107];
IW1_1 = [-2.3050296254862181478 0.26507061810132614976 0.49591148749099922544 0.53598765428189831805 -2.2367649025718012723 2.3334844658440143661 -0.61572693873544015464 0.9896161523393508519;-2.7965464624395992566 -3.2219706620965018473 -1.2674799939565746243 0.050909569886890770896 3.1099172580129241794 -1.0796850892088631912 -0.44735491957057327284 1.5447953966396834602;-0.10195697933104572175 -1.0116036656452804809 -0.34229865271022491546 0.3937223906564752407 1.8221034266107125621 -0.48703700344799266686 0.64280737709295909799 0.4143660471749717189;0.017930782025863841594 0.036498844944584156458 1.4133362036421046604 0.18659114249720765155 -2.4120235437038539672 0.6563997130912340161 0.72675579337920048761 -0.26796826538047813226;1.811747051403278741 -0.42300742211890751676 0.058931454461830311276 -3.1108181865261110843 -0.30082377666535581318 -1.2418134803320850512 6.068486179999526442 -2.5960318844245078296;-0.47343169346396701824 0.31410617980539462524 -1.1551575142611052449 0.86683227884423963427 -3.7841780594541676308 1.7384315337199696128 0.86662341124110386925 -2.923496845423569912;-1.1156510803737216264 -0.36829468961847749986 -0.010966134188826841961 0.62157999793939100464 1.257676985321505736 -0.074705832383240231009 -0.034468671992273526572 -0.58472370866273570833;0.91633808181985942909 -2.0657912697708984417 0.87517581564284385731 -1.7787537841588800802 0.33341276403883990387 -0.84467666744756353392 -2.0963607217986988118 -0.53803932858878533718;1.0795977561695195845 -0.43266950076443372319 0.067064750275923970535 0.36608256609799821879 -2.7376032176075009161 0.54496452955733043488 -1.2771713396162962528 1.228332191848187982;-0.71790898925433532796 -0.55604813111392326075 0.21549849193160869332 -0.045445840737051643177 -1.0689231329726680553 -0.73783287293859778355 -0.93308124304018791317 0.19198121233363266147;-1.4207126727533627708 -1.8002717076812480013 -2.0079455303842057567 0.16291507362446244755 1.759431439706168776 -0.76490332493866652896 -0.047660847877003399697 -1.4402266497022102598;0.17156310111466566148 -0.41101953789164508146 -1.8541003204829149542 0.51573562493161151377 2.6366486442735892481 0.32974912138242884474 0.10173269347416583064 -0.32657472816830768192;-0.56839467259868081683 2.6987597306018735743 -0.59548822955928903955 1.5278397185787107038 2.1291833387203862316 1.1640399140441264336 -1.8965295036797127359 -0.036956726748362089818;1.9734796393995575059 -1.9317881914469394733 -2.725889870320355346 0.50264811052106805889 1.1414536779983026626 0.31226522549113067129 -0.5281696621954283799 1.1090361323840940422;-0.79360812089923948864 0.1625896464921452178 -0.38910335139592433507 -0.39062937006841236265 1.5270968356248186826 -1.0709055773298805647 -1.056003749460382668 -0.16434230330122520658;-0.95412666669943357878 0.20347377114541934096 0.2697228478647428207 0.46857790141454092048 1.9942967284953905249 -0.35189120821184222132 0.030147121779240903389 -0.824813873092829275;1.5045264514223855024 -1.2454882756407550382 0.38465567524099786301 -1.108631493696416781 0.73777533324211841403 0.90619761209674376534 4.2392751243608133294 -1.1622193325076557446;1.3054128918998255671 -2.700794696024901409 0.52957367634756180585 0.28275557925760319566 -0.34039731738236889136 -0.44893704252843763536 -0.2056257483401403019 0.16919517451616905568;0.91034394109053895505 -1.9446550470818722012 -0.90973633533512565652 -0.02292063542533457482 4.4165807617414074571 -1.5263387460737940504 -1.5020880642668481553 1.0769364045226021975;0.39662798878610172615 -1.2993610771187840136 -0.036794066648577761214 -0.3929572816644577582 0.82393199127072425103 0.5508260895490884046 0.030482475655016903759 -0.76641151329439904938;0.15873735157954513286 0.17939120148105786035 -0.43858674054780388829 -0.71619529180214513975 0.99786031075538939028 -0.68427084582502928356 -0.43598418710402725429 0.039817511433723168368;0.07105435561241583986 0.54242992975758497209 -0.28923519898155031216 -1.6344766125066980678 2.3604920385913108483 -0.14723785867354782786 2.0056449330435932588 -0.011135022387120231332;0.019048596192579616421 0.049146829423833499062 1.2429267989842407438 0.84078218929126302683 -0.21401810555088018884 1.6234129436039930194 2.0679025258671797083 -0.69689313525653950876;1.3612680861935906318 -2.1315390294893008516 -0.097396636239524339573 -0.85986033209291146129 -7.9913338357293639547 -3.5216682541557546848 0.17000046642233307415 -2.3661729817501679918;0.065632223529400915418 0.043937114174562641455 0.12309566241224081939 0.33162766843500435243 2.0797478391337751447 1.1562770325834310636 0.036240153106357697599 0.40738603849313209659;-0.81707966636628015511 0.98952490993659314444 0.28394011949591435595 0.21165014060715048316 4.4051304228971810417 0.052694204169661446047 0.1085676301924894116 -0.16465860724394490444;-0.68430007436417472633 2.0068170268741365447 0.20915863611582954928 0.57508282417276634924 1.9650282078809375719 0.043884341610440476267 0.091931762963760405016 0.9415280985194858232;-0.031520092886623204964 -0.12309776418633949779 0.85858594512674069055 1.3510809347371131039 0.084073539855752915662 2.102106260192255327 3.0338240433232877002 -1.0849759712266138223;-0.89195001119433792169 1.9553429894453022708 0.94620037067206530512 -2.2180524543340180976 5.2045973621797809727 0.44175655065598751214 0.47615541426981594642 2.4296488130118918392;-0.28261441205828713441 -0.14042886862436601558 0.35069710076367738294 0.65583619695138950512 2.0798158743204631094 -0.17920529923195979705 0.16387629656209207196 0.039564835204526228873];

% Layer 2
b2 = [0.17890986259680455306;-0.76009360627677569067;-10.494058763302485104;-1.9489633260702183826;1.2474790335200300895;2.4950064342146078822;3.4148229472222948466;-1.8038940389352009497;2.3429140092606330903;-0.23560020229019307259];
LW2_1 = [-1.3633994984509811843 1.004887568936338349 1.6492373288501034345 -0.99703436111808929621 1.7969884282595360858 0.091858021403772305291 -1.7273565355337359151 -0.22433652544643645221 -0.88591339860004636453 0.082282950219483258492 -0.69935051015226235016 2.8663953570668758708 -0.72165603632754360408 0.23100035988348774363 -2.5176994950334372625 0.88116182202628468012 -0.5817141063738242579 -0.41552699137747928226 0.17240949831662519265 -0.17664324192484159859 -1.1550074914067065546 -3.6310699513667419325 0.8224175807947721184 -1.786398651989868247 1.7341312008497926644 0.45579327312490958146 1.4784185876401494397 -2.5548763509469814004 2.9955671483182535475 0.053233820896262551425;-0.85562087531282138286 0.97121283430502447498 3.115865087405471634 2.1416363459597369534 1.3190669714160756509 3.805346787193665925 -1.1952332647412213795 2.7377365943772065116 0.30649976137299200785 -0.46223403244372468501 2.5301900815175417669 -2.0441506544078529295 -0.82009378976885627388 1.4807295008987821916 -0.97420406407951865013 1.3407288944780715312 1.0020877182586716003 0.0075215232242340310936 0.51063610118363156332 2.3725514090087553853 -1.0977998029214797171 1.5739503690223015653 -2.7775906506739618251 -3.5141889596583331112 5.8129425719024805375 1.7948168961322443771 4.8831478109583761338 -0.31987980640535873533 0.70511372097625824029 -2.2052416610397291841;1.5332544673057602047 -0.37992312520995114022 -0.48021150087229486081 1.0358500679832323144 6.2970351596338813849 -1.9119813922814470342 0.28225976660707319255 -2.9407544731644676972 2.9415521102092494488 -0.12335930089176588842 0.58260473621436958958 1.613982327932407923 -0.77750796086931139861 -1.3607390256156368658 0.71089059701899581789 -0.68211565249126615562 -0.31842794558760750823 2.5834683406315663845 0.16492791584445626318 -2.1820386775095612641 -1.8244065928520256659 -1.2244299396051310236 -1.049984091834117006 18.233314383856402685 -1.5450942817567099308 0.68933575844002092925 -2.3033949846875754197 -1.1963517774339380573 -5.9364485014567716092 2.1452583546586825314;0.54318963788194629672 0.70212866326757061852 -0.26898163929684776319 0.67090294732067035088 -0.4438681420237961528 0.91400304566592083244 1.5183577204972527053 -1.3390998794528843341 -2.7629062135616804241 1.8150305099041792456 -0.099520351725746733496 0.39629348378217416071 -1.902223466420377429 0.55927611965252033155 -1.9335108088626582479 -2.5984878748478754673 0.78856188480523847772 -0.60342373493704348775 1.7707937489989553637 0.03412456909224081264 0.0051191017779521306996 -0.65318997215304752668 -1.6294386766227486607 -1.9349940245236088021 -0.59397674844190495591 4.2297199259392135318 -0.20557498475090676959 1.9181885743037740966 3.1614056095125198098 -2.1057696460148309647;-0.051036389313244674704 0.29581356001667274969 0.22723552321713863522 -0.35789878394678797946 2.9180237323035562724 0.56494892559068110582 -0.55115629622174022728 1.3843813415621610829 1.261973737861787237 -0.29987909162838033428 -0.2271540213918328277 -0.28941855527101784107 0.80138018515538678166 -1.3360428426463830665 0.56571663761564672246 -0.69845599449346817966 0.37947039974647350791 0.11286540522314288115 -0.23616338426648922355 0.337272150944850857 -0.14421825516538192002 0.13732923639170901975 -0.91066860555003859012 -0.94632203556472427319 -1.0242683428213952279 1.0362569747674492682 1.3547139847930058654 0.17154423792892553591 0.58016495075018459371 -1.0654363556585173534;1.1490645718642082951 -0.46137596072124220647 0.038290557106058997761 -1.4005540315978415311 1.9375784609469985931 -0.068525745915794328278 -0.43106431130797862039 -0.45654612438168495903 -0.56111474352922663389 0.27587248706892497641 -0.34599456612769036168 -0.74992681187338328108 0.22802768775328052797 0.49494192474184472985 0.068342953990886884386 0.71256896201692321302 0.063806561809637077354 0.32081332323053624034 -0.002114433260485168789 1.0085961998099703951 -0.45513289949571206217 0.011225600347486781178 0.37312846191066445511 0.27771694810276603693 -0.22605584799683953179 0.40293661567650151945 1.4931285666800417911 -0.053705063700711028196 -2.734551187730205779 0.4092329567775035426;-0.81316989631300518049 0.25086685740658754007 -0.13252370908838365104 2.3282207185987338072 -2.1320727369018084296 0.5457882817476056303 -2.1936356553450049311 -1.0347263740574166313 0.014604708005430289131 -0.22634292531367922185 -1.6953440518008684457 1.3880110530102278688 0.42588675820642679648 -0.012279464093343734135 -1.0161046947921665851 -7.1030146167478704911 1.308073250496931017 2.1887865959667553994 -1.1977673748396482356 2.595287260971533172 0.87737678227051008051 0.50609109225070292304 -1.4272627385814355971 -7.7727490744655804988 0.87640616375083491008 0.24986709943243151644 -0.29123051253668985616 0.014957994034655615681 -2.8136516397023303426 -0.093958952123980865601;-0.98357216867198116184 1.7559075856432992957 1.0532959717678649358 -0.66372596392544080945 -0.01502987877650463662 2.4150879666150824043 0.25899519861751696403 0.42888184062434847599 0.035275592736968966467 0.18841095887902706285 0.78836430655753575447 -1.0102616356798415076 -2.5272819952130110011 -2.7190469087448736474 -0.25062103855092915383 1.2037310037462223278 -5.1424950766459636142 -0.18336867190512554093 -2.358794121008727096 -0.078337799979195682498 -0.068693218958285784392 -0.38106073713177424089 -0.68343790870079490496 4.4848591967617261744 -0.25651696708572385708 -0.87944987634103966023 -0.54806967834129316497 0.3277542637058634134 -1.6720233973673441774 1.18647581858388107;1.1085500020306526814 -0.846608630025417086 -0.50982298751134502712 -1.7704360113325607884 -0.41330768295469311191 -0.11129923433096595986 0.099959368134053080968 0.39785067594367101718 -0.31882232057061737907 1.3356580574549978824 0.96238059711855838962 -0.29046018878113077921 0.72554854381930367069 0.8459423014634340765 -1.7778168400065135213 0.89459069003034696888 0.77035763020985359439 0.25431581144991433074 0.4629540184315180662 -0.7101461706408545016 -0.85532693115389557548 -0.50772930514741210217 0.41482575759543000027 -1.3606470776125607625 -0.29358512097228306459 2.1019824977058982896 -0.059186717717869735678 1.1250151823153455855 3.2434654722380034109 -2.3772173903993158994;-0.2926965559836062214 0.73000337303908346964 1.31755921324294345 1.7496571922866448467 0.9242633460330936801 0.064150100068807475173 0.21545244503548355297 -0.23204413207453258861 0.23858571634558628505 0.27410130106155400354 -1.406681643040642804 0.72344179411025499871 0.038967554528517071855 0.072007580837946505437 0.053434981861404558245 -0.91521183742097389491 -0.2024265384900413578 -0.58328762442382398046 0.90354955742203035829 -0.024362124637565263363 -2.3032691512305487613 0.3359191330440053358 -1.4080522043414522848 0.19354080030940779 -1.8649016706919365749 0.84601831798582183275 0.041924474086505840453 0.7297500209249877523 -0.031369618437312639381 -1.295050922078395006];

% Layer 3
b3 = 0.41856938814681682315;
LW3_2 = [0.73172993399872665599 -0.030027528644464156898 -0.0034360209817010669252 0.091788028199296961551 0.33460919822261975254 0.17409514651732205248 0.0073523169801053937789 -0.038998842703207588223 -0.13320097229330130073 -1.1820901217288914964];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.00111018595614765;
y1_step1.xoffset = 1025.5;

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