function [Y,Xf,Af] = ESPER_talk_12_Other_2(X,~,~)
%ESPER_TALK_12_OTHER_2 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 12-Jan-2021 17:03:49.
% 
% [Y] = ESPER_talk_12_Other_2(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timesteps
%   Each X{1,ts} = 6xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 1xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-0.999999984769129;-0.999540494077341;-75.751;0;29.715;-0.1386];
x1_step1.gain = [1.00000001035699;1.00022983623563;0.0140762793578401;0.000302617642608564;0.251952632905014;0.00826143816767911];
x1_step1.ymin = -1;

% Layer 1
b1 = [2.1555597747178874;3.4165778404962315;2.3833819145236053;-3.1717856787210987;-0.68173815560778561;-2.732384896632853;-1.4066220295040543;-0.38570761206287801;0.37070928394473734;0.34937021841687077;-1.2904455861293276;0.50273622509587257;-0.062335430257908303;-0.0036235947491016123;-1.6932506508976093;-0.84722570021423282;-0.61148194262369227;-1.3298567626473852;-0.86711925586193483;1.173467018712983];
IW1_1 = [-0.36614371294941861 0.26082554353366277 -0.2210716634492369 1.085312549198455 -0.62285525005686848 0.92934206433498434;-1.1123941286894019 1.03998805164216 2.3913328593064307 0.056679221991014106 0.78038737207223463 0.099686050156534298;-0.22999825473827398 -0.47038659787979081 4.7687418373917509 -0.42158227808526955 3.5906676062523832 -2.7557008705681385;-0.40368863567858942 0.22313857806002832 -1.7347403049518892 -2.2097562183604151 0.11716366699888486 -0.48780511408229815;-1.2010721254078318 0.48108772720232895 -0.67230519574051584 -0.86243534734384331 1.1190993261921491 0.61369780276426955;2.6801712647934122 -1.8439893308719748 -0.25924192550971464 0.3493532185880846 -1.1248966129298947 -0.99054000613510629;0.12203973917289554 -0.1072406726413776 -0.013445046793303654 0.071507520010120659 1.0142882347793185 -1.2319580510569872;-0.76682174489388122 0.35282392506056731 -0.88998797589779377 -0.36306569043342313 0.97155612186168705 -1.3764726080564953;0.0484661615097427 -1.8751684881099451 0.6175037924675072 -0.51812754886428036 -1.0821077942364961 -0.40802479397476071;0.022708856367357266 -0.17720994400275669 0.69721980293556796 -0.44346791032380289 -0.37736894061596737 0.046191747476978642;-0.34455668363524006 2.9962482735585807 0.13092436577539951 -1.4650842507332922 1.1452792479623786 -0.29608354491195216;-1.1893055919222491 0.20379685940297546 0.15610450539910758 0.031520775032502626 1.8325557700021491 1.0524362872187527;-0.33940986448628363 -1.1961866751312444 0.025715915904292147 -0.24841819447657276 -0.057239668511971729 -0.34997803506818559;-0.41679466128893972 0.30309715874565213 0.22896659877022771 0.33222062023783022 -1.5316319681940966 -0.10238562680897965;-0.051269762636563111 0.20957445799582064 0.29025300194280823 0.96474611899865814 1.6717212336863412 1.7095226620948019;0.20774886614373592 0.35674390790192351 2.0040709473160132 0.010119319244900593 -0.57060469072053643 -0.78901637542469072;0.095837747735445802 0.045800615405083893 1.3387970264650086 0.13202360881263348 -0.42141516840891285 -0.37649912776264699;-1.9754983724134969 -0.95735974516091005 -0.91775833951431163 -0.057554657540464939 -1.6455265204384562 1.5162139650456674;-0.47808305887492503 -0.11904817342050501 0.55528705180712656 0.38346272214791094 1.0447209618799629 -0.8715379354253433;0.45962088263054268 -0.34824806086406662 -1.0221360794283436 0.76713815879089031 2.1180704059251174 0.62526866990505492];

% Layer 2
b2 = [5.7454813487791645;-5.7796645814465348;2.1880713608549049;5.941192124090982;-0.11671639210704829;-0.89188690844592444;0.63918876962728977;4.5272385860104842;3.3348365499275441;2.44671910220905;0.19193452177360332;0.7439509728610012;-0.71123946660371684;1.88538874874616;0.58799202093260527;3.7414328700758115;-1.4517230572299866;-1.6756438888509617;-4.1791923629628194;6.1703610877285104];
LW2_1 = [6.668167201471868 -2.2908531359578048 -0.23678601895524795 -1.1266274267806806 1.5619170629458976 -1.1553281022243045 9.9257359705753387 -0.44506703503284334 -0.84488508479209246 5.7715777803364308 0.57819106048864444 2.9826326760821424 3.1586039045588179 2.0862212036203944 7.2329507989925457 -2.0780937137184576 6.844648088420711 -1.9268143479041902 -12.163052020255382 2.1887654239484973;1.2895920704948391 -0.32716772195949501 -0.18396887863378378 -0.72820891088134698 1.9113669993036826 -1.0651769915449685 -2.3054344162053386 -0.2925699030793279 -3.2497074155409593 -1.1750117826768347 -1.1629100705394686 -3.3251407548936482 3.4313155227962313 -2.3364104070447826 -2.7298831863161785 -0.38035297264335394 -3.6379138788633578 0.58909789888833586 4.6209363953814435 -2.4856479330664976;1.5236222423285355 0.34829062522210003 -0.24503024233197554 -0.29693409901975598 -0.39449193908645586 -0.09899124346287945 0.83698061612311658 -0.58586737996096316 0.96978909518272582 -2.4111308530835061 0.013535686425522428 0.9112893861625605 -0.43706605859963576 -2.7802258908783712 0.10082588042306463 2.4720834650069032 -0.81037645224411814 -0.41737452178429535 -0.039275832647470102 -1.7603587498626971;-1.0229283931620865 -0.5476288660688331 -0.028265585494592432 -1.0985370413343758 0.86917574781312601 -0.66198302905468742 -0.73276680082365031 -0.90401475092386985 -0.012748540561312683 -4.6445855167078474 0.00010507449544579607 -0.25483871286863419 0.66153609504737776 1.5554852906962267 3.1762888299371741 0.97533992600697517 -1.7528370847320491 -0.3941957181780899 0.35556786033921456 0.39892069338959724;0.22600916368431562 0.62293442106241526 -0.14644958449408998 0.0098337235911579168 -0.079615880697987668 -0.062006920272416305 -0.29695534967751636 -0.47581645531345712 -0.65402841469122952 -1.6047086391998548 -0.33101825352273506 0.065459657764780591 0.91988921642023624 -0.67501735723884815 0.15704097576551046 2.0209990071416426 -2.0347033868943405 -0.19606606754865252 0.34241586356779241 -0.19124344373432178;-1.5231244564140756 -0.39895749648107087 1.0322914027869887 0.042444342864017402 0.87565336917905956 -1.2406273006713231 -1.0047750796131634 -1.5639652000748019 1.904842317441068 -0.44755645866487531 0.24319568913230269 0.036758085502491586 -0.56608534683636402 0.18283609679226479 0.63560949706793679 4.3061388748237821 -6.0941717393229373 -1.3256508718719837 1.164224413295281 0.011575416761342628;-1.8954676065893563 -1.4314019975932359 -0.91734561194471254 -0.41624708309773545 -1.9402225134050815 1.6433691558196917 1.5819080966153798 0.0011951102107371353 -5.9694382324769011 0.9303269314966186 -1.7041887706679462 0.54309537728392365 1.2721714051128012 3.5337158173404455 -4.6172586239805931 -3.2739027815923309 2.726346165006702 -0.94944866277646034 0.20448752132591402 -1.3882884635533297;-0.32146303229712248 -0.31787001051460467 -0.25434030857703122 0.64490008653766151 0.96264310715162826 -0.88062718277236751 2.1484168568207656 0.52291276267418152 4.2584782479462735 0.032411320268626718 -2.1133772576388816 0.99167485500468966 -4.4508090492918004 1.6278564913720066 -0.27418542701162085 -2.1438658930216463 -0.25089204173839036 2.9733495235787593 -0.18487311087507144 -1.7229361087809818;-1.1780483143311544 -5.6904260763172605 0.78175223800116267 -1.1431278342197078 -0.63598511632807497 -0.68445798994447515 0.70153827175050265 -0.33978037582000975 -0.45208033919937285 3.5141927575639054 1.3671565270134667 -0.76584842816937526 -1.8305276802180359 0.13437533952017128 0.033568730425792191 -6.0295210198641893 6.3385214204474112 0.02095106807793172 2.2598837669750886 1.0176107579487688;0.33622148988517542 -2.3068114107847006 -0.85834763006050752 -1.2043288128433303 0.27617127915556511 -0.91431641509100892 4.2993087454792054 1.1765962312596867 0.3222228838189532 -0.42805373537149666 0.97948348119631845 2.3888662412756068 -1.782930305459723 0.4088815681560436 0.42212722287922583 0.57503016604052004 2.4606462136855503 2.103386979254088 -4.9905469507068059 0.49939015494952038;1.3675494754259552 -0.78205061804528564 0.37301335824549625 0.0056379698266827843 -1.2038984612635086 0.17883944693996889 3.4853904291380777 0.76271182048596486 -0.81200484488037128 2.7276369679291004 -0.13595222274898816 0.92221468404895812 1.8476033182891964 1.6552676404782638 0.72020497114776927 1.5884407505874987 -1.3963635145203093 0.22915367819135707 -3.0501990064870013 1.8252336757719612;0.61756132903537342 -3.3280943673815369 -0.79656482920619442 0.23549466072201025 -3.8914301036387293 0.77534951689855669 2.5291994729732128 0.32406549642789373 -3.2006655983372765 3.1438751363777726 2.6514328944534027 1.9918967551989588 3.0641923180743937 0.33877776474175519 0.41381173174905023 -2.2279406887965574 5.8529778913466828 -3.2463288146516898 -3.2407991435424228 2.2788113780923265;0.55372086405209509 0.42478110471252623 -0.03460830310515469 -0.18772932295284145 0.18499351053244995 -1.0136662617208303 2.1484860035965561 -0.9420826061574894 -0.34667293624148132 0.088793209753893057 0.049061348163883037 -0.082755113572142658 0.70212149862920881 0.94624618795385351 0.28267853749314115 0.36066305148119687 -3.3379409779924272 -0.56478545392310886 -0.24910517204178917 -0.77245808621725587;0.76720578416161733 -0.076902190174079621 -0.65857401533741833 -1.0058913221863841 0.65491568018954749 0.025611332820987161 1.7997931605673625 0.098542532529362015 -0.44910098411137195 -3.0177921094416127 0.45289801126196888 -0.263545610611158 5.0115055079766107 1.2443219616544501 1.1619975404603899 5.4662193086918718 -3.3739284977246573 -0.075321828268659283 -2.9125052053582059 0.52408847039341155;-1.4234644999116672 -1.7332698290566908 0.94250866316411863 -0.011116650516765265 0.052156779588373746 -0.955349493771115 -1.263941219657011 -0.47131667977883307 -0.71894266051561528 -1.5779396058268271 0.094679991629442947 0.3018377142908934 0.81878574367713697 0.8492836530229454 -0.66235882484653463 -1.3070745899211256 2.4467772816621642 -0.10187590265501126 0.71763090337426172 -0.61132876045650464;2.9259249017842741 4.8265287409273796 -0.23293853013141705 1.7066180352578688 5.9682378495121364 -0.54674030318002786 -0.51949717807825269 -2.9665602282053474 -2.640383996716507 4.7337886906241353 -3.5640206307796376 -6.5370789505440703 3.6071325693661538 1.450228049453423 2.4352008644641763 -2.6211146304812365 6.031076021555962 -1.2221384584647474 1.7933986435680747 0.96387065298677232;2.38911245753197 0.22217152689244757 -0.67171987563613134 1.0632468243028539 -0.18836772738902044 -0.020222752596826838 1.4692253469756891 -0.37045713378687406 0.54426094627653043 1.3501292732482315 -0.36080013958169244 -0.0015693814807513778 -2.0596829817626969 0.59667399812749677 0.36829260347408149 2.2771899679701382 -3.2269684958721001 0.20301358124651933 0.52142578866205302 0.54962542946230897;1.9279329064653894 -0.13251457705387665 -1.0152074413592662 0.85980343221182154 -0.46647214298073775 -0.42927668307984679 -0.076729389658364705 0.2050207043133177 0.45555328002570689 -0.74937562406962033 -0.50261834154696605 -0.18030587363111283 -0.66007024822192151 -0.77082272183688982 -1.0896783829014107 3.7044817566245345 -4.4744216686222877 -0.3951048824637915 0.72875871353536759 0.11799047768985206;0.22717609945766976 -0.26961201362026632 0.31411532096264111 -1.175546745709376 1.1173404886112128 -3.4764229998177152 -3.4077581421810348 0.55307781337455453 -0.86818432348980945 -4.7286023390999041 -0.088818335515448266 -1.1793247927333941 -2.6459532545386444 -6.0128266422412899 -1.7024640982026318 0.4107292625016073 8.5220677915026997 2.9485095740161533 2.342960998015863 1.454564635779952;0.021796816761311504 0.28415031013895653 0.15751481681506066 0.17360604692470147 -0.278678642754211 -0.63689256786316018 0.92763845167824499 0.93489208790637768 -1.5972489453292618 0.54692486701076226 -0.49802626635580249 0.62302550766025055 -0.22107421168659425 -3.3024452925311611 1.9992854033292424 0.7506656612130389 0.22524990849040083 3.53285850025655 1.4378269132242563 -1.066620910293435];

% Layer 3
b3 = -0.28026955470082393;
LW3_2 = [-0.13993760455242379 -0.090075590302382949 -0.24059785305007211 0.49196979789135786 0.54041189255686495 0.14307033144530121 -0.0402696626278327 0.085250825623926413 0.10795112511341384 0.037793953461426978 -0.19183883464865656 0.086535304075246011 -0.38006018042886297 -0.070054306021474042 -0.54832112450974557 0.024161947626476193 -0.22135961094151066 -0.31473516267381141 -0.051321038001202988 0.17737476237931235];

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
