function [Y,Xf,Af] = ESPER_DIC_4_Other_1(X,~,~)
%ESPER_DIC_4_OTHER_1 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:11.
% 
% [Y] = ESPER_DIC_4_Other_1(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-1;-0.999998793562952;-77.961;0;29.715;-2.08184650291172;-0.1386];
x1_step1.gain = [1.00000000274156;1.00000063368067;0.0138606862425759;0.0002963841138115;0.248478071810163;0.0588359962254402;0.00826143816767911];
x1_step1.ymin = -1;

% Layer 1
b1 = [-1.7066984902835728999;-0.87128077501957779116;-2.5380956438231692118;1.8285721912381858179;-0.33440063766918243227;6.1833011972819020841;16.311585203199381766;-32.537106585428148264;-2.5125410680187316714;-1.2381231301928661281;-18.174885208048618068;14.783870541250671238;1.7652260697192458938;3.8618965183339817315;-4.8768305502124880135;-7.4172110285590617806;-1.9787813064148358499;-0.56870516267930903975;-0.18445256348471902874;-6.8411308225496503255;-1.7649393517474292192;-3.0861220396670696609;-4.9265233342984942055;-6.0008832495812258756;3.4448886042797801643;-7.3082834271371464041;-10.204455834608978648;7.7076288021296797126;2.4015575851036548016;-3.4464610973420612972;17.215153682020297055;-3.5801139538480848401;-3.5481900043355643071;4.001396102305990965;-0.31193357256289966362;-13.916827977630775948;-3.0615199872990843666;-0.061632819382460797553;-9.66053159132231265;4.5950411988418906617];
IW1_1 = [0.0010956630594930807931 -0.11596033681139066018 1.0818854552700194471 0.10513019435540427171 1.4792432345548038874 -1.6628936913518881369 -0.17067713658443980917;0.086824114307294850135 -0.18848008483816922065 -0.46575435279804111577 -0.47824398338983759027 0.79936269839669160664 -1.0209874340276967608 -0.71472469499975599394;0.18260858374191490938 0.10844611662890027426 0.722232423395756995 -0.27072272440325712006 -0.091107413809829715379 -0.92278988101498016228 -1.5115184225273281537;-3.9842193496370428285 2.9647134472015008733 -5.2347438765261840032 5.5115133756110115115 -1.8170940793616672249 -3.3831391332780653691 -4.517722062972889141;0.070832396701320035626 0.14618181410519301422 -0.23281984308050771015 -0.32652015584957849725 -0.44440462431328026094 1.7562565594354053733 -0.91512543311133565105;3.6767980628925771569 -4.6394677020965762182 0.74203962603758233385 -4.295681195792229623 -22.200457044506027415 1.5957193922834065347 10.52302802680498317;-0.29256865217457450612 -1.2943315215491528747 -4.4107091814714198819 5.8676226490328851071 -27.904700011329442333 14.194785054709504024 13.22386883556398196;-2.4972114166183012607 -2.0421844054901567844 -10.153869760987253912 -26.046356818841015723 -2.8322196282332203765 4.7756533908703167057 -3.0862745464519227845;-0.86537057924172666734 0.70099342382678087038 7.467750300616472714 0.86284286559939593619 3.1632990499296282927 1.1448084076379432261 -0.96142547162237035341;-3.9508994129793300054 -0.49523744539344100035 -3.1215733190699417321 1.39502949878685234 -5.8657683900139163669 1.7103175761316031789 -0.96789619584021158172;-0.088532747049101398473 -0.17762040339582391146 -0.47840516436816582502 -17.395832645009598139 -7.2752405106782696009 -1.474409543343446316 1.5186604055501811139;-0.8432856211291681614 -3.7262890580717695777 -4.9987756166477943509 -16.41353932928250714 3.7475687258216430209 -12.791224093323606326 -2.7775928124155346843;1.434185106102941365 -0.83327119726319343052 -9.8515617968943338667 -1.4895703972741496379 3.1197872675419837663 -1.935494972120850754 -0.81516272278123302897;-2.718649355680794244 1.7163999953390498554 -1.2255286845608410751 -4.0091715160411371954 -9.1487565821451237724 2.8176586717616878275 1.4759048180366114433;-14.552524776852523658 -10.752661096356316506 -11.172075438595028274 -0.64498711953957288401 -5.1639645878177828919 -0.17926231036648762607 0.75584453841399734664;4.4139141120578972988 -0.094268280957578914903 -2.3001358980014625644 2.9574286480238032304 3.1096083678417092422 -4.4529572800862728954 -3.8803950401807951742;-0.17891549271059875581 0.068349677420616761214 0.0032316085582062741358 0.023052842482058914092 -2.7530698800960800021 -0.30331503139204274433 -2.8438710420368740195;0.045683800989093970024 0.11011235126768791071 -0.25443833313466818202 -0.42812628709767019952 -0.29973891601178809241 1.6029617117573011775 -0.98825146909958505592;-0.096089350012517560362 -0.021375781856610782256 0.99752902525399622657 0.24564599710913181951 0.35719406474691517417 0.23180115031065551268 0.72576854340492069628;-0.36875305640646233174 -1.1540735028645903171 9.7057726608600631835 -2.8129313884109077293 -30.156358034782357436 -20.631500010167307835 2.3766060078300488279;0.15949563177632208877 -0.024443271667649837248 -0.84818942999046553766 -0.51208089167031511924 1.2210181600497183041 -1.1243058456237633447 0.22711689318832409268;-1.6632864193650851448 2.9134728362827804204 -1.6103657070894272785 0.41374906022619162282 -4.0561683343571699467 -0.64309878953304289162 0.74739819565302145765;-0.41850598571260477776 -0.83688809416394371876 -2.9996932720844209364 0.069765190847709018795 -20.72363284013302831 -4.7351110086545658007 -3.8860345299968730437;-1.5571824465848782548 2.5874369209596879138 -0.71785427726988915165 -0.64648178626943841873 10.048192465468055445 -4.8142625890333041738 -0.64089760982931343047;-0.87869732131429867383 -0.56598149997615521123 5.1469299743275955805 -1.3939429006983585779 -1.666156217078093027 2.3391262893380222465 3.8612873864055958784;-0.30833611930864840422 -0.68760282342850453929 0.55317116487504780498 -0.54277350595243845088 3.025015727157287504 -5.7371805170979337518 -2.0299554279240856758;-3.1596723516053968339 0.028428055091183325465 -5.0786593943336288959 4.0848207451145652769 -18.822708969739203155 -22.584079406543832391 -14.419752153360162694;-2.4714968201376930068 -1.3877028611458130936 -7.8683813054526741126 5.7663523096344091456 -4.2852440241046094727 5.0611032697746995979 5.9367944134056500616;-0.017600197342162102188 0.011942492043374063754 -0.15639527712334114207 0.27600543867563093192 -0.15740870364606199261 0.024347059059255656693 0.75570478889454839599;0.8228139506592609953 0.35598602785384847458 1.5987830350673177815 -0.25062796941953269236 4.5570471717385281352 -3.1468785923035813568 -0.92105863069181248548;-8.3124573400920329647 5.7734174813147758343 -12.153052975197669028 -8.1163170958015466994 -9.0076260183496525968 12.761984850595510466 2.1445684066056278638;-11.721220704370809429 -0.77319397138178103113 -6.3916502666036434377 -1.1803736016949293308 -1.2276252639257088095 -4.4814187694182736621 9.3854155472453246745;-0.5361698602863963492 -3.4454836434606610496 1.8415231048964955107 -2.6971400229502346946 6.5895453658759910454 -3.5447828882861265321 3.5826440801466041819;2.149987295010563404 0.10157391875739876785 1.4136225029039892753 -0.90466005129941062801 3.4142366617683523877 5.369142876517116747 1.5116608464638552345;14.652312273886570537 -15.508818550880780762 -2.0663538085869568306 6.0036015915753706551 -5.1325604423096296003 -1.4785019298399906251 -8.4349703960252622181;-1.3992178859848876105 -7.8129699545898620983 2.3424718983319636756 -14.492169113560640525 -17.573128593949995491 -8.281353799197654908 4.8432549507611781792;1.4461198582418131053 0.56411064572324987054 -9.7247435808770124055 3.3262300056739748655 0.4935051435863925251 2.8667532350006155184 -4.9748345586538560426;2.4914772987402775684 -0.22640465225030523277 -0.73816890252029310204 1.9721178864781203632 -2.6195857160258024798 0.45779624448213868115 -1.5069701178808363462;2.1558547895804567318 0.1002482076589675114 -1.3770801579863192199 2.0044771261159279874 9.5741696830589138045 -8.431003634531126778 -9.4352124117547955251;1.1996111095306076777 -4.2010793424439505372 -0.63903331932523954606 0.16354991774244240221 0.42993705447889868365 0.4666697452863243778 0.31647609505890483828];

% Layer 2
b2 = -8.9606169224332354872;
LW2_1 = [-1.0466166901834605962 0.72356301424739177808 1.2481786002733825658 -0.017741270406654949254 -4.8736808854067223962 0.013248415826555099903 -0.0051715320180592723454 -0.02589394852332144864 0.06957244552241229063 -0.035496057654340371079 -0.11250498616795875928 -6.96333608462090492 0.064019298652715495024 -0.016594723500379087583 0.0069728723766059414954 0.031654941961070208467 -0.19879767537588932669 5.2144198462913511349 1.0409129557002747557 -0.014811074798192306115 -0.47914559176062254675 0.061596343966384706214 -0.02392004305169229128 0.024301491317887592702 0.063293714653985602703 0.26812156995856717412 0.0086245236241414124229 -0.014611531437868465963 18.200259617663046896 0.1096005702776241264 -0.011637001117463569369 0.0076382343010660347718 0.016997763968434567716 0.041743543706606579813 -0.0060620645654208913516 0.006316764891199659962 -0.028375309806645123561 -0.057451027694888648356 0.027078088574641945296 0.051173189098188097412];

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
