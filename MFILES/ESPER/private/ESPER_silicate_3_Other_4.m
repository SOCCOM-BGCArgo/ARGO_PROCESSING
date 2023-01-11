function [Y,Xf,Af] = ESPER_silicate_3_Other_4(X,~,~)
%ESPER_SILICATE_3_OTHER_4 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:39.
% 
% [Y] = ESPER_silicate_3_Other_4(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999998793562952;-77.961;0;29.715;-2.18999082656054;-138.904784718693;-0.9];
x1_step1.gain = [1.00003298043486;1.00000060321889;0.0138606862425759;0.000266418009857466;0.193068829037552;0.0582073296585987;0.0042772329208737;0.0415843642790311];
x1_step1.ymin = -1;

% Layer 1
b1 = [-0.24282358309676463359;-0.2597182394862295185;-3.9475949065005835337;-0.83844038135956278079;-2.0308559442051392629;-2.908946029460857563;-2.6530621307452011415;-1.0951115929981316555;-2.4517636933149775658;2.8754880126920827621;1.1121570033454526438;-0.59781996203443121818;0.33613888451222856402;1.1289048302917892297;-3.3348689405202383007;0.91014132599486818886;-1.5278627779799542896;0.14503010509098160363;-3.8320670454950920281;1.0666945095644637664;-1.4879177609199483179;-1.9695222706469466267;0.44875730826337700829;0.13977483458353084922;-2.2334032400066092805;-0.08371268628900824238;-2.5065891079598330293;-2.390347335940857576;-2.3234570926019562087;-7.2839401417270339678];
IW1_1 = [0.75004061164714086463 -0.73339780777891905306 1.9899783344566621679 -0.95552369210151066881 1.9082886560632246198 1.9205341422641426519 0.29989934604115026096 0.53970170852693144514;-0.44363196118169107018 0.79799865781449574875 0.63224944652632875108 -0.51384713008538496215 0.370822622644596811 -2.8691458828415292714 -0.40616931514888388044 -0.10688076389131240251;1.5978069099627756344 -1.2168960171660041514 -3.4959046764126493478 0.18029991276139128997 -3.6876938396885217131 2.0549619888802870271 2.4819697024847342526 -0.51965850477274266872;0.034340697387434941656 0.32609283152671830752 0.72163567087585123172 0.56886281933382931708 -0.043565928203282311415 1.3880187874526692671 0.42825143723154746045 -0.94212093620624282053;0.59962609806850541982 0.34209420822502889514 -2.0913735485846798312 -0.58943612805372291863 -2.8199940800426586129 -0.75288380609691418766 -0.11518109982539115332 -0.62037732880280971326;-0.079170251765137367173 -0.066395354052558877567 -0.86476117143709752 -0.33743099318048269675 0.45865430280363289617 -1.3014795474569267331 -0.55629027476060299851 -0.086491451384986067774;0.35349346635405709538 0.0059797825942146248801 0.12594238132738172498 0.4418723853402106716 3.5628009775883495713 -3.1298741250707600159 1.7501080148874144093 -0.36010715979041574286;1.3429211688934112878 -1.4752079388073950916 4.9152042449262456358 -2.9598598512575042818 -2.0318426003915912581 0.27499691317220936737 -1.8631228938541111173 -0.31682818014347458035;-0.19951972776709264723 -0.46204222592479426979 -1.4443147765148047768 0.36720222128322843602 -3.5276296473554835664 -0.83392965785453920802 0.097471118539218720711 -1.7596374259657396877;0.012642896136748373731 0.45338213038885794015 1.0372195024093973093 0.032387760730983990165 2.8228626282205664388 0.99395711007877585175 0.837744732880807641 1.5382798824358672185;-0.18211710854939572313 -0.053570808649756364306 1.3903080582961184763 0.76304465511798669208 3.7447312492152571295 -1.1999345350509127428 0.0056644964167833645324 -0.46685344758555041844;-0.11320791426075096919 -0.5893421001460839026 -1.3926225570204293636 -1.6891040396574295013 1.0375789875246028693 -1.368537964672082774 -0.73584102915477522799 -0.1278855834788407142;0.87774733252697079156 0.098053597657864879911 0.51254428877386049024 -0.1540610949867253221 2.6569350903856903123 -0.097611632910006379538 -0.64861614272162915196 -0.36928298157778033373;-0.73924768533652496227 -0.60910943383813198615 1.5013703591984264296 0.34604416923748693158 1.8015611773972226572 -0.3727473253659457364 -0.4935134110900379345 -0.012309292516304984333;-0.17732354640509098709 0.087532244536562053838 0.16261010179644297402 0.45800007899585465365 3.590443988682670895 -3.4179563353237010226 1.0478450447614691132 -0.91445518115358193079;-0.099257724244897135857 -0.063774837094870417964 0.072199229881774823481 -0.48355603461546808886 -6.1867293348590477464 -0.088848103503375172596 -1.8791624166547953312 0.26277573064380715939;-0.026174804765367067766 1.075828414343030337 -0.43945614998728593781 -0.42353042721512429436 2.4544123931399783878 -1.3671369415197360819 2.5041579416593191532 -0.40925617307156952895;-0.39529359177825734717 0.062948506162698966704 1.2159342108562869456 0.26782289073256182466 -0.050848987133611397993 -0.044377456775499680985 0.13928545391838131895 0.20274554104473202498;-7.0021969704640207866 5.1562102657973811759 0.70936808606666901245 0.50781414209072894828 -7.8573261372028859739 1.27638456866512362 -3.1291513510619499705 3.0044593352268282871;2.3677028701142188716 -0.66268378250806536656 -1.9132266026765170075 0.84025456526977959637 12.412187685914679136 1.0712114820156675155 2.5551146864675655657 -4.5881931380100251516;-0.47638713902502549669 -0.44893403922997032751 0.55000312934164863332 0.44708258271961753838 3.2114050708943278067 -2.912909425609509384 0.28506073599099113203 -0.43823959915662652476;-0.46054595516002938549 0.67771376285801743222 -0.30779381936068422876 0.37173579187404603985 2.9527456294281666516 -3.3712353716514518176 1.3933370525642074789 -1.0458524936640183345;0.15107647944339266943 0.14135836682036653777 0.40660047968269774632 -0.7689413029711208436 -0.88675878508767502773 0.51914938159295476172 0.55326772276958025021 1.3380524838093170725;0.18199881098218628761 0.10047890265099981844 -0.26482175129972213012 -1.5858206528084182185 -4.0749182279987712363 3.4269956687922160654 -0.18287275592519652845 0.49523082432800141772;-0.6415276303974076777 0.65682532855366115854 1.7033727763637955821 -1.0231241723643953456 -0.10640006102554667999 -0.37473219235282118067 -0.31710077622152343935 -0.45368298948690238825;0.035621884631382612352 -0.1461205409940147526 -1.214847165008409613 -0.7220813006929241995 2.0897855108382099587 0.87157006103161460686 1.0838494325843599508 -0.14365217998116991405;0.23849810158347162226 0.30026106529182500227 -1.281288345900599035 -2.7237998383855277496 -0.64258208374188641177 2.3261716199714355113 0.091931411043610883338 0.22552660092508119982;-2.6969880669031462972 0.10753699101084272793 -1.6091582304842926288 0.8807616103224483739 -2.2838805626372700175 -1.5236631630780788083 -0.27145726007668563584 0.48754482140579924465;-0.091145918716635299095 0.13317815273046812452 0.57706785972011698238 0.023919499969964379965 2.4246951591695120598 -2.1520842507632571561 0.3223611995590870416 -0.93712761975595837427;-3.1032652063813750409 -2.1490440123780953918 -2.8007971489153149491 -0.48851885925498383934 -0.1195960821398139029 -4.0117406598888001312 1.6771728364467513028 0.89490234659138079287];

% Layer 2
b2 = [8.795257568968496642;10.285146737152665253;1.6092688625890709808;-6.5747851552094891403;5.7260437270168891288;0.59794791024628213272;1.5150661049156901594;-7.0213689313998166597;-2.2787099090711619986;3.0911633421155917389];
LW2_1 = [-0.4649732942173379624 -0.14430280947736845443 1.289548048871685948 -1.7594261578299521531 -1.4521003399874987672 6.2737093110592105916 -0.94796299800849603834 0.44025620281685234891 2.5347878816235236421 0.69576790447806013962 0.64546096203967173732 -0.76626177191123456822 0.27779868033624249835 -2.2564429687156484228 -3.1246771650771103523 0.1994328031463212525 -0.9145178550792431027 1.0774095650339223784 0.17631530825239088989 -0.13010780068954727606 0.2586706765619045667 1.3505642705614575583 1.3168537183752198061 1.1108233298078626916 -0.61442622303114957649 0.84311490270902234023 -0.44323620643045347212 0.09316280405024562894 6.7520824137810224741 0.3853663387845549515;-0.79322400981869711778 -1.7164603095689499046 -2.9402240832524970671 -1.3140475527450681259 1.2600436786806752387 3.4975694793305067698 -0.63124832587696988462 -0.2203946047750847681 0.14298501874407806844 -1.0726448203693861316 0.56334475913315884554 1.3078439748050909586 0.38669439573973307667 -0.10962455065185990055 18.948376611302283123 1.0037718922122103038 1.4027104690213494376 -0.11302784898553280835 -0.028268554389187943582 -0.41342516770100024015 -0.2838917680864205817 -4.649923205712983787 1.6168468091094070349 -2.3623237158522836232 2.4262808202617001108 2.4466450379466317955 -1.7655863574115995274 -0.9094997095889902905 -4.9999869618033407903 -0.57363361642288979159;0.19855902403197908668 0.51611232205654888894 0.96059603802544146767 -1.4570989708118937678 -0.31665924251798477362 -0.041013663417358334207 1.4052065556769970378 0.46133434501976444686 -2.9457527214350736067 -4.0165815946425347249 1.002459037528563579 -1.785078481343744583 -1.1150483117186358317 0.079488255919804851057 3.5065183014877527867 0.55117129380657015059 -0.659956444281301291 -2.0819279880777394887 0.18546905734572960434 -0.19248906400820312013 1.1676686684739849831 1.9741954928887084364 1.4459626114364956262 2.4217082321153351465 -0.74487170366951249711 0.44960822913397779921 -0.30573921310533802576 -0.46598775321719948517 -3.4748696403368928287 0.85578864150363642693;1.4756444750703956448 0.35870260253017849372 -2.1024617936553782194 -1.7811047517419589248 1.2648761768375340608 -9.5475665066950714532 0.36725760036876914372 0.0087413535619816039313 -2.6739317343001127725 0.097043541417196177434 -2.4723103447839611668 -2.6335368529949438887 1.5864520939337201089 2.2796803378367722637 1.2616220403102762493 0.9757465895839466663 0.2129236389606873614 2.2890767275571350936 0.23741238683483464889 -0.46154966277088282922 0.17766453940005894063 0.71951543859825328742 -4.1628627245143006164 -2.1311529813317102722 -2.8221479397300801573 -0.57268001004145663302 0.37038230265913596462 2.0977420379039966036 -1.8028317517595242592 0.21478045572329956814;-0.43794922711999711984 -0.5262886790235954404 -0.79338417861351362426 -0.57006438642478496348 -1.0320136814714806484 8.6524638003250071705 0.71236424600898617943 0.073436672147634177166 0.66202087177906854265 0.44965044676583598715 0.53103843829878760729 0.18165216316055113333 -0.063877945221993193137 -1.8958101927506183504 0.6257493564158390642 0.9760178706529332171 0.018543305752631761979 1.3839690325012570238 0.10339988559024297643 0.032322525270583858659 0.92046161316131303032 -0.25947378432543394267 0.27541692185455990405 0.4770904105608999024 -0.24043819967737944943 0.40657413895821131744 -0.23371930152025030547 -0.20104253161859045318 -0.42407432699586106173 -0.016282076434404095522;0.13696592351563538492 -1.1726444748973667842 -0.89192535795403626864 -4.976564198444130227 -0.88013768948933568392 -0.94261072724858363525 -0.75434153090556455012 0.409840634855228092 3.8687744153773704348 2.7678802805825393207 -1.71221360852684934 -1.7473656342607215564 0.84649848053660281533 -0.0018041946429510182306 -4.8237448343214266444 -1.2458243651139939701 0.04512325861978325886 -0.33467833800860835591 -0.19458491460329732647 -0.23036625675006539016 0.7224476622781580426 0.39445763792992832864 -0.6572825706171127047 -0.91848856932301592426 2.6219812409417606602 1.8954900026471990238 -1.5835433470299242753 0.66677167876889698039 7.1082467543522520614 0.28266296726212397639;1.6826774734724636584 0.85110994451698629426 -1.9529310555323913245 3.8268463699965553637 0.54942283710924100237 -5.0442669068613099626 0.27192968742563711215 -0.062942008028783663631 -4.3971943459173878921 -2.5437275741818274888 -0.55436586193040016113 -1.2897302534894958814 -4.2585542699621719009 4.2754879749486685014 -8.4753664229231748806 1.3218365703047372239 0.93835820042474760161 -7.3429727315201471782 -0.13944490799670694603 0.83699692195649166848 -2.2074265267417372272 -0.73190877903808093663 -0.49836840146939620499 0.7710000516247005331 -0.75322612386883669444 -4.0849995198234605454 1.496554031991307987 0.88750501501632605983 16.768149877045775753 1.2097195550697601441;0.66713952473721616077 1.2703860686203640729 1.1859191550381056324 -3.7853043486172448517 -0.038218983993463129245 -1.6105630930348291407 -0.63957513004009991509 -0.1143835421446211531 2.7106445445098659164 1.3720796751974291272 -0.32530566113321129906 -0.22618610385892459758 0.58859149399101573508 1.3482128535768407307 6.5313005439508300398 0.4677420470994394508 -0.35765996309068215364 -0.59512091975471992367 -0.12578097673900337505 -0.046223745576714353467 -1.4220563573381534983 0.38537999227645225453 3.6718388186042956178 -1.6404543528976174294 1.0383113111008799301 0.32069418171142927143 0.9493194144005577817 0.11736750179717231335 -4.0694630074850568491 0.1403912775766939236;-0.397005657659811384 0.28725706919480309409 -0.34180669144313324814 1.6875248846721373486 0.030087728108753500345 -2.0847084146903078761 -1.9317481560836113896 -0.41650241242048735568 2.3414042993524470937 1.8352144694642804623 -1.6426686767247482557 0.26241142851195414787 -0.33857811705472262709 0.44070162194744988415 -3.5543353202625951326 -0.98979227141529746437 -0.025242453225328871647 1.5486545430574987225 0.18618890776299906742 0.38371250209687018939 2.8231183525679424129 -0.25481643263884423689 -0.79185043370823737785 2.0576498360904507123 -1.7707028576637020567 -0.17072285050913463267 -0.051264957719648930612 -1.0092172890792896478 1.5391757449607583741 0.19896443644083577906;-3.6188432442226106112 -2.9857922672139478593 -0.35469028898224058688 1.6679774531059210485 0.14042191760067967987 -0.12154610628816515938 0.14033770239415727099 0.23084084526669079729 1.5437141077804090727 2.577375124314062127 2.1743257699539815242 -1.7073261553940055624 -0.2875016012086972772 -0.1489917186179973152 -8.2522124053238652408 -4.541214330984207237 0.57801488621424013115 3.42275214979711917 -0.029361402098681425865 -0.19669116449720380757 -0.73908856091271091238 -0.77869630232389130953 -0.79887021575415673169 2.0822482137907578625 -2.5354966755427499336 0.79365954875360167353 -1.268033552161834443 -0.20293342655807999897 4.4486981096884132825 -0.0025674795377227115252];

% Layer 3
b3 = 6.5660133712271013806;
LW3_2 = [-0.15922460400792495805 -0.085967486738582898909 0.098895869696420987682 -0.066360221130070354278 0.88406582897908925212 0.08490711323279344358 -0.063459653523348091841 7.1920968475896431826 -0.12393223183961651901 0.71874464353027833763];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.00819772922900357;
y1_step1.xoffset = -2.02;

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
