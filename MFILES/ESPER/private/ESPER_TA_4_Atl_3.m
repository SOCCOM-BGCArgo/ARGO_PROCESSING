function [Y,Xf,Af] = ESPER_TA_4_Atl_3(X,~,~)
%ESPER_TA_4_ATL_3 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:01.
% 
% [Y] = ESPER_TA_4_Atl_3(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999999390765;-40;0;11.29;-1.99149636107278;-0.12];
x1_step1.gain = [1.0000000017546;1.00031699124092;0.0153846153846154;0.000310077519379845;0.0705602618858232;0.0604107502911417;0.0145655815308426];
x1_step1.ymin = -1;

% Layer 1
b1 = [1.4157089009129435464;3.0816675028776185918;1.2913702484895754896;3.7343304593008292969;4.673344840272031675;1.5556694779942823459;1.5438164884556369039;1.3194770954886612113;-2.0608586904813730456;-3.0751382642138174184;0.081662512618215082894;-0.42309831640235945338;1.9317097843913939581;-1.5959557559817014738;1.1606175323774334807;0.75501246971960989107;2.6698074619460023804;-1.4180672796579569805;-1.2284374207593489192;1.9708489104910018597;1.8486809127821397514;1.7015428874583249375;2.4795092433463850057;-4.0368319172487447943;0.93003240867287739757];
IW1_1 = [-0.088806852006714565029 0.0075311132726756171984 -1.4786989484570185471 -0.12170052777758427753 0.49850148555655293681 -0.28890337757255252216 1.1714965331707698493;-2.3563601810642338386 1.8314143208848920885 -0.315238474988045414 -0.087130005972621826782 1.8777360478001972055 -1.2390581126528690259 1.2531418988478273402;-1.487400318077032324 1.5244672873634470545 -0.42439631744876998853 -0.1113442917777090202 -0.78383186215672373542 1.5804671526825950778 -0.2747420668674748212;0.32011197848788214282 -2.3139632031080639152 0.040292235851946071878 0.2624485317379682825 -6.0436905768001532024 0.1373696381143218459 -0.99397911393056348484;-1.8881692503577582887 2.5874407258394684916 1.5697054610584171286 0.35233477752779540637 -2.2759985254081724904 -0.62491196246588176511 0.73715181805746654575;-0.12835928464221557643 -0.92005769347805843594 2.8752208617126249202 -0.25243745671303075007 -1.4732907442160945166 0.90636826931347003189 0.35873969785410364519;-0.93253742337753486868 -0.10253730095689381208 3.2210311189382778707 0.93973383219319617243 -0.74002369467014805782 -0.88934581941983414932 1.2159986624659857579;-0.4838831574381273759 0.93274545499998495757 1.6711524018438583816 -0.14794407410615847964 2.0354475784877927858 1.399112775509923301 0.35376832596988255775;1.4320325774345075143 -0.076577251201022711324 -1.8131427852003176859 0.1765988016174070141 3.3199502393781292753 0.28850980030401984777 -1.0002502920065279302;0.2214844770608856428 -0.78856066210962949548 -0.3490202836771422179 -0.22347706498897196292 2.4455937970600531983 -0.88207345926317781704 -0.69369677630623050035;-0.37083592482089594577 0.1219327362387420538 -0.092947487167511833039 -0.12787739164947467541 -0.97277120290353824572 -1.3128380598441229399 -0.90475865043308378866;-0.49141848537881421821 0.4260954639273474287 1.1685168879038794199 0.37414775378502734826 3.6229716900608717189 1.3302555048919475489 1.8457255559903036346;0.48278076352062371246 -1.8362366810471075773 -0.38828987333650510916 1.3419754705779900483 1.2656493987601540319 1.0618771500972783439 -2.1452743781292031855;0.16066864034472944089 -1.1453481326368124549 -2.5251208763210652997 -0.73508831838431709649 3.8559780160404684146 0.50476215015158021338 2.1760339243024011679;1.2620956145647972413 -0.92884477991722891321 -0.19054516655514011636 -0.020723719522557935979 2.1783208744998252548 1.7567352033991265881 1.9106731279432382919;1.2241070819819126481 -1.371280699763777422 1.8798639253736912913 0.084640578210163908146 -0.67074227924760898567 0.063383998842255179218 1.4071909712806291992;4.5786658180226291748 -0.51141350545435515507 1.6071159636185279318 1.6336772833712664976 2.1134270784874180649 -1.5455432812167089196 -1.6211573163393204755;0.70693620423049430279 0.81043725760262685931 0.4157417132397288051 1.1866428829645814247 1.8817790021709899229 0.99744769406755207086 -0.62131039667302856433;-0.15066348046662308979 0.32398592741251952543 3.6208076931637984686 -1.0098506131629232385 -4.1430736706948625425 0.66954014773859438225 -0.068782719088402249907;1.1255646035707209673 -0.81250241841408521459 0.51799448575151918472 -0.78988698668679557802 0.67575929693572822643 -1.0902986993737284749 -0.98637232759705895102;1.7964423173916050924 -0.48728347344328315494 -0.034333785211741962851 0.94898721218023540658 1.1239251654981423645 0.24433140918901008476 2.4644131376410793877;-0.32113478785403604121 0.77210329831487112262 -0.26341383457966771608 -0.27978045705650766983 0.68151864101718873101 -2.1847044475447954248 -1.1023099853159421535;0.30728861882807745332 -1.4709339561695176091 -0.33729318391787999065 -1.4253581348333641543 1.05823488508913921 -0.37457331120531522917 -0.93407425775522057165;0.4404601067521239921 -0.38625349424364163742 0.84122997108065211691 0.086116177146354541683 -1.3434083278740931444 -1.097841363982993812 -2.2479139471670053219;2.8751279798226270223 0.18064597122777048543 0.76794399733458673651 1.5060924682392975438 2.8155255114661521532 -0.94976253842728874943 -2.9376412448323079296];

% Layer 2
b2 = [2.1057567092101625583;-3.0925916144684153508;-1.6297605998675133154;-0.72534785985633820626;-1.1274791171945564017;-0.76160215499779326809;-0.2702004347325779321;0.45733335279752412861;-0.14850312011305727422;0.093040725699261320192;-0.98969223379081316949;-1.5428104363392600717;1.4248417722759416648;-1.5649579918196045991;-1.4005085603130769112];
LW2_1 = [-0.67187604559253921721 1.0831317218949729764 0.84770256281873535009 -0.33795446336198109449 0.13042899049567432956 0.28061380225362048124 1.497706407128386763 -0.43186172494434460489 0.17582000354942833864 -0.6093797139894149284 1.0588738869481126859 -0.27094410345034614096 -0.2627710670759572098 -0.99912273956642971839 0.72849017244413472039 0.013138710116135817602 1.8359382334213558696 0.70338936629608350692 0.65196584739550311749 0.90340267311259381255 -1.2698471951469547658 -0.23182143456865619835 0.032433871508066897116 0.97158303659364353422 1.2002473795095973674;0.59583662920370927818 2.3747295145239659675 1.2514011207581072771 -0.35542172275939398762 1.1394493845096338447 1.7814571322816248333 -3.4076122383977187624 -1.5987089257144655452 1.6380593141009802682 -0.33853254931489273893 1.0132308171711206235 0.25173319439304731038 0.98312873004463507964 -1.8653557261395099509 -1.6258512870342975987 1.455070310034946246 -0.036594370401252955538 -0.51189179137283424836 -0.30646517585946930762 -1.4437750187203914898 0.062888749417309278145 -1.8027851162197288115 -1.804504499707123566 -1.8628935612654646814 -0.19621876097166515707;0.71786636005353166023 2.1324885187687190324 0.40642730498900975222 -0.56543341520699885727 0.3499584009559238984 0.096524797938003281672 -0.93217083214554241621 -0.92850967139598172828 1.7429010556353643491 -1.1892905305348335787 -1.3292885918068177986 -2.0388685640104138308 0.34935594410522358144 -0.7920618342385892241 -2.9163403052930463311 0.88084675435642645347 0.5708049224408643374 -0.47395039941884187229 -2.3033489746588133862 0.17568742517308114093 -0.61926392692377218996 -0.5910631395215020234 0.18950046192442768223 -0.86704431010853821871 -2.3171014896426367358;1.264330740762521943 -1.5254875726007923742 -1.0959751497327439029 0.51594432649201171781 -1.1211290788651511185 -0.17008611705739945519 1.4486273638401025998 -0.084247428760480025378 1.6119996652328212328 -1.9387453661990750042 1.5801964513146984359 -0.075378806199436432234 -0.59117661339282467914 -0.33021196404847341643 -0.38532759961174278818 0.34501921369638810511 -0.71434313754228973004 0.28785900911604028751 0.82313175469947985619 -0.1828672453367756745 -0.86305763837374760428 1.4277865795408881677 0.19272189484455118369 2.0176577267890118961 0.13781733025179582208;0.55319469453412561943 -0.19295616095881401586 0.04405591676260552092 -0.69537434154736244007 -0.00069276308503487225854 0.062854164802322490258 0.39569225516663109055 0.35093616918843539842 2.8748171144556597945 1.2808834786297960129 -0.52128771653565253441 -0.17183451949385813995 1.7471592545966878607 0.33745833055779694831 -1.2003348213734996719 0.49081764622436896417 -0.46228956556237471576 0.68967249518986517565 0.84352261507644732408 -0.78652155850612703869 -0.66058659617204074088 -0.47168650719045146502 -0.1729993669572564563 0.41550999655577181002 0.50789836307729463449;0.96667502632106949711 -0.20639185172101004451 0.20061005337238313473 -0.14862006495257273908 -0.46616887611769780575 -0.13829169617407427451 0.55544143255859890207 0.45947092706687858721 0.82156892919133717168 1.5260015219380504981 0.065651041552008515545 0.9660233483394138343 -0.03030038055466077393 0.18761807392156715224 0.55791905184557022412 0.59720542127066367222 0.019026723023857466943 0.15362883004468216863 -0.97242188252902883594 -0.88011333724704499115 -0.92280166242642380769 -0.31511384879646697765 -0.8503789534998922317 -1.5901461530889728913 -0.6297738745515791825;0.98199293432641276969 -1.3878746313925804934 0.64728468372725833646 0.80559828178406089894 -0.76409859470589003294 0.16961402289702479096 0.76141212777374100984 0.59445652085656330854 1.4071410000812600849 -1.7331273813096319536 1.2973135566429061516 0.75850320927072845123 1.325024597615888311 -2.5662513175281942424 0.47789911106267435326 -0.76537537710475489483 0.39200446159834245741 -0.20906485756852366542 -2.0624148451044868047 -0.40320712198880909183 -0.47404734010871868444 -0.84268907900390188637 -0.035784787093589787155 -2.388058404513202948 -1.8745167219433362948;-1.1205536148342321034 -0.44204425925655776997 1.1660937805077720153 -0.24703472296957917909 -0.24605820958151555788 0.21656227018912488358 0.064641584967498111336 -0.56529152788700620658 0.77609223052274400878 0.86887026673605050053 0.33018559558312154945 0.043410850177856431797 0.8100532307433935264 0.28013957799221828049 1.5167223152127533581 1.3891185048867551455 0.68243604857879380798 0.08675215943395822582 1.4299331844319418394 0.65575971995993398433 -1.4227547075430961865 -0.12104774531502127122 0.78973627207280616691 -1.5500506935180189405 -2.2202034062411648563;-0.059775475155026333574 -0.27098637459961377738 -0.3798976322967755892 -0.14342788175309170717 1.9508145994394694434 -0.92148355300947226709 -0.42706240525340738445 1.101453805654710294 1.5411183989128229843 0.58986279784868722587 -0.2743053477549789454 0.2579702398884965997 0.29497402884694212677 -0.061279011571659491864 -0.74309172579827309502 -0.13485644851280209622 1.6060132284416641646 -0.40585840561572250618 -0.10440733955144776424 -0.023454384976626090759 -0.012837556907557473584 -0.088495366206154196798 0.35800950347243104543 -1.3007753632821477119 -1.6173616721346946168;-0.9057601340818480784 0.29851695091041791086 -0.93078598044362925723 -1.7764415217918994028 -0.11628089961961779464 -0.92807494008336699487 -0.50563645042091587278 -0.11588296300978725983 -1.3646835929039968782 0.71440096202257707958 1.5228116029311362656 0.67794675405916815514 0.84348727525654987325 1.1476981104908963172 1.2778792723883418692 -1.9046596421217403261 0.74619053850352801138 -0.43973573676982413705 0.2747245114910648045 -0.65523584220298825009 2.5762656009847670902 -1.2298208356925264262 -0.2167623075937221977 -1.4848834342435484057 0.72668058813507618865;2.4804508084724208317 1.1183913337498698759 -1.2292244217510543969 -1.5281653208720498149 -1.6286713100017788936 -1.0431927740226405454 -0.8332253604937085445 -0.2731582891953907577 0.25304329969286071256 0.033771356098294014192 -0.82040474003415242432 0.665177501031068652 -0.65681293939584439734 -1.1942318878395552506 -0.90736541918714952448 -0.1381953788351387713 -0.0073215768274777145996 2.2180342707652411782 -0.60851799360743841216 -0.86754509723183736991 -0.90205208060560038152 -0.2347437012599083439 -0.4240695766483124074 -0.27379837111314414555 -0.0054701349830033330068;-0.29246127410354505471 0.89791671974156150871 0.6653243831014816756 6.0902899093326325897 -0.16157135987579104852 0.4738893033484402606 0.558639788483222266 -1.1519480347813051413 1.2764484165973215557 -1.490739486994739682 0.72818026732401830436 -0.13926588022965544211 -0.95201445084914704253 0.69205315791133681369 1.2850437213592127428 2.2293940439299042566 -1.8280729328872491468 0.48152739906564484551 -0.60980127335797262855 -0.44971479522241192273 -1.4128070143749804632 0.11561668091130739011 -0.53435040766739350104 0.72951076152201410618 -2.2409960880649455461;2.1788597431205349864 0.30887420605005772023 0.32032320759209897965 0.89894096989038918544 -0.18473912681536525104 1.5719393983043044027 -0.40174407322556454636 -0.2621971371333886025 0.035933360970745661933 -2.2054430630889414466 0.27735064067913622887 -0.15516528739165660888 -0.14765553268434911116 1.1647171439244499158 -0.92759896785064155367 0.98995747616556006854 -0.25747380981539474964 0.023199786522257050125 0.05454189151010645914 0.73170720983922377112 0.40162335468203469446 0.13420336073148703138 0.19280622040712869225 1.8769428265907956987 -0.37227684824211743075;-0.98811775059049733461 0.97855454837846844729 0.26479434985061039987 0.44323626301126761717 -0.15564018191146022918 0.1420265317280640649 -0.26541527117584184925 -0.70219716691784161977 0.56403957850485852976 0.47361012108440614377 -1.0222337443145914637 -0.091948864376913602703 -0.15271211169551615616 0.26978253879017755068 0.90642870779693718308 -0.1430415559610775933 -0.47219563659109492493 0.44641242910614148398 0.35747660349434112748 -0.26885102522862308616 0.89095666986964350276 -0.72306789626875944688 -0.56815223032177797347 -0.13841890417258706503 0.14589871733976811363;-0.46643761239187148115 0.71962600571866619514 -0.88646926128964687042 2.4766474788242596183 0.91373300707616889049 -0.98572086255938229371 0.045820577898465941757 0.86034001246096247772 -0.92456445009308185767 1.1161338373862705886 0.78665630183588042268 0.22213398541983403467 -0.97700739500516242675 0.67241350264850230367 -1.2390350693142788696 -1.1593259787218463952 1.0205604961158845079 -0.65732734741855802785 -0.41988533862023580401 0.52900162605132305949 -1.3072306233364041628 -0.14537451075036242343 -0.88933113972089461452 1.5073895102437655513 0.62021391679627768667];

% Layer 3
b3 = 0.52756327222628618401;
LW3_2 = [0.03545661974752808282 0.30556698525879588679 -0.024069415334456577271 -0.096109603799346049469 0.16113875345159031638 0.060655176322555272306 -0.056649119743338541277 -0.029752473061087657269 0.5934317579032953871 -0.021839806423274869818 -0.029595697974191446933 -0.014915079759989544744 -0.052165759831542611591 0.27021229210463348913 -0.011848816045101073779];

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
