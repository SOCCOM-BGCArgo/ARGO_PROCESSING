function [Y,Xf,Af] = ESPER_talk_14_Other_4(X,~,~)
%ESPER_TALK_14_OTHER_4 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 12-Jan-2021 17:03:50.
% 
% [Y] = ESPER_talk_14_Other_4(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999540494077341;-75.751;0;29.715;-0.9];
x1_step1.gain = [1.00000001035699;1.00022983623563;0.0140762793578401;0.000302617642608564;0.251952632905014;0.0422119037568594];
x1_step1.ymin = -1;

% Layer 1
b1 = [1.3824571662618659;-3.2704064555290109;6.3744264525942622;2.0370469532030855;0.69905359013881085;-0.96893587258700808;-5.0691046966448035;3.3687305737634827;0.40802730615624028;0.8594234463496806;0.95276499796802783;0.29353093497924165;1.5938266622761073;0.43418081305755174;1.9450026756093926;-1.1017203198666208;-0.74340865175887294;-0.076509800071105083;0.43406839230572447;-0.87139311530457886;0.86298889734819717;0.11421534847107773;1.7476446193732418;-1.2515040988046355;2.4092905456365519;-0.11431124511748955;2.140960809232991;3.5886880715371721;-0.22791659335445238;-3.5671909423547308];
IW1_1 = [-0.6597774198052575 -0.95542360046329311 -0.068832925322848001 -0.40820769705590165 0.86053665272161195 -0.081683811210412863;1.5568584289697844 2.8907237441772713 -4.4749240277044882 -0.76888538769711423 1.8127504439884992 -0.2356763328420306;-0.36184377251979111 -0.45929306314612767 5.3762312349124244 3.1218598263323352 2.8158476491505584 0.46843536027006139;-0.032433081281063258 -0.38549939835628366 -0.87626329001942327 -0.57937504281061458 -2.308107662661893 -1.9868167939745345;-0.59549013462883982 -0.20755365556070446 0.41955182973230548 1.3056829054652483 0.59581166763140625 -1.3188339679112722;-0.91881911092348101 -1.2528342601929812 0.10669843729963477 0.23542143648651928 -0.74463783813453999 0.19567884165699015;-1.3988865991418724 0.71032504797285645 2.0557931233537934 3.0602171097662998 -0.49228417427292265 5.0194069785729711;0.015875807199452623 0.076427091559787108 -1.0622589745103397 3.0551193031267943 1.4683223601400301 0.9178509782252291;0.35179331712643308 -0.71854411810171526 0.16905666707647546 -1.4380752450089085 3.4416828283436423 -0.91624153981266743;-0.35005181165053839 0.49517908655464743 2.7455543157693398 0.75762843350466014 1.3160571383526574 0.2622724746349816;-0.24283900755223317 0.55853919707202049 -0.11856179148279655 2.0889775082946054 1.5767220517701122 -0.038167222533907745;0.30462323648026807 0.23730637597016704 0.50753265191998742 -0.17783954857892562 1.9553180289479377 0.00062144960649310188;-0.80384623267837119 -1.1684527191246643 -2.124773970073345 -0.48834529945083233 -2.8136511336383379 -1.089479824806479;0.78221653175407657 0.43207747900635529 1.9503677887661639 0.24991136375233317 -0.076388043323657909 -0.80046752067884963;-0.59086062749857859 0.59179717509309426 -0.5148844790999203 1.904333747874726 1.1836601217959681 1.1057198317646735;0.15300218001812738 0.019549768409831577 0.11933945946327751 -2.2841195994596344 0.92469964903530155 1.3329509305019529;-1.7444774314253872 -0.35622621865023973 -0.48814897217851094 -0.4836031494467713 3.031986838288995 -1.081575373496823;-1.3855278881173618 0.2564906661700051 -0.50099047980908451 -0.57647636876495589 0.95134006556959061 -0.76071880828182603;-0.34456374388050248 -0.53579013955888311 1.6945562469357971 -0.82187901891859771 -1.4333044318029668 1.3141607309063408;-0.026901208181879889 0.16349098205216028 -1.2335790900531312 -1.3912468276338843 1.0502951298744234 -0.5238652872392372;-0.32847228616024027 -0.84026282696483934 1.2303493480385304 -1.7079737989103525 -3.5460134497114444 0.20802734229505049;-0.19921681910200875 0.11392355791144156 0.9593562446176751 -0.5787771973784861 -1.7865545053609129 0.34226644101149922;1.0486298562987861 0.79073634528432313 0.098447361385041801 -1.1050022735828653 -1.5389982603643266 0.028045525670591237;-0.14573382820079397 0.028488607069568503 1.8562427482672639 -0.096023467679655541 1.8304372436174876 -0.18045080989130777;2.885260441140562 -0.70787494236595949 -1.5246400419745267 1.7179179022547733 -0.99281684470209997 1.5041240850238462;1.2803829972551326 -0.99308140340003304 1.0968265122054235 -0.89898089625950117 0.63458270114062476 -0.037672283214052027;-0.18509152598450662 -0.010833705850526341 0.95721536144234909 -0.071456223360615578 -4.096430676738084 -2.0585354021002433;-1.2236781907749055 -1.1908127834322315 -0.11748030825258277 0.7135537911572577 -0.9092170218794231 0.53573400128526871;0.029034994379495926 0.16436710997430379 0.94138643487837936 0.11808230425291294 0.62422935066440977 0.039696482888364575;-0.50309179654439451 0.0090583560605715721 -2.9956527808469753 -0.87458152551829171 -3.3393727593139642 -0.80619429066193504];

% Layer 2
b2 = [2.2978679478611497;-0.59356812530788738;0.37230189352450938;-0.41555886065967379;-2.2002672507728662;0.39390966843960096;0.79888816058832479;-0.29539000240108138;-0.20348421111477683;-4.1136212001150927];
LW2_1 = [2.6539621550449111 -0.38725489408175351 1.0065327337389938 0.26518299845924415 0.98430627094525258 -1.7788574504333585 1.5322240458576195 1.2411928520860431 0.1744917207475418 0.31725506266680287 -0.87011592075953503 -0.96901817417635372 -0.21719883885361704 0.82040770393527684 0.77345321723176108 0.65843879434030783 -0.96810632760135684 1.5420600214994511 -1.6175623592955144 -0.7907946865830997 -0.98509959305301609 1.1038084091705918 -2.3665200735329606 1.5207492857630844 0.077737150637363209 0.85767775038776339 -0.40460727721834289 1.099518391314805 -4.1748811951358764 2.0857357214307535;2.2947557474763856 0.25141997694148877 0.51384907126874335 -0.6392534924588088 -1.2932695162502421 0.63145935113815266 -0.10287136532715602 -1.2756555403814043 -2.3312386526046436 0.37471453253692683 0.26874380695075567 0.27743848772022389 0.46790987931596284 0.90664404682839883 0.24285986847808694 -0.71299261655975654 0.42565878760085735 -1.6798906594064256 -0.12086035622272356 0.18763049772680071 0.27506533507719411 0.76796238389290572 0.81027404508870293 0.71674063550870759 -0.60242163082963207 0.23838582017505966 -0.28904708538589241 -0.18298397982547429 -3.2332189615773519 0.34964395573349144;-0.50511285135150064 -0.72399183371915121 -0.91857488522560937 1.4912460938372756 -1.8843038369651086 0.4014422603859305 -0.50688145700203302 -1.4500773054669902 -1.1776263704778371 -0.90489206472578765 1.7353656721357067 1.0322813606218915 0.039221744802831346 0.078654278831842867 -1.0159752212026441 0.09595770044195101 0.69133038151662773 -0.28678194760410597 1.0152962741980158 -0.35395368440560437 0.23348038603635224 -3.0260002052482209 -0.27498565272218767 -1.0404339065487453 -0.39822252502381883 -0.54586675645276428 0.69600270328652192 0.006061981212377901 3.8083705922614906 0.58798310557095079;0.12223619211691227 0.073865269455077473 0.13730409575777308 -0.53214418624118376 0.20229057436009409 0.59986435479111355 -0.19693538143859898 0.25362082921014184 -0.21282625486364579 0.34620590194507855 0.65697264753104689 1.8918150458756346 0.12698036639816063 0.23751860756561233 -0.29101094678440043 0.084847034287679007 0.092541036537514818 -0.28168978004390277 -0.014326194760926361 0.19419004274442767 0.046924658828816489 0.46254105489688668 0.28862762146068577 0.7070548849754279 0.043997706513063159 -0.0074780985110570394 -0.48569174645774749 0.11809840421134209 -2.143498249247604 0.88496337748417697;2.2315552462164812 0.14568572544514108 0.25860945205526342 -0.62504574584351846 0.22286913130403149 1.0352629496637051 -0.51969909999378139 -0.34368620396227589 -0.30749110815729574 -1.6874383706781375 -0.29845195767025634 2.009925998358951 0.6873527251733017 -1.2320894486495957 0.39464042287050055 -0.19663949568744427 0.20131120823077303 0.52007099907284882 -0.46381395751403892 -0.90423889064759566 0.14245656652675026 0.3789706499231501 0.50976543470525015 -0.97813033291237972 0.08080001357874704 0.20396537347344962 1.0972079859709341 0.014174255084392205 4.7187339124355923 0.64013823847653195;0.029092956464177394 -0.10725477122622026 -0.48105895222676626 0.72152515030731357 0.10410070146950476 1.6124266928417772 0.021753850732168612 0.10399326633697251 -0.42128794201814751 1.1173075598091715 1.359344275807038 1.3301525500344311 -0.26474786666814026 0.90538729872171597 -0.62552005935369426 0.28761495572390683 -0.3539539484498751 -0.5063449829454052 -0.66732222131584562 0.67645376194072782 0.17974746599961022 1.5591440585184524 0.99755636961327843 1.008884567352996 -0.068191657663654454 -0.23926707394513111 -0.97519477783938735 -0.034121164095746101 -3.9824253200777688 0.26486663877574479;2.5473277765816236 -0.084947883060935286 0.33787241071284291 -0.21593713403362835 0.5543274971428318 -0.71138978928773355 0.80528416926983915 -0.89613203801916197 0.18970010464299875 0.67262629245137828 -0.07252698438121355 1.1526415695807946 -0.11528886234938002 1.1254964702145085 1.5182029062098743 0.36143011266388875 -0.87690197004175052 0.91269611735835887 -0.39362722220345669 -2.0095263987841534 0.86007235448691133 -1.3509935140710452 0.31488726364156044 1.7659395371507784 0.081089099994436489 -0.10645547639545022 0.89200282439301348 -3.3818732589932274 -7.4412006617270983 0.5847480373129369;-3.7900360264364088 0.27122828580594277 1.2345598322281464 2.9430061705055617 0.65998894058155233 -1.3187853552676134 -0.60335126489598578 -0.61514706176339828 0.67074342729183067 -1.7217518006162433 -0.56244433369224711 -1.4290577001477134 -0.80486533096518942 -2.5530358775377731 -0.12887188955112688 -0.5246838004583243 0.015363272249908806 0.12386582220341287 0.42484347250152987 -0.25768625215836299 -0.97455938528319686 0.23986586746898139 -1.4171087008424346 -0.83144445012488033 -0.16686475574835369 0.38894329224353108 0.89519974440742622 1.5910046090146559 4.2496492711575105 -0.51027791568635694;-0.05445678372922922 0.13679948350101145 1.1763036479865681 0.64979582824650628 1.2055391359616514 -0.80495777062649232 0.17226870085802296 -1.1587947586167646 1.6277053469931793 -1.8259996897862747 -1.1222451756980403 0.85254160608914031 -1.7887685035807859 -2.71741295381953 1.0736623771583884 0.23751342699732181 -0.11526450898383875 0.73075771849772508 -1.3198798450861864 -2.2841151223106073 -0.88717922174802544 2.2507464964912005 -1.6659431992021936 0.26744539031707615 -0.27327703607453746 0.38490610540859704 0.88617225191195559 -0.050015681186084661 0.02902969194725679 -0.49499806875320418;1.7088467718815894 -0.28570432899973552 -0.13652703535811911 0.40638267167569453 1.1526114121406057 0.67607712756862881 0.31939273849729061 -0.99726563590592199 2.206969920010986 -0.8029383402648379 -1.528058793598182 0.010335773390483598 0.40362313090199498 -3.1582079922990882 -0.49916442259777749 0.5091573180111939 -1.0019871352166649 0.92414125565047212 -2.0669187710982841 -1.4462920109823056 -0.43477988088355596 1.9570889892318146 0.33966989659665453 0.8746606737752356 0.26354961561298323 -0.42316641311623815 -0.491914258515529 -1.3729532633814303 -0.22610273578116119 -2.521825042402988];

% Layer 3
b3 = 2.7118605742045365;
LW3_2 = [-2.537822281117907 -0.19096813716047814 0.2615487530244795 0.7302364012992687 -0.21766786931112039 -0.29878444954508793 0.40816817918021764 -0.17994762923178842 0.14680246602722322 0.073396090793604221];

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