function [Y,Xf,Af] = ESPER_talk_13_Other_4(X,~,~)
%ESPER_TALK_13_OTHER_4 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 12-Jan-2021 17:03:50.
% 
% [Y] = ESPER_talk_13_Other_4(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999540494077341;-75.751;0;29.715;-0.9;-133.803853032615];
x1_step1.gain = [1.00000001035699;1.00022983623563;0.0140762793578401;0.000302617642608564;0.251952632905014;0.0422119037568594;0.00432440760873659];
x1_step1.ymin = -1;

% Layer 1
b1 = [-0.094899716929643288;-1.6906844803078509;2.9977103895432267;-6.6205093165435169;-8.2534843865072212;-1.4442721919944259;0.30728714191690132;-1.0806294455304877;-0.20174961522724061;4.4098026244594966;-0.065383495793968099;0.67467538034293395;-1.398818944129381;1.2597944611057359;-0.039466179288769443;0.047924187359726728;-0.20953120988504606;-0.13459989475682532;1.9409671637773294;-4.8874403345595043;3.920419889768572;-0.53459249801718789;-1.607060099583961;1.9720041743792032;1.192778747114192;-13.334731189827224;-0.02972404243248708;-1.2996451403153582;-2.2212912880759936;-1.2623005613921878];
IW1_1 = [0.39304668273566956 0.42684290516653373 0.44032238495819337 -0.35519731244072889 -4.7426939967012354 1.4536325973688726 0.36994775633109261;1.0846679735391365 0.31620289303036375 0.10986168489426244 -0.78547251864975887 -0.69139137788214977 -0.51022055293831725 -0.62031197014335182;-4.0175730120195103 -4.805850240914312 -10.4192682339516 15.211758293264282 -29.759980684951479 -18.920355295208065 2.6232697976480548;-0.48120777287025679 -0.14837432835719253 -0.66861540371696226 -5.3417783029962402 2.8287716149015067 0.73135896342737672 -0.58688493948052056;-0.34887967195117886 -0.41890725238423249 -0.48034845555123806 -6.0836025111384524 6.9455847287999495 1.9851544913630013 -3.5822350244563044;-0.79614093480860604 -0.22881425782235038 -0.85247703900436467 -1.8539882330881676 2.8943941325256519 1.5535672070320099 -0.93213790026745036;0.19037847726976623 0.04682546006843729 -1.3885310759216725 1.0058967155865586 -0.97635880411172571 -1.1414160815570058 0.046013216622688619;0.8913091186123836 0.28013551647026214 -0.01616767726926293 -1.56472225002106 -1.3880622659342827 -1.1050830542286021 0.47346855995131765;-1.1178287835048668 0.94032133180939326 0.75490317406345431 -1.0450360364541937 2.0757578688521598 -1.1029987076232266 0.32390039852320501;1.6102533426322365 -3.1821669761763745 -0.73624571211612044 0.77825710597133102 -20.05430037174073 -16.247890962266123 19.422723904112267;-0.28111295370925354 -0.039475638050978978 -0.51139040857045237 1.6227273990003968 0.20863486689124638 0.35599192662234164 0.8942842364359781;1.4757493778368773 0.30121571923271662 0.78983998446198711 -0.2242573241354279 -4.0024006101688521 1.7293174836166627 -4.792227496717107;5.707619350537156 2.1702967244877387 3.9087571001104831 -0.10114153270250405 -0.095968732962900316 4.9146565772294011 -6.9461522377163814;-2.3281979841840048 -0.05496206603709048 3.6588879450989258 4.8252241715620077 9.8984271964819612 1.5876172141630218 -2.5829655740549651;0.050472567235222872 -0.075111698448221076 -0.65064414162049222 0.53853191596200412 -1.5870456201185734 0.46859779773879917 -0.20463655814650178;-0.014720177010748643 0.051735905330271628 0.46399396603504206 -0.51516139806017669 2.1646308903506477 -0.38299389181432048 0.32010942716283147;-0.047163348726875137 -0.14249579105779631 -0.4918899840456708 -0.6135379395866275 2.4792850054911946 -0.59142222931418698 -0.40187932989990566;0.019972413829600761 -0.044650149576120415 -0.5028236890193607 0.4964324206909429 -2.6145421115366414 0.35351902680965253 -0.48624455137512468;-0.077975158339885289 0.12123441827552697 1.0087204681694977 0.93341221120877449 0.0067518233851542701 -0.81006347421286873 0.2380349870448471;0.22846612553866569 -0.084217398977201013 -4.2401338598282141 -1.3356083566172881 0.82129915183945279 2.1089990387240873 -0.54920656722786343;-0.44533892687235677 3.2781048395717027 -3.5218003434547756 -1.142809622704178 -1.322683816517394 2.5785889253325953 -0.7341369103148746;-1.2687466745582932 0.93873354328995995 0.57771291903427557 -1.535254482722948 2.2699762385800599 -0.66379054263018866 -0.19669772244037156;-0.84739925680945249 -1.1348660379303637 -0.090641780702535896 -1.9632369634120381 -2.2816538236442097 -1.4035410774750601 1.1618195256346016;0.14312687725654102 -0.030657014970394122 -0.5790572284086114 1.0896549273530376 -2.020937271669351 0.32934941964970127 -0.94395506522431261;2.5779110544553747 0.9540749956544009 -0.51481600637464209 1.0749522616450595 -5.7897209305610122 -1.2436908680048253 6.8672860847649488;-8.724947891047778 -1.4270432585954977 2.1860246342691703 -8.0017217175702733 18.91769884188643 4.442284505503804 -4.8461047468333645;-0.56041910619173085 -0.0109302974715267 0.38935479738929379 2.4549342506854566 2.6843429474165936 1.0343139676041846 -0.073032419135944993;-0.11198113075460996 0.02881712602656902 0.26422615339765387 0.41095446939360813 2.5971246986982948 -0.84333455146292613 1.1189184410347486;-0.53091580602981459 -0.090522931884351315 0.0045939765285241164 -2.7888850856495919 2.0055581218534386 0.37132467821142379 -0.34201812549076832;-0.0016231556476592837 -0.59427771838635002 -0.98838251834668323 0.80334965982954076 -2.4692410741355655 1.1121740806507825 -0.48751955662925323];

% Layer 2
b2 = [205.5957672462502;0.34337329780729958;-649.29811533009536;599.70569416895603;-6.2637360491267184;-47.447027672323387;151.89220365561874;-121.13109480596815;1.9856774435523139;-879.52362132729115];
LW2_1 = [-135.00940065579877 -9.0238343797313458 -89.300051741904994 -141.23132893406813 200.25529933432182 -81.996219034144801 -90.159859934875456 -264.75630298279304 -38.096224299629071 -20.171844930553895 -124.6633486266073 142.4869057179136 -41.646825534709819 77.224111372830677 103.55019889828253 191.80813772978428 -36.560526052890523 71.327330284271028 -176.19503020403027 428.48881809671826 111.22255378954156 -33.336425311370164 83.056025702828578 27.814671510028354 -40.550754762311961 72.988786053431809 141.38050953544072 256.22286954377131 89.340928599176451 -229.3564252275649;0.22146794314162296 -0.24083813259038761 -0.019905950952852546 -0.2770610160634388 0.080867837471959025 -0.17038119231988236 0.43075354188102921 -0.11986253851401579 -0.33844758971561179 0.0050615000354023503 0.47117151518567024 0.017162881371681749 0.024701673124068634 0.028418877139360109 -8.4755664056425459 -17.950221792859956 1.312189361998896 -9.0698341642626765 -1.5030826177831964 -0.19801455769036455 0.10194261822587013 0.25290992061362128 -0.05385148221356878 -1.0760809591811409 0.018961857839462061 -0.01021206657247971 -0.24868553372209046 -0.87202225538781353 -0.28733615442210475 -0.40480247158278948;25.877684227801719 21.405149847994927 -1002.6073225709824 -18.285912584495954 58.59325272401162 -115.27423757359148 -24.843446693875055 -1.9392271004656239 11.163662515575275 0.025868540169739669 3.4007374020864822 -37.530091789180076 -7.7539209001753147 10.608774878455026 -38.992075931807925 -5.9831578380115849 -105.50699648641988 -82.201234985095169 54.947345519841186 62.346396192735249 -124.89303486683571 3.3286163542291276 -86.011804204869236 -10.617086110457068 58.400632150541803 59.931285464677991 23.957859223875737 62.373056939278406 -45.602357577180634 -30.418502957605604;225.10655459055511 550.48166756643707 -785.3507833563749 749.78900171532155 404.61505506245828 -69.32148265568182 103.09472427262943 521.65847274371572 -526.4432471125283 361.56610111815075 1240.2125566818456 -1435.4019591745046 -277.7025235191158 -872.50482869368307 487.50233006520222 641.47610675974965 -455.58429645852453 -278.74746694326325 1184.4552216487441 -592.27474824121566 451.54536625844827 743.28063783731488 212.6646117040514 540.35222397357177 -560.91514237015542 137.21461864348393 903.221097866778 -343.23042350285851 -58.099960917981093 1812.8235735452984;-2.2930127476543585 -0.50137195317125727 -44.540660805109596 -0.79068216853776119 0.38520602939124088 19.814784144652993 12.605838019533392 24.183068364822777 5.4018248030577674 -4.0176309946709257 17.02747052760845 26.263536033341342 4.2868349664602974 12.606314062322618 23.26871594586887 5.9533580551180254 -0.19418273419326376 27.993527527957028 -7.6695045639578927 1.3112986249686496 21.240845933953992 -2.0505359572955393 3.7070267354426507 1.4491616346609135 4.690417506403298 -1.8377184882069355 25.705367157586782 -3.0521660337517571 25.209039129487728 -32.030694824626053;-76.92381652188709 -26.46056229727472 2.3533494229930225 33.962565374014837 244.73281693828434 -80.785855758087578 77.458848072038691 -39.732326462496715 -65.859593670851964 45.650427360105184 78.780930106024385 -208.50575038201234 64.146640734726077 -53.208911581777748 -93.760028573850704 18.859210072912695 7.7213874470036643 30.814458086316318 25.997798497702398 105.52949153362439 20.955156605985987 1.3248595758979111 -71.378918162744412 4.0923041381070187 39.102903849502226 7.755129470374337 26.897362003887157 17.026170907141868 -22.614625784552452 -69.158884949444655;153.0817754259854 39.016567229358586 -157.27532531634881 127.72013558002416 -59.060303125263147 -8.3852843999710238 64.744428860575368 -62.523488913094589 114.32871432145299 255.43687824321822 -2.6916585097954093 35.952701467890698 21.13150142153847 10.521179281349795 38.501983738502886 38.203299563150296 11.383985984315308 60.005119090586568 -21.948073967153594 33.441478709181261 -126.65401578536405 -7.7523848991960653 -35.935484490339675 -34.146107267663801 180.72798303047423 -18.470771148411828 92.168221041731314 -68.734752752163033 80.395649000763228 23.347005391366142;40.968664123108653 -60.780660276304992 -191.78170938393953 -25.768293555969333 18.604874078758588 8.6668411176611873 6.8933734652570147 2.6676226787886508 -14.832503196851929 -1.5507164895232903 -29.094754674134538 20.760818388591112 -27.254429073908138 -24.069829282541651 -17.966323244851974 -44.485331772270385 -34.146668095100225 4.5741355859497173 -8.024793213196034 -6.3956715040088827 -72.219151855164199 15.106542518261465 -34.05656117440477 17.784430187753856 -12.636303293863454 -7.1721528019205891 -12.821744150272032 -59.753576781728277 48.625982947962342 -15.703022664569639;1.2821532161036244 23.127240692326655 -61.113230918585764 36.331587162604748 -5.6931787240166658 67.485432689632802 41.024836208111374 60.684954976306081 -48.544316707486438 -128.32683825315416 3.3497625765539234 53.126216453561639 -36.544221151008529 -40.952452364735748 -120.26649128955773 101.1619964476694 80.512339574399434 54.831381085859071 -18.719682154908863 -51.085471436324283 -78.32398236552639 -56.175022499050577 -118.4148819274467 176.63157458423609 -39.947919274533163 -9.0036734379272669 -58.80219785624309 87.608885437718115 31.369798142191673 -14.569391495266;124.25375651757317 -581.05324757958169 -625.24168152744039 57.041970436502638 247.74421067941918 74.898342777211184 -30.164279520258486 -22.337468115236767 25.72215359829876 102.29544151297279 192.55901019390794 88.162663836539593 -9.3049237901655708 -16.872855606616419 -377.82532417685769 -6.8622990901225185 52.241881892598229 -228.33876202149887 -21.777072159843733 101.15427494914135 -18.485683510879785 -41.325669523930706 -165.55524785306935 18.923276127526229 94.352688496640752 -7.3846280208819852 -48.637235284667014 373.44776810466664 76.121545263266256 -207.36164325411485];

% Layer 3
b3 = 0.19903267216895679;
LW3_2 = [-0.00038908274256316321 -0.95840391146146098 0.0034728269158658624 0.0032511781378825296 0.013010708495979665 -0.0019541611395090409 0.0034205346294077375 -0.004425495215763727 0.0040789802077916193 0.0044036775958292114];

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
