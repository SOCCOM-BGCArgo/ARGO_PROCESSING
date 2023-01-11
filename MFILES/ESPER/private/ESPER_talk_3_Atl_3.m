function [Y,Xf,Af] = ESPER_talk_3_Atl_3(X,~,~)
%ESPER_TALK_3_ATL_3 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 12-Jan-2021 17:03:44.
% 
% [Y] = ESPER_talk_3_Atl_3(X,~,~) takes these arguments:
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
b1 = [0.36198815282881208;-1.3794304649258144;-3.2676299733624803;3.8236001753863045;3.9012504661751324;0.42210443976597695;-4.5089540942849284;1.5847847644436919;3.9204440271634362;2.7275246554168318;0.14497298901510891;-6.0273195007003748;-2.7901690931467398;-2.7544169897267832;-3.2489734533604433;4.734097676905991;6.6842520961616705;-490.81162389510723;7.1779672999466468;-0.55836857833136055;-1.9066297750261778;-3.3369122918314478;3.2405257388817219;4.6164427005666901;-0.93606512467203518];
IW1_1 = [-0.11275611884794877 -0.063349264413691944 0.046653078898576623 -0.0338372343317761 -0.41293027750718037 0.066881437492071594 -0.070623254698669632 0.12817136927607056;0.83126621065827533 -0.73242726081596166 0.11002222231643277 0.19199174485090423 -0.10422587016949528 1.18973774516947 0.49084946386084594 -0.10106630230976216;1.2048046474900294 -1.0501569160013797 -0.31232831261234995 1.0713029738652224 -1.9334748440354155 -4.7821909403148064 0.54621906690460942 -1.1372917123875093;2.0104094777705979 -4.849775588129063 -0.63545995042869219 -0.14735791432238887 1.4614940115482431 -4.6144286858695587 0.53217983215093778 -5.0452083804326078;-0.53929592313015939 0.95115334217260472 1.3326516910478809 0.53988475960348214 -5.8064361329028902 -1.0447460686264112 0.23046206026928706 -0.092986329736373871;-0.034043768614715877 0.0019735805110528354 0.0015260296761974875 0.23373737464790451 -1.7363580112563868 -0.37479351292405988 -0.085146825849418309 0.037465467707094496;0.04593380490848218 0.17686026987100595 -0.098233863932005103 -0.58634827385827215 6.458908231452039 0.33098815251382829 0.18680896816876574 -0.36521570376402679;0.75633943363239187 -0.93030567412880971 -3.8437541536680171 0.078273887901426684 -8.9177509430409536 -1.8264136666494475 1.3264872893313882 -2.6730633197850682;-1.2330236270254535 3.1923173044019877 -2.4957599409901632 0.50571962490174904 3.5322783752642763 -0.10893992483595932 -0.88922535286176241 3.7692433895426625;-0.40014205297616257 0.45703342731215008 0.079849337761728284 0.10830866390361207 -1.9235835198855042 -0.73598067560435332 -0.26872234241410281 0.022006000694091293;-0.020289417212983284 -0.011666603348476952 0.0089866202083964657 -0.013482203005198132 -0.031549121454893858 0.023537267715914169 -0.010562192261531931 0.024101941728525544;-2.2993326015710891 -0.26600281016788224 -0.16895889824020985 -2.6903092530860011 15.15878301345635 -1.9478499434589032 -0.21547598955445696 3.5917695849649376;0.60751217078172026 -0.76606773174125153 -0.10816561457199368 -0.028885017285782633 1.6015806015046521 0.79174690961370797 0.29244908076672771 -0.050736685760686275;0.047242288099593438 -0.72927314026734591 -0.63371050339002477 -0.21463891586214134 2.8204984631629064 -1.3909964341698609 -0.7035506182404897 -0.19540410975329617;0.75718418763981199 -1.3765126941919636 1.8025787624194765 -1.9165072553736644 2.3361806203654112 -1.6084759464469385 0.90192575477559611 1.187260782517332;0.35896066835264268 -1.4995592824634094 -1.9765018555471259 -3.5305973594308306 1.33416983548198 -2.8874742453720827 -2.761203672900399 9.2457676442743857;2.5630471000874917 0.13747270057712649 2.2233741105615303 0.19341553519577825 0.21084426238326279 1.4663807075981998 0.50990610609380715 -0.23330950784772739;-84.84879473460964 -41.546763663467082 -111.64302621170276 -178.17230183683671 25.839588728223241 -221.34035822146262 -68.246388280743361 1.2707849648907985;-1.5732192457060206 2.8608436463455647 -1.5965075658573507 2.1978353026243376 -11.858566257176925 -0.21930769509525763 1.915644621616073 -2.1510370314548259;2.287356009740471 -4.0510428622638184 -0.99562230970490062 -0.14156985156814145 -5.6652107527272078 -4.5746229655352515 -1.3365007286050208 -2.9483721257164759;0.85440844372212443 1.3529308026734128 1.3022754605500089 -0.70085771872304148 0.80706726443639498 3.0774169385836818 2.0625563931719832 -0.91761840801193528;0.2624461617252063 0.090874766886716529 -0.54867962919392421 -0.32215677508142326 4.0477016144817615 -0.979686842754529 -0.13662428489992304 -0.21253346611485138;0.023399103219661932 0.021994899678376513 -0.021767358014224879 -0.84245154100457997 -3.9131047694150949 1.8436135545957211 0.077652304426664071 0.1032486962684013;0.038426025768483127 4.0239570781326597 1.4380654255793928 1.3178072922194142 -9.0077743200905029 3.0784700239867275 1.7130414568728618 -2.5760363193104321;0.77470803922805875 -0.48321846316440836 0.2073421724300003 0.24616147514054001 -0.29914514642712658 1.3727846190642332 0.60499940322969048 -0.090427361231936315];

% Layer 2
b2 = [-221.47458877723815;-210.39212163466692;-78.338096307503221;-94.694530680516039;-60.237132678122997;-75.95504941701904;-86.718009021974794;-57.339994248615419;244.51347133722152;0.37133489991527396;34.098463900288436;147.51653953523751;5.2532939395276887;56.79390385283503;134.2906930479989];
LW2_1 = [392.32287621134981 -27.855365838531302 14.410805930737013 55.905388074330858 18.974130517931581 -30.335642220596995 -29.211349625484416 -25.788449623815545 48.862273679170833 -77.501736476289096 -74.41293914772136 6.9369016733856022 6.0000420536247034 12.571582112103489 -52.01443567532516 -4.0558067768177306 30.610535524444067 -140.20524662654711 -56.476954002061511 -11.566921185375836 -1.2730686978352457 0.07329908025724001 -22.063246504830435 45.420475434730129 5.4342353930409342;9.0530190513003994 6.1439783841534075 0.63947870361955406 -6.824651078356549 -11.444898580873371 -3.0316771139591063 -14.950521521437068 -0.038946756383194447 6.6997928323364535 5.6029157316498273 9.5925148932373094 3.0866903217158761 -2.7962797748733106 -7.4613808492440246 27.392238190660812 3.4145091123793758 127.02104536746603 89.286454810587799 2.4960426654402941 0.35458558435888343 4.9045321726002795 5.8774972232376381 0.088217030440445376 -7.3490291979551712 8.4413518277334916;-8.2906644698701051 66.427414598783542 -7.1942465981069468 55.144521248349626 -59.414833275507014 10.88845892134966 -139.59383714480421 30.220964011915544 9.3117171110861072 -2.3475172724678544 67.3578189050302 20.644955296413428 -12.15298356759512 86.095430554742705 -117.77684655572244 -44.528945385634245 -190.03014537679692 176.1560675722254 -58.684422021016481 -3.8463490099840651 -4.7298949517479887 87.352000948599354 -40.721253866089015 -1.2588149405575708 22.87573576439932;-10.827015730507071 3.5547974762636469 -3.1765536797031486 7.9338473059278227 -5.0183730419095074 -10.962592061474552 -1.8308348629441153 1.6020159970314181 -2.4405903777520725 -2.4225920508694285 -3.7366807733650678 -4.7825767512013035 13.393105333593871 4.5415036514463125 -3.4434661912395477 -0.26725630301944281 28.125762357901582 64.877665026725026 0.32401421775117129 -3.569977635001067 -0.48826722660340949 5.2634203557215589 4.2942296097103956 1.722699140975857 -0.65526359629259057;-34.71136119808699 22.384560890076312 -12.047364694343898 55.571843175024384 21.188700650359095 38.921840382638962 49.352069649688076 -16.345071601833162 -39.09222725832813 8.9663306322952732 0.9532707840821476 9.1278745665936327 6.3325254831754689 -24.077299631759576 -32.137267162138059 8.9598990408484784 -0.61070892762829665 -8.5144213916557305 -0.67984948214885521 17.801284997268574 -8.0506469035984622 -17.470409872060962 32.178676670665773 37.428084138671842 10.336427172421004;-5.1683652985244057 -5.2264584696321696 0.52644829246686742 -5.7305976445023594 -2.3655425336877185 -10.255376770990249 3.5011028721267432 -0.023793559304170042 0.52858824790846848 -0.11603560510127364 -0.52651362470337615 -0.059889156852553527 -2.4203201432870016 -4.3896193100523542 7.767086717279696 4.1055440866792221 24.727456597018932 48.997644268217073 2.783068304468713 1.8156875925845979 0.13665403915993402 0.30349559529304981 0.88711169497928477 -2.1065995367323627 0.015831162066277216;324.57008378988678 134.90686579216913 61.76991913396364 173.14904465905619 70.962018184524226 409.23473548226337 137.75903411328983 13.870566618455868 78.159384518999516 110.73028002875418 134.33991114146636 1.7123711005087716 158.35788544243681 -18.173473508339985 -258.54855474736144 -240.73364331452879 163.2332086904967 -121.29463339090537 64.296466133317622 125.83624078611845 52.827175591648704 9.4021777619559082 60.030308194649265 -18.624746906481285 -69.266894327386453;-9.0534124246578518 82.777283655447235 -67.170468866993829 -8.0452663819774006 0.81882728051730291 -55.920962285901517 66.199517519798619 19.020749318788127 -22.537966258420049 -40.286571519309597 -103.88085674750722 77.769426961439237 -30.64874384069137 0.42052980512825011 -51.827575328212369 -21.175033780757772 -113.59643056935914 -35.222288403424777 -11.004261304567855 20.594801096837696 5.7595056238222986 124.91689708655996 34.780662759442983 -71.244968215790337 -30.463193439269986;-3.3471570080646798 -18.35933006940488 6.4952060572956896 -2.2929789114209633 24.279226293427424 -13.724477699459239 4.1657293014605239 2.0308797750057086 -11.151238369085563 -8.4872421241863787 20.953209592668561 31.635505924299125 -10.839654644817964 2.8899825136481381 53.713652016322065 11.100134525644368 -17.668802472558742 -236.80242650504582 17.297350436760574 63.980141986443883 3.3772925951174426 -13.720699752624718 0.9796767228167349 28.487744580786888 8.5980845817000375;-36.851247631099511 -41.308265589969039 -10.518134295875987 34.006304875299364 -16.900250556562099 36.341391762826241 21.050202478463302 23.927316536318681 7.2802392930034401 -45.479237610767306 -20.26241502526576 1.3929445988509055 -51.062356284771738 0.52735549809685767 -66.209341179767151 6.1328666583631106 12.751064037955544 -14.479752950641908 -6.5780392977015385 -10.034716340509304 2.8621322068781581 -66.988724725649647 -9.3257327454849666 31.321428407253183 -61.904163420797637;179.90543285949227 -11.383802221876861 135.38220327336262 -22.122777631405956 47.746595151608197 215.61028668845918 104.02945706682951 0.16808005280638125 20.632797981434059 -42.796020347332664 36.514332630651943 93.775840758869251 99.540958028396375 45.383136271119305 -93.948759554012483 -15.801595707359173 -14.906040247310077 72.174189728480755 -1.4665558614500311 -30.624956921821489 -11.082770261114728 -39.690322759903673 -67.027322419265545 19.70797965028914 -23.615533316944546;-240.93419666456637 -23.481647315281261 4.6811668531479675 116.71516592955386 -13.085135531033655 -37.423943192482241 -65.551354105324194 22.194127550091824 134.78876860274599 51.596378007401967 -193.34263798810488 117.92186249294338 -193.18920377922055 46.60952333706576 276.72530124453431 -78.155742435343782 -202.97140883197042 193.17375383888378 -22.683987087555725 91.929441184606503 98.528512523147313 -33.817535316042317 34.755323894932758 141.54008914044502 107.95557132263389;-8.0226628635151727 -1.7910346341074992 -0.0089706777772862716 31.106610822631566 -0.020359731109702888 0.97156967180691112 -0.069505048286540472 -0.0032527456341351316 0.0057346374772510988 2.3746653224867802 40.40251935491203 0.010613825474837222 2.2253083298316558 0.10731455294278834 0.0033003600721321111 0.0011362939208224088 -19.325944992276423 22.024195930984465 0.0057564358487702891 -0.0075490804733668974 -0.018744945526960655 -0.12030985907286267 -0.051400598167125215 0.010611841956238792 1.3513890573322427;43.122409686257882 11.048430617589311 46.016732336265754 0.90421829478286109 -46.099231581979303 -78.060012170043137 -11.476085009195709 -25.829303644093507 -69.641081792069016 -45.214493411160305 -39.700396308646646 50.298864622501767 -62.002679390970059 -5.3277502539051609 73.431046737270876 23.049906323052021 -41.679311075389755 -41.339328614685762 0.074450722049100734 -2.9627257594036882 -14.872032698211243 45.943456766330833 10.461917934006355 25.886636002039118 51.695641060029502;60.789859046392131 -55.843817197965677 -9.4917413302427942 -15.712748326363039 7.5883680310827488 10.393307275471388 -95.371509425945661 59.534275731717422 -17.232032665162073 38.760558859170935 -34.952367284615001 42.622823418406583 -107.98255007718282 7.334166032848124 -73.688700425560555 0.55347142778591774 -270.2605644184672 82.627435067222024 25.54065572653256 -17.489133346884081 -29.735086413562499 -27.174247091828725 24.198135450073057 -24.961973639658581 -57.215560598764291];

% Layer 3
b3 = 0.91886933982820762;
LW3_2 = [-0.0016273064215390984 -0.50849359998255528 -0.42269949856276445 -0.54923214090701789 -0.00076745553420473763 -0.46470971340219708 0.00022223373464927258 -0.00029339478776208483 0.59943623374919308 -0.0020014979490707453 0.0013675738530719752 -0.00054722475305581293 4.6318039467128731 0.0012815132865144901 0.0014795383777496238];

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
