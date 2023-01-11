function [Y,Xf,Af] = ESPER_talk_5_Atl_1(X,~,~)
%ESPER_TALK_5_ATL_1 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 12-Jan-2021 17:03:45.
% 
% [Y] = ESPER_talk_5_Atl_1(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999999390765;-40;0;11.29;-1.99149636107278;-0.2178;-143.241523403244];
x1_step1.gain = [1.0000000017546;1.00031699124092;0.0153846153846154;0.000310077519379845;0.0705602618858232;0.0604107502911417;0.0479365703301392;0.00458552991639557];
x1_step1.ymin = -1;

% Layer 1
b1 = [3.0224102412764928;-3.7959924710676916;-1.1740033870496636;-3.9960264679624418;2.3051741600669344;-0.6979911818859883;-3.1942584542982582;-0.049292636930522028;-2.2018447511820787;3.4750117287765661;0.90820366084834336;1.9278263557714024;3.1859479212467345;-1.3310067280944191;-0.16730287959126761;2.6693331379440135;2.3297806404809687;1.5647778891137514;1.6160834253406335;1.5407816173293123;1.6346469281831255;3.9179989312415335;1.1461425286027107;-1.7284075225039148;2.3693853444697499;1.0870549645615364;3.7733149200435543;1.3024557534026475;-1.6340410058430423;-2.4037770429421395;-0.5945850142269008;0.9397462075783799;1.816385513008995;1.6818448312571119;-3.1531380167835561;7.8821134208421402;-3.3199599576546981;3.6919045279206042;-2.7609462070269721;-2.1040730496958786];
IW1_1 = [-0.8974297885259227 0.78841935649511719 -0.81623382944378575 -2.3364060545261398 0.42433083377997038 0.78158743704098788 2.0917234286837458 0.95603078468868874;0.15939504990163006 -0.25339605133958842 -0.67964380970161031 -0.95828631881151127 2.2311671784553497 -2.4122651854009187 0.13484647182340842 -1.0199550923611911;0.06542296660657515 -0.63614252523158199 0.067696119753686654 -0.17040096419824508 -1.1351563240602935 -0.011483118453128225 0.23693288751505742 -0.10320411503175825;-0.13381423352035762 -0.84165337901539516 -1.2151709494228835 -1.8899699499198925 7.1145576761393912 1.2731558865843486 -0.54728594688183563 1.3921474382770367;-0.47318176965173991 -2.8297235931928335 -1.5363046239792431 0.75007324970784683 3.0453297041101446 -2.4347304426004635 2.7546574334261966 -2.5377081603019831;0.34179600061593024 -0.53268230941879557 -1.7164641345628255 1.6134849852654845 -1.5436571149258358 0.38486016526044919 -3.6773231896231153 -0.42471329092371829;-1.3179758685104319 0.0068064690532133789 -0.085564738606714583 -0.057947904607323568 0.98539471793348454 -2.8663378909284987 2.4495496113328743 -1.0084199953174855;-0.099549314471073305 -0.26448258131425945 0.22710151690431132 -0.10314865221450244 0.52563840296361053 0.061310990559189482 0.12288859040787954 0.07068611823035785;0.97653603441236136 -0.54705196249764798 -1.7809118196160243 1.7765334409725488 0.29862644373306768 0.16662511007176928 1.2738983920069427 -1.6451151518136262;-0.98324803969807106 -2.0382781879208407 -1.103345238823972 0.4514171307373806 -2.1759076675408902 0.15628260362177004 0.073233636440602934 1.4082101885566427;1.541764703063047 -0.063612106879698704 0.25154952882076426 -1.1050894343726436 -3.6323818461335358 5.4138232434762701 0.21750625079228567 -3.3157075538206895;-1.062695934421112 -1.4022467096230695 3.4373966810950094 1.6651341977021705 -3.1814103451679974 -1.064439764133873 2.6837693979305945 1.7174162540569815;-0.90539685581503715 2.2260706520678943 0.33598874373793969 -0.85926683675752924 -2.5164239222053517 2.1539434247647908 0.71851235983905126 -2.1012620639230426;-0.73901487765290208 0.33499781799563372 0.089305488533725774 -0.32565181180321429 1.0758278866929558 -0.29631327893687504 -0.30508139621502728 -0.044255211426662971;0.63146783526864214 -0.8478077491646161 -0.15096108578951922 -0.25881608288572128 0.52371470503935025 1.1561480149832655 -3.4927202041524135 -1.4477843897747387;0.25113412920733669 3.0511606263678601 2.0505205199946239 0.13358750927305157 -2.4764157794852548 -0.43603443221447935 -1.9742420512467926 0.90855676142130837;1.6505254540447787 -0.31862009122042778 1.227218118446523 0.43596790775373401 -1.3732159444832872 -0.73883287871595038 0.48606578464288847 1.8452457673966884;0.63001808051021568 -0.12804871513512828 1.5267997471200307 -0.29814483256714563 -4.6380016408376443 0.41708712320390434 0.20908786812827712 0.81592731259590701;-0.4648916892124213 -2.4221046347302857 0.75792839153472225 0.5438652121315245 -2.0090922390801849 0.62795204352908973 -2.3883649006441456 -2.3415256605551011;3.2294614324718744 -0.84024397918481741 0.94903661346734369 -0.36097331136639155 -4.7319702003706121 3.1067519717447181 0.96489640612296101 -1.435765552419342;-0.29154047528462079 -0.6168224863004238 0.35725530220421492 -0.26019157193939374 -2.2767818999622533 -0.37176395432713261 0.045766001896042206 -0.053871578827066072;0.46572605570552911 1.8417122857986945 -0.18258682186877356 -0.98411511479395652 -3.0632381613895387 1.4469553706834508 -3.0583831991409092 0.46859375824528948;-1.3335391069236797 1.34972827444736 -1.0167633050434715 -0.090418778852120185 1.3882452224752238 2.1311483448471944 0.60079060696588227 -2.7198560394256184;-0.32931549076512595 -0.23050618117023181 1.1821661252227318 -1.4585463044655473 5.1234567815689314 0.96641741592603947 1.2319916008121541 -0.10548696357192427;0.11041056968069468 0.51062923867068888 -0.92082097552064368 -0.1896608726161631 -0.84537912512547109 3.4333969522367225 -0.55718610076063435 0.17625386621361239;2.7120911861618948 2.7301981994485018 -0.44746469503611958 0.8600406562351024 -0.14434093100065623 1.9957812283463185 -0.4350528681442869 0.28946347991765059;-2.9926483692033767 -0.51402406358672958 0.26239422040060856 1.6686059186862832 -4.0666192579101939 -0.066089379194044229 1.9607361772472025 0.61302690151174966;0.021491606276338716 0.51919986883370706 -0.38911609376545298 0.29031741323060728 -2.4840820465670603 -0.25278514159895582 -0.89620070384401018 0.55738476862350184;-0.8594059838023641 4.1619842360029713 -0.56222424734614052 -0.97029553924451273 4.3317590377208406 -0.63513306601061603 -1.5115293246148713 1.6366499500894383;-0.65411866690626397 0.67452836196038735 1.4169120281780576 -2.0778192122273236 1.8022246647408737 0.43200341780579549 0.65819865551528878 -1.1911555495792103;0.46820257838875684 3.0645752326443918 -1.0592186499462311 0.11078041514069294 1.4406025950241796 -0.86612584258239167 0.54514759214275654 -1.438994258002029;-2.1776302372909444 2.7815470581493114 -0.82866199629158366 -1.8872363023463381 0.77686393993003722 -1.3112684066487306 -2.4804894163121154 -0.28360193217068874;-0.16736126928644238 0.088104991803859778 1.1198085077420683 0.50827980128964767 -3.2458437523282302 -2.691749666400284 -0.36182711787101263 -2.8536433129093055;0.58783762141842122 -0.56418836959099883 0.4763997522532446 -0.4790749090068831 -2.7802943218661649 -2.5523522853045022 -1.5467946626096083 -0.86640656387662462;0.37890254861352968 0.81429317125764555 -0.027563258361466031 -1.259787860940339 -0.29712546201141227 -0.47534157784123077 -1.7651971298858995 4.4788104820678702;2.1668810380101933 0.35460707486918652 -2.3115174923702351 0.66673452496738583 -4.1634822324937417 6.2467431078469424 -1.9351003116453052 1.8732928290169448;-1.5204832590493451 1.9690903022250228 0.086710062374102406 0.44458148428817645 0.7440222529188093 1.6741487658260903 -0.60762781573459346 2.4506515303231651;2.4473185084209383 3.0702522533687548 -2.3089030172722569 0.60607482407110891 1.7588122805403543 1.7670973182100571 1.8789432686556311 2.2195902078008927;-2.8773963191773237 0.38916497077922957 -1.4385830451524322 0.59922147723991837 -1.9805407137136413 3.4922543477643471 3.0610381474516193 -5.9795354780808561;0.4365945717696964 0.20868134346402498 0.45511712147888744 -0.24414742824432861 0.77210842193723561 2.0086296813013038 1.2230013984433392 -0.45598173350656934];

% Layer 2
b2 = 0.5821210466589456;
LW2_1 = [0.018522405654705029 0.094875773647925224 -0.89278282308318235 -0.013556233972230113 -0.1260591853270899 -0.0094138048696519958 0.042724745631945768 1.1728323726168139 0.010589974782766895 -0.033999614194065816 -0.013833793271325681 0.0070496121957619503 -0.021731305228557332 -0.15630377134669057 -0.013114716933194411 0.0005517333791321302 -0.042076372773349628 -0.041693680572227484 0.010403239173565683 -0.0033201906097816629 -0.3686627066377845 -0.010463254131241963 0.030403152144161259 0.018437720210058461 0.041978856578879314 -0.0024432900014424246 0.0033082459143718764 0.22183062601627027 -0.01990806482790982 0.025456661056049754 0.0031066083758222699 -0.017601499303928902 0.027301218183809855 0.02923185738150258 -0.024030212064447345 0.022121746502672821 0.32924989282686817 -0.0092726539116743754 0.91508885795066308 -0.19981154447739302];

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
