function [Y,Xf,Af] = ESPER_TA_1_Atl_1(X,~,~)
%ESPER_TA_1_ATL_1 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:29:59.
% 
% [Y] = ESPER_TA_1_Atl_1(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timesteps
%   Each X{1,ts} = 9xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 1xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-0.999999996490807;-0.999999999390765;-40;0;11.29;-1.99149636107278;-0.2178;-143.241523403244;-0.12];
x1_step1.gain = [1.0000000017546;1.00031699124092;0.0153846153846154;0.000310077519379845;0.0705602618858232;0.0604107502911417;0.0479365703301392;0.00458552991639557;0.0145655815308426];
x1_step1.ymin = -1;

% Layer 1
b1 = [-2.2647198201873766088;-2.9694056797807384385;4.6561882156268126209;-7.1188486643585866531;-2.4530563550216109014;-0.37451032038195497531;-0.2348614456577598486;-1.4099766643914053432;-3.0850981925593292665;-5.5158011440825474025;-0.37813741364225256225;-1.3152534746764501516;0.59748161765687701763;-0.11420522354794709896;1.1722049753302743458;0.21262155013155523142;-0.28303996724123758888;3.5127216779569692129;-1.4775092052788307839;1.0508866007981940704;-4.1879624981155938457;-1.031086761328531054;0.27715889117295017252;-0.71796509678788444919;3.151646438592178967;2.0231019276424691;-0.95854091812624386826;-0.8977401250201892946;1.5817269965608458104;-3.2426298555933201762;6.5979925758857280726;-0.90527994340597661704;-1.702732089873073118;-0.55404121091219815654;-2.8958506144425726347;5.4360247590102899196;1.793775011998949287;-2.5888064713897738223;1.9080653163166527708;-0.89092962136412023888];
IW1_1 = [0.37342093937898174216 -0.7515574341019077842 0.78192266782337571662 -0.11728579280079490244 0.15656402478667558187 1.0604992956600680643 -0.20680073024481931965 -0.20824121423627189875 -1.466044006820303558;0.18578211526329507786 -0.2837160309766019739 1.8154531140835648007 0.27931707438739777905 1.3263112923610718141 -1.1925007990736125052 -0.29636496320458161069 1.0562432080559085446 0.34512333809331868784;-0.87804089950349617499 3.1562957925510732871 -3.2146827451355077798 -0.2921534551824291337 -2.3499627992874563276 0.22617256224461412817 0.19403685572605439136 -1.2300813476904353116 -0.30606666835520762726;6.2120960349492237995 5.6011798067554101621 1.9796229602285599913 0.38017884491957992799 -0.404234428334025786 4.0005518604891161516 -1.4681567421401504614 1.9689324372338501412 1.9436695109043693463;1.1139083116031822662 -3.4092760188339621408 0.8541904018829208356 0.36675916926565815812 0.95158755316619758613 0.8867270243570450905 0.12149896753143053785 0.78022958846722589854 0.62804134004365297361;-0.36233750914541706933 -0.6968296567538120323 0.26447472935512661429 -0.14130579644491406133 -1.6066848700635529834 0.39134759982407257839 -0.14118554805230648652 0.20571367565445586312 0.24401792147232703734;-0.98714377494805283231 -0.35756923716086647103 0.3633758586456631301 -0.8091267910003128172 3.8837267720166943619 0.25985305572633948801 -0.23935968612672167688 -0.026107299579121980015 0.076882667504614005116;0.13249552156299063133 1.510011455538585734 0.79938086062288327138 -0.10341024179717911435 2.1103489036840263893 -0.18779208532890642736 -0.52507354356726598787 0.27406125846093115639 0.39700951239757087396;1.2186135567290148352 -0.353244684667479103 -3.001139349971841952 -0.097582048101805721196 2.3283509724285509357 0.13591223988179573778 -0.14619606706902710247 0.10602516059702292439 -0.82062269493588346769;4.9853086300596665126 5.3440232968950311943 1.7898346490522738517 0.42809762350557567512 0.01026617951540641227 3.998406567637705411 -1.4650732794693470673 1.922559540419604307 2.6642922657031768985;-0.089701688331618262273 -1.0455172646238897194 0.94850064211314377705 0.1554364396590302011 0.7040472823565996352 -0.082936334581681606104 -0.19782995609971765205 0.34631639548375658944 0.11769647176258556109;1.9611290711991156144 1.0494146077498847092 0.33518477815598690661 0.98649211397035596605 0.40811827325333738248 -0.050086078450100157688 -1.0878257678304203981 0.421300136140165149 1.0369049960629099782;1.313772894850601336 0.75149566917929666587 -1.7848439531986146633 0.29818655006536842755 1.4655326819653513581 0.14153202572344353305 0.067736510132623994895 0.1016434510851516454 -0.97071744868354137292;0.13068544461629988396 -0.90310240095407834371 -0.33949902283505628553 -0.22804093821877802206 -2.3899846724867495062 -1.7925242551188498297 -0.2691899281657164078 -0.4851055692050463275 0.24368469863916064622;0.80024010166433856561 0.97765067876092370813 0.3445298083470330619 -0.35393106124212569963 -1.1408778701495123631 3.089447767173044479 0.4891467137829844658 1.2348746783646085134 -2.1077043516790210376;-0.2790573914374314124 -1.1700482464061996435 1.0118926154451610699 0.18948397555757806598 0.16106905072688243452 0.04251368286524997292 -0.17427284112723903653 0.37015571545851000401 0.17862320731995831591;-0.051487171722101189586 2.1727665167348191311 -0.32237003536117359026 0.56166061497545372028 0.72261950287119525171 0.47011566653567471441 0.39560812323907906807 0.82153719143761105848 -0.81135440985728701779;4.3343500347164143349 2.7420870044999205994 0.67610541836670001725 0.39399079328906133624 2.086774081314400231 -1.1959538026898743901 0.57183840901284344049 -0.38735110455090637771 -1.5703928910809830111;0.4572415305578780842 0.32672548113772992284 -0.70797766478970769466 0.38742977708402631665 3.7352487431034866461 -0.030365607011187074893 -0.019494808765154292979 0.23816220136996735235 -0.12482928428502620966;0.35300894536994553707 -0.50932096900065282341 -0.071974115725058240844 -0.071930066450495114272 -0.81610513444613674849 -0.53555961138587315773 0.011058117170701231916 -0.093193107116112061661 -0.18000317362287884815;1.6030995354982824352 -1.0088708028541830686 -3.0067571622882134363 -0.043385045794770024652 2.878389140560651871 0.21049382847918601813 -0.2101233668698007373 0.39568334351135486804 -0.83186201165251683065;-1.8306184920627099899 -0.59495495697749012365 -0.57201698471697781656 -0.67944378944206185178 -0.53808736877266305498 -0.020265500165979409564 -0.32134799190018614734 0.49624075605171996672 2.8292300215224481796;0.29456233361404787763 0.22544541471433024449 -0.52573210177187523762 -0.78237968161234716646 1.1568747356582655961 0.32830940872058267388 0.60793362337768153125 -1.4139883045985361942 -0.95840301280815287654;1.1392022696381083335 -1.3928477477638074866 -0.87733184852754120264 -1.0095016903956215426 -2.7887938321439844636 -1.6701637680109608919 0.75882204769522365595 -1.43416257257916957 -0.14724597061642932805;1.2556310604163325451 -0.06898160773245400057 -1.7555218132704009104 -0.13329551165428307757 -1.0803839215459998258 2.9221634365828630386 0.082237294694242429882 0.64694036715801950255 -0.18483648235590816045;0.99034909978566976196 1.1420853744928598505 0.5420904873694969428 -0.55895844212204326418 -3.0302609392979649172 2.7757892532033956456 0.4871776954956837935 1.0346671702479197652 -2.2312682454074561633;1.4917652936762721971 0.72404150519375543027 -0.68598781535077368865 0.47295215805126195052 3.1568458174734987942 0.2744465847046024809 -0.074213873297478472058 0.9653780949292356528 -0.59016029991019369749;-0.56004044731457824913 -0.11608375587105165816 -0.44620248188928662847 -0.69595444394414796019 -0.15254108795618298577 -2.6967875628668651622 -0.04644335981416430198 -1.3626043679471508963 0.43908520655629329488;-1.3839934912792406507 -0.53145325083721495218 0.33364565066312146246 -1.1944598477273891124 1.7750097193393516992 0.15761350412119043418 -0.38792816489668702218 0.42913032038202753027 -0.063846067854810062414;-2.9818345229808715935 1.2190376055811800171 -1.5601355211926331013 -1.0999038808320009686 -0.066149910545606149803 1.0544315057288027315 -0.092173162200882133299 0.018839130456424466736 0.22068104071213601181;0.4244827023519175957 0.88582600343065043269 0.70368579500133954241 1.4670458029191077998 2.5559662026497247389 1.3729836870614620103 -0.3469213573608594392 0.40723337252389807395 2.8167738367492454188;-0.35307761061949050063 -0.40125424133381365488 -0.104779184049259802 -0.16409721829147547179 -0.95634833195968016462 0.019819143343671973695 -0.061907645398103700041 -0.08229627647878172414 0.025296913223628741257;-1.8191237619750526555 1.0039343109439458068 0.28321078779291047889 -0.35554641555445348766 -0.68586494113601836897 0.44334601878585028212 0.36483496522331337886 0.28707188100253494545 0.060104366550625662469;-1.0187212633991609145 -0.80636465677185675816 2.5845056095703036192 -0.21063729689130086808 -1.8910246975118085722 0.10030359263893759425 -0.034135992913941572147 0.083335911356124850458 1.1632769467185022361;-1.7264497730479264614 -0.3804485216074774967 -0.76787203324446196007 -0.3893996990484873133 0.49650166847541399306 -1.8434178016781452847 -0.33369156487866774574 0.35448071777358913792 2.0369169792148826836;-2.040462927431640594 1.7811252318938131101 2.9165314319440946988 -0.0076698849384695772999 -3.5665842561724954152 -0.29926148140885433513 0.28361809752770800186 -0.73481839041219942921 0.75165112122576838161;2.6860418897151645545 -1.5788779580712082318 0.30132238199884531227 0.51915279721893958786 -0.95696176870091498312 -1.4489379368708155127 0.54806720955840249054 -0.12957852784679566027 -0.78947144456358653297;-1.1353038500630605601 0.077459105910713574317 1.7677965828011721339 0.15158713715644547837 0.38785774920133453403 -2.6534894557867212406 -0.061914876015082026539 -0.42299725557119854935 0.12131797615301251336;0.90082701082865024311 -0.48111763633894172321 -0.72679591658705888779 0.14667829866259587956 0.22250788221032788416 0.30651497614322953389 -0.35224889984059160053 -0.44747095591663155201 -0.1065061495478320519;-0.052279934224165831858 -0.33033298982343711758 0.7048599030478420735 -0.32062597899204048346 0.085528428262788211422 0.83000680193850806532 -0.30174225534457921327 0.081276433595518310371 -0.41720182916929654526];

% Layer 2
b2 = 0.24142216338106570195;
LW2_1 = [0.14909764897560628949 0.17751488011393565247 -0.090894988189561765068 1.7414917429151361983 0.097952782432132026735 -1.0619229332387583042 -0.82002635115899713547 0.10516108679073238319 2.2797080691990387891 -1.7953216269219109957 1.3882202890844543841 -0.045252469254146186761 -0.68785324521976809464 -0.23780063923176378293 -0.10776861522635372559 -1.2758417878754977348 -0.051166255130031111642 -0.35576903055470937565 -0.46204644052707255497 -1.1771736155851062477 -3.5504317499878883524 -0.14416576629822644429 0.041259305752846954507 -0.036379678631571266134 -0.95117863586478312854 0.076877341884281302775 0.29682978571312801375 0.099183279315599254256 0.50237720841895239232 -0.77707192280888537361 1.9167682684920772651 2.3012861683088132736 0.82167600687981012797 -0.4940261123401430976 0.217865933678293211 -1.461188033983269019 0.19593116164303922555 -1.1215095464639426304 1.0432028883551625054 -0.3596477479289235113];

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