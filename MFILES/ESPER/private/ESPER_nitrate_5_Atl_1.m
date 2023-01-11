function [Y,Xf,Af] = ESPER_nitrate_5_Atl_1(X,~,~)
%ESPER_NITRATE_5_ATL_1 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:33.
% 
% [Y] = ESPER_nitrate_5_Atl_1(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999994516886;-40;0;9.1619;-1.99149636107278;-0.03;-177.232103709389];
x1_step1.gain = [1.0000000017546;1.00003626070655;0.0153846153846154;0.000309885342423303;0.0656325902039552;0.0604107502911417;0.604229607250755;0.00425400504498441];
x1_step1.ymin = -1;

% Layer 1
b1 = [-2.9418172723185751494;0.31690158141255314606;-3.2275712183306470848;-1.2593113924539047765;1.5931937314911177594;-1.4079824467714996405;-1.7347087452619334158;-3.6142524829577706846;-0.75533834081998718624;0.56411658388008267817;0.90215736514268207724;1.4941404040603005754;1.4752451899557075254;-0.78983021506624750963;2.5249438662612346818;-1.9185522382961095467;-1.5728149543955682876;-0.97867265399800507097;1.2638745531881021567;-0.55014204785826570365;-0.48748162527815702694;-2.0708509301084054854;0.38855474901805037158;-0.88925090968780595713;-0.19946400819149845418;0.029266242956865626268;3.5915421019717888029;-0.17685042316119056749;-4.2011228422666793136;-0.44093087991417140303;1.1491692275738727425;-3.919917722173194008;1.7834586768799682854;-2.1715720358306445092;1.4264502169248711727;-2.1177164763104650014;5.6983419570384352681;3.9497572626590704559;-1.788782675915271847;-0.67229460281113850861];
IW1_1 = [0.24000279509297808844 0.24686520012912008259 0.35351777748179025185 -2.0856945633057177858 0.17997924290564890071 0.31079745781156231654 1.0134231147410088791 -0.77204671752025411102;0.69318275979893151373 1.2339654694344319719 -1.008949396516065411 0.92720607777288477358 4.1067092921309544806 0.45017473869048796331 0.78037021285890539879 1.0845770184752305809;0.49778811148075796655 0.63363032953984521622 -0.34612874764545181749 1.16424004178500895 -2.7049605178079065837 2.9337292240658867648 -4.7625852685873617531 -2.8538317422939885049;0.1776160402324866161 -0.031605002896580455807 -0.90597809912543914024 0.62335971732954986724 -1.3084345587388934717 -1.1176484795463985655 0.42811657062564673071 0.76386529050717899558;0.073760229214482456173 -0.18241105220223374772 -0.14248338915217714784 1.4173698290240703823 1.9330377526358664131 0.89999227067047626871 1.6845608314095392721 0.084497924165135629426;0.65481844748666107048 -1.012584917961728781 0.91022458609495549009 1.1140150286443302186 -0.32921510535953285581 -1.1148624627814420762 2.6118257585903816853 0.29846545888606107466;0.27759838778349221755 1.0716124939527855719 0.55161765853390964942 -0.0089417732556496352242 0.49392126975215377982 -1.5357101920746003554 -0.80350297366668721466 -0.41640077572607575584;-1.0259748747383519163 2.3920988354606800286 0.2955726278241271876 -1.2097009876147952756 0.90215473458968387099 -4.1940333802535070262 -4.0420989389480936893 2.8687741587428070567;-0.30286084414615582761 0.24779694327507686968 0.88703157545492505864 -0.84634728038449524856 1.630531216179268883 -0.31495343635873918053 -1.3355242259573110175 -0.79891669946752430498;-2.2565386174600488012 1.9667985695970093829 2.4045995555872363703 -0.3223142778634337513 0.14388483075230182329 0.3930241204722514925 1.4987804616573567884 1.7209835756086335756;-0.26833099376792185486 -0.51658896612247884317 -0.77643326765418385627 -0.092633396281314006604 1.0231278471490108384 -0.20134222260268402893 0.83339985824889484256 0.44228601310114068479;-0.45332304731992228986 0.21548988147318415276 -0.93481148142727266137 -0.11203576270152709737 -0.20532289585378737939 -0.43283633981943459013 0.24315389327449588031 0.94737925825190238971;-0.082784652975940870023 -0.45307649853652981475 -1.1251441247216891206 0.10101923910454066458 -1.0356789576152372057 0.12543945853109428334 1.6320402394787472122 -2.3280797630882421601;-0.22825936804753538256 -0.057773634536083286195 0.99440750099951558116 -1.0574955209862890726 -0.78166574440004676916 -0.65917163969076930652 -1.9415681615193289833 -1.0921368287154031851;0.13480113565138290843 0.16383179049815510542 0.036611496742124521875 1.9720011383410727124 1.6822151255491120114 -0.26217297093548180831 1.9024306101720038775 4.3050583642998940803;0.82773896532763791001 0.37939135842984372804 0.55457233755348622761 0.59734655284139015485 -0.060202501414780097699 -2.7675552984298783521 -2.1117260314014671962 1.862584096582972526;0.38228608031602562267 -0.69063975404468680441 0.33362254532631263615 -0.18549082490531509926 0.35129464995317094544 -0.19047046351571722544 -0.069109210266203813244 -0.43440355471906993978;-1.7437716659536866004 -1.5224997701608871115 -1.1295958527553997541 1.0450356075289599644 2.6452605498740724421 -0.044195181256444532325 -0.25303318435790755903 1.0973858181770714815;-0.73151743483176212557 -0.3012750158445371973 0.64367702143749061516 0.16936384792853415027 -4.22912373920469431 0.8487344392095049761 -1.1116776472882223814 -0.16882073466006297657;-0.86900460126805145755 -0.27033065038683223857 0.079731642877540240266 0.19585673866423949918 -0.84700774695127700742 -0.22497693239371843732 0.28063650311157378736 -0.32269507192431129416;-0.99978853423621649732 -1.4077071355857266877 0.62535816976776348319 -1.8037408192772468318 -1.0167229663807686446 -0.64163244081149706233 1.0939546706145821808 -1.5549437248868149819;-0.046426865393402995397 1.126850491410111621 2.9375444332231612776 0.33802086470045611399 3.4946057280788864752 -0.10151954753746721127 0.70505060883734804733 1.3067669304978639655;1.8450153043139072562 -0.12164807131275383323 -0.86613854084929620125 2.5969249960988354253 0.16654238244304239713 -1.7257090407017943168 -3.1001905501971820378 0.05786751692952421533;0.83419713350326452073 0.1691174903430750498 -0.074037032143792633199 2.5764624934441937576 5.0345231603628439743 0.65034614216720498803 -0.51298979371255160054 -0.14106228514471502722;-1.1127891787074515317 -0.077411762575451709134 0.84023260063053262225 -0.83220202533169962411 -1.7219615947061854833 -1.0838485750810149533 2.2856083761658814169 0.0531847744198047645;-0.49706513477793057287 0.43125034040728127405 -0.71190570265978148079 -0.61733587078394569492 -1.6411434463162093422 -1.5998802927105257865 -0.81666507253447895298 0.20764687644469251304;-0.075346953708827296703 0.035519039277560532364 0.01808375009035631531 -0.35359803730802091826 -2.5755219398036062195 1.9852641253232068408 -0.55272239926890409656 -0.087092911195387159862;0.25640641199188352051 -0.21832645347540585723 0.023895622172970192032 0.7205023835004701338 -1.4096462644946359255 -0.029025579792790325712 1.7855651335009137259 -0.27111443691796777156;-1.1083086447560706045 1.7417070838291048851 0.6306056570112114068 -1.0203474305473170247 2.8743970168635701867 -0.14580492079825180252 -2.7943869038077209233 0.3850823704470600628;0.65183699625469104788 1.5681314467695919301 -1.2130339368057625471 1.1293965065588029706 5.2174834289569149703 0.012614418274855302218 0.77863361881960979627 1.0966335317432025498;0.93328049528347878372 1.2831794167615981905 -0.60377668056851940825 1.7015614194465438214 0.010907789306962543258 0.61061100966143122548 -0.89029277415431351539 1.3740600371064124463;0.057384445497183793339 0.30528414182013463574 0.39922007789246116705 -3.5254790614114801528 -1.1507882655507370995 0.40499171466248184981 0.81314105119400559296 0.33414308858090668863;0.29439595263851642271 0.62032551311402284178 -0.77716636056704180291 0.89472158837376369878 1.4494243514898559333 0.60245941392451241025 -0.67550918650087798945 0.53712764780747956905;-0.15712002483780412865 0.13192714870560232066 0.26333365113633694365 0.064508656943800840589 -0.21725130598790839964 -1.1149446376004117365 -1.1461800891043141259 -0.46323071031376916729;1.1888611404261180127 0.20847707534181053601 -0.63891183070883839079 1.0910752025205838756 0.63474082048256108557 1.5821281845421661938 -1.846564214512010027 -0.04599199073736497112;-0.8534374138154848799 -2.1307015898159531453 0.15956347799186040848 -0.95663178908342172413 0.32008715647316188235 -0.0064793545245706695573 0.58818844032159389723 -1.6061855484669500527;0.66363907997339610656 -0.088609000828309247844 -0.20556893775308426875 2.0134728317271841647 -3.5641188696702479355 1.0399129754949705262 -0.20924362394300244361 -1.1760477697665707719;-0.14623858521708738012 -0.32095995432844287576 -0.12619811722793644604 -0.27614295052173343237 -1.3338173173602427823 2.8690779556531200711 0.043664624584986391764 0.2918813769007599368;-1.0537855730227272399 0.87430917945119779322 -1.47243289387356735 -0.43864777409789051932 -0.42516919321113966213 -0.85263534160857179511 -2.1487532621484199957 -0.83828295628095228498;-0.020999214615720566968 -0.085447880078360308187 -0.85012822434224599721 0.63684903489644195673 -0.098246614578575447441 -0.0048531930032516303053 0.41910978359837380847 0.014976006449056399744];

% Layer 2
b2 = 1.8367861213695855227;
LW2_1 = [1.6928078752230524273 -0.35234740212442772211 0.091386792821619683425 1.6808709909600556553 0.78092294141332652124 0.15491475500045845393 0.37453649513549858163 -0.14521885556535316097 -0.74417157339587924891 0.046477299487933838862 -1.2930481173675689721 1.317509040008909782 -0.26740447616654311469 -0.37611960223216139854 0.079990098248885310728 0.16024261766684005925 2.3399549761910467538 -0.099677304350978226943 0.3246889186662104132 1.7012006297008734901 0.81567036075985255295 -0.10071469547750311335 -0.085207183591642934872 0.30019208754248083437 -0.76430081407037997465 0.45262336635842087551 1.8210073626549194703 0.76190572051229865647 -0.10115256505005765231 0.33655418053140390144 1.0456959528714875329 -2.0203028559929285812 0.84604494231561644391 -1.4903292459101487033 -0.76898277305149276284 0.19142738230073555461 -0.71728690492394864275 -1.0433579977529776617 0.1950233504108448146 -1.6415260959796595941];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.0470499670650231;
y1_step1.xoffset = -0.28;

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