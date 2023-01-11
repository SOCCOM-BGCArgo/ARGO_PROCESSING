function [Y,Xf,Af] = ESPER_oxygen_6_Other_1(X,~,~)
%ESPER_OXYGEN_6_OTHER_1 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:46.
% 
% [Y] = ESPER_oxygen_6_Other_1(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999998793562952;-77.961;0;29.715;-2.18999082656054;-0.04979];
x1_step1.gain = [1.00003298043486;1.00000060321889;0.0138606862425759;0.000266418009857466;0.193068829037552;0.0582073296585987;0.568038194888224];
x1_step1.ymin = -1;

% Layer 1
b1 = [-0.55052960362305380926;-1.4321780831682702217;-14.378248709875991906;2.8691821908684218556;1.785325865247140209;-7.7455541041001785274;0.55263717857208127793;0.50601346308488337478;2.273562237028825983;-4.1941745706160951102;-3.6474435833672678164;-0.83811822559852866554;2.7301003107712071838;-0.56771389760389312507;-1.2879232792330670421;3.3592031756636258066;0.0077126777907828587139;1.9561182948547095606;-4.1330996530672887701;0.38118434076870810756;-0.10078827880981318377;-0.17843861165734495322;0.77634166861440456753;-0.79719386148782334089;-0.099540318269005376162;0.99794371630848244248;1.2443119666474899976;-0.91282709697478325861;-2.8442218176338394109;0.5412967377617315945;0.88422220775248305813;29.990666358372337186;0.099178446378909476677;3.3193807094996645723;-2.100463866489281628;1.3203900512893589347;-0.47966368832241007558;-3.7249316564177954447;-0.36915546194575160621;3.5542985674988711864];
IW1_1 = [0.10366157438805685687 -0.010856861120604265059 0.40688002156930852538 -2.1653099898981693627 1.1562299408942640078 -0.49085454907240327893 -2.8326543108096515766;-0.61885822659413003421 0.20807536635272566583 0.42157012157432899313 -0.053430252972567586489 1.3428393624000547479 -1.0593963952058667033 0.21982302084041108192;-4.5811107188588060524 1.4516955497214416582 -0.12238146421280816922 -13.282925576076552332 10.706717287601067667 3.5192800434434614587 2.7863733562363814933;1.4560604052313583789 1.2408038836135624194 -1.0090614780703364151 -0.85911240918792863841 3.7687326174872119999 0.58662061420821476343 0.87831419845705382965;0.016158479382050682654 -0.5734329782370870543 0.0073671701137597422199 1.6377206853338501791 1.0643721517690494238 -0.85623594181453521745 0.87847438244076514469;0.028293274576422014488 -1.6604895248708695377 1.3832646332628333319 -7.042129369618201018 -0.25959420561594520604 3.5326721790483879282 0.63563627576048498113;-0.10307345805711741471 -0.13444745143725547742 -0.23393578907863757377 0.10942697735713023677 -3.2579338882972508351 0.55511908084443295497 0.0027842397668792402167;-0.73015926295234878651 -0.34759877121377524833 -0.48789493341289780215 0.30314688765665936954 0.86607477624739059863 0.21478754340203956574 1.4597903770613902541;-14.476416734319938584 -4.3947728750620589366 -2.4978359020344074182 2.5785643617494922175 -3.0362155391015557271 2.0401951014491577219 3.3842535464777663279;2.0360450550189530539 -2.2401764926455900451 0.47949230586131946774 -0.79413541297291379006 6.0260306578646227749 -2.3706617638902476308 -0.87736083648847584815;4.1859332855083124514 1.0851883015615737538 1.6279224413511603764 -2.3814315714040170668 3.6759628257894401138 -1.9034788019738924092 -1.7733895407928019772;0.59639467612924845774 -0.26093143598815826678 -0.69546447528208199529 -1.5027881663415953994 -4.238719069070611134 -1.4160605097504899152 -0.32085050697044165124;-0.31907506744929226006 -0.17818100693822891412 0.54582918868477314778 -0.24477087818865139734 -1.7978575885224585829 2.7160454084255505336 0.61108828537476156395;0.45616800897458698039 -3.4198189212361742584 0.29957953655837871176 -1.345178330202557504 -1.1385548587975018453 2.3369980560150094284 -1.1017476551239873128;-0.50270676562012062139 0.14460532325660663266 0.26591864645633195297 -0.47753407130680147752 1.6576402420549383177 2.5880568823775287335 2.9994385437302568675;-0.9911911428122176515 0.60883530602706414925 0.76855180671681522231 1.160158280265203512 1.35075904603794128 -1.14355889998559479 1.754325640730717728;-0.34545893884632067339 -0.27566812609350505037 -6.0878678358114823865 0.50031090578278403225 9.2959140706966483236 -3.8436741257020017315 1.2375110476837094531;-0.17604007101514773659 -0.10857515999749932367 -0.94399514250090943435 1.3164496352425394399 -0.43278907155719398281 1.7960581562833062286 0.061167550429120649824;0.35931098657169774357 -4.1206558109720088368 -2.3427914601081893942 -4.2745357160787094131 10.708724615631787458 0.67753555488242822591 0.68290940555283652369;0.42717843428756907542 0.11059366984571700199 0.08243601708249233273 0.59258046839541100859 1.9463627646663750337 -0.77472835855765309621 -0.79415547569282374241;0.25182645604811332296 -0.57696131548292772795 0.38565627376662192649 0.39683573831747481764 -1.0036118016501289762 0.11701072653319517791 0.40004901203412429611;-0.32609002324572272657 -0.11211091353068562693 -0.13886544748079451828 -0.37366722274699459838 -2.0990147477325593606 0.61137887782230149192 0.46106145868999720561;0.36919333595622644273 0.24093304092553155527 -0.67408696124022771112 -0.21182355049626130916 4.0108433390732374235 0.73171049928146270069 0.89788542461518194848;-0.017697782047047423692 0.12585287858564853103 0.27745929063524737979 -0.3945120468240952416 3.7469421196698102783 -0.71980701907450295174 0.12028585881580079742;-1.0992558874106772571 0.067043490100690203426 -2.9145920744291373694 0.021290359055116180736 -1.8875180162079268076 0.37930215046853760219 -0.75665922342408020729;0.43204208399478422065 -0.17147122447625262609 -0.70184104481534936859 -0.60261980095707989857 0.03919324829359005502 0.83759285048264398021 -0.40860749580064564812;-0.055026838527206460572 0.054782222095036144094 -0.12425253404748462516 -1.0462857557194293889 -4.7186506180266984956 2.5063826982262051857 -0.42551885036801145867;0.1041338567314948893 -0.0060685165026279899378 -0.091075023443426716963 -0.098085186978370403343 -3.1111237817181973675 -1.2255201962986301822 -1.6403835198982559263;0.26701681336034144287 0.15448694956691622071 -0.68055176719144694353 0.24532637760738304489 1.6525534940892909752 -2.8022892080756229838 -0.56382812797914338976;-0.94255482199309448266 -0.91437273597563673011 1.1678781278352001749 -1.3535050409883675471 5.2915131519613423094 -0.0017870210336896766648 0.74250372863869773532;0.65564338680045397556 0.043598528304250046272 1.2643579544941025983 0.73295912528249473894 1.7246828607345594531 0.92892838445366276989 0.68390024058789244954;-0.026173102570099512909 0.33796052769402590288 -0.37358489264822208886 29.805128838570489336 0.40247561598357489698 -0.65256215932588146345 1.2540314567980253013;-0.24975698236845944589 0.58095674040010114592 -0.36398930954633523793 -0.44503284312311930382 0.94195152743541532558 -0.14240441435349493471 -0.48859345886204386122;0.15326565380804194061 0.71897731678582377235 4.3987695582617885037 1.4048623362435326278 -1.0715870766035664463 -2.6726409042233205682 0.97086412733749538617;-0.40353470702995952735 -0.61631326409061371052 0.22952483483043870094 0.23773293094619843413 1.4809427043629719289 -1.2521982941713907245 -0.22770399283011402702;-1.4148458216355208616 2.1884052697204889881 1.6765493976927270392 1.3881562390562152132 1.2324731719213599668 -3.2362308530283310226 2.9152829104454665021;0.79954784785326027396 0.39051627062945359903 0.55866015243221744946 -0.27313369785757191099 -1.0892333979344770434 -0.19943511786793585716 -1.5000626790704401081;-0.3615140293815146677 -4.842209323575359825 -2.1927243298442298247 -6.4743689945679330222 5.556882803103412094 0.10952685281966524633 -5.0618290398295613386;-0.50726476269872955438 -0.27595184166890179833 -0.34717087628521187126 -0.85181686658852895366 7.0447152954890812637 -0.51585062370718670621 1.6936575553657891025;0.17749979721351283324 0.052405875565030804464 0.70939050946385762142 -0.27811091982050828442 -7.4937422629672250096 2.7752903349559430168 -2.3745471559425643804];

% Layer 2
b2 = 3.2890830000002249101;
LW2_1 = [-0.19681674529107875959 2.2091730592916878351 -0.025679970088693420543 0.44360740971314321923 -0.48765099998881750176 0.10188438756017956233 -7.9079144945873194317 4.5057125865186300473 0.016033163132924223293 0.15545924931642507438 -0.027305943724738202666 -0.26613258158690844546 -6.1114609359065150684 0.047622708823378302745 -0.28189760727711887789 0.34011242941702563014 0.065340015792523922777 -0.47809892161764250273 -0.0227017632820424968 3.0376781130857315461 7.0456104866823094923 5.0367412996177556295 0.5077202230900792479 -5.8947446757973551712 0.14639193975095074474 1.1216832711117445953 0.90208945076553614939 0.93065068596812450252 -5.3794015666106513862 -0.11232159823848129887 0.29346676532713561469 -0.092822003483966589177 6.7338368871317122455 -0.079574696882464934777 1.9614631338679113792 0.04342674708176351428 3.8070900038795327802 0.015484639414615742475 0.24523414179743541208 -0.64714595858110812721];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.00454442172233583;
y1_step1.xoffset = 0.3;

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
