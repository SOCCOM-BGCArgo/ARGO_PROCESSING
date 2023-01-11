function [Y,Xf,Af] = ESPER_oxygen_14_Atl_1(X,~,~)
%ESPER_OXYGEN_14_ATL_1 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:48.
% 
% [Y] = ESPER_oxygen_14_Atl_1(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999994516886;-40;0;9.1619;-0.03];
x1_step1.gain = [1.0000000017546;1.00003626070655;0.0153846153846154;0.000309885342423303;0.0656325902039552;0.604229607250755];
x1_step1.ymin = -1;

% Layer 1
b1 = [3.8843792058181345617;7.0024939947386899419;-1.6737257796215458949;8.4228308347092255559;-0.72977719179670119676;0.72937366768888323154;-1.0459034392366177713;4.0453554372902544856;-3.8198686127292482695;1.3301285543848533344;0.26619712474212553488;-0.22809637201177479504;3.1563082708076044547;0.8926080863985271785;6.0146109181248057496;0.66459915778236644535;1.1375758125020394296;-1.6610893637201373796;-0.90354526067298401859;1.6704447959299946191;-3.2024382822070243115;0.38893446424485589308;-4.2932944273298634386;0.6201656545208868998;1.3259921055778987053;-0.27986779288539387345;2.4967831542336398343;1.0315500107069512925;-0.84501005390585937072;1.5919715790713562242;-2.0957743230346506103;5.2276393505670037953;-3.5656202051002439113;3.7180094396107445931;-14.305436584837398328;-3.5457431926336582251;5.9660765652182812246;9.9814127673392114559;-1.0597167402031322236;-6.7644583452520947731];
IW1_1 = [-0.45288297335273558497 -0.22293305051039680187 -0.17414476211950347118 1.7578006588290528978 -2.9775017380781712006 0.0040606366392567216222;-1.4514278723128324344 0.62938785507112671169 1.13556020508372435 -0.60847000284102426004 -9.0373400808103721005 0.094172477022558542803;0.19327257772290049043 -0.78444985930304589328 -1.1841978030357815932 -0.70483141502221080543 0.23101910052129279172 0.41837136240367700513;0.80833615562664040688 0.10744586137727030195 -2.932470353253886941 2.217872118832216799 -6.9398899040454935516 0.94281017658081800459;-0.39652565618824625604 -0.18166216402371632155 -0.9397548273059910251 -0.74652656562042751798 0.025591036068675399734 0.33639004069905742433;1.9560743058637350966 -0.68870489975355753121 -0.98340107542847265965 12.207134634161416642 12.773900084551243594 2.8391344582648372352;0.78250020223983129952 -0.6592263083196786555 -1.9842891983304342673 -1.5320606083033578848 1.4164467226343657735 2.5755260899755607618;0.14870538987457335156 0.26114026725916533733 -0.66781122048326191365 2.0570439099429056462 -2.3277012899307303329 -0.54051684702366420421;-0.21058058131717671824 -0.17208374943551008052 0.60594195357579672123 -2.0966790878320615832 1.9473906153299240795 0.39910092698594473726;1.1537898591282615968 1.038589650222085492 -3.894906081367413897 0.89643705111808635877 -0.38056850575107575851 -1.591742709710731063;0.91563263371964054294 -0.20622598591250701494 -1.4175896116408830405 5.218622641963124309 8.1679654602618132486 0.85296763037709810718;-0.27702964594554657651 0.87841811664318059094 -1.3317702471256567254 0.35872291920689752809 2.3868428087660662484 0.94421822504719710167;0.35907134325199852043 -0.31289111385045681457 1.0506633895456665151 2.3709259621046836486 -3.9100136860146683659 -1.114340311560818364;0.90209077145589111613 -1.7061797478731657307 0.84123748508738238616 2.0002252002524580909 -1.0894990279512435372 0.3619123869459328624;-0.18989115253342733758 0.83262465339176305879 0.092357702983855663459 1.2878472345700959778 -8.262376159130097264 2.7846737749510555027;3.2079606630697177216 -0.88630648726670868687 1.473827007384223009 -0.017980988225710040768 -4.4126056775040343183 -0.045902865185831857875;-0.69238016873998997625 0.10193151566316724366 1.0123534566496190656 -4.1664220464287460288 -8.3729799457753433245 -0.73982884820089511724;-0.62346734825168226291 -0.54389417287394581368 1.9946558349006995314 0.13411254584888140928 0.2527435469149601488 -0.12781505916359436892;-0.85785354471251806796 0.91182119076009249259 0.25876770671617060904 0.97147589743184892264 2.8865322331246394683 2.5328631895216893177;0.75419897991440199014 1.6177750712888867568 0.52818909450794604421 0.68415563588211203427 -0.060099379237479225169 -0.038547573580354009548;1.4143234804101065816 0.19561990634416609458 -3.0190465024447776798 0.73798671617344113738 5.8729687522804452371 -0.89113425181130367925;0.087635484921741377273 -0.73719403682048967319 1.2187385260294651967 -0.038408345760390635859 -2.3063058836708538735 -1.0422309457127907351;0.054582277192338553296 2.6445494552433497226 1.4526332879833181444 -0.43046258732734582475 3.4010544138401050951 -2.1000122164448415951;2.7881528151453465192 -0.042950659405993206474 -0.34917746226659357456 0.62937352970387527851 -3.1051223749995253698 0.51879099292432473423;0.8014307757464691484 -1.505841962629380637 0.65624217972799103737 1.8910768780912965426 -1.6894942050642254294 0.16264279783451648931;0.17209600998960325469 0.97219700738190950595 -1.3045524646215271858 -0.075377833137975608357 -1.4793266804246867352 -0.07489756097314577421;-1.2858360986031771578 -0.27309500011932819064 2.8052033998535947923 -0.71988781285582237324 -4.7351340750658827972 0.85096877906695389537;2.7342238406898320413 -2.894272989024359255 -2.5400844217954339221 1.7477992439974676664 -5.5767547561790271615 0.71820944423058419304;-0.21331687401920668146 -1.2124532548000519494 -1.2282625902219161951 0.655413566112613033 4.0723329549954936368 0.21159539128222779047;1.1325837634555595557 -1.6137415083025226092 -1.268053639910087016 0.1159579221790741077 -0.36371078024507591531 0.43550731982076623661;-0.61361240921287896199 -1.953397823188987914 -0.25610643244990455303 -0.65955710526842048491 0.25210715811197648151 -0.11163128905520472434;-2.0644480952882222624 0.33471553133017661885 2.6583084947679651577 -0.60140863315136439127 -0.95267545520265284331 0.31783348785047110807;1.2723093327542325248 -0.61168820011842239825 -3.4383171150805096339 0.039129108807772365608 9.1834033398688141858 -0.16081214268510632404;3.4055628545312162103 2.4928250403049418971 1.6623680932353068762 1.2981559805633324789 2.6284290756476966777 -0.45301045596106159774;-5.4153771539219919262 -7.906730134577700575 2.0893256331909157986 1.5693929087154439017 8.0163810148091663166 0.38841077678333535417;-3.2811882871357083147 -2.5809724465818040251 -1.3970686345149967433 -0.9501382665304503039 -2.5475234961697106328 0.30257911635439660447;1.0321392324444349509 -0.29488871037533614938 -0.76296143735600052693 1.8692662713598333646 -5.1867435574772491336 -2.2278853808076872767;2.1166637333851787695 -4.2596554572762084234 -5.4738357782362205839 -0.031043621034360570354 -4.0100216669161286021 1.6312983812957526641;-2.349665742965366011 0.50697089886209911747 -1.8195749000360663672 -1.1585496032013069545 -2.0229775674385681583 0.43278392108152002882;-5.3385574286719501558 1.7600951301412148808 8.0881871550282564698 0.48673589479041495798 -2.6821829655759197308 4.6303392163072754073];

% Layer 2
b2 = 1.2509167006804779998;
LW2_1 = [0.59179692541170592079 0.92767277467247943878 3.5226244351814446354 -0.2581580908065181923 -5.6794551128009551633 -0.04371174546989160109 0.29224007102454357954 8.2449837848811551311 9.1931109050337358468 -0.49626507839903655217 -1.0172592791332342177 2.5266171125064214031 -0.27063380268990977484 1.6784560245456976446 0.27112650561339840527 0.23862533703358365655 -1.6499722447650266854 -1.0804908116378248994 -0.40491446010613679762 -2.9933229868423349984 -3.8345211136564500976 2.9483364840575498889 -0.19900735345772221918 -0.34226524763503912707 -2.2287307437388763631 1.4569767870967698098 -4.3191027842175619966 -0.14002698465483140455 0.47225554384998896795 0.99244676303357659553 -2.4177168219691225204 0.42123701517158573804 -0.47660264665902812853 7.4170137947918055588 -0.087656316534258316153 7.6228125145699205945 -0.59539474086574040168 -0.047545368519273614738 1.2387604766322672933 0.055707589905215855464];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.00368391969055075;
y1_step1.xoffset = 13.5;

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
