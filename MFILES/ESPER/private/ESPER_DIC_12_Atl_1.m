function [Y,Xf,Af] = ESPER_DIC_12_Atl_1(X,~,~)
%ESPER_DIC_12_ATL_1 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:14.
% 
% [Y] = ESPER_DIC_12_Atl_1(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999999390765;-40;0;9.1619;-0.12];
x1_step1.gain = [1.0000000017546;1.00004398834998;0.0153846153846154;0.000309885342423303;0.0656325902039552;0.0145655815308426];
x1_step1.ymin = -1;

% Layer 1
b1 = [0.84798203003128069799;-2.2439987248762420791;-11.514249598135727126;11.797870681351321664;-2.5330179250085027576;6.1078183826831518033;-3.5573413285870083556;3.2346543422997253714;5.5489509183444063467;3.1259043164192905628;-6.9741070083650189559;-2.1535986265583493271;-5.1165145419016830886;6.4669105761025633683;-3.1643829748784124689;-5.4946796119382383949;-11.272566475955583698;-2.7474936154920799325;-0.62795064542386391793;-1.7512328027530361663;2.5269494038145277059;3.228886837297389345;6.3504089783488968379;2.4313239707777039023;-1.859663508355940742;4.2025464023695517923;8.5485469468488322065;8.4324768917709089777;4.1214094022798537509;0.87997660971048541345;-2.6409040836743504777;10.678972137866443504;-3.4304787859284231466;-20.101523617615626449;-0.96896686530578446384;-3.3835719622715756394;-3.4405673665130578343;-2.8301575129329812341;6.5363848739433381141;1.0567873097033451124];
IW1_1 = [-5.6782515476482489092 0.76254987214155212083 -0.86151859486727566662 -0.55838496083471300668 1.9224403584488085617 1.2478916890179168053;-0.48651538019823670211 -2.4363860445412641731 6.7785681013963774433 0.087678338342613293399 -8.3609590916125444693 -3.5243504487594403685;-0.95067125813575337201 0.58350947287550247733 -2.7047931571184236077 -4.0954118205479641546 -2.3050739068147882982 -8.9269901481229645412;0.90226197987383027765 -0.53714966691105181518 2.6020839871643182128 4.3558723254673878955 2.1162593604284460902 8.8726093230189810868;-0.22041015905050628798 0.93613827221776335641 1.7063397457813023195 -1.0672472557491656797 4.188641376226645896 -0.021885522755234934622;-2.2014062505623774868 1.4699364170061084423 2.9897624273726375321 0.14520052681105230508 -7.644164100158006292 -2.9372627743276611412;1.420353298383455698 0.37080731659928262145 -0.35930445113815984159 0.24319649659375905393 2.2021743254471046747 -0.29951324819728980531;-0.80466120558284526254 1.4976915798187167006 0.23962530854684968573 -0.536079472275206248 -5.3371265150656297394 -2.0670383444639939263;-2.8173140296152410222 -0.34379143025263186884 -0.38946893008453947749 -2.5477414210209494705 0.28315731612356864755 3.5103245326190863906;0.65409067864085224109 -1.2103805371965874471 -2.9343370050855113185 -0.28258527959857554501 -0.22202845034849658101 1.2106354323037773479;0.15922192019920047845 -0.15971891389030584696 0.4820705704708655448 -3.5592132731672014856 1.4999825324981117536 -2.4774269548166079069;7.2604045922186877249 -1.9645791751967545125 0.78106787389283249823 0.3982555874469578816 -0.79805081059893057116 -2.908470401641578551;0.22165142475008059453 -0.061740523378751008265 0.4253038876854736694 -2.8529578268962709764 1.7152679669208974467 -0.72789972535657054831;-1.3038139061034537125 -2.0707489112999009073 0.75465347980456198851 -0.80918368978395327495 -6.5817778049382056338 1.3605488049249692128;1.3130640457449473946 0.2736488317936173198 -0.20649945065898064889 0.1798440400241740067 2.0693336256452994348 -0.41992277940882472009;2.098979925558093651 -1.1618381112154685386 -2.9678799290467265948 -0.12300262922759618023 7.1656482023072918253 2.8577960801282054071;-1.3110963855916357712 1.2911198267477157486 -4.94861804385526316 -1.4043783120989374869 -1.6606135784925353249 -9.5787145890764051614;0.23736709560873150981 -1.240826470577173879 2.6638434112601241388 -0.19358698606065846692 0.85210105480154008095 -0.75621639763615799978;-2.9665077039405023385 0.034108882555484422439 -1.0285618541710621621 -0.35679756976230536658 1.1106047069346858258 1.9943960837059937496;0.15225638181653333714 -1.7218025489391268756 2.8171804364148358424 -0.18245211826285692713 -0.97891702013988013409 -1.1030961085059418636;-0.21106976057815288184 1.32649757580740979 -2.6823551371015721401 0.17299408127972551452 -0.50721438093580872852 0.79442890866694404473;5.006511104799163725 2.5320381146230759306 -0.25097423703156723773 0.80300133597340872527 -4.7444673418723430913 -3.2108087321575755624;5.5437786154860599908 4.0609057801521499798 -2.7904747853203364372 0.64595294386632573502 -1.911879613648195253 0.44703866047516122029;3.4946213958300238467 0.16746819000659504684 -3.2689644221193594475 0.9800775122085234603 -6.557062764152038703 -2.8010045687583868812;-0.83393575051568569201 1.3802732026539654697 -0.27811886095022886334 -0.39029540346251689886 1.615014650575439914 -3.1827508811699707358;-0.55385753785683000672 -1.0324449260478583135 0.10459452056133934639 0.7302701785103145049 -1.5010483063162072259 3.0369128784086210082;2.6279365107890653164 -0.8860646655404120775 0.81534160919071307916 2.9066532324737544002 -3.2149339710497173428 6.2347451203430699351;2.3387595425084999334 -3.1514461815700842706 0.79869021441810406081 2.965669504730453987 -3.3583232160756515405 8.6807226965124115736;2.183145680661572996 4.146524440441635484 -1.7231974216132703859 0.20133389647966840053 -4.5912766448547746023 0.30093874264807923824;3.1484868091706208482 0.0023719157538832080778 1.0495142546798015282 0.41300387826122220147 -1.4994252099108180687 -2.1436290259182220197;-0.67287508330229572895 1.7190568313206198159 1.1072065431567188476 0.55569856579224219395 1.88441292557114215 -1.5402688211396011742;1.5667388589540822696 -0.77386204219460152931 0.99339300094050742018 -0.46841178248421172503 -16.413270879918965051 -0.17093501065999380795;-3.3837245472915160605 1.9596747674222230184 1.3813119881967579339 -0.2482724214142890129 -4.8644217526372184324 -3.7483215112804844082;0.19883142551387186781 -0.16352535948122354226 0.011391806153203712718 -15.416521308980733806 -2.1349185296564545666 -4.8588444359857385635;0.95028958934740437314 -0.83295778229602057685 -1.7263857978153087647 -0.037990934736473400135 -2.5699758286988267919 -4.1540458687134540483;-2.3610385256653656505 0.90310402914063925284 2.7260030316493337743 2.6420741885960650208 1.7339690896972483891 1.4695337632926179783;1.1180063417430470274 0.23132350545567845201 -0.049030145616631701233 0.099544466957747296254 2.5523958287790939892 -0.6787785201556743031;-0.6637232890684718889 1.3151259872691349884 2.5940647426714535939 0.31784184393396752721 0.14387841427595102206 -1.3242304514353993028;2.9375957432525976465 5.4314363991477581095 0.17996997518518706438 0.06748238961899118149 3.273514917065583063 3.9725349098037425755;0.94895561219559354527 -1.9154661598622066965 0.25224044431238490382 0.040062307116792490735 2.239878099498314068 0.84606977590342724582];

% Layer 2
b2 = -0.25630373941511258584;
LW2_1 = [-0.059770736011896648254 0.046506161019057364048 3.3608590525343662314 3.1950799473233884029 0.18693847378032205087 -1.3412253128453945905 -2.222946566361916787 -0.27471515260710172024 0.057104044848680174384 -3.0147915442134496367 -0.29532959363834865307 -0.04276216280146449833 0.46710465672367884604 -0.15796553006959654364 3.5259202937279692236 -1.4229932016549986518 -0.17107652338359671806 -3.6702890454862666658 1.5662615997176465577 -0.96336175322818262678 -4.4929150271605919897 -0.069343074418355341026 0.059285691291991078622 0.051196721801828153076 0.10480805357928581201 -0.22758521654283878077 -0.025443555501656448664 0.027363986946088044577 -0.024730644603197494968 1.4886798796844813442 0.62132408058112975624 0.053075182409580728971 0.040325364310902923892 -0.46123954203621664316 0.030942301137350913809 0.024964396322684410251 -2.16108678223301931 -3.6181346517902408522 -0.088566808997544529958 0.25347415967061170949];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.00115315130997576;
y1_step1.xoffset = 678.286531932822;

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
