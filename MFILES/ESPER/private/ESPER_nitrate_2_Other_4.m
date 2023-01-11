function [Y,Xf,Af] = ESPER_nitrate_2_Other_4(X,~,~)
%ESPER_NITRATE_2_OTHER_4 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:32.
% 
% [Y] = ESPER_nitrate_2_Other_4(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999998793562952;-77.961;0;29.715;-2.18999082656054;-0.04979;-2.02];
x1_step1.gain = [1.00003298043486;1.00000060321889;0.0138606862425759;0.000266418009857466;0.193068829037552;0.0582073296585987;0.568038194888224;0.00819772922900357];
x1_step1.ymin = -1;

% Layer 1
b1 = [-8.8903349321058477983;-3.2187048024244493405;0.73390685251559839575;1.5173969680550114525;-1.4519790772945699864;-0.96258007089436903314;-1.267936271412100746;-1.3917754091460861066;-0.70301390038887667799;-1.2221479327142186211;0.70936541663209962127;1.6476414319658765351;-0.9358655002104541154;-0.78319949997472992909;-0.14162795043195547384;-0.0083627195769413707749;0.65535813779512019916;1.262360257970404076;-1.2369961588379052575;0.82152646242415972022;-0.46375826030360695063;-1.5695603332272698438;3.4157617584016515266;4.9899328371870010912;-0.70948813942211952366;-4.1854567403895392985;1.3051695879416433943;3.3446067536524801689;-0.81326695302570450874;3.8233929425026773607];
IW1_1 = [0.056362908346392702474 -0.57608118838045241361 2.4840600787649971259 -3.2964671820843225447 -3.1782049888447030384 1.8952912834369799455 4.3190189709598465484 -2.7563746133400832683;0.097158280779560832197 -0.48014797066600184428 -1.5473579254656613102 -0.097150464601380645546 5.8249723352942943677 3.8003840311943326569 4.245966238961997874 -2.6544622736994671897;-0.75790962553681862079 -0.75212219445215733327 1.7150643676452264064 -0.79931600785741385717 0.73033292266123783776 -1.398958878374027659 -0.17019907952051174793 0.10477825606967920702;-0.94017330881906191475 0.097559696211218682738 3.2746179132088770736 0.6012754563476320202 -3.2387847693770281943 -3.2754435953977161233 -0.9815524038342318125 1.5278666947825298994;0.29778954720488887142 0.093585230719523099241 -1.0946515314346996206 0.57338760536026611714 -1.0416244474652207597 -0.45018108077723006932 1.2885134517499774454 -1.0376430145808266392;-0.33119738888336236471 -0.16538558057523142231 1.1635336740386026388 -0.68130082791086732374 0.30155253945861759668 1.4372248902255118086 -0.82246912254883575155 1.471773701177893745;-0.1919772421729232792 -0.15250012089066203314 0.60809930785291077537 -1.2520431232816435774 0.24121887346633733018 1.2013498236817039544 -0.92074306182128873299 1.6425360704672815615;0.089871728246635071424 -0.17969618322868299853 0.1098145583799106928 0.081471983338306616651 3.6015504144962364741 -1.3582086075078390053 -0.63507093225863664454 1.8406512413137252615;0.035638776525584892063 0.22689894440442284762 -0.46631589937992651285 0.067841586296138775558 0.86673737779562831651 -1.8658977502804168225 -0.26258723742904266363 -0.10564534315051697655;0.64542017860605571666 -0.037524730790488382282 -0.62624214755045404246 0.59734636406577223955 0.507524267510226057 -1.1986389669005117486 -0.23236433787595958678 -1.0113842674022222035;0.054554787847300931225 0.21035239298891003634 0.35331420288305548372 0.39890537131878023303 0.48537845785272620835 0.61582772078392322701 1.0898709877232657384 -0.45901047793875421332;-0.17193360765280352176 -0.2260044756216735562 0.64854949511735859868 -0.43506367395834188994 -0.40337261353446984824 0.35984622040998703874 0.54841838067943593504 2.67946072427561921;-0.27010565913786543923 -3.1698597659281344541 0.27285374027033942879 -0.68149236637275689699 0.55968302065877340024 3.2548630329364791614 -1.0945527780142307872 1.3385955870918082589;0.084666309963604857725 0.36046695269605488754 0.54154405954075102159 -0.74255276194278019286 0.11031343508326238445 -0.24940330327425580537 -1.2935688143710144526 0.56079260827377663823;-0.109559294384860062 -0.037466560404419653463 -0.28516657567159148234 -0.1116605537321450764 1.1266612965105866628 -1.8790300794398799855 -0.060905882707103627116 0.53165633644815646353;0.05915323858254636219 -0.044743455251726439725 -1.1523589931009345921 1.5169244237758532368 0.48187311162451973434 5.1971060072122119422 6.0194015269254466105 -3.3815097654328396537;-0.13774469569221831078 0.23843656964378026619 1.4644680475135458675 -0.30855676684543953403 -0.12444326545191643218 2.0946630376348003288 -0.11727941557779983373 0.54140612227722850758;-0.56259618363322305701 0.37042797951690181746 -0.76152000883488968341 0.15634824439926142681 0.69872989608180535726 0.96606617326289356829 0.41133228113428549344 0.70503508025371641743;0.22823935726289246162 0.16703163282642249921 -0.84951838812910407395 0.4948172754820003294 -0.26139224760974100192 -0.69489870136483478102 -1.0270723374749970436 -2.2845711135134112446;0.34658855783898023795 -0.32426123519637500836 -1.0026031415872138641 0.7153068084121321224 1.0041067625017507936 -2.1160389742613316955 -0.073993079395578115198 0.72152765583310130015;-0.11443370052677738491 -0.16602883370252752293 -0.6604140488360391581 -0.066126244895806943203 0.11111823995326643166 -1.5192681689275879808 -0.005693544021515216727 -0.022178020123371545297;-1.873490990202116846 1.3197682999770588008 -1.1582450613597705669 0.76453229350406992637 1.2819078602760545227 -0.39898848824487281872 -0.74910289285554643168 0.46641343397043344421;0.41803348278285834549 -0.33226698509287688754 -5.2229696400561964609 3.4037786871174411374 11.918226000328866476 3.5718131692414893941 0.77387886671688965734 -1.7011822505575469044;-0.0048206256604794075343 -0.1733932518473830231 2.7036303053179366884 -0.060257316149992717902 -1.4904581232957243575 1.4568569226230172475 0.97040868208316632959 3.3533961551329865003;0.16031780591366120992 -0.11672197553612453924 -0.18351564657982119555 0.011899968945582643787 -5.3122803192240581893 -0.500122228199940011 -0.9795816732705537877 1.4529075781401257217;-1.7672331603591737714 0.659676480518110675 -0.094933566147087533604 -1.4332943886087636809 -0.18481429548899369619 -1.0306007769831004861 -1.6759141615067787434 0.21421486749567841823;0.13118709316870663661 0.063353497284261234279 0.15536063237705724505 0.39811135425651239528 -0.57586345023332830007 -0.57899572455977499352 0.39023738776613020596 0.47650457097973786258;0.22042331474306819938 0.08530465127874797826 0.28295972132892993622 0.3706519171784657396 -0.15747099050309398827 1.9599953945776766862 0.70882547239050452959 2.3059573724276347662;-0.12251465025534076514 0.16937902989993933178 1.0938768972362560294 -0.77638627231466861289 -0.14619843115425529279 0.27658440046272497748 -1.396670758455969219 0.66508883871226387452;0.63456394441633667824 0.20833854970102935744 0.39675507347368338396 -0.021517243916131544024 2.6011519586185851693 1.0173893687962345389 0.43962949830781933303 3.0248483216178110311];

% Layer 2
b2 = [-2.1345725756647069282;2.585459777908783785;1.0946839140008501889;-3.9725200277756362865;-5.5245538537677667534;0.7487811241461408418;-2.9532169263543064375;5.3840126011578322007;3.5713014371629387966;4.5860580299331976661];
LW2_1 = [-0.1260695693863627953 -0.29446537193699434676 0.19802461765912282576 -0.11394088543713878914 -1.2408525003292092403 2.8211470480254061677 -6.0470374637683281804 0.5317612754252499796 0.66013586288524250811 1.9977605900079624668 0.63130652848438328917 1.1267163408837246319 -0.031418626318778314355 -1.3405671152147340397 0.033716227839733235727 0.3144621134883288649 -1.1763696536996417485 1.4196739765371171327 0.4836285913541231718 -1.3506043359321608932 -0.62226957392102943967 0.16125448193523217455 -0.38675917570811496615 0.92060793030483967225 -0.095086746446950060063 -0.26280604968989601522 -0.85501421774000130327 -0.27284865273860375856 0.79668514133048629056 0.0032817465092558782344;-0.8978705133754647516 0.90517658036404557986 -0.054492826814719078077 0.0087475056821867645807 -0.15636863089060676346 -3.5659168143328532352 1.8714385078353172887 -0.24225833414638264784 0.4415992502205683401 1.3194298486008475546 -2.3038476136166332608 -1.4356292724288015616 -0.24524718309053336607 0.18694625697270825238 0.54256692263906958207 0.0022565654328518774778 1.4964429140548034525 0.052694540206209762589 -1.8420398017416599323 1.1593388142153946418 -0.88336280715993564794 0.55997977777149343659 0.24026604185392336221 -2.0369444636768689882 -0.13241500938076380955 -0.8338672875101490467 0.073339613833801831078 1.0312881607136306084 -0.39418970426085075953 -0.66125166720463934844;-1.8764001593811372448 0.033226179790641158229 -0.051690596060321288008 0.05131537543014774716 -1.1250734336798899893 -0.26687032376527591859 1.89562815132235718 0.071607863217619929785 0.47490680590574269049 1.0230886951739581736 -0.054952305682999254277 0.42758942922693737865 0.24497607058148257453 0.57583967282266190768 -3.7896995619759472618 0.3575959016930621015 1.8391903729648573229 0.51975851013867557793 0.022643805233766647889 0.38286521344001289702 3.6970159415421872673 -0.014234538748885682338 -0.17466140730966284167 0.27347289466341873077 0.32634368058870549811 0.43025840578034846295 3.1817658132263542115 -0.80596638920111707272 0.45054593963026107195 0.054605605545596222483;-0.43499156680271661157 0.74573345073832941754 -1.1257239749316796384 -1.1099302306010245811 -0.61103922046295355841 0.74229266492335232996 -2.2152285535895162027 0.61095121902743043663 -3.6059852067703670464 2.0855720622015980759 -0.37533789909658776818 0.57204943478028658088 0.6143390733249504132 3.9433445319885902514 0.82204697377843316541 0.079228808284157459241 1.0110663019301802112 1.6470445804856321814 -2.1921686589682618518 -0.11208209605950353316 2.0784538669801242428 0.54267158631366396726 -0.85959341717147952355 2.0533391835545464765 -0.38470827124818690601 0.15198199665902947797 3.6855615660489569052 -4.2111505764993957257 -3.0465074037862178713 0.14605305395005194202;-7.92794975405612945 1.2578047253842432607 0.034048685412170535958 0.14213856072779082473 0.4172952468143828475 -5.3018882970017990885 2.5853447722849098867 0.24013971158525257432 -1.0758001699615085123 0.30444884392886240043 -1.0395532183635831647 -0.060216054336111476952 -0.15963206967253010249 -0.29634118049110103055 2.1247569808782960088 -0.18283235485763668637 -1.8043034125331172124 0.29336116365454256716 -0.2388706840384142438 -0.052405376329076075459 -3.4982182510677297671 0.37862077395146209735 0.28715797001603732275 -2.2881984882680295179 -0.14100470459899649978 -0.49460306124780423209 -0.70837897809547323291 0.83751437203182821367 1.1882676396395224216 -0.73124145471622936654;2.1768072081994769285 -1.0024546839845540624 0.1812251398399291058 0.078718393154567092851 2.7263014462083354772 2.3800334024520317655 -4.4908647347052506049 0.048723031260240072526 -2.0433778790735330766 -0.20330462616659181196 -2.6964338670253127361 -3.0023243065630356163 0.12678937974242593389 0.7318627993440509627 5.1230419122002990306 -0.27983055550585689275 -2.7326985809019341822 0.089344512166447176416 -2.8248571018886261541 -0.92319178883582575512 -6.0143073414949537181 0.52757200618380806212 0.040429756652352252799 1.1351411922246426034 -0.30356393615473342695 -0.75771559584388115116 -0.35583198090921958068 0.55405044727601837629 -0.98335821855545757497 -0.56776916776198849135;1.0493442541995270378 -0.23988023021118845546 -0.17317950128115103259 0.047543622010965611024 -0.20073838564529514206 1.9391674767596918105 -4.2138777808958778337 -1.0580863521538641514 1.6417846761243217735 -0.50555887594522130879 -0.27456676239082244972 -0.38681086780077655884 -0.0070001066506299091871 -0.74505721598126495309 -0.7487625177701821011 0.036832834904382348584 -0.67503764793211662454 -0.6441342582681973683 -0.0098378542846520103393 0.015973300488306661771 -0.077935844379089125322 -0.037990858710416124111 -0.12194489657628794499 0.65300824675191826429 -0.0033911766542152150736 0.084902427859504026042 -2.6807160899858493153 0.65567053759365612198 0.71680742869345692725 0.054050991955182307869;3.7110434176146909202 -0.059677011137796935614 0.31816010289988227999 0.35011627599488348572 -1.140643354970398704 -2.6035263622705784314 5.1697542541566603447 0.31105857243917056643 1.0599792638760610064 -0.91384988165602454302 -0.59296962092078542295 3.1556255228718317341 -0.028235790635429580675 -3.063393043971189833 3.1956532213765687089 -0.23817135167339734125 -3.9809070314386456602 -0.34729206265961404254 5.1003865141591688825 -1.9491986034993389332 -5.8223355606847935206 -0.092352810975721261189 -0.20301546084060093089 -1.5170755243450562144 -1.0791370785349991834 0.59377319297807218401 3.6327201709065364277 2.9753824144602067925 3.0456800636889886924 -1.5883498750410551814;0.93276147637128314916 -0.78306840659488230472 -2.2851176039365790338 0.27780983575052431434 0.637782583071066278 1.0803639857985194439 -2.6041603691426025158 -0.04169561220073581681 2.5897971211682540016 -1.3381657992337339103 -0.6123603423891617048 0.14836320397958591144 -0.17196817348880713561 -2.913884079688457529 -2.4919303376001771255 -0.15914359127723934351 -0.90398813966137958609 -0.83432108919269598601 1.1493880268013503265 0.33415182516622637943 -1.7626114304037185931 -0.36999989804963512308 -0.06923449260030786756 -1.3617872394216807841 -0.074823873788535361085 -0.39041804367239280049 -1.1545172359530204709 2.1056639286626190888 2.5654133743816842639 0.14261039269734010482;1.4628439767835599294 -1.2591279968967254987 0.18918711431335161199 -0.06287270106062446684 -0.52279854754252663707 -3.3501876754026023164 6.9675771204780065915 -2.7967475236436745156 -1.0894806175658209302 0.038371393639180899349 0.71904645671218525571 0.85193611335708485655 -0.23336630060309093704 1.1836664428336858279 -3.3024452709138141238 0.29010934458539155978 1.6464444522269428006 0.60648081669248332659 0.81375033285666087401 2.481237277536649799 1.7754852149885909096 -0.05669143304166035513 -0.037786089945305668947 0.73386037059795650173 4.159104756789817614 1.265800445271823893 1.3489994168101042682 -2.1865673241072953736 -1.6062088550367747963 1.3882059761694918176];

% Layer 3
b3 = -3.805208537344134001;
LW3_2 = [0.26468371178673627542 0.35290241730553562061 -0.70415729121757542064 -1.9307105395092298838 -0.32259768529227017542 0.33815145744711255782 -2.820207764596196931 -0.54736340547220407871 -0.68485540059428873239 0.67614142515220243546];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.0415843642790311;
y1_step1.xoffset = -0.9;

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