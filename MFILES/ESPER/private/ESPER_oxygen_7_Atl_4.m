function [Y,Xf,Af] = ESPER_oxygen_7_Atl_4(X,~,~)
%ESPER_OXYGEN_7_ATL_4 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:46.
% 
% [Y] = ESPER_oxygen_7_Atl_4(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999994516886;-40;0;9.1619;-1.99149636107278;-0.28];
x1_step1.gain = [1.0000000017546;1.00003626070655;0.0153846153846154;0.000309885342423303;0.0656325902039552;0.0604107502911417;0.0470499670650231];
x1_step1.ymin = -1;

% Layer 1
b1 = [2.1854451345088001801;-3.054446111409474085;-2.3603653421679280378;-5.1462767602841905656;1.3719925407514346194;3.6427420852108860494;-1.8894526342739466074;-0.4657928666303359222;-1.9583291072986794834;0.98631688300204700681;4.5970549928805377604;2.1017564196236220297;-0.43775681894106321934;1.3099277373438589223;0.89424240203983029751;-3.0366439084184619901;-0.7402511794739795592;1.2622502112357862902;0.092931418821123337737;-1.5178825045264496829;-2.9130948449404745482;1.9202443410771352639;0.66213787208362195891;-1.4860468597526634849;2.0748706904253646499;-2.7064340561498587689;-1.8209824260325648826;3.0571034772908287636;-0.23407795367203884651;-2.0560654596266103589];
IW1_1 = [-1.5361119535157219396 1.4040410716071727304 -0.27859945637084471137 -0.99498218057187737262 -1.0077045800996822233 0.085212337525968406826 0.84245184476600332157;0.88611034139065747439 0.49757107177348081084 0.022087350661811024771 -3.343944207688728909 -1.3278145853358904738 2.90695526724236597 -0.15611448250529430526;0.77456482149224137324 0.54949428630028795961 2.2829634624110664021 -0.99376256418302977291 -0.32310058496635324365 0.31683047914369039111 -2.0714760263751945679;-0.75608558600348096945 -0.35742104648743655559 1.37936544082584156 1.0629709146722237101 2.9497934510263426056 -1.8759774807242552974 -1.4017214675889086983;-1.0975742842820956913 2.5361870977146958595 -1.3949801041954468239 -0.69358705906583661971 0.68353385744009631519 0.6355927259833104026 -0.03885136630040956357;-0.463928033741547452 -0.356101822908488852 -2.2646829935975456571 0.57589685936801937594 -0.92136536982713179711 1.4582350346158547971 0.8516411507142200854;-0.092216351356733963773 0.92309369285863340782 -0.27852491474241480862 -0.48180728045418430572 1.7520123291911728902 0.43270561351741426703 0.91392929381145127099;0.37013606620741029696 0.45605154756846638664 0.10601952683293236479 -1.0890726109710300307 -2.5371603442412138385 -0.88067249452981655189 -0.13343242423683066011;1.455246360267897332 -0.07314930161765260952 -0.38318638963750972781 0.034353223649698623465 1.2381036172765789161 1.0177059924794258361 1.5877845321554417968;0.033430575417266517046 0.14498515460743810479 -0.95735014288838415908 1.2738978509719867116 -2.4627493771281789847 -1.6093435617484574962 -0.018177817896613383125;-1.0707595829159168854 0.24066153405956366984 0.00064781430688714593258 3.0756830010609439441 -0.27442957427403408266 -2.759487846902823005 0.33973418973656005493;0.4002972837886111801 0.74008855785311933317 0.39958191983765106725 1.2828885489087735738 -1.1530291324417729637 0.30589399107906728448 0.92335193773112778093;0.57247125117871533462 -1.3400972564377739982 -1.4123498510998673261 0.70228813630246200717 2.0208040556053319392 0.76465769448684173248 2.0540380159937581261;-1.3630724248224259654 -1.09642296323406363 0.67508996156172118575 0.11531477711290610833 -1.4656252493527519576 -0.41428633259438274461 -1.1130830091874861498;0.18309358289173480183 1.4170402481879380563 -3.4436181651645818569 0.6935062328061353254 0.7575767913570106149 -2.3136993270437171688 2.4245692628846997607;-2.1120157885118611141 2.0049015839083730306 0.53036558346900275485 -1.1495340776395690519 5.1000890613029339349 -0.20561012083221907876 -1.7086534007757656628;-0.66757989058556466144 -0.80584136131643047118 1.82905258396543835 0.53385081875355600012 -2.4830967282591931422 -1.6933791534795747946 0.91459673114699424623;0.064257299786757662852 -0.51305493384572575799 -0.47858728596831479063 -1.955449462192449861 -4.859514118418564621 -0.51771306961109453404 -1.2434428303033160379;1.6242621751730210722 1.7143849094794905152 1.2106161162475641557 -1.1944669254289228544 0.049124335699661562937 -1.6319565130522541985 -2.0851907044880180209;-0.12105251882951290066 -1.1829440820317771443 1.6425485488803894807 1.3517630212884654828 -0.43819053593237133892 2.042823313048296896 -1.550064994299857668;-0.92042824211089901976 -0.66861802347061594087 -0.072421726726673865682 -1.1636394470996889172 0.75486808609640276835 -1.444723610162373495 -0.90354629158645360931;0.58752053018343730617 -0.376199070434428684 0.76295063682309227637 -0.45745047339544786524 -1.7197184989311815606 2.4319265187543268603 0.0062425883068585319435;-1.0521910812858297213 -1.4865094357695636251 -1.2289915604069949762 1.7904965885598549402 -0.095967615874905265461 -0.17937605804947523436 -0.085481227215346497483;-2.04279277089040745 0.91024407540850804654 1.7235761353775826965 -1.7658932849716668656 0.23105216397919584037 0.82338804115537123085 -2.103789676851908208;0.55140164177665429435 -2.9531239435578608976 -0.37208562088782898281 -1.7521321008335226654 0.70745535949239668483 0.090188172151595163983 0.41928818926542116285;1.0358494441751409099 0.67245594714889722621 1.4774411280145585135 -1.4882391389884661592 -1.141518339745225985 0.40805913946383531865 -0.31568124917143786101;-1.3104400161830815907 1.0048088245265927032 -0.40921022079488705847 -3.0455259778151142491 -1.6104466763639644356 0.96141466822115162483 0.73819824800595312819;0.61355081014195989031 0.13137798191294786099 -2.8204800599436357089 -0.71600494020330507361 -2.3880703794284277031 -1.3916741925132345425 1.2415709076319281223;1.2812261585617858994 -0.38044338347209649687 0.96593547812287317011 1.188665958238288356 -0.43978043183252185644 -0.85973408153840247259 -0.96541640413617391125;-0.28544899482644175981 0.88577151304959134137 -0.98202501359153515637 -0.42230262222476577794 -0.093027774761792261793 0.026311270166644099011 2.0737364212731130841];

% Layer 2
b2 = [1.8189559138702204866;1.1554523584482869758;-1.5414421778266000906;-0.96177988969054872825;-1.0472294799047232594;-1.2513266375045155598;1.5155810180786843944;-1.0300391992428303833;-2.5241990209457125793;-0.87375726733587932848];
LW2_1 = [-0.65719826343161857896 0.28057160577662859646 0.88609133911987914001 0.20185769149042906423 -0.34943440466757064167 0.78256939394701596768 1.2002426981152694019 -1.6694102041712424178 0.0065742443867086355208 -0.47964343838142181609 0.93932746410872702825 -1.3656486824733422569 0.31119798278712074158 -0.19979732517144951909 -0.20632814209709263165 0.24736841843422521792 0.39040967350697108884 -1.5675184438296736467 0.39951721001438350589 -0.45449876248483711549 0.57130224029861564894 -1.4104870194151279961 -0.39870928172852337923 -0.17233421568343695873 -0.34595300545670126802 2.3595191159037809392 0.71822415279908424779 -0.16195589174262120524 0.070893893699885385251 -0.75542540247523681796;-2.4100711675183870142 0.24873797055558230507 0.31937145991510029752 -2.1342893143469519046 -0.25558987036949221139 0.09871898466972166275 1.1943407795092659818 -0.88499768854371330118 -2.8628701470512418759 -0.77120369124359822166 1.8894446621344138659 -0.68038274279040900083 0.78954131223664569816 -2.5954685015372658263 1.3536797565445768399 -0.35927580744892528841 1.1085003925491154764 0.65259574966891931158 1.3909735105834177826 0.17890926511372917673 -1.0873342220281831949 0.055540459162856922748 2.9754028643781946784 -0.18128509831539779684 0.14395436902200398066 0.22552408749018404421 1.9615389188296992184 -0.47815311442804292463 -0.97390567475069089731 -0.086992655061136417438;-1.3693341243151320175 -0.62942560163412786256 0.067129193413949961311 -1.8497840779367156205 -0.42687208130016057739 -0.023335767098316403945 -1.8447975512436474155 -0.099673478213421207172 -0.5352749016667773807 0.16211720097286896891 0.75570979629924994736 1.7106810488658008573 0.9777251621262518233 -0.30678011025263146605 0.13178673626230608118 0.78095461899882911982 1.3081080531987450133 0.58381244806900112199 0.073327179165154257601 -0.34182358463802492432 0.71780021537509297591 -1.672942243573284582 -0.24036689386830678572 0.490642849549984994 -1.3013375821979187297 0.7400292458703785492 -0.50061971097234725381 -0.82838278516444396882 -0.95246797559736706695 0.28577394838866410076;0.67943375772128977719 -0.66700517179543983293 -0.60164522923930985332 0.58277249717041279542 0.1422879698589276376 0.55706049043015759548 0.55726263250281804496 -1.1319518113578890262 0.058419100727560964448 -0.73309791658166756356 -0.65469938673251115713 -0.45693935753182080006 -0.39114202841044398706 0.44843540087663097404 0.019406707127496702903 -0.013891051509450995025 -0.55558290954398137362 -0.34874907618513945051 0.064891940800612654683 0.7526118954764495772 -0.85263898668659865443 -0.31747145315548513933 -0.18585857947530526335 -0.56298329373268185538 0.32400202515039133733 1.529466619007441075 0.021542422887738471582 0.16348150053599405651 0.29036073782435933843 -0.19932342289988749795;0.8515726879390232984 0.12788536417511381371 -0.68040232839573988599 -1.1397565558613149683 -1.1564713568508795394 -0.99890205817028954183 1.6225254356151739277 -2.1711924459691447353 0.58209557261679767404 0.93747596594994875296 1.6906279183694215007 -0.067118176119592434792 1.0561752950452607536 0.29027036389720345388 0.30653759865842328258 0.48511712438875326248 0.66846376008414398751 0.86860185519880772631 0.76199779369772169257 -1.4602572313597437237 -2.5514954816166683926 -0.31121803748314119664 -0.61024768328066159029 -0.18132760250275875569 1.175816017646279521 -1.2526595237327249333 -0.52082814540213717702 -0.77554170417019585582 -0.936898860193650318 0.62133949461106829926;-0.2048179992375192271 -0.0042639043400775050036 -0.71668218847679787675 -1.9824202186289210381 -0.077066871543745085749 -0.89151155386796010482 1.4455636724557503481 0.51608039421151252757 -1.170909404207300275 0.40883196753890588271 0.17971599814316960297 -0.54393273747978621024 -0.68464444354427411721 -1.3623847166603058056 -0.013231249912193439322 -0.35904767256014930021 -1.9801312596800002552 0.19449847884664103748 -0.23463188053931063881 1.0897278447090130715 -0.87721566760514280148 0.88206282416539305569 0.47715737073373937216 -0.78446997022967535074 -0.10823857506838693854 -0.8554454376890123779 -1.1775773739068420287 -0.32058712686419799098 -0.82388771937865779016 2.2535542131859100934;-0.46963609056173538647 -0.51659887151328109933 -0.6661347338783846217 0.054188371013915409558 1.3808398700685025062 -0.23103554453716412809 -0.37095138546643219302 2.8300972368009240654 -0.13677700645321355655 -0.75804912550268610527 -0.34196198293541280044 0.17879402359993018523 -0.40582404091224366649 -0.670469453087202516 -0.30942926468432663256 -1.4968266967855212712 -1.582893366209269459 -0.56318710841032748604 -0.082525673564037865804 -1.7521313516476777483 0.31145297771383123608 -0.059853035660791670258 0.38631995091082466098 1.8691536302219242671 -0.78672451337797866255 -0.87133560846401014732 0.27750713552662725769 -1.8285701421135680533 -0.73289658026383719402 -1.3667879589977280919;-0.28767337039410767474 1.2204738394038083804 -0.40960261034124079149 -0.31980509156086911471 -1.2505236340751164636 -0.734186319360627615 0.90879626494486365207 -0.29022879899088527278 -0.29967356923106747146 0.30370915951271809652 -0.40563670257934064667 1.2988347402477773418 1.3010737549175230843 0.12432966863384314893 0.17510001897575500074 0.18401551301540991501 0.55474321998102138487 -1.4327169694626253893 -0.10421247876376979469 -0.073528570486313768462 1.3988535524908480845 0.50340003049306736305 0.017083494749658791056 -0.058814929037930993583 0.038334153220272490503 0.54267702830577990714 0.23777631424975684271 0.21531388494824943436 -0.57429441712741402526 -0.28132121077412836074;0.58684502199973687109 -1.1783985260301461651 -0.08322312494696774976 -0.1590292373637688661 0.023494248194537603602 -1.1684231452136233376 0.29338564894872459776 -1.4646944792554437154 0.87904375275024215775 -0.95612256413931628884 -0.33676940163284646568 0.58593011305909747932 0.46873833615982757328 0.64111736375317984749 -0.28140787783921983323 -0.074193866802529928606 0.27149827086525679665 0.30895907434471664432 0.11494611840557519833 -0.46880277092606037481 0.72381384933315606744 -0.12091643134727894648 0.17156868869780586562 -0.34208763925046720145 -0.51134274463389794985 0.60171973307328796565 0.4983431000213571771 0.89924287400311098128 0.34509696733383571354 0.48024698871641480213;0.59842188644212801485 -0.68323424028093326932 -0.78383473561747363245 0.14810912546158200298 0.30239569937414528811 -0.058215856955579445287 -0.98019089011168536452 1.3029889859277985487 -0.62822719660122383623 0.25079631138532598733 -1.2219366554484394349 1.7869542990822906425 0.31240148582688392898 0.47560630221964489417 0.22672366543694011387 0.17471672081026776158 0.21431525383605204427 1.2923663913518173008 -0.29090713700344678072 1.1064334155019399653 0.23257700360279739971 1.0872878998864479172 -0.035192114338799776463 -0.18478141158607913619 -0.005184257977672157916 -1.3929960402211380988 -0.5906809999659135002 0.77782400865519019106 0.1263709375995986095 0.48889549244321089949];

% Layer 3
b3 = -1.553826014166711067;
LW3_2 = [0.48424200683007356805 -0.28370289415783861431 -0.22182046215214476503 -0.62843251045297165991 -0.15064837384702042811 0.3092678499089080435 -0.32391014954418773097 -0.58758057719367451366 -1.213066762758342465 0.44339363342021187453];

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
