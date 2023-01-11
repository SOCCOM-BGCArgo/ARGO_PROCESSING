function [Y,Xf,Af] = ESPER_talk_14_Atl_3(X,~,~)
%ESPER_TALK_14_ATL_3 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 12-Jan-2021 17:03:50.
% 
% [Y] = ESPER_talk_14_Atl_3(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999999390765;-40;0;11.29;-0.2178];
x1_step1.gain = [1.0000000017546;1.00031699124092;0.0153846153846154;0.000310077519379845;0.0705602618858232;0.0479365703301392];
x1_step1.ymin = -1;

% Layer 1
b1 = [3.0078365418079631;0.42389890632421434;3.2592429743960314;-3.3367953276615765;-2.8653156148454517;1.8884592758531613;3.9369378981999419;4.2631951738839531;-0.15766768522775274;-1.5811938349830437;1.4940130739784105;-2.5162675997237991;1.6046288474696828;0.068728626496701112;-0.29373833025739021;5.2412673028669232;-2.8271122664146078;-0.99050480507562633;-2.3511426825762101;1.7104254098017784;1.2091857148110972;-2.3322888237867585;-1.2996969274972681;-2.5498480400123009;4.9108371026851847];
IW1_1 = [-2.3245179424251168 -2.8266530151189926 -2.4616902615402765 1.3573259906230175 1.4236297384906005 -0.74762761038532954;0.80654507309940493 1.0997890367250018 -0.16924083980581911 -0.15133935384729469 -1.7648794252148994 0.050342100488949139;2.1627000191968433 -0.74333507236015828 -2.7119999579780152 0.29443460219091466 -4.5589477328513821 -0.34237506051273647;1.7426738396397499 2.0975411563613058 -2.2215603128804822 -0.3970021163937974 1.9111264168893973 0.31796488532673817;1.32463705809447 -0.67141452056117756 -0.4228968072937469 0.64666574608388605 0.29168486237910674 0.95462360885508313;-0.59099578654959484 0.84067941239987209 1.0401703271234832 -0.64977137354201231 -4.2966157497289075 -0.58854570721310295;-0.9002795452779776 0.96587128918925158 -0.013957006434437119 -0.87916284037108539 -4.1839188425901757 0.051791060788398333;-2.1409658532908131 1.1697252577354478 0.26621457958829653 0.031177630241356304 -1.6955817548517766 -0.83942803508710662;0.25531790581656771 0.7149952512675456 1.0400253536059987 0.41083124377033597 0.92629168741615131 0.87197756976191321;1.408210565742356 -2.0878129155252845 0.075165236546491779 0.18693990444352415 -2.8560202385457636 -0.45170632202334871;2.1632440444345464 -1.414611500566433 0.83954568824092213 0.79550670284945901 -3.8216158338160247 -0.44625729803995556;1.9696096614316951 -0.89095550410247848 1.4382561064310995 0.93637091088720104 -5.7513548518077178 -0.65466554545147004;-0.47120922509119101 1.0053041164993977 1.2806319910304496 0.031871494930990091 -4.1136663915212166 0.26091924960825552;2.1300864335866159 -1.594705942786931 2.5507926365657747 2.4741875662281854 4.609365803312012 0.0038185838339399435;2.4482593169764058 -0.34102657141086562 -1.3007921487395868 0.56130817471371419 1.0694009768697801 0.21146180439036227;-2.6094810227299594 2.5527500464729052 -5.0601895390460623 -1.0165776217958287 -2.1438332147147294 0.17883103353498536;0.20287863367194151 0.045004470828263379 -1.1371259520628234 0.40482077388230403 4.0837700486776134 0.42578689700976291;-2.6593373854153426 0.23066576429177152 -0.95318091419667006 0.47753065239235887 0.030439209224498984 -0.63924305964763273;0.50286919176439548 1.431588915042693 0.63906316612763037 0.012252919154759418 0.44038418371080212 0.27295802461959234;-0.72950982105190088 0.95520938151388091 -3.3704102192580434 0.19956567822159535 4.1242339259092784 -0.70378806268977545;0.74039986684949177 1.7862936265099487 0.5583291654876551 0.48138224567914423 1.0132441080043644 -0.26983397469795645;0.23544966012473559 0.28299142474319261 0.25285965659677573 2.0465571335386006 -0.55145942880740484 -0.55015296145590498;1.3569839859540649 0.012060464083056988 -0.12204490866425467 0.16094643862439581 -0.3894676068964586 0.72161315198110854;-0.733583331173844 -2.2916332165875311 -0.56413032732471924 -0.58400614211805257 0.070524325098208496 0.39116941866491645;2.9086801616548574 -1.6988716996307047 -5.5440127641153705 0.54557819456173007 0.36700607062081458 0.7105402299096133];

% Layer 2
b2 = [-5.8146749017383206;6.0562546878442962;0.55497331370805703;1.4772020270344794;0.56512167180741346;-0.20520491553411288;1.6229091457497726;0.57634643198539204;0.94021333663286188;-4.6670863806872882;-11.851583517306088;-4.5143085994081753;-2.8977842157093403;0.57323388347025239;1.1668360539010385];
LW2_1 = [6.8150371696847953 1.7638078679865905 -0.69295128494993363 -1.6097901559421568 -0.42493770670169051 -2.5249576597336354 1.1452032860874473 1.5459666208559495 -1.5592702824540725 2.0627712973905781 0.87429303705838857 0.64386940319237229 0.61808892296251627 1.5786515378956956 0.94872993231349501 0.22598636589260909 2.0209036186535876 2.5314582094479738 -0.24023225948441185 3.1719289272130831 -1.2371352380928455 -1.5390559434975293 4.5254237817000407 -0.14742213319011524 -2.7170909759116268;1.6906163395265887 1.7393466888601732 0.64838976679006277 -0.14356996528513261 0.25818474435561001 -1.2159545798049214 -0.33489903707338164 0.34961576208316597 -0.60593823599432062 0.68819776704022761 -0.10289895364074786 -0.71920401952223534 -1.4423406948135848 -1.1175637141792729 -3.5262286868206085 -1.1542205789702102 0.24714230191496109 -0.19416812645532583 -0.93439782702532637 -2.1855203429741836 -0.31153308500365789 2.768999911993713 1.3928141669048291 -0.37457508308552884 0.34025860470997915;0.57628951511071069 -3.8197332347858279 2.2822137618342868 2.5714543709097137 -1.433755116107803 0.2876739343792995 -1.8465762844783158 1.1942631587360597 1.5497301872165157 4.3975179965950781 -1.7530533192479911 1.0367547532288368 -1.8045065938621048 -1.3471240368781054 -3.4514403383303609 0.20912348271290171 -3.9089187763762596 -1.2629157014363743 -0.89334097089956466 1.8731419566953642 -0.52428533453559734 -0.91345511248964617 4.2724569905197436 -2.3988391345874414 -0.74314084050440521;0.19611656479392661 -0.48733463202380278 1.6943343633390904 5.7917196769846084 -2.0963750502907783 0.46708611540842637 1.1053873169948951 2.3729584674087558 0.7803820239220195 1.1294486325632891 0.53904066090441982 -1.1069792874488227 1.5265113710649818 -0.80826750454753893 -0.43185544720463781 1.4766215924418797 -0.022720479840106673 -0.16887759073039815 -2.4971980416739088 -0.756265768827869 1.2914152892552906 -0.38153413475452369 2.7912077420105974 0.26229986931081745 0.9147612562014753;-7.4017962715678447 0.87562867217681395 0.36722152262831775 -1.8925720976699509 -0.29684704598745304 2.5244254232818588 0.86872088011440207 -0.049563749450261915 -0.17944871117599864 0.20707463302242832 0.53923115754947459 -2.0138814593256185 -0.57529867992805228 2.3011281758071913 -1.5817997400521284 -0.88248922037536959 0.39662024702966098 0.55171809829655993 -3.4115878014984489 -1.5928825743886239 -0.027054920205292987 -1.7119190398507007 -1.9761913534173485 0.65014651052338979 0.014991907105945376;0.080404104775369911 -0.020860713834155897 0.89842425262931014 0.28928132644759724 -3.0937650959723966 -0.5199963425707097 1.5554458287183199 -0.19351881901178261 -0.16452969673995493 2.7038458574757116 0.26405143740422532 -2.0083664313885072 3.4093005520005852 -1.1231808243355974 -0.55385316765849413 -1.0983363767433516 1.6834834126242595 -1.8232865912074563 -0.9906364606635526 -4.2587964098560862 -2.6388618535029038 -1.9046253588871931 -0.52390118834255794 -2.2902495665353295 0.056372881748134339;-3.1045327226997927 0.73834799095752368 0.74434684336862955 0.17941044339568751 -1.2338725464202445 -0.88049310270138315 -0.4643398511631156 -2.2801558727239786 0.52153144213665148 2.2490319044633109 0.12426118601231967 -1.3100971266079526 0.30611801340050976 -1.4749363323868634 -2.1728221008864215 0.71354147662438849 -1.1388791540635124 -3.1469741331621397 -2.7170770826809236 0.97073776502912656 3.1189953939415549 0.72364892691350735 -0.50450797557890559 3.244270479590873 -1.3479461553757857;1.9763223549781885 -0.64136648831290433 -0.63966994057141102 4.9526033621236092 -1.247864711819181 0.69788273030617198 0.71473460397334942 -0.0042619022805310894 -0.25351549300466553 -0.24446646286394857 0.74165027001870065 -1.9404609686341481 0.84061219065295367 -0.83894428187891124 1.2144007663032468 0.4934027833205834 1.0757639683490243 0.23152447715531457 -0.65785287702178774 -0.5698056713305425 3.4251870003451748 -0.076427089917452534 1.3684127138341327 3.1018441360262079 1.765229390347141;2.1382579450901442 -1.0987738525861885 -2.1195183202540475 -1.7484622866241644 1.0712444316907874 2.319746025176094 2.3651727283306765 2.1213803632733965 0.17902289087578907 0.13993836493570572 -0.6652971148480481 2.4618242906160832 -3.2002564761623105 0.19611499667044802 2.9411964258805807 -0.50145495542561891 1.3201204087805778 -2.7453638687148896 0.94274523250190667 -2.4482068968715511 2.5983283904067447 -1.1898509017184258 2.1941621947078058 1.6636466656652056 -3.1990122004930068;2.497107731804689 1.1082271424856394 -0.99235110047116282 -0.85694468859362549 -0.24234138835064234 0.75761041738360602 -0.62497273886568183 0.94573995118715382 -0.77086849423129922 2.4920360561924633 -0.054525004879786995 1.3575010267085574 1.7667630366746931 1.9231489342243642 3.2380224348422932 0.32436180954435007 -0.9739766475588999 0.10548722771203498 -2.8235705105414688 2.4669502482695571 -1.0544485934275412 0.9474969268010911 -0.59806564476950907 -0.00054113497107723986 -0.17741504122747559;5.9031562495371785 -0.59391217839310007 0.5449976746327917 1.1790156700138994 -3.0388640036869714 -3.7234999422352297 2.1140720652601037 0.05879401290662345 2.1601061578148184 0.53997248336686376 0.24742209154418052 -2.450072723144149 1.1827719671343577 -4.7451274503313234 -1.0305424534899188 -5.3583416902303815 -3.2027440406584398 -3.8731822945335153 -13.825350678198454 1.214109447297687 -6.4515601014193491 -0.10276994549446364 1.6224756222883856 -0.14415705915827334 0.46793944682971089;2.1158784450671915 -0.93103273390688013 2.5427006902585529 1.6181367396841873 0.18830214455234034 -0.1074885519448204 -1.0704903521236444 -1.4456036483302257 0.16223246164590077 -1.6756807360948867 0.10941508197520006 0.92183659983184563 -1.1141306629799903 -1.7240496905224374 -3.1857207974988957 -0.016766742967567382 -0.013835555785117163 0.88280498921814021 0.87565796729888123 3.3447699514065556 -0.24665586259836211 0.57270399155046248 -2.908859448128692 -0.77333280275692995 0.72548233165334208;-1.6118815896867713 0.66764045698473695 0.52752940477721455 0.6774896654841186 -0.62357893134840048 -0.70840788987958248 -0.78661477619980391 -0.59374821368247088 2.3593882562342694 -1.8866114385212553 2.2903248903907309 0.22803948530195514 0.34829458891073817 -1.5423767647306046 -0.4358685869318415 2.173104287275085 0.062868888495732772 -3.0713800648317853 1.0048470865883017 -4.6048922131570311 2.0431318105260274 1.4622847840586159 -3.1454283509940888 -0.86103847147053658 1.7355628390801765;0.45829213714151612 0.97043795543492839 -1.1900653045798144 -0.76085361354139669 -2.4940252028293295 5.5814303572861039 -2.0596588816115284 -1.0867438913434662 -0.094365820873680853 -3.7696335361074147 4.3535265339328211 -1.0220594671058147 8.7358552003562391 -5.164930915443291 3.8810370097292357 -1.2006122825472185 1.4745531795433735 -4.9647764884798793 -1.4591047705840687 -0.93109458897496677 0.54003028234930761 1.534966035562529 -0.79985218503604649 0.93261777074730423 3.5046503844110974;2.1935682944155519 -0.041176036698399428 -0.71101092813868028 0.51142164677360435 -0.28745018358330221 0.65236441872701612 0.60297862389088797 0.17904042732486303 -0.19398593947134976 -0.15854592471802048 0.6865469703188527 2.8767079804625286 0.1530847204856354 -0.87858774400281137 1.1513574423511934 0.33948824880342238 1.0932258559317005 -0.071756666759519833 -1.7339175366748092 -0.7358986491725602 2.3062641742834957 -0.12938865410302444 1.3134954689661311 2.0947283952968916 0.77251501559680424];

% Layer 3
b3 = 0.43986616294878517;
LW3_2 = [-0.066778529716170584 -0.070907007546424192 0.029256765949388742 -0.03427633883525795 -0.088441515300224288 -0.038502389847894759 -0.088260208986621957 0.54532113948211292 0.084308525114299745 -0.036135691431038688 0.011338455643289162 0.079468570981286082 0.013543180681420621 0.0072043841319880744 -0.53685438167597976];

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
