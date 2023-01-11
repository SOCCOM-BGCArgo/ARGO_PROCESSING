function [Y,Xf,Af] = ESPER_silicate_2_Atl_4(X,~,~)
%ESPER_SILICATE_2_ATL_4 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:38.
% 
% [Y] = ESPER_silicate_2_Atl_4(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999994516886;-40;0;9.1619;-1.99149636107278;-0.03;-0.28];
x1_step1.gain = [1.0000000017546;1.00003626070655;0.0153846153846154;0.000309885342423303;0.0656325902039552;0.0604107502911417;0.604229607250755;0.0470499670650231];
x1_step1.ymin = -1;

% Layer 1
b1 = [7.0014197369842445795;-7.0782551855721536072;0.42926080756951906503;-0.40236990541294709312;6.0528011723216126327;3.8949178267897588412;-0.6820589996184158732;-1.4003931842733761126;-0.44131785747912261053;-1.5263087928874534693;0.77549622462238487763;1.6833942311174221462;2.2726796367977861735;-6.2922228392126537955;5.9550074232742042923;-1.0536178660335155932;4.7029354629600508275;3.5165711196878191025;-1.8083805238974122176;5.2146551794176536632;-1.4369480547989821329;-0.15253994667918796546;-3.9193371104245264824;0.93476568539613036357;-1.6906588658557966109;0.49651233774612357763;0.0087727188295085179598;1.4779597774853709957;3.937347531320318339;-8.131322879568610773];
IW1_1 = [0.80535963574427082978 -0.93731847026669312672 1.3843369452138771791 -3.4521088800482724324 -10.341127066279092617 -0.9237862235506932862 -7.585242330571279723 8.0332956748560562232;1.2316330334121086132 -1.1272249912665432969 -1.2525521770398180266 1.4292944280585395678 6.1676141059995082117 0.86258340747209338861 -0.29024619239382659641 -0.37438229497517810929;0.63348415323616191763 -1.5389146575511856163 2.106009707653492935 0.48159716397956048306 -4.2433872095315035011 1.2625575813939444725 0.75627356299551928842 4.6231761799684978342;-3.1449558536895185767 1.8794025671939538213 -2.4620278410556113435 -0.24894362705596870255 5.3983193603143133643 -0.38797033353731863459 3.081086014138620488 -4.2787658278177023163;-0.032735010512377558933 -0.71015816737214587651 -1.8374391364173838781 -0.46114033548315602751 -5.0393618853886019693 1.0119414082667812504 1.2288844630440218797 -0.8427566743592169729;1.2643797273255956704 -0.017325292622556304201 -2.479869450398109354 -0.27053877932600095546 -3.5895039029773947448 1.6528407778472715695 -3.0014399567250857892 0.82517877401903061418;-0.4597719820726554607 -0.0081228417876122900598 -1.1041846946615501146 -1.3151877700733589638 0.012790472726043677215 0.92502629911277356456 0.059909293795345881783 -1.7005041628036463575;0.84905318119158790235 -0.38207679240048492142 1.6817071182253293671 -0.16397223546758973867 -2.7100546039891892924 -0.36046432350287982427 -0.3972433223235658506 -2.1724046408411190612;0.35776870284146450318 0.029301186606472698348 1.9004721465174636563 -0.53906914031559793887 -2.3282496693192120674 0.661197555386711322 0.00032847790516032579489 -0.74410240656329407205;1.3843744524834660226 -0.67151829288328968026 -2.5475525322765100178 -0.56144256925577340489 -1.6057975650868798745 -0.31841996788279686159 -7.5869327778974362531 2.3442804971454407337;1.9751777384151025796 -0.62863956386504615814 -1.8772135233569975377 1.046186113796489181 -0.33813902036906545723 -0.24696205844796748075 -0.99926680846714821538 1.659732542551727974;0.78357168850480196642 -0.25083353442042433112 -0.037418591938331219326 0.49949497028080053207 -2.045246135114837216 0.74552749097432124969 -0.91010662218140259139 -1.4862203627189545241;-0.35077802952860998031 0.15666506793673082298 0.199113930832505448 -0.14906860445078820576 -1.9011078888787642516 1.1931401514584540902 0.62076877932613627742 0.82397372388218526051;-0.015867846574785036962 0.033782894331827642564 0.93390405384198726946 0.38055028002324214897 10.282713932858186467 0.98680359485369217865 0.62721263717767705703 2.6827594721924508114;2.1971061422786308803 -0.084816431193109387277 -4.3344886134037281522 -3.6790637625536413147 1.4395384009393037417 4.0064193684094426828 -4.7422212739838585094 3.3262896980397034241;0.628831062709992894 -0.027886083623121019071 0.1161073069592155893 0.3907687066690490596 2.2366154206915647151 0.54828704216142609695 -0.9794404758931470889 -1.0651870652398571782;-0.34048839307434303203 -1.9659198689427435092 1.9041611996520566574 3.4428767856771540501 -3.0668062869115839675 5.6736007023323660192 2.9150983477422807155 -2.432771014725383818;0.82316141319537361465 -3.4675107048420663602 2.1367539130915718104 -2.0971641354127577017 -6.5838412979656935065 -0.91663145467362783236 -2.2431123250693545401 1.9708958387466071205;-0.57141439645507929868 -0.04403869286773711561 -1.1683605935323271297 -1.3188956831147831128 1.6835857505728808459 0.94847182109396732663 0.31027032139294991042 -1.3422006562779666972;-0.096532721613892813473 -0.30300475538627846817 0.88975332951121210989 -0.46146922547355712085 -4.9315158324247390098 3.3065417658979350257 0.73923587879576058146 0.15458718427081366564;1.1081315783571084044 -0.52546547752378114371 1.5194886077070841601 -2.0000673688926116078 1.4260183554548566676 -1.0101909188218891611 -1.2177877579198324032 0.61671541734085677522;0.7756648843212858857 0.11963355763733078685 0.6049671390204895971 -1.1939794478969527614 1.2927018053954113519 1.970625624974380985 -1.0814957182540523295 -0.35477751165883980589;1.0080733450220047409 -0.82688916692051983581 4.5370476205855858609 0.32366035327059033966 2.5765751981892033307 -0.629626133486481665 0.4532865020220240293 -0.37104581025640759329;-1.6283012302561747919 0.45403068168861321352 1.147096295918197395 0.89453528918861358044 1.5165107197367195813 0.17616074978555346098 -0.50797396301155006704 -2.8025891974319048217;-0.6017329424841457719 0.8113736529100485928 -1.38356165923726393 1.420348370318450959 2.0050631620499221874 0.12187358105083431759 -2.0529003684706310651 -0.25313516867810809119;-0.56403046575459958589 0.032670360986471989884 -1.9099424490007732569 -0.00055136153166961045929 -2.0909516490640007724 -1.0964478864356905774 0.7034829459630108639 0.040972342864707647636;-0.25360985518087197921 0.1435613975225170369 1.3934906837300413684 -1.1121599497489449426 -2.1883506976907658448 0.91483556548260658303 0.79906418236412546374 -1.2532347227072622164;0.21402941940894115724 0.076965276483141303876 -0.31490157836176868589 -0.33381740757897254701 1.2844324438972578584 2.2929919728857357519 -0.31567726826718822686 1.1944992882382972077;0.44992034229876209928 0.038725387329820080973 1.5946364768347656149 -0.19434139165016658946 -5.2941658771736666012 0.19554342180591516298 -0.60279289629494281488 -0.17502387785351786187;5.4752558366050720906 4.6982843733770387473 -0.2675903701642589616 0.51629420601494102616 -5.5019067896423594632 -3.5641377585105886716 -7.4065346358243671077 7.230278556359313491];

% Layer 2
b2 = [0.41163694742106399227;-2.4909205540787846545;1.1289428230578812951;2.2175278976916117379;0.56806139976869551855;-2.368805972716687247;-0.89582309545729743583;-6.0358960745389778069;2.0861685175191189145;4.1361466120408412195];
LW2_1 = [0.27641545140716800022 1.470485073709558721 0.48690660226106635688 -0.61321823893643223613 -1.6904045724822076568 0.62627579320377968486 -1.0363962921466443046 1.2777307914260089827 -0.53870431449163225235 -0.26604780756031326527 -0.21686610787181076487 1.9236827742879805125 -0.038801703878718479812 0.21078605786945284195 -0.27724222085392269399 -1.6198505068569442056 0.030785033828984786908 0.27749730033118541472 2.692882867313001416 0.513671804401141463 0.67285108636120938286 0.67323919568738677111 0.057047656259676324253 -0.16095350594972662739 0.33577660803200776174 0.24281581452316136627 -0.70129210307284617354 -0.055035873001684054717 0.6457320491074769464 -0.89783539960114289169;1.2042059532790143628 -2.8745603734475562163 -1.7243414261551239797 0.085962133686400571397 1.494777361621100642 0.39061362847102193419 -1.7685804821808384091 -0.62794446011416105158 -2.3673921019474057914 -0.036712391386692022621 0.96642050033539306231 2.868213190099620391 0.83033152661476161693 -0.6819655649511879858 0.35731725647716683358 -4.7107907064932312835 0.15267716279099477483 -1.4618708426781341636 0.46675505075187978532 -0.89277105985176818947 -3.7183170893330950157 3.3235588971501042188 0.76850506422665965101 -0.23150614050647863595 -0.64040916599191821668 -1.6016868495732303757 3.4046914153111864643 -0.81874799091277183916 -1.3549299699666246255 3.2096044921825472862;0.46632620601708818509 0.23446611108961951109 -0.053384550543468460571 0.47670199828854908919 0.056570370649457600576 0.68731219972410706465 0.80002196726475993671 0.89717910465255756858 1.2025111392858054415 0.0097040947282188566581 0.42950378812302864029 -2.6553269902424543325 -0.33639479711765651881 0.34758192637668883229 -1.1206962838353251577 2.680178220388129251 0.126075952274541353 -0.60876350465579243121 -0.67717029571185616277 1.6750441754589462384 -0.73702842898777121139 -0.19288592212298738549 -0.49556875869846389593 -0.092107035724607919103 -0.37837334805980954444 -1.4336351735077732528 -1.5050342505663398818 -0.5409315449537829279 -1.2682216353307091961 -1.499935675388083256;-0.49984067715558239131 -0.46736676130113391503 1.9311656360865629978 0.041656463629324555487 1.7386882638871092865 1.6971514934571851896 -1.6794408258864694794 1.3302875944205803105 -2.5137824066263330103 -0.46784001787925649785 -0.82035648738082556086 -2.0123439183391007568 -1.7397121240903636874 2.2005998774907573612 -1.2577838486603234447 4.0858447345843744714 0.10826988532319659075 -0.66786919533374766722 0.054744078183161334894 0.46412250395751847298 -0.24768363445315530469 1.90975366876719721 -0.69522060255688089558 0.825152669846977882 -0.40486095399050103794 3.22252253244470932 2.2269752575283847662 -2.2227490871304578235 1.9965973933721441469 -0.14000864297050052976;-0.14808399272891689669 -0.50872833317427890165 0.50367771738898259315 0.010494489795332226367 3.0273728256903393863 -0.81225252727113239626 -4.4109747580694200764 -3.6340486336286090818 4.9875168842881176801 -0.1436696202256394217 -0.53789194141131324578 3.0813201332380204533 -4.6375482861060621786 1.2756105716512282644 -0.43436492626306755227 -4.9406574181286986303 0.003682298858290541585 -0.041404947471868222553 4.3702299586157202782 -1.8809870676652218346 0.11709523226673047847 0.46058512556031644403 0.07134045930702276328 0.23699747768391307701 0.19300733810328107865 4.2081822733582212592 -0.58645615870332357122 2.1310349672566970547 4.9562378893130407675 4.2493654313181012938;-0.35884807013340580539 -2.011487307520646084 -0.4259879202789996766 0.58488813152159602193 0.94286364307978420829 0.17360442170088516423 2.4958942438392819341 -0.41891344635019162057 1.9798352789216349112 0.44266351586222463244 -0.71840606971490472699 -0.581845384263768306 0.65259429801707935503 -0.8898305398317731818 0.45922927317065603825 0.049437861245075265249 0.24572773791979182767 0.061379293865763018268 -3.036164233108692212 1.6313225211849586849 1.5858240200058837388 0.15859664481486174248 1.2703903281031536654 -0.12466444089556312502 3.5600224381617064928 -1.2932516528667035161 -2.1146182215573183782 -0.40291232613857425049 -1.6436954267579153566 -1.4203940795560880961;0.58619226885710118413 -0.095678677463479522292 -1.2078542071952824699 -0.75033531179261780952 -1.0658377368164373422 1.0457102066813226404 -2.7577194499206081346 1.8525031401303961864 -4.5907856026375517544 0.36381236914069914334 -0.90255032613259267293 2.0008064698727530129 -2.4495314584693090865 2.7261834953338963672 0.74859294532593900762 -5.0757049005371142059 0.15741966472870713711 -0.35008858419163824127 2.9848890505440035525 -0.15018783548735512023 -0.23899614344050693071 -1.0366150278323504885 -0.29807444631283147585 1.7453118794561290894 -0.95959106970457330732 0.22644857044011210134 3.4188497927346741534 2.1549434127307010023 2.2523523414409418741 1.4323748477741784502;-0.60947607991067098698 -1.2354221637929621025 -0.42046134463708256801 -0.12364649433188740213 -0.76781957095245345624 -0.91780857827561279638 -0.31466640472046730759 -1.0652883256152261282 0.13057682758432065384 0.10100038248296529242 -0.17518388415760521859 0.97395194324619194237 -1.9106374723498749102 0.29361361946828934411 0.65470353726269125971 -2.0120837140602390036 0.53672970336646297351 1.166925714601948938 -0.53693780556485515287 -0.55427888913986955721 6.1460771798108906339 1.0755142172679885704 -0.64973941823947689578 4.1518369505983532264 0.25767231779663657409 2.2762071951995745955 1.4532682736820281555 -0.098300248819650687393 0.17877943039860588748 4.8046434503552761797;-1.4142432350227263083 -3.1831630793885627995 -1.701573141741551165 -1.3788129982092149906 -2.9017785236413948802 -0.90542932432737532888 -4.7830004984980005744 1.8218084273236647697 0.082487323409108206929 1.3080224379259199896 0.20898780440142117332 -0.13339584349605865254 -3.0883250251961271715 -2.2630541643737287494 1.3912492613670057739 -2.4056382114479144541 1.3736934984340596344 0.41471990784841511468 5.3550097435232624221 4.3338806490248424907 0.1589574493346337658 -0.1858234103082617994 -0.61845568544148110668 -0.53105287328413508519 -0.1046313365743061774 -0.12348656617415017245 -0.39573447532419292472 0.33601719960738346549 -0.89837492113024763896 -0.0062296847311391380375;0.24523716629156278035 0.37807275836344494957 2.3864525408645627103 0.32236307342498615736 0.26762252766689226258 0.75420135670200605382 1.4075722394434115881 1.2437601065932788469 0.54626970606329583457 -0.023117950499964276534 0.35297752890462563702 -2.4119970848137080033 1.3753590682483063379 0.25416599867796485501 -1.005513991556501896 3.2896170263159456582 -0.11087751292405180481 -0.78874425775397583038 -0.87867190074204160055 1.0812461800048855931 0.17539296585343686075 -0.69742887906136419307 -0.079599397092109142049 -1.9023785338292868996 -0.5400755681896263205 -1.4869898656712570162 -1.4800767314328417967 -0.33085705468652959471 -1.7389459858018272431 -2.2583232444649827819];

% Layer 3
b3 = 0.36774929542106599145;
LW3_2 = [-0.079450594143737160446 0.12067150982945305804 0.38986562321634377071 -0.054195928775147297896 1.0540450597605082184 -0.12805858588093210759 0.04504211136867870835 -0.15905603716478330933 -0.047144551122860492431 -0.47258496788985254744];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.0145232735458572;
y1_step1.xoffset = -0.52;

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
