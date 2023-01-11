function [Y,Xf,Af] = ESPER_silicate_5_Atl_2(X,~,~)
%ESPER_SILICATE_5_ATL_2 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:39.
% 
% [Y] = ESPER_silicate_5_Atl_2(X,~,~) takes these arguments:
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
b1 = [-3.0245712975509340126;-6.4072759077948937545;-0.19368667575569153971;-2.4569977727540623924;1.9972974959274394369;0.068342917251925122679;-0.22226129555337204025;-4.5710628147894443174;-5.0768523636142095512;0.30381319169214415776;-0.99468676027297164843;2.4270050406998855408;0.70469923589406624487;-0.85152797913983846456;1.395406995650910309;0.24521475048219448789;5.1447421906413230985;1.430690583458037235;-3.305161512656864975;-0.8849193172614039371];
IW1_1 = [0.66418710052499041474 -1.2223830131430688528 2.1915872027903202657 -1.5783879980130501863 -2.0592413769624475783 0.41767004083570258288 -1.1672952629863875273 -1.1606617158883689722;0.24847064135607360713 0.076941813458491672972 0.55419334651094320776 0.29558324236223915138 2.0243233851518431798 -4.0701123759054604889 -1.7465699074025806681 -1.1882322145224275545;0.99979523330092168187 -0.51903304894310742235 0.63468598033028333916 1.1301008295584478169 0.33006602103900795075 -3.0006439206896615701 0.32580764377669246823 2.13492966164161313;1.2652739222292024301 -0.38405304689711017296 0.89090336847719686642 0.2739015478218334021 1.4452327908560311798 -0.86263647725434111635 0.3635120867000022904 0.64505565280799892491;-0.69814073477020111014 0.26974752489626796503 0.54101176700494990612 -0.36000169599423526012 -0.48819634271744444254 3.366205233919275841 -1.273861377307756948 0.32064937335124665818;0.23639403202147635241 -1.1726302139761375365 0.90436145612706675401 -1.1269644083055949668 -3.0070752965028475501 0.17899115978481094924 -0.64583760890222152717 0.16546163108711253731;-0.62165059325087457243 0.032788943586134253194 -1.2492844086258632075 -0.0031555827129659012528 -0.20615517416662126893 -0.90885248364294879408 -0.0447392010024404152 0.009014781492484367012;1.329041006170476491 -0.19840992788624209342 2.5654178461482928242 0.22892624331558975026 4.2342612681289182319 -0.57573198763668576028 0.40389501998067828259 2.4161177201588355246;-0.28673423952492649303 0.11609580206213325482 2.4832596358412968307 -0.93417143015090275959 3.7126678249059810177 -0.022173852554213332317 -2.2520566124344085779 -0.17340804901098078039;0.11290829406454247485 -0.006413436459691361155 0.15077605801074026037 -0.11535843419613918404 0.2023979649341799969 0.13746118674466217824 0.13342283111171279741 -0.10366572899742960978;1.9654972269024537646 -0.5357687577262263412 1.8741389836019088921 0.34301448447047705459 2.8508490092827774021 0.051953613517567234781 -0.23520005684219624786 -0.72071478902116781207;0.75887792574297086023 -0.63755075566037211487 -2.2622537442310002298 -0.18953406772106198841 -3.5770806885947719067 2.0933263317642274437 -0.58884100358322732927 1.2599894137373937042;0.8855089292037437998 -0.16495141082743919325 0.82156475887504665234 0.47777991693413507379 0.19422580248352844245 1.852449640647594542 -1.2796971294548942222 1.4897366905189193265;1.21920945376992651 -0.25311400568072128792 1.1747198202761515162 0.10663589512466334353 0.51881936244694470961 0.13458040046273997414 0.14400702733706641223 -0.020649328332849915812;1.0103489096533562197 -0.095269684727702766414 0.5646549876816155944 0.091832031617370260546 0.023794345479039476221 0.081826921662080484565 1.0725647645573106459 0.61691581346230006666;-0.11225198956356227664 0.036984439969021995487 -0.05397886265453594179 0.030974746891857460329 -0.60379034654526064774 -0.30328689636451994227 0.015994190791914504601 -0.18473525799773665712;0.025104505909700161603 0.57724831319822522691 0.064473779147981247206 -0.81733641716276805855 -5.7955857333360594197 0.98770109010364681712 0.64321819745143382629 -1.8865228097378183847;0.11391629916701555281 -0.16556531846606895031 -1.4821462133751279744 0.38286670403953859676 -0.19065899077506093251 0.49279651286939207822 0.91311484129131914589 1.0576779754603429229;-10.265635885969370378 18.102120852112083327 -18.990324367609950684 36.176454717709823683 -21.681274116203603342 9.4791510958254097119 6.403140637173962979 -5.5428157366234547254;-0.25607373290657664544 0.79485360368168567202 0.86139334044279469715 -0.83528856619640434733 -1.8896211125628334937 -3.5237108331919211324 -1.7749456100967868633 2.1048795050760027436];

% Layer 2
b2 = [-1.6678346841907907283;3.9058244925216465404;-5.5533748474516198002;-0.96940725393485005501;3.9904329775133571268;17.087710183311632051;41.025170136819994582;-6.6913302771067701968;-1.9992018703348775066;2.5054847859870403859;-3.1601117505623683712;16.655794387277165214;10.213241098711145938;-4.8211800445137376059;-3.2837466578891598701;4.7267656969775115883;-13.680444911627997584;5.7547525440064193347;-4.2833386378863878008;-22.95776414487069772];
LW2_1 = [0.25055798858389288775 1.2652500029323827668 0.89485977660426430091 -2.507557053636006561 3.7841997756877243653 0.51356146034303018499 -7.1004961402901178857 0.11769448454439894458 -1.1128118030142313266 1.4832146139287960551 1.2217237860725962406 1.5794539529326632188 0.3815439493085951983 -2.2721362504346136468 1.1249582872620711438 -1.5502041121768337373 0.75578185892430860271 1.449836052724545965 -4.856522236013070426 -1.9414469238641998139;0.39928286829857162488 4.2129551809856060629 1.5440605180523203632 -0.73803413665477246575 0.77303438918970679428 0.28053508982557257001 -0.19157222917843844123 0.43493875985448171573 -1.197815453124416063 5.8586965238061559091 0.23145894352177462383 0.090374338912740950769 1.7275474042179350675 -1.7221768470938019213 -2.5635944420366265639 -1.0910592946984740692 0.68185383022744328052 -0.50672217270959785118 0.04058070783310371854 -1.0474704204862368417;-2.5729142219384186774 -15.304693026619418106 5.6907577861268405783 -0.45407402393877094315 -7.0956725528465991459 0.10778758274228944503 10.073482959877171794 -0.43124516083433400704 -1.9915539230940750848 7.3070244972049129828 4.8624429231915389238 3.3202687460066413472 3.8714919663308147157 2.2503848689789762716 2.8282271650018029341 19.834065998519449181 -6.6218668258844228802 -3.5990978635123296314 2.8585329859713679568 -0.21100205313678424557;-0.8647088028379555924 7.9431079896300227361 5.1112825894848992192 3.9827169009120946264 7.8124304872976368941 3.2787764358265558684 2.980502685540676211 -3.0114129548548818072 -0.6525466114431788256 4.4665117515098016199 -4.4192614882743450622 -5.9169713983632696852 -0.81945384724380032626 0.6491930583916234454 5.971002599622948992 3.7124698848329558132 -1.178890712602910229 0.92554129467052492775 -5.0586300425922523516 -8.2576011079771198808;-6.9185810296357193039 6.4002595844585465912 -2.2532940346916481289 -1.8584990464714785663 3.0023943717411674115 2.8202923618916870829 5.9355821558859345188 2.0358736830311308097 2.1744705233986603155 -0.11582653189442509989 -2.0064224789856237585 0.11985062003329816083 -0.72019641390354172117 5.2153974023509759306 -1.1996063738089899342 7.5038741173446439348 0.062299607897658655642 -1.8378531058526603736 -10.104315771566870552 -1.6886521058258898709;1.2494766654915743231 7.0169315520956319077 2.4219567523701446277 -21.275968020192866703 7.7267330521031434287 -13.21794974836160641 -6.5936865245538873737 10.421247779522774479 3.8355052764045631797 -1.0075442846424293997 17.327408745973951198 -0.92859402618285091791 -20.954470049622443639 2.6513144279256701807 -5.4472601004401655089 -1.8147959599116352258 1.6560412901315118006 9.9245342221569146091 22.825121039136831769 -1.2214300168625435994;1.5007393207692873993 -0.6018401162061658427 -0.54643481049115971704 5.0603356897563944727 -0.95241074746093357106 -0.57789355801084130171 -4.8073341202517232418 -0.54017570132463144095 -0.56094974936920805053 -17.496391088842464967 0.64258780499433487687 -0.36117713068869355952 -0.032056981770340838012 -7.4165344240619477389 5.0663380890594424955 -3.0991687887565877446 2.2751346731974591009 1.2316156662822281831 40.586900599970711312 1.0292433915022514146;0.12830278357129126787 -0.55813004484989992005 2.5970192846189448233 0.59964345284521747281 -0.74955777053676586608 -0.36646115087024627721 -5.2556082503368246606 0.061145573970743209713 4.1054936358187354983 3.0084002208901701181 -1.2128233896214974941 0.67544103423109735918 1.0046340678431060489 -2.9255485236997116338 0.48058788290061937865 12.750089239102809557 0.31432456076982057125 6.7508207347493831563 -2.0256082670244430766 0.33238384563478501343;-1.9911255036933737106 1.1931220670998159417 1.7986431417544228584 -2.40657337879339428 -0.62896860874613347292 0.48715253554013510096 2.1252297655391987341 -0.68965874595200826747 4.4354664462182498141 2.3350265650972090725 -0.28783850059273491384 0.36780038504480283557 -1.2581248493494117113 7.1248134740322566216 -2.3963647843487150269 8.7213449699734582055 -3.8952236143032958893 6.3053865557573258727 1.6484492333094558081 1.6473520065443911786;2.0789291754350198715 2.0156143278536085006 -1.7410648803580746069 2.1533350065933065665 1.8211144177501337182 -0.7279012995724607471 -1.5292869559025652482 -1.3359191289110599055 0.34479109755410580762 24.835422335732069143 -0.98549453636712769367 3.1711360945183213467 -3.7776559944667922331 -0.40237670998791563726 -1.533370007490167275 -10.392988687974291651 -2.8733107486081559401 5.155890263397449047 -0.95280392703845995239 1.3196559037777455536;3.7702373538875511016 -2.9189570880646900086 -4.49094769969828711 2.7044643855446515168 -1.2658480212406753029 -8.2612851349505884713 9.0009286945187660933 -2.13432403386142866 6.1225636739465869596 5.7637457137622831738 1.4766088261546137517 -6.2185968027039031369 -3.0925002164702872776 -0.24206438555085374453 2.3849600611291816854 -0.8901043406982935835 -3.4629985126019100505 -2.829251156503162612 -2.8143976395314824757 -1.5370968526671608956;0.11119626601356372153 0.67086813637116737841 -0.8508369848493145815 1.2104052424623150674 1.4399405200759756163 0.14555429656691834528 0.87355194191827911787 0.14815117126755183996 -0.90908301673713032276 -2.8861016666679049081 1.8183573898669262991 1.3870847846503961964 -0.17709885399782976378 0.42904246752637698048 -3.1340956354829270225 -2.3290658673162698555 0.19014633800894839566 -1.9939099489462694947 13.376604017599349916 1.0465178166533353199;1.9110141239725741968 -0.63125913347879236692 1.0828538303743799887 6.2784519500003668213 -3.37783493123943801 3.5819918103623669126 5.5064076047941590275 -6.3512607048348090899 3.1531272299169756934 0.31116273587458603833 5.3847750286198161618 -0.93715642921076236505 -2.4577275877598716569 -6.8986773729665440413 0.032337783418676831204 -8.3361069803100349418 0.11401193461550281072 2.5964011370008335255 -3.2276624888935492308 4.2398212651666451478;-4.5787857154315041441 15.382275659236160692 28.221632108876825384 5.2056411417002914632 30.645706332517413273 6.3478648287612946532 6.658122639234120399 12.246521882036216056 -0.9825077559697950047 4.0631798711832871618 -0.2825878361891404067 16.032027945049794226 0.69628919321423854427 4.1465916089014314139 -0.93783223755159161872 7.4611033216554716319 -0.79837061201694814461 -2.5577034326194367786 -17.816580660342197007 -1.7753294430284680683;-0.4050905750066581823 -0.11334553966852034534 0.58874976754452812422 -6.6102915295227981929 0.90145054019397785616 -0.57557508995068495938 2.981246170175094079 1.5815026928612243573 0.57438035266108877774 5.8334450882645256797 -1.5197258046990818059 -0.65592005060195690547 1.6916472200731051956 2.4084493870464904752 -0.77690482963948748196 -2.8011072594280910231 -0.88920086187219993246 -2.3433439599069689763 0.97319422221817641905 -1.3795311669472631877;-1.9271059857308763252 -0.3376923346038214202 3.3635045796959261999 -1.211416870197916662 -1.2539161953102775371 0.72561094005534965135 3.467919320322504273 1.9071092358098458153 -2.2472022194343770352 -20.399071927916796909 3.271910867646638188 3.3581366951030404699 0.59349015568329221981 1.6787437545048902177 4.5856127709387939362 8.4567218293405854723 4.5234470903438319667 -7.0776485056773674742 4.9295894091675913629 -1.8062784463449899075;-1.4791768389187676114 -6.8182902663935003318 -3.9613684335711050188 6.6626901387605501625 5.8034022056529899203 6.5070059857399638048 -8.2804940634632053786 -1.8375526531007848607 -0.53177884691830001884 38.904045658602747437 3.6019280294114741636 1.4416684802399490906 -1.933120076765673323 -2.7726685926806617921 -2.1124339385754282716 -6.1820847926936242089 -2.4081820642142552913 12.674318057169490004 -7.5254390988906116888 4.1779238331754848801;-1.4025327040212203222 -1.4632285678217984426 0.82081330082456915864 -0.59999772029662368222 1.3152944642725765689 -0.46295888573156923629 4.6513205135426796133 2.4949873778727686968 0.32492460668607364704 -19.374644721584676432 1.5809955497162071492 -1.5527093585090654404 0.79479916096898317601 3.3632028616716542757 2.7128873589035054081 18.390631400325350597 4.0709682513946887639 2.8869088498383015384 5.8702226353225626099 -0.39798557891094538119;1.3428943990326349667 0.67966702138561185809 -0.57951733431287455467 5.0764221027888902782 -3.4619563922760212549 1.4786234260331525547 -4.4707001632862093388 -4.2247486546222301484 2.9195515952900841228 3.3256099002819681232 -1.4819864320656732648 -1.1847963862154780035 3.7448094090650370092 -3.7728123643575437463 0.1530937257647022598 1.5357404572704587142 3.6532741007392992572 6.7971613883345947116 -8.1403112044733276775 -2.4311231960004247021;-1.6144151398889841076 2.3494273198914670253 -2.4988121458530763164 4.3790621310400874577 2.191442597713809004 1.7585249147772605838 0.6244256051862271395 -2.3590779426659609364 -0.6687154338651140062 2.2214676672110900135 3.8265539800803187731 -1.6314317258605230965 -5.0926611954878158173 1.7137941044641480381 -1.8885781712070135452 -6.7466534750150248811 -2.8811144319350732168 -3.9392079056731716946 -23.404493328453575884 0.65299661702655509288];

% Layer 3
b3 = 0.14392867332151348947;
LW3_2 = [-0.38513816530369954227 -0.64946283397269111148 0.0017589518133932656365 0.009226572257850380529 0.15122817572677119125 -0.0047783571288466386573 0.17870575109997891006 0.081715671492873079296 -0.067902813463462954879 -0.067503871050544311805 0.010762575601557461749 -0.050317180862798689678 0.012641934827806839206 -0.0054966639008387808973 -0.024184929057614230202 -0.017513762912987829862 0.0062659595377633379651 -0.084987694363525836949 -0.035835893701420366841 -0.020613108191875423314];

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
