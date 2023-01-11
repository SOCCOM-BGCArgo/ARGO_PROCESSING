function [Y,Xf,Af] = ESPER_nitrate_10_Atl_3(X,~,~)
%ESPER_NITRATE_10_ATL_3 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:35.
% 
% [Y] = ESPER_nitrate_10_Atl_3(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999994516886;-40;0;9.1619;-0.03;-0.52];
x1_step1.gain = [1.0000000017546;1.00003626070655;0.0153846153846154;0.000309885342423303;0.0656325902039552;0.604229607250755;0.0145232735458572];
x1_step1.ymin = -1;

% Layer 1
b1 = [-1.240300160602192836;-12.258064857016105265;1.6726429277261973816;-2.4178553608858468493;1.0231239857289855522;-3.6870490993466562557;-4.2666251295906034713;-0.12211930143036418828;0.81992974471657253943;0.30610645209000958422;-2.6198187991519694684;0.38581149283564153407;-0.60639968683612099909;-1.1943775798527302978;1.6808630746807853296;-0.94554066302993966531;1.5138642968892923601;-3.0079073392823421074;3.9778091286904522228;-3.5585694733725956596;1.3357196139726117412;2.506869621607865195;-0.47358792646795677683;-7.1222640000523975701;4.0769422529883145856];
IW1_1 = [2.1588824833888353716 -1.743089993645057012 -2.8026561357802184382 -0.36342000887163755651 -5.6064920385789758228 -0.30218443522015953073 1.4737494267961386285;2.8207380050670436411 0.15036817071734151496 2.4048289440545480566 -0.47183731969732267686 -0.4187905967126638318 -4.9516311706759426059 -4.8534083017287317219;-0.23042825561688015656 0.76050646301468827204 -0.51401621037558198335 1.9967293223658355572 2.2907183948733766954 -0.49996411923192229931 1.4726504318233675228;1.0317851471194146029 0.22821584954183954386 -0.61418989044276350864 -1.8505530997912491564 3.7480836763179068427 0.05561140644768993585 3.0655758270896353856;0.61186548459061307348 -0.38058936224518113978 -2.2886905687415088906 0.8410557575243590156 -1.3836254342985465637 1.119781745153255903 0.22020284432570075284;0.19361689723453448519 -0.099455376357841379642 0.0779392630588855162 -0.59451781373462919689 -0.71849649310748964215 -0.56333226653891410596 -2.1531877123743483615;0.97719865882234169252 -0.13118276206523515892 -0.89879764005259332738 0.024342441849351512118 7.749586168668095354 -0.93069865670131524116 0.71849067964076673842;-0.42100279101377890978 0.093846894548212878018 -0.42919223205204692917 0.58606091365046319996 3.2256611457218942185 0.086974103505911851131 0.97419676402892563249;2.0036996487147749058 -1.3345556161456673117 -1.3516356509833580279 1.0045828587144305377 -2.1162910479859178281 -3.3976897213641992224 3.1767871718417954874;-0.016070921217200433151 -2.798977408043629822 1.847900848182217004 0.48303339283160640072 0.95265322121306184044 1.378835450065065471 -1.7258093507211837991;0.050087775916692819922 -0.5448480738415372171 -1.888989742340700273 -2.7474006176083305419 4.0313210386508924188 2.2946463495660864851 0.2259458709573786972;-0.60528821202093596732 0.24071478337330892705 1.102080336387948778 1.6468303116349001236 -0.1505535954612207683 -1.2736753736930421788 -0.14036973447808429039;-0.49952957546321080651 -0.14812115142854220817 -1.7197187298620746354 0.1908672940534035245 2.2149559364752784418 0.3019076382418331983 -0.87219808840942070205;-0.25626137152478223324 -1.129776731924158506 -1.0952177946604089165 1.1986451414607677268 2.7848381957751064952 -0.24517654262138735555 0.1198885364252259278;0.62370659690749974402 0.20635639361973892592 -0.13661402366875841086 2.3679000765122495409 3.9906165825100137212 2.7028187197884698811 -0.015525511088714454705;-0.8507258903541613515 0.33004196730511320412 1.9454232564479128076 -1.0687377785329637003 -1.0764848464951530804 -0.52662025580099580679 2.5097023044837722949;1.5642124395972636375 0.49474228959209420875 2.1176940338989540535 0.41058428686467041135 -1.4451149661610953601 -1.1253820791419757441 0.21639541727785469027;0.039693754434769232264 0.21007040261570822381 -0.14895919607360363912 0.17785280313083814963 4.0875355226680278875 -0.42141921549804034486 0.3317405234535703018;-0.76092605701868920587 0.43043142138177753475 3.122326880228134538 1.025347377822737549 -3.3138456432663563334 0.61513208280304543241 -1.9589235159091211713;0.16688415863085234725 0.34279704569948565851 1.4639497688955140919 1.4529443972196460777 4.6958951354847515702 2.7190689603345106207 -2.892508602830884179;0.36418524873887736515 -0.22957721927989341304 0.41212308755839738295 -0.28320178722709338226 -0.89038461885864317313 1.3082482199027163361 -0.81826564206761509546;0.81824738656926698788 1.0007670941746491522 1.2022881386444239649 1.1700872310436170043 -2.704696116549703877 0.17840721960470987129 -1.1900597577855718257;-0.77419342888013154003 -0.21705646986274271248 0.44869921702546788911 -0.0086506074556290286914 0.17421745034851895961 -0.35664819825577603885 0.95996207378963072099;1.1192518770378310222 -1.8947055209911769502 2.5158152800103765756 -0.60217088826494413745 -1.5891660727124099495 -0.52550192969160924861 -5.5565438862122169539;1.0580097168066469049 -0.2470670903838582122 -0.77432809472648844729 0.29293714507716517303 -2.706175992023256871 -3.5278771162818509666 3.1704817306554287093];

% Layer 2
b2 = [7.3988914699717254209;-7.1451270122119971973;3.6611381015864177257;-4.0538913582752469367;-9.6340640952244260831;5.1606930487782243233;2.8357799691458187574;5.0075003601987457813;2.5495461979706126598;0.91884858710863115405;-4.2885305486108160622;1.3491292228248839002;1.6652615715016312326;-14.310051703854325567;-0.37202969440112171329];
LW2_1 = [3.8856205480484060644 -0.49906112603138069739 -0.94472153852796603157 0.054134146925622239621 1.4020089503427555222 3.2060492738215695674 -3.1410235330012143962 5.475615475775135188 1.4640697217102434102 -0.091885232869900843045 -1.1421935191816297017 0.50581768617722022707 -7.3485035216354157583 -2.0676140667172933085 0.96024516539321680408 -5.5061203783824286262 -0.33406208988804714233 -1.0289022898165036679 1.2320375038263449152 1.0005048556344193766 0.24852481993915320757 -0.55263694508668503858 -0.50500792006064099127 -0.36941407466984294006 3.1346624199499113494;0.57246052359346821792 -1.5244619655218911713 0.25938613122118581167 0.073585621661645742031 0.69367325352951914041 -2.8474887272433568874 0.9762204105252293429 0.14537694995836575318 -0.41781789286401666006 0.053212232998451125532 -0.27092562978101852877 -0.037792391107611268686 -0.176055144251424045 -0.057840637710127727755 1.4657203303549619644 0.66839099281145675224 -0.16507470902956208891 0.01701584110928918564 -0.14423409932297442948 -0.18818888056194874037 0.67797411218879244377 0.29581432859610440333 -0.002353798667525355999 -0.48853835332115863599 -0.34857786818306973897;0.50160943216724507288 0.09517425991717455569 0.91705104514216262412 0.17590363781310555669 0.65371898882585122159 0.73868368738191902967 -0.50258178957872601256 -0.44050016008858139394 -0.25685162986180282152 -0.13832773440783124874 -0.14271192009907257559 -0.6682649194998793174 -0.44854089711109795813 0.48167540850990409051 -0.47442016202170200767 0.3669029049683163457 -0.8412978492032943123 0.59006053695581794916 0.70745802658222489701 0.12431080949078465725 0.66153834175256887029 0.56265803252867363504 -0.15990271140051329213 0.37372594445457124612 -0.53399361048522542283;0.59513640274150425569 0.3422471452879696785 0.39635789027146722807 0.42513059215556853188 -1.3474995307597585903 -1.7531284301826279837 2.6469577135346451513 0.46262972524247691908 0.2145175405544411662 -0.31462241935717921715 0.66210744276222821547 2.4896258968981963555 2.99570076965469001 -1.6164018726802789061 0.092494549255809313526 -0.2929862369332226768 1.4764532136239412186 -3.992484062032910952 -1.6461532775124814165 0.36006897837179480115 -2.5488312537757744636 -0.90712445992849732868 -1.4126957738677361487 0.45281376935475964318 -0.6442165703612964256;-1.8965194481774516611 0.97178394632576503565 1.7150779882355540185 -1.1863110950930741971 -1.3030404601179683688 -0.79501257857521745454 2.4487156257367979073 -9.8869370468620640935 -0.20678647918459755473 -0.34255083008200665207 1.9761441926722715934 2.7149825658520487792 3.0477808254648781627 0.40699491090024386697 1.4421563381581701258 -9.6193343815668850283 -0.90056775254819121379 2.5929195373587208273 -2.3988469755724519672 0.41005215389980198548 -0.39353089013725761625 -1.1467323328973073604 -1.7370699679826837603 0.89674984866786733306 -2.9621980456302248186;-0.016302207546023260443 0.87746194244142594609 -0.17436153158026640453 -0.55194364926097161383 -0.81001259800073988071 1.5906048484183148428 -0.73461189666134318887 2.0261548879898776399 1.0250847763204573582 -0.41027581511601729503 0.48201634781902652493 -1.7134039059499484026 -2.5104207365565831545 0.67710933612068846532 -1.0582954462372264892 -1.0892733836818790927 0.96717801825118687731 0.70193264231244578699 -0.98557845990190462437 2.0407669742983673977 -2.1436310177477815309 0.47923832039138830607 3.4568471974395040824 -0.57661057859969100381 0.144132860671766172;0.61901691542645587152 0.046582089930518703891 0.69101720573739267017 -0.033925919085738727443 -0.47321959575851063606 1.6347593293188524832 -0.20664665896808423473 -0.38186628404394235003 0.027982221324064256163 -0.25813741479598534267 0.15682876038052975809 -0.32007885652422152223 0.071633610663207378244 0.32613247948506607354 -0.36268891834077277014 0.26973352490999014108 -0.67324206358221416746 -0.16654038531870468676 -0.21897645825445993806 0.064241763866253215332 1.0911523730259400367 0.8029048536633943689 0.44855087699244750032 0.23108905489123196531 -0.30506642306583342528;-0.1100800768361043791 -0.36481729276536650763 0.14749367895673728968 0.15474998093130495902 -0.10548589969898329588 4.4507104891068145491 -1.218514063538301162 0.71288415895332324368 0.4622239407721255855 -0.061422864773251997783 0.27623081544569799695 0.16554155709921755668 -3.3178531704584117712 -0.33781016111785316935 1.0200955886548823681 -3.0971763941964405298 0.61427455300866118382 -1.7071591377052908189 0.026878308553822768956 0.53277322553368244851 -2.0142974125612642666 -0.52506235064238793164 0.51995260851904556709 -0.57975399539802974225 0.59020803395011389725;-3.0121157835343601583 0.46971978280330117794 1.8821907333950249619 2.2725786453000429432 1.0514068396540219297 15.435382731839808912 -0.13717647107337868628 6.3927830488892434246 -0.84145390277236609133 -1.2999979162057926363 0.45593720187885156125 -1.1350045169526890021 -0.27056434047285493616 -1.0967047219821448589 0.73244437062076139799 -3.4558104406329599456 4.9740404303143126441 -8.915224380172192653 2.4566006548510337382 1.5861536827009408146 -9.5424488597047520955 -1.2835728875131957594 0.31700774605086734503 2.5998235950131198813 1.3520203132012742486;2.2394420739852711755 -0.31126731071970770159 3.536600460416258862 -1.1503465130092171975 -3.2401143845824234013 -5.8818102530769440506 2.7035931462837052841 -3.4596555684591496771 0.34739649990179233274 -0.55121032226828559075 -0.59877792527776085851 -0.30339243286461392568 -0.099715672265362043092 6.7869118728278365893 0.30164117000442602912 -2.9182766905261190438 0.26095938872107810047 -2.9255806224325588971 -1.3731099674701747748 1.2569430480125041161 -3.2832922373804978022 3.7163197514456682846 5.8299276591820712312 0.23624201258859134356 -1.7580999801188648135;-0.69611698303291547063 1.6535274124911192306 1.1320022294488936154 2.179209054397710954 -1.3821507795974365962 -7.3263244293729066214 -1.7743125199788329915 -1.3389999939565904175 0.92709056289196578593 0.40169449000828205687 1.6565462452700556728 2.3920546398472994909 -1.874911067910343565 0.83700471501517625939 -1.1435233214231572685 -1.0747567743093247028 1.8004068592648001523 -2.7229621383903674925 -1.1225173572018143364 2.5833660903401582942 -2.8680361017448245953 -1.4124242226184804405 -2.1599455514123775401 -0.0069893807045705675859 -0.012056798086258101871;0.49010231222579420107 1.9558113566182848686 -0.12309554065049321814 -4.1959566385454181159 4.1893421568004587385 -7.4902407406640056919 0.58055774600307041844 1.2355358440092276684 -1.5395134360976441279 -4.0237860510775398382 0.658537193588812797 -1.8088031776309203558 -0.76831608826732877837 2.1779922611310937874 0.70239789404058494693 0.55221960539144443469 1.5849281408743449706 3.1296000255862765904 -1.4178785893818699471 -2.7428399915620720328 0.40889189959681959685 3.3804585413071990274 1.4979316180775774647 0.35989656112076967576 -1.4628942511234366286;0.36491775422993866229 0.66116334847262236973 -0.2947686279604092574 -0.16568503522836533493 0.33426591594502663707 2.2659361002741760238 0.071066106083660707249 1.6194837470836858095 -0.28852248634715921272 -0.10014188222834029263 0.13415915033266384571 0.21965363868885703913 0.009099804580656546682 0.3344760072245926974 -0.29578577378385201291 0.54244753904691267454 -0.44605765861923452054 0.70610501950541659788 -0.1661413113568057931 0.13056420712816738106 0.43747238505212687754 -0.092781785455277057673 -0.5932656540100230691 -0.20032473741496586994 -0.34423741324646933704;2.3358105027866216119 1.0246022057747019574 0.56233703317098215901 1.2631631789515127551 1.5547564001862397287 -6.976480244154571686 2.930036802094969417 -1.577504130052755027 -0.21604324876251695953 -1.6011561932846938294 -0.65895209667663778852 -0.78597392768128016005 3.9822076790192570961 -1.2993335357650401818 -2.1724440467922052633 -4.9759042453255686311 2.8209463052814367501 0.70291166080582134779 0.0031602993885430762794 1.6997383967509289704 9.9257536412660147818 0.17237381079263136563 9.0891504839699877039 3.057577449939721248 2.2671133516490864501;-0.35121981177835859267 -0.32219104996772623961 -0.073570098322996588247 -0.12281885005821446821 -0.65611646271199419722 -0.88064853035940637849 -0.20081503866274186265 -0.98392316630296172697 0.2088856382865462713 -0.08466693163169783376 0.19084237416048269176 -0.37412688503110419491 -0.13330535147137875795 0.27292199621341162175 0.13044269238586753135 -0.073563085919607348884 -0.34399516815098013511 0.024932952771296375272 -0.21044360000862316618 -0.10824981953525221468 1.4262509877522375756 0.66888217460776533407 1.1058241914318929222 0.26454526538428690952 -0.1470250187784349627];

% Layer 3
b3 = -0.57332015727190055454;
LW3_2 = [-0.32810294969874886917 1.198420127028953841 1.046235329146268267 0.23886278015617129489 -0.22708957747268176952 0.16310219661194078067 -1.2127998654428755465 1.1994342150892487453 0.39447251615624007171 0.075551726214732342624 -0.10520694226230965784 0.10782487797425220621 0.81919672104992913297 -0.23079005914487801809 1.2929899885381388458];

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
