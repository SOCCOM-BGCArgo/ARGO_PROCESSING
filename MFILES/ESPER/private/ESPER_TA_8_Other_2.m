function [Y,Xf,Af] = ESPER_TA_8_Other_2(X,~,~)
%ESPER_TA_8_OTHER_2 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:02.
% 
% [Y] = ESPER_TA_8_Other_2(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999540494077341;-75.751;0;29.715;-2.05132667284751];
x1_step1.gain = [1.00000001035699;1.00022983623563;0.0140762793578401;0.000302617642608564;0.251952632905014;0.0584431850247228];
x1_step1.ymin = -1;

% Layer 1
b1 = [-3.3138251026386265075;-8.2577172528411537655;-0.68626739865953045872;1.3538662970000636676;-0.82854053165318641838;-0.84387923146907839467;5.3470703603987885799;-0.81899475710910940585;0.88922217278754323022;0.56672179504299557973;2.222994700293464021;-0.96595524966034074232;-0.73947074417617364439;-1.7470240945181385595;1.6263332558304848252;-2.2639600239755033328;3.2090263537899428314;2.1583041824441506762;-0.5263054787540993873;10.947994515873752164];
IW1_1 = [1.0397362043449229763 -0.15622314036865192488 1.1765216155351911897 -3.7967878926938638529 -2.3787339045972171192 0.99094541972652272488;-0.99538515420430162539 -0.080524084325955258556 -2.2699625892876045263 1.2644939421934324741 6.0547066200225518884 -6.4197547280007833947;1.0742476017772659702 -2.2881005610519027371 0.63959943935164154105 1.8564759195196056663 -1.1789833664895710275 0.4239146342233787168;-2.2158647658304557737 -0.51292764640209176363 -3.8414576754046603213 -2.1321475928977076109 0.74766418796885758535 -1.7045922072429668948;0.88164467071021757061 -0.44314049051773640731 0.51516092956883707821 -0.12166065914558356342 0.34552210269648520047 -1.4895598863683854329;0.095041902797332419861 -0.12550172697416533962 -0.28057017578634046995 -0.13389546059879320339 1.6900404348115785069 -0.39901152437806530759;-1.4665634240628535601 -1.9455985048707984397 0.17740774184459998097 -1.0726644940534155737 -1.2963477167277386304 5.0733792364531211661;0.66338730248377730803 -0.25663123452794223134 0.22323511175097751624 -0.080911621166607633482 0.90333915703050937651 -1.656446785833022517;-2.3875308575554523571 -1.6648405441722753384 0.1279147113944968106 -1.0142432151622950975 -0.74734494132051965121 1.2128968974455294205;-0.12246525236714853047 -0.11569251177881906267 -0.21694716401132388417 -0.71189862745721588233 -2.7574597459218521678 1.5686950295642800857;-0.23474962418228212391 -0.021696883768582357621 0.42590621479249068937 1.925213516108652767 -1.0528562799255141158 0.30509826590308919414;0.063083038365393459546 -0.0089069396969409953685 0.64377479850966368158 0.96126595970911821176 4.0690357317063483578 -0.31966060183362404734;-0.46223901616158963312 1.041047347274599133 -3.0422506552824826187 0.81569886982620010141 -3.4894654906137678729 -0.89224602031943922587;0.37074480113253982516 -2.6744356731691061313 -0.63783435404391941592 0.71455238883498506386 0.21377166269175729152 0.8635081823768785414;-0.059883976748758900588 -0.40758845102043911313 0.76791794529674373138 0.84091934120684752596 0.63622388889837433457 -2.1521129496093895916;-0.89033067337993165857 0.24234880918676837691 1.3503341568783335802 -0.46452573544776692982 2.7439349603788105725 -2.6780870164657630816;-0.23537257430961577853 -0.033876349695212019941 0.99655403635781825145 -2.6958558140861117458 2.1125627014875894005 3.4116879117291407653;0.036402737661234078892 -0.079705880133258230624 -0.82150325938654145119 -1.0734860365560923157 -1.7642639419595895767 2.0679119552899294376;-0.11402871402410083568 0.13591149100221588952 0.18434503080209010961 0.11708167839325410442 0.14555794462325258309 0.65616339653615241012;-0.070964711480878600591 0.028010142730611954642 1.2561038924850147502 -1.4552243519702166008 -18.476925765740954688 7.3408154050556051473];

% Layer 2
b2 = [-1.4846272137367533261;3.5892328191972144857;-0.12883381308661012721;16.44292120068395846;8.1023310390844933693;0.23567708954723370907;28.342774515753632159;1.4155671693193423444;2.7280938490213757497;1.1810128923912837084;1.1501045553066839844;-5.4994787146119676891;1.3492665370175056516;5.3445377898468731104;1.7115830411918959708;-1.265979523718434363;-5.9675168558018070186;4.3921356059409504269;3.1917819055011333873;7.4938350183501043489];
LW2_1 = [1.314743916941041979 0.38554091247921895569 -0.14323733018607190393 -0.36569387110730755541 2.41238288243869059 -2.4644701792615135716 0.5727539022633844068 -3.7231106313240847072 0.91995972451825713812 -1.9580484098853274677 0.35827079727867544579 1.697984198532414668 0.71335621885256961239 -0.22629314270085829075 -1.3540931183647852976 0.20427831081651823264 -1.2066597329453707665 1.0905821903083787738 -3.6580002483614930675 2.3817117693285592139;-0.10757952426206754404 0.1660422315107346336 0.14194447730994508294 -0.26570744634555870656 0.17593374194236002794 0.16511843317532762243 0.056051870104727410937 0.080144591144931512372 0.018643510279062185109 -0.64283717969854548002 -1.0567664476048128641 0.011638283283037431301 -0.062653518473095276553 -0.086630277487831003347 0.33053769114337488988 0.0498242610488829063 0.44336263842294676429 -0.34507158136532783965 0.79541175397638042011 0.57186697549679221808;0.78323643417969157632 -1.0881003203738646157 -0.059617508239259188463 -0.42418415371059975394 -1.5482576997217125658 3.0268032326326053116 0.39057851309432856102 2.2043294704797795625 0.076109604598199404046 -1.0543294410197092681 0.97401623935975634172 -5.2980289652187275706 0.07410843137822736848 0.31644141327196428914 -3.3104076475894133047 1.5408191490400342527 -0.58367722144727185452 -4.7514987629809004588 -1.9997660566892410117 -0.18547321571720642996;-1.5798101370403669907 -3.4402617728296718269 3.6860607622950567475 -5.3052346504932934579 -0.88017156590774892067 -15.469366116076370687 0.86263651360983006899 0.23794702901981643062 1.4236441291534798204 -0.37815304166097130079 10.172296252703764807 10.74196918324427763 8.0626999571865312078 -2.4413490406348126172 -12.699763388444122114 -8.5195580956060439348 -14.199222063645235536 11.248610628912986797 -24.820870882339093555 -3.7399954254159522016;0.50506896212963059423 0.32933005284612315577 -0.074664587431923029603 0.14491554324033287249 0.68124508041742626308 2.936049419331227206 0.061999786455217566206 2.3709837576014130001 0.74061028564186781598 2.3004075970749786606 1.1723016017887231133 -1.2075374151727995553 0.15383383580950923952 -0.95115130795039248834 -1.995238558782487015 -0.37806138887974166662 -0.0060749197352429392532 -1.8609108491782482275 7.7254196426985206614 1.1371789138835799449;-0.3823521510166770776 0.14710047197424888643 0.13123144649860241095 -0.18974146696034846493 0.092042665441667836923 -1.4234553705219648023 0.018563232561293267764 0.11465176356382139611 -0.0018170812474951069954 0.022914607017468998501 -0.53485045030550837364 0.12793683453023738106 -0.044413911437169960317 -0.054604794351306953915 -0.17302945982433245997 -0.094221942164258637442 0.71428744323163950725 0.019831733621337780532 -0.79529702008069380614 -0.7666304703102631013;0.093633654473992841694 3.4261288369086773109 -0.40706116967420080677 0.37569300862145155184 4.3545850492729014647 3.4894001734327009245 -0.42925708711002374951 -7.6308825082040785404 0.51239474509456850448 3.0921318472854437687 0.52814862889195945961 -3.912924985874631556 0.42413815270319904505 -0.52907300890390107639 4.3173340104165403019 0.60751434405739557221 2.2803695396954255514 -5.3841837733065789351 29.613385560858990431 1.0007394409528616031;-0.28373079003016565824 -0.78380819635036524939 0.16420906893431672091 0.16444944337525346789 3.1723763698846703107 2.3518185568589142775 1.5338624807002616546 -4.9391910058442221043 -1.1954726219708062818 0.78750151256759615137 -3.0763689589051166351 0.53593650938651893334 0.10839140625403979878 -1.7576458822948037852 -2.1207429536717450347 -0.73084119031321359561 0.28456602454161095128 -0.38043780761282691705 -5.5050786233726825003 -1.5924685398313775941;-0.23328149461758043715 -0.34181240890687203393 2.0705087002894622117 0.98443154084435857687 -2.0422238496164184518 -3.5035250281471812528 -0.75320370217139831492 2.5959376967619771648 -0.98890894696780329287 1.4934425420963253206 -0.7337203542381990351 1.9496403092576144456 -0.86698749846136813346 0.73341841918462391536 1.4889238007246552442 0.91340715176657361063 1.2060286388295708981 -2.2631976400609712385 4.8817113133239011802 1.7529826539482902259;0.046684966830508295432 -0.079590063242267319876 -1.9611744218356952363 1.8109778727111403018 -4.4451869550536082798 1.2664283281945627468 1.0778761359095501149 3.5426258441139850675 -0.71186333286482073035 3.1052620118881155697 -0.01155829320392332861 -0.89365322304232430728 -0.12301339395800765797 0.77471297867408039917 -1.7218958332549940682 0.66440183720610468399 0.12684676554198473508 -3.0020468684389780556 -1.9790143705109490835 -0.0085472531099046265141;-0.505141931762788432 -0.08147104789709674888 -0.71882037138293453093 -0.84775591373293146535 -0.21568042147735738689 4.198720002701816334 -0.36433118707715694828 -0.42368988001047042102 -0.32237003063493946398 2.7008989192583783279 2.0373301809735604451 1.7070087187087419345 -0.6695012580622766718 -0.21862481630731855908 -2.1359199853697679039 0.51689630927909602232 0.20530169815709828351 0.10557451424735070367 -5.3061185087039728359 0.22598575647225754848;0.59200967441619800624 0.6042194868117700679 -0.24647467937764891421 -0.078123139640798991423 -3.1890199173242921304 -2.0263720096392772163 -1.2105064029230998557 4.6930948754732888517 1.0025015435262296659 -0.35147725037934185455 3.2051519548126132797 0.38801343627403045433 0.062636131745267484283 1.3005723952252929987 4.5282687064290243484 0.50441112564263157481 -0.22589674149427363403 1.573849228582860249 3.9119968804896592118 1.3799619499038651149;-0.2885098812505653143 -0.30965754060216710286 0.38211819207778829899 0.10085994663311562136 -0.15321666675770057031 -0.68693567550306622316 0.28286606169700823088 0.29101558210390543469 0.26153927233462842539 0.28381065825742779474 -0.62137355503731928064 3.4551060652812748408 -0.11631560720673367326 0.079815155225946218809 -2.024333580981295011 -0.34774430639672815335 0.40417536105973522753 2.3111659878820094427 -1.3324637273223927103 -0.60391059904297061944;-0.68886836032767972959 -0.85228025623050140958 1.103043947666104474 -3.575938607256272217 -6.9126834301700723628 -0.088072217486456327862 -1.3284342099291013284 10.917463417885253207 0.62516420635637981462 0.050832041962312370098 -1.3968656958181542294 -4.4932785171742963826 -1.3801346407724612853 -2.3119693992698930884 7.3941937225214457996 -0.49583806708575445255 0.67006114900473745521 -4.6049209170068383656 13.885995040159389902 -1.8225926969485395368;-0.41626342874365507285 0.037122527465667590207 -0.15626070810889675089 -0.71993884412699515885 0.67410611528029384942 -1.0781972993928199411 0.00047078870348327540949 -0.78507363686882258857 -0.086593358718578497535 -1.0848059869766519814 -0.57330192612401154051 -0.15833020352364118866 -0.16354028454772520651 0.30564382420268787488 -0.92679398349516972644 0.44172443385889537115 -1.1189745720717105826 -0.64052496438183659233 -2.0505998933737306089 -0.77923621598783421316;0.33737327577781145438 0.26694599018967790638 -0.38698721242428907319 -0.20817604816761620112 -0.057885059409196679137 0.70499260078803871465 -0.27796585393195877289 -0.088556977042663242994 -0.22740296191774828083 -0.29433026606977502615 0.60323093710149922053 -3.2913126627877726804 0.12024847296673091268 -0.12728623432086891287 1.894356494343541808 0.3480416672205748907 -0.41357607083239189549 -2.1425103992126115493 1.2669140541992010007 0.60236039661421403313;0.72189455947035785499 -0.13906212716110988303 -0.11232346177389063158 0.54290529724327196703 -0.73326627092088625393 -1.7358081413738370991 0.16444380970814390097 1.0235130806032999295 0.18528470772692795521 -1.317736555554587019 -0.97328371687862680339 -1.3644773389271205311 0.46798892531736963063 0.65734787089183543962 3.5107090070590856179 -0.2768498844770757028 1.2719922377284182957 -1.0560131742155658419 2.6000135810415523352 0.20212084946812616804;1.1984602199516558585 -0.24047581691607100574 -1.1115248167933893253 -0.1030701156033655419 -1.6275126912785384814 3.2555822770728730653 0.56259047615615409565 1.5352565217657263563 -0.38790476126270362345 -0.2694206781132455486 1.1399812411543748247 -0.86751572054406356216 -0.1254836070901928502 0.55533955607511842345 1.0747364680560125283 0.34522496782726014297 -1.4705732369295831941 -3.0515392544143526976 2.6208681109664047071 0.86326725964271378011;-0.05834864751408466027 2.0225257071706539413 -0.04731402075743391944 0.1874756661211613995 -1.6222751099322043622 2.4485831339344685809 -0.42669185346138943649 2.6107235388532736842 1.0112412161373565134 -0.043576795413592599171 0.02277941317717557812 -3.220877574049230585 0.24834483502724177817 -0.69187578764668455911 2.6163564279104951638 0.73384626709499478636 2.0539266795020254541 -5.6315108732875458486 5.8094437612730427389 0.46749859023046635986;0.2848616130987237649 -0.01242707559266680549 0.89606254741856083967 0.3298077566345574696 -0.50118621720221345761 4.8327508475873397487 0.39964820423959285245 0.1214286559463556292 -0.13560021493505242107 1.3476645270287335165 1.467529986500386574 -1.6586369644328267103 -0.061930483816533869745 -1.7605916834413137817 -3.142224841798437307 -0.33939957045343943642 -0.42773491567277771574 -0.64476213182034447513 5.075115324633875602 -0.28069612045763647767];

% Layer 3
b3 = 2.4383153502953769376;
LW3_2 = [-0.085579818957453743233 10.691611906078945182 -0.53373644418507215637 0.16578397872129199686 0.22716707283415663032 -1.119361351728970444 0.15395753897718639225 -0.38526930238690976216 0.086606400565191643537 0.083640659690737598209 0.79771263314425888691 -0.39128548932494994217 3.7013128419385190604 -0.072648279320433065132 0.99290843281882845694 3.7977461297224266445 13.941769921631683715 -0.23405900558455913774 -0.19451360818258595087 0.270874822283060257];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.00492610837438424;
y1_step1.xoffset = 2075.6;

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
