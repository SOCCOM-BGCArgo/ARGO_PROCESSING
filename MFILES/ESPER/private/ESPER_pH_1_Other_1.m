function [Y,Xf,Af] = ESPER_pH_1_Other_1(X,~,~)
%ESPER_PH_1_OTHER_1 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:20.
% 
% [Y] = ESPER_pH_1_Other_1(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timesteps
%   Each X{1,ts} = 9xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 1xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-0.999999984769129;-0.999540494077341;-75.751;0;29.715;-2.05132667284751;-0.9;-133.803853032615;-0.1386];
x1_step1.gain = [1.0000770825255;1.00022983623563;0.0140762793578401;0.000302617642608564;0.251952632905014;0.0588888685539285;0.0422119037568594;0.00432440760873659;0.00826143816767911];
x1_step1.ymin = -1;

% Layer 1
b1 = [3.9952876352811999006;2.3204218251670609696;-3.543898526990486797;1.2986451817678708043;-0.87782007825017338565;-0.73785999158602522119;-1.7248884132730688545;1.928405672786378755;0.28177862306842094142;-1.7177477540860817129;1.3414926452894999631;2.189401779275183646;0.87356106385604626041;-1.4661488271055105326;-0.52162520501995568445;3.8313207494670780129;-0.80683289607649866237;0.44014085058248464399;1.2843162680460986635;0.60918103570963444415;-0.10118293669760955455;-0.13068139146774815362;0.63958400181110175442;1.6764983717297039512;0.8538421456040369284;-0.22220422861993327501;0.25904060982437016047;-1.2522694390262794073;1.4857887176343194113;-1.4586748442597741171;-2.7124224965417700339;-1.7481161042843942344;3.5941931122789054776;1.4202832723201461729;-2.3962437578091115853;3.0245684921700415337;-0.48158975166934064882;-3.130894468580280865;5.2595690192324289214;-2.0941780698607339595];
IW1_1 = [-3.5123241506489444319 -1.1474393494444308228 -0.15969704026511546968 -0.1791284295386659231 -0.74351617030522576535 -1.4137872557842876553 -1.0634981005849428737 0.98912939909468244881 -0.82865544429014104111;-0.43915116557479083603 -0.44885567428293204895 0.90721440383470375313 0.74912108133399457621 -0.93349876431451161984 0.63503672408088618351 -0.34466306980096117485 0.15405385255594714056 1.2484661128879659397;0.21056516005145456405 -1.8601754695924839211 2.5177212870072773754 0.24197255555316768416 2.4790452385924055534 -1.9586839203102224261 -1.4600740873227928684 -5.6056194488200130621 3.5635782235374433569;-0.40487063697860686817 1.3515579737587193332 1.4828837889750721146 -0.26425652570041496814 -0.016196020256848128044 -0.6602768754377552618 0.13870890463218596422 -0.13631098542273853425 -0.48514927630623377874;1.082596054169487898 -0.42419199034486004818 0.31449852539281419883 -0.037592584885232725733 -0.39827761699163960474 -0.57966865660371091895 -0.5413784772746313978 0.54072106409293008156 0.019805059012426223103;0.88030543359750601518 -0.37089536200882639028 0.51698389032239022711 0.15484595430067504829 -0.68345505672060602453 -0.55809206934645205767 -0.48981980911322009486 1.1795952625570880201 -0.25734622965739617495;0.28859263232199211524 -1.0839950190896185944 -1.3814163348437000778 0.34568806293296305832 0.53299928151390396636 -0.27147775302511639728 -0.69884302333972458943 0.26372725204308611335 0.51228326413033464792;1.2601610604288027151 0.29075417088800870058 -0.44552080968603441358 -0.33224811284749256046 -0.64166922979956098061 0.77262140930899392277 0.51702383474203228442 -0.58934468212047919078 -0.24537360110498335208;1.6175343760670184423 -1.9761051742642226081 -0.59214416242047740457 -0.47740419070214168418 -1.3541202155144773656 -1.145685990632513418 -3.320296952130682655 0.16531554640836162995 0.55076702002471700226;1.0206303089069461709 0.72064474347771467766 -1.6119772187267320529 0.94742287318942419727 3.1865085082397475169 -0.17481667097474581918 -0.22620819639538483625 -0.057973477609663259802 0.19797495046031957089;0.47857147221344420451 -0.55107297409635747076 -0.51702986489086211197 -0.32608790760998673797 0.499381352551323765 0.83021868601977777402 0.93698037358414099085 0.41200449688940588011 1.1820802157182408454;-0.68303056457990041217 -0.079725040146717801592 0.49975003341337059792 0.3346647955731503199 -3.7112722449715112027 1.3786786268551425305 -1.0720957419062284721 -0.060738626710088949867 0.26023947061910063239;0.33893168233964882097 0.46331100721257584274 0.1560286439328843755 -0.42189650715819415616 0.78329661420557072926 1.6445803773277578586 0.19413515755651028494 -0.69109491804177669572 0.42346875292925106882;0.97062583349418618539 -0.82492737103543345434 0.95639770187804418189 -1.1270973947976064 1.9326232584809799331 -0.58484905903524508908 -3.1059766872042442465 2.4571733078347786083 1.0551887801182082693;0.24159637483917970058 0.13590045965635724756 -0.97258592742158922473 0.57987936155460406695 1.2948262485288206669 0.12296779294766128909 -1.919260115335248873 2.8127241609673787615 0.11680544986391265538;0.76769271042432440932 -0.67736901472005484059 -0.77188476152255236507 -0.61481258278191652789 2.2557119049839284841 0.76961399529480278847 -1.2256788487109959718 0.8269987559399643029 3.0375340078372503605;0.39908541670191599593 0.10106874496695526788 -0.38564846784046774575 -0.41298442568743437553 -0.52697838651182171699 -1.2862787938028921175 -0.45084706953682085473 -0.065033419513865153849 -0.82051230158271348802;-0.53890341742769265476 -0.26574395390012856 0.55713311436883516947 0.58334330807848433498 0.31758201034606703406 0.37190757121362466586 -0.049035287821716906642 0.2558848620124454798 1.2973081758025801768;0.48892976964309448062 -0.52417747257764180091 -0.58502792709419382255 -0.33721584899047468697 0.38725444429792549306 0.59497843925940807619 0.84577396984641417532 0.52883710758117863371 1.1689700659923651926;-0.018698119894352990128 -1.0231769867853008904 -1.5031698234233172862 0.34340340471977670589 -0.77171082212713837389 1.6001315557769313358 0.78369826535743170037 0.5330571338039820084 -0.31966382590888803295;-0.31453028604207483587 0.45559438887463382661 0.16730119814979296788 0.022135227598591152715 -1.1585660352358566616 -1.3886845298755370326 -0.8577503594414398469 1.0990282491265443987 -0.18902088853934501111;-0.47087804642219455165 0.4296053411430081348 0.065643712485747976459 0.05325533202208331679 -1.087196135942963382 -1.1808209827471274345 -0.49946707788306587972 0.91999952467651857635 -0.40991536750279866652;0.24773883301317745076 0.62713221635852678393 0.21271461977341771243 -0.47138520002501282935 0.46448621302466069904 1.2334314130885168925 0.21710836997426716244 -0.47260535967718736439 0.26140571748012852993;1.5314960679913722785 0.65939345455065112578 -0.26408284415712129922 -0.43617183323101232917 0.094058407763255122425 0.62265989764573403953 0.61361928453116276749 -0.39012830611487930854 -0.52087936375132659617;-0.29649903332625399122 0.12299467314736688139 0.28371658982897335344 0.070370625649467208329 0.99801653209876117767 -0.99402337715957878039 -0.32572836432402957341 1.066759376784204516 -0.097614623328865265939;-0.62243670831524200437 -0.71530290479860192399 -1.5611890053358381536 0.44500589533461343272 -1.1340776783324473254 1.5190329746521160637 1.058354327473337797 0.37299582476254805607 -1.0966469888893339846;-0.047528386018969873139 -0.12899590831218343823 -0.27195981426278342807 0.27881148614679368558 -0.2282513473900513179 0.45841400249604108419 1.0503117969174935453 -1.4175186096388250334 -0.32613932362766395734;-2.4026467302139304927 1.521952071084492486 -4.1928380593580936164 0.31448692162397551275 0.077416961597755615254 -0.28600202977660293957 -0.91908092993351064859 2.6800037989790004289 -0.30233305027582529911;1.6353567074368680068 -0.7786796079588382824 -0.17825060661017610997 -0.56390247785718849993 1.0788040070697220152 -0.21756888408424546189 0.12059772666229834348 -1.0365384523508296422 -0.011886244131905466387;-1.5340085812957702238 -0.81705241356614821679 0.3196668083761444179 0.0018964340035272520685 -1.3286765609232487062 -1.0073682423918226547 -0.53714446948443239638 -0.52291454261454106689 1.3669424980625934918;-0.79120441640168270858 -0.48351013897351141635 1.1058177120712591179 0.61673323057042772621 -1.1753939034345668446 3.4482693568343738555 2.2544556176860290897 -3.4910382098526033801 -0.60193122109064667313;-1.8009232150006231343 1.0610110624005375524 0.40164611746599465647 0.52017329910499743306 -0.59100183649195536795 0.020114437070399587837 -0.54398917783026201622 1.1087051048333445014 0.41569385537870912062;0.42185229975092503052 0.97936245780487696244 0.257241146260716258 0.078276235341549496582 -0.81095965081380028128 -1.5336513299927994058 -0.33539522477960470592 1.3393697415856815613 1.311749751909275119;0.32475023246071127447 -0.058167747359514548888 -0.50537000905605145018 0.35575428132596997877 -0.36642683410775567232 0.39597529934408032259 0.61167542632186133744 -2.1171777500790831539 0.20942614943727380639;-2.1271461439555610085 1.105822901679641479 0.66432454380336269928 0.34865193252233633858 -0.64129781166205668264 -0.40624185570513582366 -0.71933702584578640415 1.186147382465909228 0.10091623870280468089;1.1785554180537998192 0.10137372010282064139 -0.59082343092139022467 -0.069462640665406194529 -0.33210164497723199117 2.0353971941130555301 0.5263594689120066894 -0.25880818683922229839 -0.95987200545814355124;-0.51422397385608265008 -0.1966300174122563027 0.22792751806199371134 -0.54317858382309269771 0.82021920215846577662 2.2511032656531786422 1.2437062905854792128 -0.43053426737462768825 0.32741269660378452544;0.33791225293142929109 0.32434838955159017715 -0.73885120262989223416 -0.69180848515993209613 2.6461595232949362 -1.362826289889683018 0.51913558803066395253 0.32367689700425239518 -1.3714209867804769694;0.32115503364559355726 -2.0096228955295756258 -0.048909622041770281498 1.6645453269891239412 -1.2106439846158809637 1.1711236686470587465 0.18530005487802653219 1.7723156920761582445 0.72530741867603609752;-2.0226635635215162168 0.91891871251993484382 0.4111878470752547976 0.33305727286476805027 -1.0214898087570625229 -0.11685858562900998392 -0.29226599718076362588 1.0342065099704884634 -0.089491423668364600275];

% Layer 2
b2 = 0.77123646679933977044;
LW2_1 = [0.36819643098752169452 -0.56577599365698527656 0.20153122597286701745 0.90231312910644800862 0.93328027112689782019 -0.84022629531479442466 0.89844468958081691312 -1.3747762234032008255 0.033565411627515612247 0.081006815287856137608 -1.2656993963869205277 -0.24950636280716761983 0.72366094697119409673 0.072697960021876020797 -0.21485240600465541072 -0.29567650083789592808 -1.0349808163653433368 -0.72851571161937722643 1.2725583626867675857 -0.23808806391778677969 -2.0851901396705683212 2.2139116217440597545 -1.026708821064529209 1.0376490933963751218 -0.78276207603110936617 0.20889786785448671846 -1.0013479489022254487 -0.037655641514711964135 2.1964104693572301663 0.40791931735231740097 -0.22896178063597197316 1.4340645160167249106 -0.71551763545882141049 0.43317818843735605983 -1.5795407267098369175 0.75589978141669478262 -0.47683235414567648469 -0.53182062962403320316 -0.21095347496955743627 2.1285788860850183291];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 1.90239167587468;
y1_step1.xoffset = 7.36300175998845;

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
