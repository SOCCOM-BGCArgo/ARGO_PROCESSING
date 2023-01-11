function [Y,Xf,Af] = ESPER_phosphate_10_Atl_1(X,~,~)
%ESPER_PHOSPHATE_10_ATL_1 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:29.
% 
% [Y] = ESPER_phosphate_10_Atl_1(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999994516886;-40;0;9.1619;-0.28;-0.52];
x1_step1.gain = [1.0000000017546;1.00003626070655;0.0153846153846154;0.000309885342423303;0.0656325902039552;0.0470499670650231;0.0145232735458572];
x1_step1.ymin = -1;

% Layer 1
b1 = [2.5020298954439863515;2.0050164135188706105;1.3118089280121809637;-1.7655893108716913531;0.62357172203422728263;1.4851417189350666348;-5.0088966984874252475;-1.4627799273445167838;0.74650383260748542913;-2.3697820772167834669;0.18953212115101744639;-1.7905020721276647677;-0.67286988660976743137;1.9880651816560257483;1.6160395813960064526;-1.8206394246821804295;1.8023350298823472837;-2.0416299571878684382;-1.7815549373562646895;0.021299494905689358992;-0.34349028289597804964;-1.3897135881575353178;0.29804658270786943408;-1.6688680353415983326;-1.3645807691881686186;3.3610075366604963421;0.28041857183662394215;-0.92667062866052929415;-0.92510281051970189026;-0.92448735139550175521;-2.204344158845670254;1.7759419727241165798;-0.82798863456758531321;-0.10204129479595214214;-1.6406844801970490444;2.1634579291459412609;-0.4583013662435320712;0.060037673982119829641;-0.96562952604226903741;-0.088004415487801052054];
IW1_1 = [-0.79494771602765468277 0.45765242441841902377 0.8444191554081652562 1.6422242936152684489 0.98177415688829050566 -1.1350738570475660261 1.0939309750229659901;-0.04477732988489440541 0.11757560019897347892 -1.4677994416121964427 0.12929764347639038924 -1.2607423636285790547 1.6604606959534740085 -1.5541683599300522012;1.1889370139015671413 -0.13832680524916871057 0.82509602003606141718 -1.2516123857744929637 0.65408477340380022369 -0.9863734209692981425 -0.42269900484089845127;0.2036690764836553702 -1.2044208166235588209 -0.33507331733974454746 -0.76814438466275702933 0.84717034855701645046 -0.51826648312618384118 0.4194415145566408043;-0.0099312656399388118955 0.085606330185659079701 0.66516943855831323074 0.29595569972594415464 -0.30775341010543710851 1.3805190066163202101 -0.70358040444211289444;-0.26213081264085541511 -0.13129814784288099117 0.55306358771118879769 0.84426050689022791307 1.3572888543071250922 -0.2086618114626607734 2.2545655734140490623;-0.54376882049632291327 0.13837609684168164503 -0.24082413744477718098 -4.1107965218582869227 0.29491614302853363938 0.064038472174439101114 1.7214546982414538068;0.040509341396264067492 -0.2579538668944560964 -0.55721843768340839631 -0.55712379890098329138 0.99897689494263897103 -1.0117326736105651808 0.57474531298659803369;-0.13356924131664729694 2.2069619622503582335 -1.5316349267077662422 0.17861322020608666472 -1.116674638069128811 0.27315066429095735678 0.0966062888690913113;0.16646966333901283552 -0.31422679126451419362 0.48885589853854505282 -0.59508346008549672934 2.2900271930797080344 1.6186042511495939422 1.0547213818075016345;-0.75028290780918971414 0.35057576345203611545 0.31291473975009387321 1.0639128629009506621 0.44142863407835786216 0.41615002635180409207 -0.94587455713096890708;-0.88212401838786835828 0.27988466680055013036 -0.40447744384585915256 1.4076006377701248962 2.879709571474712515 2.6204529551441568813 -0.89285376476452305017;-0.085886851205040165569 -0.061045228340591423111 1.6435493170825923315 0.42347494224417753239 -1.6971686553768352557 -1.1473048428223178163 -0.48386909134247907405;0.11249458409855765673 -0.69635229106272344968 1.0668524368297145344 0.59182195700011108741 -1.4178127936798474718 1.2550560294339820899 -0.12545167656269459733;-0.37385264166089077209 1.0957911607627017769 -1.0019925890939720414 -0.10046840923546480895 -1.742587970035303746 -1.6648081350317489591 -0.74992157268950720983;0.090350449687013625888 -0.11067159553256490645 -0.5175075892060442273 0.68456654644075676241 3.110251042988032566 -1.4333155255727421107 -1.3073590245705775104;-0.39052584459844535969 0.010026633981645060451 0.75739891464531661214 0.66965155012233124054 0.97410711566699448305 -0.35096910225372518477 1.8771300565065176524;1.9192843722785248861 -0.2984719981703480296 3.9711039978767890268 -0.20284634139883270243 2.5561826915735741395 -0.26382880670825348268 -0.45750143588744890399;-1.2482520581590883513 -0.9064141274266160897 -1.7344120176328474692 0.46217810371023726468 -0.15948345329544108551 1.419750247887914929 -2.9607539886046825295;-1.1694069930285804215 -0.45148656485343080513 -0.34849008609213483201 -0.4925486058680084156 -2.5475523318558570551 -1.7808185633082704857 -0.94167581068972083802;-1.5907993664090407382 -0.20302391463422389539 0.51359944467829399262 -0.83078983160194197488 0.94021946741307382833 2.0263650322701449724 0.12866893558063816072;0.40058879289463644735 0.30553247589165638542 -0.93157753419185140853 -0.086438536666637574846 3.2985371280012549811 -1.8044771217475719105 0.30232909966510829181;-0.27604050444084488047 -0.3211288155964774349 0.94545023938939110053 0.26593275627231205416 -1.4863475059437025738 1.8751055936676193081 -0.077588235062105978534;-0.025182472866664584732 -0.0028355971212688600019 1.8912713661026108181 0.53261942268136597267 0.15960664422033057908 -1.377797627665939828 0.25490224807427935838;1.8827466864805737501 -0.27349461380559758616 4.1112103244703401828 -0.144475011340375048 1.6380393349861988028 -0.2550511385685343968 -0.60442030808157842792;-0.14870519400301149271 -0.23342131238666713466 1.1241244463634376682 0.86966691782912075581 -2.8198679667397681392 0.84721904103575207401 -0.3392200621967367824;0.22095820915353323044 -0.032414868305112899682 -0.48011898197616226236 1.09218211563756884 1.5354722631905248775 -0.52296976413609996825 -1.3124825236717190258;-0.65989217547718204138 -0.54232400627155330497 -0.81254761538423814216 -2.6758723789322798048 -4.7220456658991771803 -2.779319697117414556 2.1646199271558104194;-0.42424516930215089294 0.33184895892568189169 -0.43944189011052525728 -0.014259186590526953808 0.76694385306785750611 0.0088422589425698368271 -1.22632295043665418;0.86307773405861532012 0.54935289908793683011 0.73596971187629045552 0.046907323050993994518 4.1311566625439493095 1.8757516561824609802 -0.80053861909899048932;3.0473414326963346177 0.8871629770308272489 -1.3366805984590519607 -0.0036816928501217523716 0.98530802014143037404 0.41775276149559648076 0.41675362849062153048;1.1012686426891888924 -0.42249157530138498062 -1.3035070473048160355 0.017479239902562056924 -2.0121613926286538998 -0.039768584501817902388 1.2740649714774017198;-0.19462159219953406586 -0.3027813652793893695 3.9830861326528768096 0.66458609365884180331 -3.0788749451933332679 0.027564026122583122624 0.017216551583911138268;0.15290840469001623636 0.32811280101844214618 -0.16314318924075571493 1.5970221891282150306 1.8124184646808059185 0.68352797857798353043 -1.9261465496019936072;-0.77579886730973568465 0.44666413301280866177 3.378085943345908948 0.096640490058389968975 0.8313108206836602454 1.0732696268008115048 -0.78273745353211365838;0.68942344697900415085 0.51641207585835846583 -0.69075500063164918529 2.9131026637272556812 8.9122784038326123834 -1.9804528212552199307 4.6487935312647366004;-0.86805550979220846841 -0.17821843988262026093 -0.92853781810965541066 -0.89529720947508750761 -1.5018612935651931561 0.93090124153118369854 0.20899736564080856338;0.33165691771385424902 0.24689617384282824819 -1.5134400796520774968 -0.12066960749596479519 3.2195879670219658131 0.12035327030708649343 1.2820826176335229363;0.029836950322045485001 -0.26356890615876599204 1.0676042088209876102 -0.48710574413317347942 -0.095444809603513397489 -1.7409960991010624554 1.7888401949308030403;-1.0215661809052465347 -0.03812583132163292754 -0.68678515670538575755 -0.64895664798841734555 -1.4890421344731523678 1.2452447316091059815 -0.02388397805538244531];

% Layer 2
b2 = -0.021887051126089046571;
LW2_1 = [-0.33119897925624319779 -2.1952332726788523765 0.31626078627865811255 0.72631830569983302226 -2.7298201505562875724 -0.80093345088115663888 2.7214018960044907836 -5.0074580722143524625 0.14813058231567410861 -1.3820696387894053636 -0.5631105048485955944 -0.14099502368089325555 2.1339086539472051918 1.592288507131362163 -0.41644690517023569276 -0.889521225084517142 1.6042592549202123742 0.81452892936804843771 -0.05484552940840577151 0.19887137603808852093 0.13031441285344200209 1.5814490780410870308 1.9160498524580897861 -2.0634870184897908807 -0.76455567149692893469 -1.9924247847838287839 1.8502015057258294206 -0.81605615116241336793 1.1468438739084840794 -1.1603935922782016021 0.1087953284023027295 0.34532480308391777513 0.30402464439308563016 -1.2188274175362194818 0.20215464627012796894 0.24622281745014099053 -1.4250509197474627854 1.4077948401968121139 -1.8550749772901466983 1.525800511498928369];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.604229607250755;
y1_step1.xoffset = -0.03;

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
