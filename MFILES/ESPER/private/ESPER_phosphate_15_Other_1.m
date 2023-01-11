function [Y,Xf,Af] = ESPER_phosphate_15_Other_1(X,~,~)
%ESPER_PHOSPHATE_15_OTHER_1 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:31.
% 
% [Y] = ESPER_phosphate_15_Other_1(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999998793562952;-77.961;0;29.715;-138.904784718693];
x1_step1.gain = [1.00003298043486;1.00000060321889;0.0138606862425759;0.000266418009857466;0.193068829037552;0.0042772329208737];
x1_step1.ymin = -1;

% Layer 1
b1 = [-18.918104938798869341;1.1197305927249521229;-12.068224800962651599;-3.5675061149753930501;6.4286650528041393571;7.0776012240944794485;-12.524873718209038387;-2.9890712143760609365;-1.3295789067048391363;-11.222388898573408156;1.8440591240326209554;-1.9997019233350432632;0.39670990137867434555;-6.9945063785931260725;-4.5965895076433254118;-4.8112202137537734004;0.84027454492431896682;-0.32766822634908215894;-1.3286038002569211525;1.3230772491818023617;4.2398147356516799888;-2.065444029253248992;0.35928269131506224143;-1.2916484234673457632;-2.9003085442115557058;-5.4660588224383630518;6.3172980180292039876;-3.5505490278818743199;2.6068860648658933599;2.0133627066501764702;0.3246258486677451871;3.2428356307574861539;1.4931277992562244084;-6.3032427696675972228;-3.8438110060771197851;1.3538906649140960425;-1.4263113814291921155;2.4631151806508899682;3.4738022049550290049;-0.92293655770445803821];
IW1_1 = [1.1311464780701949717 1.3885377096891029503 3.4323669251429924643 -14.069980641819892497 -3.8504851871083567438 -1.1305213510762743834;-0.034355815008317795722 0.024181356197255341606 -0.60730513240639061312 -0.14583890649529279715 -1.4736516300038404736 -0.9557667841908622508;4.6139254535958000858 0.309004629979065637 6.9739581841419102659 -4.6074618971160195002 -10.630628444770726659 2.5015647302866965696;2.942946943057912712 -2.746639218996255849 1.4167539757678082601 -1.1798525759826097303 -3.8146106823833920707 -0.7769608746266360999;0.27382012561875135326 -0.016633805691684086514 -1.7171102806894102955 6.0523261207039995568 3.0070714242684886131 1.242834262141316648;-0.17782463414943514013 -0.13873573074881412337 3.3179929189866346206 4.3908925174653408163 5.8445492270161487625 -0.24324431920024455756;4.5665784067797892121 -1.6963514501026932813 10.490546589884377937 -6.1599819476831765286 -6.7873798112787371295 -1.3406635871814169825;-0.15240702857310811646 0.68557840950074588893 -2.6673103705184111334 -0.54779272448900184767 2.2818177991106414204 0.056090963483394651323;1.3591342882166108641 0.26856236696944152387 -1.0636988588495170927 -0.78186366410994390108 0.48624780369083953069 -1.0479884084798700705;4.7611412150194993842 -1.7818563632878554603 10.59218285636138468 -4.4343465714881258322 -8.1739470589349689789 -1.0446094330963555663;-0.71215007534169705306 0.10694136709837541444 0.14463460251107831955 0.0265330412976550363 -2.3434665931658837934 -0.34709068600982040387;-0.071254936933229540141 0.047880931695692276195 0.29867339600633824315 -0.35557269967614324457 -1.3258596444450596596 -0.64406393462338229483;-0.39160405082678101385 1.0953556574624865316 9.840974431794059285 2.84537869228119078 1.5187031138136743991 11.85212708189302333;-0.2996792255097983082 0.64454639748386544085 0.56918692399659920333 -3.8061277207279178292 -8.4380001332360947686 1.3309216290920025116;3.585030700144399507 1.0436423888308294927 -1.3733823845700443567 -1.3415498253171560528 3.4107585749864917091 -1.1196661867763342268;0.26275033530498553835 0.034777556212565764737 -6.462171536878679845 -3.8472450673909559526 4.7959863323503739707 -0.020948811323992858058;-0.015509383842362054959 -0.12605471274522855873 1.3680145882783190103 1.5421722371942114815 1.189911902715179659 1.109014820492910891;1.6623916173946444719 1.4641701041800705418 1.3035366677999751239 1.1368332567630032859 1.0896772846756512809 4.5729031045024237656;-0.19351974980362451895 0.50228273328530959052 0.26072432570459641576 0.62659688127979851213 0.12827321360837240505 1.670006151305174269;-0.11757760120418504768 0.66252880759406418409 1.2593164732112804849 0.87109945448937353696 -0.51493701663218383047 1.1282805958774779764;0.1772918878260077713 -0.032739732729531692557 0.090406762518899896897 3.3588853062744687961 4.3178951637486813198 -0.43755354567540183375;-1.970485394777631738 -2.1619811867916935988 2.298059826707025799 -0.047018599837527284835 -3.0535036418887835374 1.0934615257489548323;-0.26636308425093596641 0.24288818829429892343 -0.18161878025204208909 1.1359870681192281161 -5.9008576346440753113 0.60721418587664144351;-2.172307816774026179 -0.27464169621857786874 -0.69164239543656247378 0.33144410488005110471 -4.6605675743125294375 0.34290561077649778232;0.26962842253065133491 -0.51962843239652167782 -3.5206787912176644717 -0.47587111906780626969 -10.38962486518552808 2.3072696609817828772;-0.11445601742463135109 -0.37565275713576251615 3.5901245553586216275 -6.7247231607973141365 0.57765573182684082365 -1.6857708597851832621;6.0491144746434279966 3.4083260502082701748 -0.58872581757095499011 5.9748839749938094101 0.52904157809064145113 -4.0849190447877266053;-1.3629035566393159495 0.13499565769798660408 1.6558834588662718623 -2.7407790326531662473 8.8796773046824561959 -3.0738076473531559252;-0.29547485669000284059 -0.039007048227920601557 5.0211063071456045748 2.9214140934457355492 -0.17467112363417061638 -0.76790608336096044528;1.9198689201568215346 2.1312912848191190029 -2.2845828682952844524 0.076137775015820557956 3.0763212280192604453 -1.0333354652848154576;0.17561075967896086492 -0.43540894661119644571 1.2376934268856312205 0.68463853592856660502 0.033919525968840937014 1.1066265726253541324;0.25249378056277504978 0.053940432722437266222 -0.93884946330851148932 2.570248799979220955 1.7324607168078249853 0.041681477240194185485;0.11672294574168982917 -0.56972979152553293591 -0.24219281327254255087 -0.68123609218884240146 -0.01144942163790056526 -2.0998259001309911298;-6.0239138820219064741 -3.3906610315378129528 0.56692843338765563921 -5.9553036179428024965 -0.34734795828907039095 4.1085898573577059878;-0.17425578283168924321 -0.072712063653335948232 1.4317119764995458464 -2.7525915670538876334 -2.0521280720532715236 0.52262104378119356163;2.1766389927892544343 0.24646561782866427404 0.87240980437974191553 -0.31973221693528725007 4.5224593765487162145 -0.43629735578983286359;-2.2150039245176693647 -0.23218363985703233121 -1.0625834470396027243 0.31307452373801608614 -4.4798746585659667829 0.53524206535872620449;0.067866390691584441153 -0.2099788391609809235 -0.74218564558379829244 -0.013719500880748285859 1.6529696291133835206 1.0566958249790494584;0.22134817118312680173 0.059313279362155017782 -1.1690195143242951836 2.5993071629075554796 1.7309736534129418484 -0.26110060655978700739;-0.30997913111601416158 0.47175627417476623249 -0.17625608969526090908 0.27305032556459507553 0.23331741808898839863 0.60670581717658633725];

% Layer 2
b2 = -1.3527667380980954448;
LW2_1 = [-1.2800396982973556437 -2.417077038303697023 0.35000680811722190144 -0.093823800018412339963 0.45132879993707702981 -0.28334494355359102213 1.8200737877310295687 -0.70977367144500780327 0.44998056652140372424 -1.7345678376466946435 1.8217674966994603647 -10.153080963988470486 0.033366057812018422302 0.53731778733829693717 0.10854419182926085263 -0.13317656804834557271 -0.76039140725886111838 0.054437290285290526604 -6.9608189436584311238 -0.59509102257818236037 1.02919806992529006 -1.7916857081632131532 0.41502573233541228959 -5.0145236527215368838 0.18065180377571254944 -0.23370417618290045825 -1.2535189617208861712 -0.058241442592858301552 0.14899414109720726662 -1.8339927411406045099 1.0689219572784349754 -9.5565342980445109333 -3.7363509262863723315 -1.2549065881777814457 6.7420895141595913458 -9.5120413808242059162 -4.4728096714895455577 -8.383059287394926784 14.644363177961652767 3.6954514270325771186];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.568038194888224;
y1_step1.xoffset = -0.04979;

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
