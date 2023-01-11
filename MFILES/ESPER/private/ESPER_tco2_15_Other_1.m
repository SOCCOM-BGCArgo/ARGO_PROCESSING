function [Y,Xf,Af] = ESPER_tco2_15_Other_1(X,~,~)
%ESPER_TCO2_15_OTHER_1 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 12-Jan-2021 17:03:59.
% 
% [Y] = ESPER_tco2_15_Other_1(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-1;-0.999998793562952;-77.961;0;29.715;-133.803853032615];
x1_step1.gain = [1.00000000274156;1.00000063368067;0.0138606862425759;0.0002963841138115;0.248478071810163;0.00432440760873659];
x1_step1.ymin = -1;

% Layer 1
b1 = [-7.2144612889682174;0.23619579453121281;-1.1739571245525908;-0.25560037502495891;-6.9055020712468815;-0.77429483054410031;-0.314925765567935;-0.76461058392913617;3.2263809373126073;-25.601032566262479;3.5013378433101607;2.2518473245640274;-2.8334859365846792;-1.8521887414406195;-1.6450950518745582;0.99838731257141056;0.066483358297638434;-0.2028309169179825;-1.0680635764546911;1.9345677678040789;-0.3091768607062868;0.85595818310223493;-0.47638387240896057;-0.867519254073298;-0.88366746698823984;1.4684754566574822;-3.5015836489510042;-1.7716368919460175;1.2834366734112159;1.1609679702142197;-1.0964739837975541;0.26510491129387859;4.676878321078779;5.6355455146245097;0.77954592497262221;-2.480340793500845;-2.7746009210166536;1.8706295986620332;4.3868342694523275;0.97765582593771405];
IW1_1 = [-0.25227762974898849 -9.8613414524924732 -2.9168628478223293 -1.7016473672327128 -6.1248857255335984 3.6357140651855402;-1.0082005985757416 -2.3342619595979714 1.9532042134654173 -0.26752356442708053 -1.0159030147298074 2.5149579738961991;0.26807452972749568 0.08972877596919987 1.6707993195419339 -0.90170242993423755 -0.79928975375934042 -0.70127582058950577;0.13044562677966004 -0.21784068842412974 2.7939775855471503 -1.200135878490652 -1.2154022890530001 -0.34221098064162908;-4.6214708275054424 -0.024540605976310033 -6.3820214428703972 -0.87752067510921461 0.95539039458460129 -2.1890763714939645;-3.9267438454749946 3.814714127556972 -4.3929822620329668 -1.4295026862883233 -3.1977232209059596 -1.8031798339783525;0.55418204419856865 -1.1354858851997967 -0.49742532849864518 -1.6831130569010904 7.1553256308326656 -1.3773661659181344;0.12872773206879187 -0.18453318921529668 2.2575981250202579 -0.94784608768989986 -1.3152272896664574 -0.26518271816082056;0.64201650817020628 -0.89803518094019852 -1.1316591850077555 0.12894250169122684 3.346333217883128 -4.1888459201640122;-0.31533803410010103 1.0909081798150277 1.7120880596721046 -22.656134391286276 -1.258058302449772 -0.54791734440608209;0.14076727295786909 0.46771610763927579 -4.0259882279275763 4.6227485223739277 1.1141346296479531 2.2154540284568687;0.33234267181891702 -0.11921899731083589 -1.9365212969041097 2.5652793225288715 -2.6874643727733596 2.5866083544850782;-4.7390487674030091 1.0561522506716747 -5.2422987192160289 -0.10460955374564752 9.018378610490835 0.80896024259731181;0.72235781417322398 0.42886007689673533 -0.16061641237198748 -0.90805490081393558 2.4258184926627124 -0.22670340357218435;0.96904005344970723 0.93614261325811543 1.2761468392850777 -0.99598451322159365 -1.0462839009106317 -0.18772343937994743;3.4223223782852492 1.6426830459532089 4.6583534040922236 -0.34484512055433925 0.90670753597796105 1.4264337809235601;-0.37961531352978178 0.094673086098949449 1.3657702092157293 -2.1586864008382687 2.7399496280028757 -2.0712316065259135;1.6164201253270734 -2.5727234251326982 -2.6911110775376601 0.89522446842972558 -4.5427065172163443 2.9983812176216458;-0.14870859489013183 0.090492608212582676 0.59621179647539047 -1.6061045587952094 4.9091941815536879 -1.0199383348408115;-2.3276789705828382 -2.5869256758091819 9.7899502844566371 -3.7736761719176473 12.765211811708919 -3.6119721381248002;0.051057422960194297 -0.16958917405683624 2.403456807919151 -0.32411780349584485 5.2977728900328351 -2.5687224796546837;-0.52618924546661927 0.28150732805402356 2.8322556668018057 1.5842232996255321 2.2432401871074203 0.86345992278923256;-1.3701539027542302 0.21296235696955829 1.7797915340843011 -0.33433555716808305 3.3463734412625485 -3.0641249943501387;0.47424662981761967 -0.25916842666622625 -2.7580520424692199 -1.5517959761024542 -2.0174304311703901 -0.8978115319472435;-0.58981401707164061 0.52599849503971452 0.28108410940670542 0.12021003425361389 1.3852782909321739 -2.6160804747159045;0.31309928290555078 0.51107781371013483 -3.5704416842729598 0.3772746156461832 1.2980449532172669 1.1323756402785099;-4.3562714722167915 1.5823540245502983 -9.6796218415518673 -0.59237630436409849 -3.6149479433003218 0.8742321781689415;0.041361344567278288 -0.0016310520361834324 -2.6128410499602315 -0.31333862308485178 -1.431810559831429 1.2656706169296457;2.1038988236827558 -0.27020770577511993 -1.8550142386927055 0.32232248160817595 -5.0000982303049319 3.2317113507202651;8.615284205328523 4.8345082308098473 6.8602395949136898 0.7454006251944495 -6.5090148688545089 7.9505812344063669;0.15135248276810265 -0.030289100925745552 -0.6139271934316527 -0.96962021478447291 -0.63183541470152693 0.41030655619004947;-0.98127698332833235 -2.2611270147948637 1.9146808156293671 -0.26965767581636035 -1.1138799646398476 2.4390095497021691;0.018983575175866914 -0.10854988216448413 2.9317389843095492 -1.3400529056809272 -6.8369325566794972 -2.4736694895975582;3.8707482006231912 -3.974162616096212 -4.8177705121116245 1.6725260820864414 -2.9846923594238604 0.47322198049480096;0.64919472169190939 -0.62512242735990653 -0.3653322740343824 -0.23629936696388992 -1.2593433373079348 2.8514828354655086;0.69239002944451078 -0.84488561556763697 5.904493253626387 -3.8677026160456425 7.8503137951092503 -4.8738660150250608;-11.059976434467302 2.7470164035810964 -9.689469882353384 -0.30587324378659231 9.3690583996860592 0.036648601902442055;-0.79308784107102237 -0.46143067999614645 0.10791273473883982 0.89455112927012415 -2.4232412557562037 0.20451823561429103;8.2999020476997263 0.10516787077169887 3.7073200799687238 0.39046043496683341 9.4493442166472263 1.7862727249698296;-0.1389480318689702 0.20285092725766143 -2.5845678267875605 1.1330045715177197 1.5245517037899965 0.26253556461721456];

% Layer 2
b2 = -0.40157832773202234;
LW2_1 = [-0.021417793882459551 -1.4580149235695328 -1.345945530934018 -1.6278375881929081 -0.026400996529933381 -0.023364751467720088 -0.086577671706384918 11.670301740959459 -0.08680477808844414 -0.43189175608502034 0.18162567996722653 -0.21535854115297934 0.043879394762896186 2.916206858193624 0.33416597983798962 0.045880379081958148 -0.14678855281800865 0.034930588926907244 0.41062559126640297 -0.017397982808942062 -0.25773517322397171 -1.6564866691314659 -0.1520490006326134 -1.8699904187903564 -0.92312539717264386 0.28892685924726141 -0.028542662622211876 0.56190333233847145 -0.11669573976709803 -0.017463306664281458 -0.90444444904694365 1.526097235896765 0.23876262819649893 0.030858819048163358 -0.73893495371383067 0.060409666873663954 -0.028978083851680252 2.7673399793711249 -0.029272453988205667 8.6358965065580726];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.00314743155299998;
y1_step1.xoffset = 1774.26765577007;

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