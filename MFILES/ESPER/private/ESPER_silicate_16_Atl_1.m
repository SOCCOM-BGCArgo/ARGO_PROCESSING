function [Y,Xf,Af] = ESPER_silicate_16_Atl_1(X,~,~)
%ESPER_SILICATE_16_ATL_1 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:43.
% 
% [Y] = ESPER_silicate_16_Atl_1(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timesteps
%   Each X{1,ts} = 5xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 1xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-0.999999996490807;-0.999999994516886;-40;0;9.1619];
x1_step1.gain = [1.0000000017546;1.00003626070655;0.0153846153846154;0.000309885342423303;0.0656325902039552];
x1_step1.ymin = -1;

% Layer 1
b1 = [-6.3009731036510689961;-3.3201079626619902641;3.0295070164974489124;2.5301059109203705155;0.10498356112610070234;4.7584837422622117131;2.3757385902398717725;-2.8487576460336825335;1.8464284002980071442;-0.5767863115954245945;6.1814327378058058216;1.6287708348443372675;0.9876076375712523614;-1.9714147448452488121;-6.6671676197361513516;1.8683959699743397653;0.61045649972314575127;-0.14151343884167857934;5.0896680918402834592;2.0571301725762722867;-2.1461303393249129279;2.4413545574682244954;-2.6811995292525292456;1.2894059701788092021;3.4925412417908159313;-2.3755828892408836772;5.3322185088513629481;1.0613788664530665606;0.69739025786341890534;4.0152850162020738978;-4.0971610508481983715;-6.1986108382988360788;-3.045196587387182241;-1.2332854233489158879;-4.7342321785100711651;1.8993035371754483087;1.8453618599833472302;3.9159492032769338543;-6.4638309500933122109;-0.14595352332482497415];
IW1_1 = [-0.19118525893791327164 2.8244645905586751589 -0.65059902027281824033 -1.5169839006491570021 4.8633240908566124361;1.4100472952879985478 -3.2507639295191665951 2.3921142378369411041 1.1242363067881291361 -1.8096859282327359963;-1.0212585787368968493 -0.35196138224876960532 -2.054494852300035479 0.22640271412586107647 -0.0058299434519672437699;-1.7610618236353878796 -0.0097508459812351549328 1.1707366770308955672 0.42950464112033365227 -0.94541987099673430173;2.717788624249985574 1.1496484449915613801 -3.2118627899821516891 -0.99864778355228289719 1.3435984877210285227;-1.2066056674656757686 0.72612418265751654811 -1.4573013095853462051 4.8682067819824021981 -3.4785617906859864767;-4.1326945418276004318 -0.10738947641500794228 -4.3078202663478215584 0.93473100265587505397 1.9365228347281717713;1.421528477711474503 2.5962167681794769969 -0.3402575521617951515 3.8759430825749263683 -0.81800065774505958949;-0.14614912602081203685 -0.3449670568917753144 -2.2394304653542649319 -6.5467914230829284961 -12.607420661471808288;0.17268728649336484438 0.098953899854038990247 -3.7728270009637596694 0.44635970440507694024 0.9442297238377277413;-0.30144339275094245156 0.90017770694643473028 5.5605641888661487826 8.4729766604082623616 5.6066146932140057402;-1.257317266364115671 -1.0021485611438782914 1.4943051580574373549 0.9932831091811071822 -0.53686116100730596479;-2.174178110228164762 5.1325058079307330772 3.9200242642090254996 1.4006077751449972357 3.4969781030678497302;1.6234799130545589652 0.091611749621189675352 -0.92783806981025684912 -0.14837909645703756256 0.82720241781595160724;0.55736842587397983362 0.1976311311817872951 2.0065576249941416975 1.3084728924371125114 11.748201631597732231;-1.7441521708854017536 -1.2777730237507245459 -0.55345489842196016195 -0.71471076953889578487 0.42832891096014791277;2.1745155144480778908 -0.19850539319931545523 -1.2607485187390292669 2.5932950997886243805 -1.1924336873217715738;0.18429871947968776147 -0.58837463829632230539 3.7248262802464804366 -0.40904801364437254341 -0.59541342508466421624;-2.181037573736254398 -3.2161753927201965908 -6.4731928653575359078 -0.12594518213834046683 -4.7047870167282992426;-0.59866936967581763085 0.39259985159759136497 -1.4137843985606430852 -0.84882179709865823725 -0.029769648826448780188;0.3913742215690662607 0.4005001943142447729 0.77979858839300675299 -0.45600763197946042604 4.6987920626491135323;0.87095143642977967957 0.62659143922412752215 -0.040059603097216456291 1.3366484310499309185 -1.4318337506145748517;1.1400034847705167174 2.2664517153461285481 -1.4557853062753445261 -2.6587868311101487961 4.2167111049258245714;1.038195678036331282 -1.5413206173333824367 -4.0005675922568304514 0.95418616349440732449 2.2711846019604622349;-0.63822558750869262667 2.2868118540617357048 -2.8384135881567189585 -0.12304945135732159311 -1.1250719041002239518;-1.6831072958852644206 -1.6592130337739838808 3.2003612755845032289 -3.1216051689025223759 3.6944895965779389613;1.4804177734293504365 0.3219011625081021899 0.71085031238841656087 -1.6704413667169459856 -8.7649490741044697728;2.1472699664984249424 6.4112088750804048232 1.1367750196975863197 0.43623862721904671513 -0.82707489271530187924;-1.7357857245745156316 -0.98400714525453691905 -0.96339546667989228723 -1.7163063885874563219 -0.79290846400862913246;-0.86359189613768561067 2.0994180365142840472 1.5190131690783517104 0.89221090069089836749 -1.604292869133454591;-4.9189448428584530149 -0.31726085059666742083 -3.1668319980113892598 -0.13405919724753431743 6.8357142451532695304;-2.2797511360704953987 0.12363939702090270822 0.73968605886278127048 3.3647044732191226935 12.691004607255493397;-1.2367113307739985295 2.7525766436084615485 -0.84407339623027799469 -2.4167293305384669466 3.5738929683977525009;0.074960932540813179914 0.87872633976993652638 0.9353807812638161856 0.37474784210007427987 0.23937305769060687743;0.26392268920763567452 1.3543723462225853993 -2.0958029364433650166 4.4502599369755895253 0.72000862214702188169;-0.36936418634059015043 -0.4269883907575811044 -0.80612898319899306543 0.13774300326931229743 -4.8530323270844739625;-0.23743398871841908671 -0.22986968518912934134 -0.58032316132597960845 0.37765990215145878173 -0.50057304625854059932;-2.4783761689721286459 3.1722230837757385835 1.5561648731628778819 -2.5889260111607859116 2.4685636745452748109;2.4258655052124553642 -3.1077712566233923752 -1.5878297578731952466 2.5918683967902165222 1.3256654380113228608;0.3253151028835437053 -0.18986935190309292554 1.6713977900260270637 -1.5456243404887637372 0.62135599790234308237];

% Layer 2
b2 = -6.7933692738160775804;
LW2_1 = [0.48837649188970055913 0.80885297711350889749 -3.5484322756022770307 5.6731722096217804818 0.24301846340649918554 0.11805674678724110971 -0.16238801231034458161 0.6088933555745751347 -0.42297171738271160946 2.6211425555762031259 -0.15425726807908171634 -1.1059899976357614459 0.20281108582638698867 6.0496923934424904346 -0.42016146232772982883 1.4405130608921694613 -0.19945728516217892112 2.6949318954866683384 -0.12098649364347137392 2.3925470876040990653 -9.2191334468191623586 -1.3652776121361676065 0.21045029503485673805 0.35379358205072580779 -0.44130447948850304307 0.090160581407042819646 0.25298646727274204382 -0.06261392933732569388 0.47120590102183268488 0.86437569146742698756 -0.091384812025680065006 0.1647097734347127973 0.15935297789320093664 2.6677845893586642489 -0.74912790502142934113 -8.7601302259955460272 9.8584192878385028536 -17.60845956271154833 -18.226768814456310253 0.61191206640270334738];

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
