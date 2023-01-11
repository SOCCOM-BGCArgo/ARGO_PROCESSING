function [Y,Xf,Af] = ESPER_talk_14_Other_1(X,~,~)
%ESPER_TALK_14_OTHER_1 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 12-Jan-2021 17:03:50.
% 
% [Y] = ESPER_talk_14_Other_1(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999540494077341;-75.751;0;29.715;-0.9];
x1_step1.gain = [1.00000001035699;1.00022983623563;0.0140762793578401;0.000302617642608564;0.251952632905014;0.0422119037568594];
x1_step1.ymin = -1;

% Layer 1
b1 = [10.151137173231668;-1.7121332337263275;15.973741035913944;-3.8092522317167283;-1.8729538632661067;-1.2730374274697989;-1.0660968230835808;1.0566658380091158;-2.4979677526665425;1.0717390524462005;-1.3905910027242698;0.31797119712476779;-1.6386795837760904;-0.11895037744468905;2.6776381758508898;5.2038345555651642;-1.0559674692771548;1.5572897575208373;-4.66423011523284;6.3167186383843816;-0.41242692580029627;-0.37611143190115803;-0.049181895136873918;0.46242378048116872;-0.31121248356648668;1.4485598394805346;-4.9738459225872784;0.1054010626978234;0.77077268316840475;0.18292328152831139;-0.68999238302440347;0.09820723679528659;-0.025549202301986595;0.41231042706173265;-39.47555330747479;0.6487536972551009;-2.6671014526526324;0.38801165966734763;1.8446520901302705;1.6928625511196631];
IW1_1 = [0.11998442135291178 0.76792924134917506 -3.3164930350016641 -1.9177797697932086 1.9652251418487667 -0.79571532183027183;-2.880966021253772 -0.014471832836813201 1.7244372187676771 4.809957581928729 9.8578575618557185 -2.5661333328366416;0.8403129326786517 -0.75953340308500972 -4.2799344300451194 10.656833959917499 6.4806225283506409 -2.4466566083666654;4.0793134160356015 0.0044387258845075897 0.64890261151090145 2.8902556799612507 5.4708449678612787 0.59981952007604677;-0.68668461221925514 0.69587351279176746 0.75596422727505896 -0.44108430827891226 3.0855066521474592 -0.8159463524563445;0.20967032054701332 0.03348214262437 -2.1703574757967079 -0.044747674143688708 -0.5792055461074771 -0.074232230140633265;-0.0019694433158565709 0.060365983858592237 -1.4293485990691395 -2.1114290948062093 1.5190873913049621 0.22326037775693156;0.0043932507341258665 -0.056486821653023406 1.3634872117636729 2.0166199294535265 -1.4469655134824979 -0.21074595089809636;0.53525426107565899 0.381777241868061 7.4612692276635917 -0.82133195360306832 3.2416960318960615 1.7503731879314826;-0.2020142203290258 -0.12079604455793722 1.4927030893021984 0.0083115315761354402 -0.019742237686476243 -0.82117133604268344;0.22793108879632851 0.15327992731552675 -1.8819473762184191 -0.009566889067545022 0.36170897420017023 0.98078783309385875;-0.23141878130492322 -0.20575103437454931 1.6528442845263198 -0.69217771811810425 -0.017926602295007869 -0.38439044538595052;0.077193827327414713 0.074404131894323156 -1.7533575766543026 -2.2554251135747556 1.306813249193354 -0.22460397660168557;0.48433985214740183 0.52148847317251312 0.64912942205178392 -0.60324566562116433 -1.2598755097812673 -0.21052211754433808;-0.38779021829560223 0.071791567491889072 0.21626533461214162 1.3511765694954507 -4.9803482292039227 -0.58604076324586896;-14.661202816512198 0.56278547383115041 -12.793103450169518 8.2038057072071222 -9.7391469500890278 2.052945242745388;-3.7563300462779896 4.7511975472295127 3.4262331643197483 0.71810496673733093 9.1825768396319614 5.3909684928285477;1.627675667816231 -0.08418802108886056 1.6416085500995961 0.81862977271204473 -0.054808903287005759 0.1142361859830809;-2.5590988068367366 -1.2703844339179355 -5.687166000323451 -2.7401830614613334 4.3746713425806094 9.2368165865789233;-1.7713209367638119 -1.7606619146064022 -5.2215141601136379 3.2874708951279739 1.2958740199193512 2.2005714531770186;0.068603913884052778 0.50791772604153484 -0.073842882868732432 -0.37567425318374104 -1.6422644047641748 -0.86491744870876319;0.24006921969825867 0.22416236457270478 -1.7811977907916245 0.79235921244712282 0.28062741495687904 0.56377059616722114;0.22312990837765312 0.12289276900862653 0.084739848050482444 0.032622991568779139 2.0396441021889089 -0.53664365931906177;-8.7916996778669159 0.99595406465198388 -1.0791187073051181 5.0151868946087541 -15.834306564050873 -2.1568960376791195;0.22524145791812178 0.19401390402184274 -1.579884166079178 0.64103349874919024 -0.094655974455614883 0.20293494065055537;2.0446372422180845 -0.57343797146858955 1.4056186859242434 0.86226889357196734 -0.92670932797980632 -0.50350595943039911;-0.18356130240260371 0.29760926366804991 3.1449817644073534 -0.3113778419976988 -0.62373390421437414 2.8595506668381314;-0.23436033690860233 -0.10946494316408101 -0.094035483893610508 -0.07501754176513252 -2.1568180639795491 0.68782503545536944;1.7773702017427426 -0.61975639275034788 1.2183705154205908 1.0567635439372167 1.0816914663269197 -0.4893522940159522;0.013712276872009148 0.48470774572737235 -0.44804887032741564 0.41126343630623829 -1.8556195968367799 -0.86413657203925842;0.45728322352220341 -1.61352093758741 -0.8634649352693925 2.8138824336511847 4.6904499602656813 -0.44169171042885563;0.56955902895390598 0.059840749720538856 0.34616409592454239 0.87565831024978491 2.8051623054638193 0.33083201717610283;-0.48312448816543718 -0.058429328967875734 -0.24513426671559344 -0.88288746276691532 -2.955038543312043 -0.33539554593714688;-0.049132872277787523 -0.58836108937846943 0.039126986073443719 0.38935749074821391 1.756764049203243 0.9924682227348266;6.3484410333603707 -28.082196755039128 -11.62520335871435 -1.3512245936490921 -3.0826202509205678 -0.5110298690541355;0.48188972093963989 -0.031272248961093201 0.80373217915472384 -2.0429154477655773 -1.2115450107554819 -1.3590771464795381;0.31552425414973856 1.2281783433873825 -1.0643172187849317 1.6899430616094249 4.7404170983045386 -0.050252961910534501;0.46245410080351806 -0.087587671140344894 -1.8468845935587761 -2.0074458294561275 2.1681000857316768 -1.750787770237479;0.68885080455087766 -0.69576788227585684 -0.73188430858036324 0.42966557212208151 -3.0829738856198099 0.7815887337432359;-0.32952847177527306 -0.10607014729709371 2.7895069152462408 0.21927195364823546 0.10385961664348159 -0.56804110908520389];

% Layer 2
b2 = 0.63080315538591458;
LW2_1 = [-0.57496667740157925 -0.014900466446470769 0.057671988700326457 0.04251361996260325 3.9084538922737067 1.2416035124815354 -5.2606194760380696 -6.0700164422050298 -0.035116098366249016 5.1608980844865711 4.5726877592303561 -24.927087739163667 -0.6382880601546308 -0.46637991935960499 0.31440720764314711 -0.010540770928791944 -0.0075830960479506045 0.20564898538318721 -0.010372597479472943 0.049589379311572729 3.9504753837381399 -10.606713801481744 3.5353317406481342 0.012475628612591815 -14.414500780903333 -0.17798503259191867 0.14681504950449809 2.8510767291372114 0.18158174384809933 -0.67644539357939448 -0.037527442689139365 -4.0313642662134512 -4.55571013167964 2.5851940235544006 0.15439444946031441 0.13154776213527533 0.18957955931999942 -0.22218913690650111 3.991116292157487 1.0235989082360148];

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