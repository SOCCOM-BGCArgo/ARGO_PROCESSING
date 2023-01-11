function [Y,Xf,Af] = ESPER_pH_7_Atl_1(X,~,~)
%ESPER_PH_7_ATL_1 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:22.
% 
% [Y] = ESPER_pH_7_Atl_1(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999025224415;-0.999999994516886;-40;0;11.29;-1.99149636107278;-101.945317545701];
x1_step1.gain = [1.00000048738803;1.0003169936794;0.0153846153846154;0.000310077519379845;0.0705602618858232;0.0604107502911417;0.0054210015363163];
x1_step1.ymin = -1;

% Layer 1
b1 = [2.7583780204935992586;-2.7670968767848078684;1.6500425320173461863;-1.6736829786813771648;-2.4981172932347548432;-0.24438863075245118028;2.9494022954374377044;1.8941003945629291128;2.0314683268758892964;-1.5310032991719924045;-1.5311558990013729442;-1.8688245833758447834;2.5101683104132495039;0.43468270644767426081;1.5254390236586943264;-0.96705286029291415151;0.2763784417858101139;1.0839164805782364187;0.66789897726009872514;0.80912419546464997566;0.15453506151340945518;-0.71064610664556648878;-1.022845336013378903;0.60386780646812709961;-2.9237590787203640019;0.65662822085095473579;-3.0459391760703429242;1.4227690923909854792;-2.2009787250004362313;-0.6739262097933564899;0.99881297878094910114;0.74596497793048111014;-1.0541179039181787225;-1.8293086518242951399;0.025248310823091566629;-0.69872119059329118596;2.8716971893408502758;2.8651210466558918455;-4.3944237335023492719;-2.2834330702127560997];
IW1_1 = [-1.8040328274557526633 1.6903484837204680868 -2.7733595160054202999 -0.042318992299166101334 4.1511922570955306355 0.68196585751174609857 0.40122517658328560319;-0.25234187680837039114 -0.10243622531565685363 0.94242216515027377177 0.23206194056975074091 0.31076107140044706423 -1.5230516505451130449 1.8781584963722506298;-1.2835475028165941325 -0.74702051694972471108 0.50148970800599712572 -2.2764176075862594395 -0.37486494720044782891 1.3935864496238370425 0.020208747889479038162;0.40221274862582356846 -1.6481356358520060645 0.33060785610548121261 -2.1023666434750469811 -4.4920981798288686093 -1.2059864219289417875 -0.62650727215111179369;1.3737635769595915036 0.15155198145624834694 -0.53627203859894667648 -0.50946424962156378058 -0.27456563759610980169 -0.13900428355625360499 -0.40382643620264363493;0.9392415137513889789 -0.76003196418112672905 -0.11025213090452495868 0.58231992194899906412 -1.2643878444018961105 0.1300250453530614192 1.085576515608677628;-0.78118036054478001695 0.31358857008072615891 0.28021933079361321939 0.10783945407843487096 -2.4092874990208610342 0.047532958475454266389 0.16020139462257748653;-2.4828112481898583219 0.50610816716503959078 -1.5150163269000713751 0.0198880603148052712 3.2603069293833986109 -0.26620502087936448454 -1.3672768596811817332;-1.7163606898288050751 1.49997177402357873 -2.1322372636190776163 -0.38987479875326563894 4.234413431812536821 0.65562854250737911954 0.3417020489546591655;1.3672931595735346999 0.67302263932780781452 0.51894180671798562088 1.3926926140832094347 -1.4571015015676755944 -0.89627038263186153966 -0.54321911504047737385;-0.71317611280776749272 0.38415012912633050846 -1.9841430519622684781 -1.2365473005004907758 1.2081911610324225226 -0.24400083025795582525 0.90877794370078424091;0.16126876675715684506 -0.046984696398005766638 1.0525823128048024113 -1.0520069441187036574 1.3735841197149245563 -1.0307914587329112255 0.45265853085703000591;0.27119896261202480758 1.4919298074603053461 -0.078581126849476920504 0.49847738764222659924 -0.89244271641859718169 1.2839117177763981203 0.61289201530832504439;-1.1077596151796123003 0.43083003935085867342 0.14136509198805086163 -0.5724761677306803298 0.30920247025063601898 0.10675284442257694517 -1.3134183154158398654;-1.8636498375299175301 -0.62699090859522688124 0.27080501105921228611 0.02141749515520788294 0.19599693544611737628 -1.2980917124505293447 -0.70726148479468686769;1.8067977494980185238 -1.0586184502768705151 1.2831596090032415525 0.77769746153775465292 -4.4870449110732462117 -0.24178681307140900691 0.16123183555010309198;0.86252812184359728409 1.1773961603834373069 -0.43450690901519667397 0.048141899476108619693 -0.044664420270658339307 -0.81921518870998244477 -1.3521412887440387163;-0.89833752282572176817 0.25099820501524089389 -0.44255699910295476229 -0.044973657028602773145 -1.2140799488084195179 0.19320362930600132323 -0.35709292667872855054;0.20495562065121861184 -0.65029440829443230232 0.3187586590793221597 1.4933668150520358342 2.5056674042879900632 -0.41104610974776251142 1.5865492615407856825;0.23997023096883940418 0.94835675858129320215 1.1321436932215882543 1.8483354485897514063 0.37312415888697386679 1.687274418030515255 -0.48678292398199590751;0.36996680424863975789 -1.1349459764213982726 -0.28440926356021462018 -0.52194925223452992569 2.0857026645324694414 -1.3835004055075517115 3.8418232356575687803;0.1864867953731774719 0.71690347431520740074 1.692776383585607558 -1.2981560904894549058 0.61806453291487883916 -1.24565329795153934 0.047133330030912405639;-1.8896455474756264081 0.3718684580057041722 1.1417630339341402479 -0.22912417908420346091 0.67771152446987825702 -0.70683844273631879496 2.4547020655930769095;-0.80431900179115600746 -0.82022359035297420782 -0.26054277035630285519 -0.74670369338824293948 -1.7080594421877210998 0.66789458810188828686 0.35540298728142327711;0.17992006288943226955 0.28983305534600656284 0.70432929860325210747 -0.60203365761499560982 3.4130097432531396606 -0.48244113098414243002 0.30539061260924871277;-0.95927457523042181098 0.35408353822535770794 -0.56829592868266198824 -0.28148685712376597667 0.8686818833375550275 2.313673906158792537 0.20189195588705555195;-0.14872548068570978774 -1.3629446012939230304 0.023760620770693706427 -0.8929041731614454358 1.6074083394881863995 -1.3224563515743825715 -0.69192042281013554916;0.086503589795781862803 -1.9299411006788644851 0.75666076912827073819 -1.0795095500473592676 2.284724668104082923 0.52601587217206224878 0.91249704782804541559;-1.083737812346400764 2.5998539372551467963 2.2215045048543826667 0.24147790930994614267 -1.1138063867185832301 -1.5191468448817646575 1.1352374551506192457;-0.95344009339710766859 3.2186459615439035886 -0.74079657420733369122 0.62523486595353794382 1.5021705768777220413 0.6356307677467674111 1.9314700074591337753;0.48421517928444263168 0.82173525453421036868 1.9731877153491967913 -1.4127666546159198813 0.10501672820948158849 0.51674687663002694116 -0.66723406330659262853;0.75996751602373824319 -0.098572864172220539536 0.45008648999846007088 0.37508857920805915853 -0.92823568708913839487 -0.067623762920562199841 0.80019919545146089845;-1.0059359031215351354 -0.18815503686726642307 -0.77725585529807816876 -0.90052784635441562866 -1.4909809650223557753 0.87715879936669527428 1.7505286779071380643;-1.5748511456431428801 0.51852869057017036969 -1.046037216364984701 -0.01664169100879439589 1.6387037450811035022 0.93783733042692707382 2.0165194222830873372;-0.48153849491592465304 1.0621709491562201233 -2.5908229257496895315 0.22317476654861181928 1.5705324877334341238 -0.38157181865827544121 -1.408049912522702618;-0.65925604134958570857 -0.33427624637339575786 -0.2650605424626585882 -1.2473572149486391591 -2.7615403452624334157 0.93473159648553327816 1.5831132942796652952;-0.17643087907522556068 0.0097969066215340475101 1.7190148859822935723 0.090850648354882879154 -1.8216674664483045731 0.44406232061177436377 0.63759746816379914147;-0.46693743851061808581 2.5249257931430237889 -2.3712883886143529821 -0.17614167189502366595 2.1520437360125677095 -0.11993096102048256268 -1.7700883790181436428;-0.42320442027896643111 -0.31010656004487680004 1.0130651335536804414 0.59067946231215351371 0.21731507420586657076 -3.0848481738927078233 3.2436848219366090085;-0.9592271095100070033 -0.89925241404729083783 1.0525221799434185677 1.4051890256355019648 -0.65157915688447565206 0.75706998518960078037 0.97574260473662555171];

% Layer 2
b2 = -0.69653510523026274104;
LW2_1 = [-1.4406389540088697654 1.3934499020497319588 -0.086684607075579095392 0.2350708413902827254 -1.7376842062604624939 -1.0949144495942000344 1.7699641472638290018 1.3424192545964208634 2.2594466276613687228 -0.50219618829149781369 -0.21729139959224222367 0.85840623440745278394 1.3980470163756373303 -0.78954079698405055243 -0.52705798995850927646 2.1797694527516036089 -0.33194135272709412998 -1.7560538569074608528 -0.77576342265883346805 0.15883251955186891702 0.36650631333835953507 -0.37018320273848787094 0.090185339366710873255 -0.72108797046807859843 -0.66880081172356276031 0.25984738782775113064 1.3056514644609054088 -1.3071091604135394171 -0.17681983366357820198 -0.23259615137187278355 0.08367354586626288282 -1.3202343970412613317 -1.8404682340184297384 0.58415477120891390328 0.2756410797614654018 1.5054738252895487616 -0.82607121004634809225 0.7599374072708948713 -0.85033370339997849108 0.47603436137819155061];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 2.03830178109291;
y1_step1.xoffset = 7.47181919175367;

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
