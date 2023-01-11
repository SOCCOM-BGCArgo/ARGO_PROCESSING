function [Y,Xf,Af] = ESPER_talk_10_Atl_4(X,~,~)
%ESPER_TALK_10_ATL_4 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 12-Jan-2021 17:03:48.
% 
% [Y] = ESPER_talk_10_Atl_4(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999999390765;-40;0;11.29;-0.2178;-0.12];
x1_step1.gain = [1.0000000017546;1.00031699124092;0.0153846153846154;0.000310077519379845;0.0705602618858232;0.0479365703301392;0.0145655815308426];
x1_step1.ymin = -1;

% Layer 1
b1 = [-2.7213076168642525;-0.19524193309142909;3.53256816085432;-10.435084486554231;0.26793119570155866;-3.1177026483907802;-4.5041500981396547;3.4113642767952852;-0.65947950062998706;-3.3315697544143106;-5.8960102409085779;0.57124663139830956;-5.0996249563624936;2.7805926833331207;1.3380722743073445;1.8954621936580502;0.024916437662429614;6.1033031328158298;0.69126594661190199;-6.8071765034028093;-4.0887509997913458;-2.8849933978915465;7.4373383842089789;0.33745791072886744;4.8275519588587326;-4.2532697512553224;12.685289213552334;3.5984930612395294;-5.4821790181210472;-15.44948709189056];
IW1_1 = [-0.59985394150301174 -0.60458351808190636 -0.16512122362933399 0.3598077168797213 3.9415153713569322 0.62708140491325604 -0.36252956461935637;-0.96225981906786562 -0.28158458964983818 -1.2776475910073153 -0.0024621647534891033 1.7610395452987138 -0.063615091925858167 0.10293688278868957;2.1592548842238837 -3.5941060442567565 1.6877071715602634 1.4984594358030967 -10.597319873390077 1.1576288842511697 -1.3228940883989124;-0.59360538917167405 1.3271024091600581 1.3228732868187456 -0.75500240777646532 17.992198879034472 0.42864344486879857 1.099734573791801;-2.5045859066744351 2.8763433991424523 -0.2742587149540629 -0.48088798891894691 5.6006094566426912 -0.48672158233664431 1.0627864199678509;-0.59512580772690005 0.27239355920546537 -0.20511636783733153 -0.30250686480295996 8.085578219010479 -0.42972145237289178 -0.71463322692224107;0.0076650051307427134 0.038185604801617222 0.94203009693252426 -0.12161569087522728 3.02431391221702 0.45049269044486256 -0.68898780015756167;-0.75760862781088723 -0.91787221109271988 -0.57730119088051379 0.77235773735242974 2.6539435090458117 1.6683972657385817 3.2757040260947079;0.79325425986233 -0.05783667440785293 0.94597945869553246 -0.14805220719173037 0.41897643662726569 0.18804738786562597 0.091026401375148322;0.79207333282977588 -2.0001429347234958 1.3576365532712706 -0.28299569768660754 1.4617154386915228 2.5908457729582532 -3.4635345674152815;-2.0964458656989318 1.3837799497301595 -0.42188204877178814 1.233516349153871 13.75728119319491 0.37236930067789875 0.66282097666962791;0.86744578628373881 -0.40196493568286634 0.50442861637580727 0.075184964099683227 -1.4107393170247453 -0.11624327487635522 -0.1830768899902038;-0.5304905772095353 -0.023585703489589934 0.48880461966405991 0.183966680515294 8.1531538135239856 0.35110875337626535 -0.12995512996896652;-0.4982420054669186 0.090683656281306563 -0.6946655054237143 0.054725097153132937 -0.56874573086515656 -0.19174172777330742 0.063966339454967519;2.2685268659270466 1.7069041065022741 -1.007846909752746 0.48590550176250574 -7.2750394472766313 -0.19354402266278029 -0.19403343262796713;2.1430675256254754 -2.8083614451514896 1.1901159014787068 0.11512432280204771 -2.4206766855800677 0.45707189385640667 -0.64922335036777923;0.25990289764176794 -0.099902983796326783 0.39732386545349124 -0.0031037620682735365 -0.97899251300356627 0.0060307270772764964 -0.12253529151053057;-1.7718436451270845 -2.9890892559437683 2.3703056213630425 -2.5756644461253488 -8.4867398430220344 -1.2633776204063891 2.8783321794087322;0.05114591172706482 -0.029477882160390836 -0.15353958846778545 0.026080902533106437 0.35273384120501589 -0.0058312473084146336 0.38688224849951858;-0.92159226693703955 -0.054015565139866939 0.82215837898217936 0.46826300102935881 13.517231092443607 0.49302450086847743 0.085891337938128143;-4.4991285999961885 2.0530835551347133 0.56269080832746798 2.3426112428571533 11.093408382593488 4.1781441421190006 -3.130599807751226;0.13534158834350046 -0.011076011001529554 0.639402958334779 -0.067588728764530023 1.3846504749204345 0.28157185905723769 -0.32301647177648468;-0.78971438934444782 4.6927371887171088 -2.1601883820655989 0.45520995927415486 -5.7215917555844991 -1.7553091579551294 -1.6036863666406018;0.029055791475755474 -0.0049863230106927227 -0.24793991789514916 0.029992844738951332 0.69121328967834506 -0.015702868172670085 0.41726773864558936;0.51389423774113252 0.017708536275613786 -0.46269745418514119 -0.17459201742841801 -7.7413040315352601 -0.32632312571239414 0.13484272449085111;0.63795842425478277 1.205111523410076 2.1252573831410442 -1.332956335779369 15.529498189346956 -0.34224746156605318 7.4927984533947152;0.64484794602302808 -2.6524832148383322 -0.31246405124866017 2.3067368708793072 2.5566220087734992 0.8544730667851077 12.258904562736868;0.24784012647916287 -2.0081939816011216 0.38394381960257151 -0.55603563850690529 -4.6531578726513239 3.9357636175444415 -1.2362898403316045;-3.0730077855250437 -0.99476716533499787 -2.6289848692525202 1.8103469251452422 6.4686105270887699 -0.25837850802688594 -1.2701679831842017;-7.4416638259404975 0.91544179461132691 -3.2468757967074184 3.173313576767069 27.642808637880758 1.5410396951187839 -8.7676819198153755];

% Layer 2
b2 = [-12.443151917940849;-0.35928562950347565;-12.224636357721495;-93.097480094277188;-110.96797265077342;-21.112417905901811;-22.62966462659606;0.16477864930139785;37.973771209801122;22.047638343172117];
LW2_1 = [-7.906528779236754 8.7695948241200714 -1.102299793001039 0.19444828393382851 -3.1662712762036476 3.2423485520666189 -11.771240939946171 -20.430240693326663 -0.083786951024526199 1.6055561390002981 1.7089385662609777 -3.7156321637378396 0.51836364936524681 -22.730997444888921 4.6577780828453124 47.137905987105476 12.883655335950898 0.84658607577760847 10.526248920771407 1.6669658526372866 1.0916748591836227 26.007061813492896 -0.26445783569290615 2.3062878343855817 0.047726412771715861 -0.14261704287088112 22.3839082957805 -7.5326687048088816 0.32423150924630573 -0.15401484917340166;-32.019662329067636 12.409292973156386 8.872342460125779 23.771805910195006 -28.223561977509011 19.730491966020704 11.831125306206498 -5.8665649667936721 -68.58723287795695 6.6551404084224508 0.61690714335146979 31.812400391790721 -23.54664553068233 -2.9562082712205457 14.374083363747211 -58.595545997205754 29.582831638700846 9.6719067134933621 -36.983787360629144 -2.2533721034520986 3.2313539461757035 19.825183784219988 -10.821982245453963 13.167939670499955 19.179727586811328 -13.384160998978805 33.141963136671912 38.79693759071943 3.592193138107755 -1.5193670280924843;1.7797186813752235 -0.18394900770066455 6.4076730933217796 3.2657761229811659 -2.9572576304426512 -11.491199273360628 -12.75983974029514 -8.4173483101548001 -79.324001480815966 13.683314930374761 24.817084897213853 34.454952987076318 5.6914774043402554 -47.708275950836601 26.99946814038497 -27.988291860091373 -2.0682354713269784 -33.900919948822448 14.598292994647766 13.015644731370692 107.47411771552008 15.403159320294355 -8.2919535383012004 34.670579009970261 -61.374745425596387 -2.0604848386446499 -53.373455922780273 27.825277515160639 36.415729935671635 -36.025520899051841;46.049359524417333 37.1403666089805 61.542961469700828 -33.58995366403736 22.872406678789464 -19.106072296389353 57.139629974844986 -83.613624310605104 11.407022067751964 23.224381444849644 -59.89574571735055 22.838038995809114 39.698697928877706 75.554816934848773 24.036314504660034 29.214109109389188 22.846169068335598 -2.9148066700754764 18.926838390888729 -13.188243613946952 -14.721655587129002 70.112682964839593 -0.051798885421009543 66.256142033428844 119.35333450825861 5.8951020532346812 57.643318952151951 -37.339505348273292 -22.792189983779252 -3.1889143972598482;-34.317593240707637 -16.26698015959159 27.922238230375005 -79.001059603328116 33.80598795557308 -9.6520591228259605 -21.583156110306501 16.605518660642673 2.9049172951958417 2.2066462241949378 -15.627359154345713 64.404328807893464 85.745562352801699 -36.85721614680282 -94.167035079950566 20.695368225339912 -23.568139480577258 -2.1224335901739604 35.402669791636981 15.345969274907681 14.232147170817338 56.599567701522197 18.876660757185288 -10.146272344119172 37.056203124965819 -4.0088206939251316 -19.361959610387085 15.424773471290482 7.6819432362293716 -2.9305789320797242;5.798726574440658 12.925103955799141 -16.800879558595035 24.969429896753194 28.072360298333884 11.512094214525581 -21.000645063120594 8.3185420025305934 13.126430424643813 15.202128205465129 -45.250576061681521 15.800892559690229 -21.393527141554674 -13.930420647836069 38.934529958547259 -33.993454268001408 3.1949778691987367 -6.7429059489166328 -49.998153419728752 -12.310670667816064 -21.727065445390515 -53.852568194401691 -26.171901654277793 -27.357550098189023 31.091621488701922 6.2323701587179396 -17.974911012955285 0.028544760475656822 -23.26079782367562 -7.6046856141821904;37.614318813747701 42.560607879928426 25.495915049032121 72.488561268547272 -86.037411842002584 27.298799070792974 83.652515860480463 13.182664234787667 49.341428778032018 -13.643034273676832 6.2780083122192529 -10.3148783814966 -41.900468503370327 23.563907113467259 69.631949394041754 59.386871769703269 -38.824144545394468 19.826902483076601 76.908053103943715 -26.847355238412522 -9.2623265714188712 81.726740832588192 39.670159702475729 123.62444610493856 -64.669370458989917 -13.630144517316687 -0.52901720903471749 -12.584273012517132 -17.272467235584124 -5.5971145723506215;-0.050014504075717843 0.25574565163728702 -0.011084219279494848 -0.0057468909949865352 0.012886437081119898 -0.15071991682212027 -2.8901855550902429 0.020039028027787325 -0.31937930989077906 0.017598856612415181 0.010735167029079701 -0.42612543481443516 4.9220314383423878 8.1068172196561754 -0.019225940004135127 -0.028921680588736818 3.5586572301496968 0.0078712830539394586 -9.0316661264930147 0.15382090954367167 0.0037046965325789479 8.7493398080964742 0.011163705679393567 8.1250965205175483 5.4101522474045591 -0.007717216520453706 0.0028728600891528002 -0.0074840785825465517 0.017277208711372609 -0.0030988702291825216;-53.487815030419711 -46.313198228335786 -45.046974221991441 -33.60152873669626 -54.754201107773682 1.1290899897901971 -35.315101624358448 -4.8085875345635545 -40.724396545415338 19.27017452803403 -4.1299693032792586 -18.687327949344056 16.010892616579465 -16.307342080114527 -30.330687372034767 16.575249501366585 21.446227087041887 4.512974056427411 -38.812081762934568 -10.323800875261149 6.2143170892882189 47.353093720217103 -25.370127947014446 9.0747017433949271 11.395534274413674 22.518569428009243 -6.0089112598365126 -4.8838908655820061 -12.249988533529814 -9.5630261010545947;-9.5924622637130117 4.4445389915100826 -2.6852978703922878 -18.850733881258449 2.9095418527931005 -8.3227487254026862 3.8328534563479888 -27.258421382177893 -30.420681520793003 32.653986138094908 -32.533347154893434 -2.0145784715746426 1.3293278617481967 -17.017999066557746 -63.778529223060694 -15.161677956165882 -6.2719118859964684 16.500039411131262 13.749680654994402 5.7324993353417444 -14.166748729378602 -28.27263112594401 10.41357153424627 0.73229083060799516 -1.2188132288516047 -7.2117975977928568 8.5196923140210625 -10.452033001399499 1.5181490314127117 -1.8350231126750918];

% Layer 3
b3 = 0.11315309930419923;
LW3_2 = [0.018576892171498573 -0.002075268411027504 0.00070274002811569715 -0.00041069210815829308 -0.00077504630692927349 -3.5081975872209923e-05 0.0013402036302874955 -1.277194797051046 -0.00063605565928800489 0.0019756531588893799];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.00111018595614765;
y1_step1.xoffset = 1025.5;

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
