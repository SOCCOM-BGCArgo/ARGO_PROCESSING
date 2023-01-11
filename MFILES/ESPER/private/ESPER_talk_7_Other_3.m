function [Y,Xf,Af] = ESPER_talk_7_Other_3(X,~,~)
%ESPER_TALK_7_OTHER_3 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 12-Jan-2021 17:03:47.
% 
% [Y] = ESPER_talk_7_Other_3(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999540494077341;-75.751;0;29.715;-2.05132667284751;-133.803853032615];
x1_step1.gain = [1.00000001035699;1.00022983623563;0.0140762793578401;0.000302617642608564;0.251952632905014;0.0584431850247228;0.00432440760873659];
x1_step1.ymin = -1;

% Layer 1
b1 = [-4.7866977563659701;2.2587217696664958;3.8900325981600425;-2.481332828261781;-1.1861072971582851;1.2306148897285683;-2.5831855935365731;0.2780690428858153;7.3366301351105028;-3.2964893392075538;0.36695738903780339;-1.2829794733440416;2.1839259670741247;2.3325228267468905;-1.5458493384234775;0.69608378361612189;0.084585974500337668;0.066993574833213984;-0.67180702264266889;1.5612901852834435;-3.781830145741909;-0.49271365673289219;4.7763764691003043;0.70187951506314927;-0.21622922822295329];
IW1_1 = [-0.50770547578629954 -0.4447333349003405 2.8089109837677122 -2.0879654481602348 -3.5455657974555828 0.33850128140579849 1.0684270603634267;-0.092218170695238411 0.18739221197666614 0.12083178864423391 -0.74924319375460713 -2.4736331856858853 2.621181235694114 -0.12915801087019771;-0.058241546093747966 -1.8888177115602622 2.0025217554316015 -1.5697249682018248 -4.4274199744257121 0.0032243495493209013 -0.1813641688189874;-0.021373192797321826 -0.096139652956670846 -4.2979220984355244 -0.1614169034880045 -2.6679772282425285 1.5655804595533878 -0.34248015790841824;-0.36337971554101306 5.9508825945573811 2.3342584251774174 2.8583611383840428 7.4409681359952495 -0.20584522904805772 -9.3896709697781198;0.3228978249105402 0.14575471510927199 -2.1853199831061478 0.67288938689940603 -0.30613223609484846 -1.0883845523810636 0.19351010104999944;0.12428593534213225 0.28510491624789935 0.69307002782505545 -3.0632561415125199 -3.0068095740042584 -1.4839054533275411 -1.5297283650094624;-0.1656341786660504 0.16308767646711475 0.037425702643840555 -0.16036663619308514 0.68828751648276254 1.320874081700752 3.327957670633706;-0.35498783092502584 -0.7943083335680009 -2.5787827324602519 3.7572356046866786 -4.3018162321239384 6.905293068800832 3.3030733513867871;0.74230243941454843 0.81917933578720115 1.808559570766832 -0.29998615025490394 1.9868097755247998 -6.2862722894757219 -4.0845462038213816;-0.36720111815071743 0.52736200689519619 -0.90748210287717745 1.1055253697586376 0.84559645616146695 0.22169582579525957 0.61220653998265173;0.60378991705346996 0.077945114729170706 -0.98738767807861072 -1.0428722440509399 5.5144965657580229 -1.4107938102507596 1.615154215366633;-0.041120164610168104 -0.42528966584401839 2.1110934906155903 -1.2753694767909802 -4.8856169315551421 4.806795020141247 -1.2756480527158691;-0.061503366601561277 0.044263798503585722 0.04997510214496103 -1.0559837306093034 -2.9594023287632867 3.0358169521202112 -0.19469763501735976;0.13513910094891057 -0.089021338198651934 -0.13546481412543129 -1.6996776730298224 -0.86890172848401648 1.297109146574182 0.39531176717904704;-0.25749580298791497 0.22129252174067746 0.044723789191459966 0.21491750905108326 0.84438959372057498 1.5663946932486856 3.7950396913245474;0.33663077691985704 0.67306189673139705 -0.48251448122039137 0.65194754098783347 -0.79931046177514098 -1.0145303454017804 -0.69209743262104451;-0.13322894194109469 0.0088018467284928289 0.97033795971528747 -0.17444920122939869 -1.5178180498070226 0.68214355671654459 -0.13365728140595909;-0.19454196431373302 -0.11144057276379131 -0.85099526843686268 -1.4060401540542931 -4.1891146844800256 -0.70019046597995449 0.89924472628744223;0.022108755479184887 -0.063365353549535472 -0.087627075130497084 -1.1348667680498736 -1.4652174911363387 1.191389214661984 -1.8723241165297573;0.79913452989757927 0.89165859457101881 -0.16369302670790714 -1.5143915319900907 3.7499747640077663 -2.4825107767708339 -1.6883084782876818;-0.18542699373528629 -0.5201986617126434 0.35058704075892844 0.8077797193246492 1.6443324386976974 0.63324343013998052 -0.45936852094308805;1.335684045437705 0.64727598038292555 0.22957842734338338 0.24900954374062673 -1.5012409759699756 5.0609368407421931 -0.11033184301555503;-0.007735679308500931 -0.0071107247381800688 0.59569145753615016 0.84262426870915919 0.44102053309528949 -1.1739430611872068 -1.0226886527666645;0.12125103492491127 0.058378032650666087 0.54108831477613661 -0.86029913156512605 1.8782963459657436 0.006987558582960786 -0.069291900448021901];

% Layer 2
b2 = [140.89568677710702;42.471715757340398;-68.29585543486219;-53.570295582706372;-139.36873992893644;-365.46661389217127;-1.4055391585685326;230.65274723623;155.28671446760154;-501.24008985843739;30.201918072671752;-147.30397822546882;9.8582729676403975;-50.239177028994625;-924.89616625871929];
LW2_1 = [239.24437256227912 21.02832643595195 -19.69789038902378 -3.8520565302544014 -19.605599729819073 0.0086049050395501751 -15.355471659540337 -83.713176112932871 -24.140795024954059 6.2691905345296339 -44.5806417053982 26.197468572696323 -11.352515285946724 68.690616595095491 -16.289080558744448 15.887350401280804 23.761751533616263 66.247458218630541 -13.090650541699082 74.451957154251758 -47.388072049304526 11.29198264218245 -5.5651785452548657 5.2151467123073267 24.448851601653892;-0.54832398503665081 12.197980562634346 6.5115275432439779 10.522483283416337 -1.0898818919730364 -22.712879639233158 2.4323985448802095 -17.969539809807898 0.20640938138588408 14.785746867581668 -26.924454269363945 -13.20996638438762 -4.5124250338531979 10.985758786598504 39.910788736676984 8.563132928331191 -95.456849931044005 -83.376957596116853 20.99481415180923 -1.1763920347458676 30.999393249978358 38.777148415891133 -13.691904513857928 17.970408122621755 -51.328741719262439;39.686746994115168 -13.435021871456659 -16.245137160017826 9.4955740954614516 2.8357658663531415 -46.162618645087555 32.792164518748656 -55.427997413876689 13.903092672738902 61.430723227817623 102.99989745549387 -61.742103125536559 -73.512533589354405 22.302098758364895 19.44767807852309 34.848988834879975 -226.10433675577389 -154.62058016519245 -46.617976308617621 44.116534815887363 9.8425865415486964 -19.289575314478224 10.007868647548982 -51.085929476392742 167.70490038172832;31.271573642563588 -24.53322264736012 -25.145082816533264 -54.724150006944093 -29.215958105334071 -38.93488881312625 -5.7107004282278462 -186.21749327891206 -24.546684123618636 -91.989019608798912 -43.538884151871258 112.61225333711795 169.34732234421912 69.669798485084343 -86.953114307334019 -3.0722382287731995 10.451505725146264 26.763275448797614 47.303972025093636 49.511056492918961 -50.821340323570865 9.1823658406603652 -35.458472444174781 -86.975894570648137 13.397794204951099;-255.60266625897927 83.338321762881648 -134.98126585857761 154.30549735255403 64.071470080124072 99.617156064904677 25.5658242296116 111.00021621651997 27.779015552628792 -15.704331079517212 60.636938290130118 70.097584751757338 228.15430233873749 -129.64270425433844 9.6496508819243108 -123.05431487762208 298.39481236602995 -194.55316759489216 -222.03814656031355 -151.79224665942587 -247.15886973128238 219.99889068192508 -7.1151536976737813 44.08303816846233 194.247810833942;118.1314991170063 13.791618353838032 369.25366810243037 145.0022561771994 95.414685596234875 39.805464018123075 88.166641578048228 163.23070740370372 -56.428149431285419 27.426114931762921 192.78814993692757 -22.846674805670443 7.6380586250077789 137.52154502684823 -66.873623050513416 -21.711914325688227 -26.497879134451363 124.45329002999223 -184.35136982908401 -79.67690413955242 -31.029718720981478 -110.69623703886951 76.363626591292388 -84.569117674333583 105.66816322428558;-0.02361308227904876 -0.21895716702325627 -0.02093076732592912 0.027281583740262894 -0.0017736485760639418 -0.045363075686626661 -0.033992668577530818 0.15635304399937894 -0.016639144276127116 -0.0069297246102274631 0.05299986563963216 -0.018874288141694479 0.020436220542934067 0.19878711629338516 0.14171041356117503 -0.11832120190167524 0.034868332177078097 -0.15255436204957604 0.03082340221729021 -0.084896156418132335 -0.017092500161711082 0.039983706640567211 -0.012635642686472293 0.19272059610459363 0.19009850661281522;197.44442113562079 -46.877585391139853 73.529478649349429 -30.317812839072698 39.302825039686212 61.287945351721071 51.62120356909962 190.51044152773005 -26.933605776645862 -2.3612146555161586 222.23445714530351 -53.688538134815929 -57.233462686020971 -145.17628514626813 40.721345785892609 45.089756646239138 187.09251201283604 281.88145409346345 51.662577223120458 0.28802023448609626 115.07025376218195 -17.248975788631423 15.629754786763804 -26.225127287406419 62.353094987058782;4.0293411694077514 -67.804497028999762 11.49696293553488 28.011105916892632 14.313065132960052 51.428338270145773 -12.154437679969851 -45.238468102090636 16.253475591559649 -25.750766424533701 -101.52736984373935 -45.607569522840137 29.616164088808667 -128.84605807895699 8.2935741580131452 12.747257476638719 49.167696832190771 83.838927051128167 95.833096441084905 5.1154726979417875 13.085313698833087 -81.20480931985098 27.576392209263346 -155.00594680284669 -11.995921757507331;-70.636117010593537 -21.815673371004081 51.85408982657539 157.36664168920419 2.3202408094473173 -36.265915137413835 1.7033171386771955 -65.140214783454454 -29.765658401105199 -45.529613304044197 171.56429247239805 -1.2288643118708833 -36.750219870246681 213.46993870206111 -58.558153949951461 98.278139056733366 134.04829084521515 47.235838501034557 -114.85825692284331 34.852350575715199 64.219696814692966 -26.130980248435442 -124.53322888411181 416.87458757568839 63.307783226524904;123.27012810550235 54.630066432933759 -44.232526231635546 17.278436669120524 -32.812934229299152 87.734286357340864 -10.702131458933387 104.90373517186764 23.947039227211906 51.693830482925016 -62.486348131276635 -21.838159897068994 63.372650579709408 -60.52332084861218 -26.203732478207336 9.8274541584062458 58.446842471545111 94.073704529841208 36.518502517484798 42.07813709636136 2.5207151483340948 -100.18104958417889 58.614289439351616 -212.97974102871814 125.73318400677273;-88.576568214943961 -1.9473964633404564 -16.539463836706251 -35.811376357615956 1.1202150251742364 -8.9247916510899543 1.4801114812786207 226.34899605696691 -11.227649840376557 19.566377373988029 78.952818597401617 -4.4064652686733758 -28.907737923232318 -3.5383138105276961 -29.267069176164277 32.783206743130421 103.37323099064926 -10.401184668093258 13.098611578021732 116.20609787693749 9.3749765695970879 20.625264727681817 25.360481911089661 76.355858789556095 41.709321711389478;343.40645367516839 103.30168672675886 15.322849331935158 67.531834087066841 -104.98530119919774 16.584908862225042 48.824482274585272 788.48319249329143 242.09391064846031 140.96029179186365 295.4226789838159 -94.264279956566142 -81.369373251348946 -10.677943666537013 -216.92554211853638 213.60080633564621 -293.12633015080183 76.503871852654328 77.881592108764238 -335.03525526987818 216.41053118351081 -182.56503237210865 88.746994284397758 310.5698316433099 200.39123834886618;-11.568919411694068 -0.87157996337822352 -4.4820476752378546 -0.65015514637292982 -10.991177824777001 -7.2968019182486898 -11.396694104818794 -129.10810920329064 -6.5466350514726424 4.9668811966902968 -44.004428725221096 9.1004988827275124 11.62204870955313 11.671467379686977 52.361476892868104 14.066127026859851 -161.45519270620758 -136.32472879859574 60.920825117723851 58.837950104269993 11.143112320001155 -14.352822447457346 -7.3469704121757164 67.190023395885447 -85.646801133731003;-389.05122042065631 152.87042314052721 76.647244742045444 -2.6204269827558284 55.687212469628498 -58.501447497741921 27.02638074266649 -192.11983912998704 55.635113880247886 -63.317347614794777 -0.89693729254189636 -66.610910194858846 4.3638828794866624 -186.10251878334901 -12.520632213671423 112.94333012190923 354.56760557581794 288.87941811671618 285.35909812803112 144.3477586731789 71.321094155113371 70.454761306320037 210.67131544486426 37.942556168814392 447.92666747281515];

% Layer 3
b3 = 16.952425513186462;
LW3_2 = [0.00097695367679216205 -0.0074862918138553503 -0.0040781726321698899 0.0038075290348015258 -0.00026354445733808602 -0.0018868940178936933 18.872469640737549 0.0016759920907340935 -0.0042438291796818208 -0.0049274512528376183 -0.003070058345050209 0.001379319415652431 0.00095243174350165044 -0.0042772218181910083 -0.0019778438447705031];

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
