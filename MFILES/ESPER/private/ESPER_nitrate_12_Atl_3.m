function [Y,Xf,Af] = ESPER_nitrate_12_Atl_3(X,~,~)
%ESPER_NITRATE_12_ATL_3 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:36.
% 
% [Y] = ESPER_nitrate_12_Atl_3(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999994516886;-40;0;9.1619;-0.52];
x1_step1.gain = [1.0000000017546;1.00003626070655;0.0153846153846154;0.000309885342423303;0.0656325902039552;0.0145232735458572];
x1_step1.ymin = -1;

% Layer 1
b1 = [-33.013703764377694938;-1.8730609672095699381;4.2699508860465629212;14.184743883000757947;4.6129907534048717466;20.222317905712195341;1.9538273149383615923;13.746963291205844371;-1.1676179339072960683;-4.0252004176424778237;-5.236446129259467952;10.527668734103373893;0.50132225549189524116;17.989944005481248723;-2.4990511996643864912;0.0097854274157865555467;2.5453692310158668732;1.7995457648500199355;-1.3034054709633275948;-0.1962811272364390347;-3.4952834932116090272;2.2200713878994009498;10.230603519810744118;6.8233970198206606028;-12.620729076457660511];
IW1_1 = [-50.264833645849890331 -25.19262545813617038 -15.97576872055657482 9.2341579964592739316 -4.4460565476269460206 2.2895157533184775644;-1.1977000194634082053 0.77578451394587311096 1.9062487692434189057 -6.0515164706240334169 -7.4791499284183844409 -3.2538663632464910336;-0.18255848462997487092 0.057457117415605662991 0.11364390751340339158 -0.25960346222451297082 -6.0988963426662481382 0.51064454799539882313;-4.0518848210891835748 2.646558923366526539 -3.7179338115962679012 1.6947409994105457898 -18.685700078904407206 -0.55974006933650399809;0.87371685196800552564 -0.36712619005978308628 -0.93211913292995107749 -3.2485816607896960129 -4.358631323439630556 2.9346588456443498671;-9.817766408326587424 -3.3016558546288132803 3.06079246090803192 -4.4492692053103048977 -15.844705817498004663 11.007776613309674829;-0.52936846598467901526 0.23032719770638226953 1.0425216963972749884 -4.3176797084583036224 -3.0748660396199687384 3.9260718415181545815;-1.949838994679289117 -3.4291974723852058915 -2.9224681061159918372 -1.4225901111335113036 -15.074972872493274778 -2.9456902601537628961;0.4761712204668522852 -0.053746131060469937157 1.118074714248189272 -6.1655771486674453286 -0.74661605309090506299 6.002729004504795185;0.13737688778682430324 -0.063874962506100807502 -0.060151859568970650338 0.24437492332969912301 5.7287981406231152093 -0.45056673253806844848;1.3074809920888594839 -0.39160843094858999169 3.5210344423289199334 -2.9448812401934425509 8.2537180489611241541 0.93266764363024468043;2.882687242247156334 -2.0349721180958528777 -1.8742326204818935409 -3.1064890139448011475 10.212555837967913774 20.577285497904401979;-0.17675083087890142153 0.22125888020090789254 -0.40682085136712026197 0.14684604236725054549 -1.0624728023591996084 -0.032316616295257935554;-11.629133423475797926 -4.3653953066163886376 2.6534121787086770716 -2.7767713197447942974 -10.292682123040666653 12.438472281661750785;-0.41178465197508995477 -0.29409255801316613876 -0.39538020042887522809 -0.34470022033287855967 2.6019863889613215058 -1.0042823398197331386;-0.77806934072084366072 0.84618587430166802843 -1.4299409698954586023 -5.4236226592779814837 -5.326339630730124064 -0.82222430930220957723;-0.5242137891912246328 0.33965859126950959102 -0.41761421392350245307 0.68879768593961820944 2.3765345267505435523 3.0590619554409581582;2.9309611587677313338 2.1027245142219097573 -3.2525433343558050048 -2.4046739717152143534 3.7994178598638814748 -0.081031506621127805112;-0.75678685584126204056 -0.17560261156259457382 -0.29939089463723261852 -0.069660051124728158967 1.5012727208599601081 -0.81766662106289100809;-0.11514181498748306265 0.1020908486839576107 -0.20114859079812164033 0.090807445272223660804 -0.38686417678897927486 -0.013421888961213091634;-2.1381554234643731505 -1.1902239163832366575 0.15707098377096082031 -1.4937910599098389586 3.9036304991582184698 -0.68120523841647073748;-0.21583032778229857551 1.3438778346345827686 -1.0709754335005117198 0.11523399778056080123 -3.2848865135271485549 -0.2845898132850420148;-8.1204754236778402543 -15.69265544398163037 1.9344526872522631589 -4.3117685429676706477 -7.8068311678313282442 5.3872428523990842919;-0.36971176637482633387 0.60327857633153014749 -2.3131658945004600625 0.71871222266695533598 -3.3892983957645910209 1.7265188107466520329;-1.3240792744822504901 1.2589971853663999291 -0.13180273457011432159 -3.7119055367846431892 -2.3291993292740955646 -12.693975324193321086];

% Layer 2
b2 = [16.336861733285491738;0.57375754733692818643;-19.834638924819319783;-2.3218413366893124739;-22.2951953645475065;-11.584562843924885556;48.028662469635435173;57.983143789204355301;57.876378932729622306;-15.76348029048757482;-55.586274465578547677;19.128254894480580361;3.1248048796627405643;-23.180082281867658622;12.404120012843765863];
LW2_1 = [-7.2581231053028156808 -71.46302458856648343 30.97541801078699919 -29.102869069842117966 26.067320263384527124 -18.934587318131711697 -36.851031417708313143 -26.66970189510930922 127.21604724589924729 86.189210846440360569 -37.862993738545974054 30.21097619323907324 -96.153335907623301182 19.924008225730403865 -126.61714525483840532 14.965355112645175595 44.769297347595419012 -21.566944874534531351 -38.463145872066888842 84.488353062160342688 24.2925144196249434 -3.7039704634666805738 -23.416673955659085493 62.318605129717475677 27.477251963095095988;3.7469192498032652949 -8.6812798612592434466 -75.04315584219618529 -19.077075921770877187 10.631308269533295885 -101.80592255600298301 2.527463778482852419 -12.873457277391864295 29.100912315489587456 -37.03748127695537562 -35.209071635067445527 20.516791058842336781 32.258826961398590072 14.239975748443470138 58.588671222672303429 37.84476551742659467 60.564418035690231079 147.17446347057466483 -70.087720946155869228 -5.7602474271333248268 -24.577897349410772421 -78.387810865424327744 -104.49866134842741872 2.3527703709257381526 -17.364496067315212002;40.143456499954595529 44.522850040334844834 49.739618669549926722 56.862745951055551075 41.729039472294218172 -29.615131733933008462 -77.943326865570369932 -12.651469918588059471 -11.849962888316664333 28.319364640299031777 101.8812908916107034 2.8354120849735484633 79.212185005967711504 4.2991221035710891485 15.128892841873131658 -24.612179051556218212 13.339845584991641303 55.616631931529944666 23.618792767193156124 34.296318681526344108 -44.712180355554700384 -8.4600527117898973728 21.647922369326600744 -4.1457387529112752489 46.019689121611008886;13.528915258937610488 -5.4207416755127200148 62.498703313821927452 84.190068083298911006 -31.62318347941194574 -21.714274324240523839 13.203436359983241388 6.2366570303831529998 -8.8853899006500842717 -35.274363633854036948 71.684795391415718768 7.4281696809063655351 268.44815317197060267 -29.250064006029209196 -11.590660597837562307 56.912753404194745599 39.756825606395764794 -16.404470904310176849 -6.7403723536872357514 -79.046474550396538916 47.62755286668340915 -145.71249916730246809 -74.635815789039313017 -97.939205215772801694 113.63987257085072713;10.372693984151654334 -156.79787400157374577 -95.017345709887450766 30.353703235689074802 166.32325368030495838 -20.406749793350520861 -56.420982403161552554 -62.685229378185880478 -108.68505956689878644 -17.030457416846019214 -12.198124781580771625 -70.999044541182385615 43.537300169164140584 35.604823785373675094 49.466813908692436996 177.33376969650385035 -61.910053126319866124 48.044408871096536018 233.22095310791340239 103.34062960885820814 -33.706601028236470086 -133.08061022840760756 -91.313509356836846109 29.715661552619099695 -85.948870498593109346;-5.5521987345207932663 45.849536893490416389 21.834822421348402344 17.652660969874503394 34.944661109425233292 -21.657336071656530407 46.237189595512099061 -4.1279686384054610215 -6.0225028096399739752 -0.21513770631843870396 -43.338367778945446673 14.315594873930345798 3.6795501033715565242 -28.151222326383543049 33.174490628951275539 -59.092341111776384821 120.42398259961525753 -29.771045227943805145 -52.983169645704421669 -2.7357797861354007374 28.533102544058095873 15.314418526359370176 -23.654828398757363317 54.833564757356796804 -72.810129278321440438;-47.268304137382919805 49.329626899670323326 -135.25633789718332878 -12.95290175078207362 100.31124844523294826 140.39487032452953486 251.82141453483143323 31.384367597670440375 78.845727187944163461 30.821305474499837374 1.6441267792325302555 -51.696994020676640957 -54.880379089873230214 -0.96001593768327375678 -190.09118549309474133 -76.45994389752213749 271.84351231883505307 51.142000774754450276 98.516482055914707416 -33.965884795147651687 -103.69988064917238546 -9.7693337079005964085 161.26627548840343707 -35.137334352347522781 20.699858995026847452;-61.20457090649008336 24.308399650980376805 -55.508403551159268829 -32.805772220058557309 -192.7417616075163096 -161.35611928170038709 31.753786819539854491 2.3508110491823610566 33.087462550818273144 101.06659897829794659 106.08026466769216256 -188.81933528176119808 91.734204950063414685 47.971466781512319244 48.018049864357529088 149.47826356070979159 243.3221242078103046 86.472145694990459219 -213.11947874378066103 129.16788707213797238 -374.2040292157106478 108.82109846000417974 -116.67544297160112876 15.656507580909503119 -18.866660578286559513;-61.511268471689895421 21.656737939469934418 56.836240287168784846 -9.6008101563716330418 -52.534565095036064974 -29.813232470536902952 98.741023097419116539 22.246219015195681123 -138.52305240078086968 -5.2815750389792111008 -29.921410601935338747 56.025022906941551071 -243.01563800592586517 12.034688112194112719 55.912652915670591369 -66.839523993382599087 -38.051541907816968546 -121.00921023127469311 36.356891340518181721 113.89368457890290642 -57.08736879169094891 -9.0686760461581403803 -28.357484799195095349 24.219213424709355564 73.301036576798267674;0.090587533709025891415 0.18873467971790461317 14.343503707722586427 -0.13095642434295681755 -0.039932025941233371957 -0.021543910473790055249 -0.29703976360140499358 -0.22090638838092160889 0.17141249034333674861 15.77471215846354724 0.093358293498429445467 0.15407664725249070603 26.411480065961725927 0.047962969901471544376 -1.1833526038182009543 -0.12083857345438131503 -0.48600396073585072632 -0.083844088119661028258 1.6516367109070133612 -52.512612647592526116 -0.043488106091404105102 -0.7665466982378476013 -0.0046037591656591303133 0.61771079789252891779 0.43007799945143532216;63.262122574740530467 -13.409788220667328673 -26.673666098390597767 -6.1908013604455520706 29.260694205464606199 -0.59197715595468736183 -19.178750573668107648 -57.332518492223996986 -49.958232455780581915 74.400670356389809967 155.24863684666075869 89.815539201067366548 6.5287882270110166871 -51.402077240465594343 79.734714528331352312 -57.403548724572218021 6.6247228018677066785 -0.60085803674276527797 -79.889343146515230387 -1.4687006636653483493 122.26832126149527369 -147.01836378565454311 -73.689390380562144856 -67.88410272902630993 -245.41484612672886101;-24.926193694528954836 -29.517702254407172546 12.513814584678499386 30.448555783270393249 -30.623721643477793464 -14.472739581538176523 -17.840935582645247592 28.939232773459782067 66.536959487862262108 -0.082708483567938390069 -18.224736953804441697 -83.535712867653899139 -226.36760236718231454 70.66877579278985877 92.010104642928936869 -29.297438818810334737 -6.3224118245041189823 -22.109198560956365043 -8.817040518325342191 178.49490321637995294 42.474443834459485458 11.182658078619679642 -16.312153511443693077 39.645065771153618073 -165.95404838097081779;8.9345722501424145179 -7.3495309572022966549 -0.24483299670537458925 -10.158688469428698298 -6.3656997379171711415 19.108169352683145092 15.332315567289363401 -29.286870262170666024 -23.717030409497382948 15.179905896389852771 -52.388448116804354981 -46.495003746451367022 32.047286719090188001 3.661050001706132484 -56.927239464092316723 -3.0915703847482012279 24.867574158498090497 67.029395474821811263 -48.157723538222860782 42.006020379928990849 59.638857087345364505 -98.976047379884974475 -23.741857616735384084 8.4630927315489135054 -1.2252317771760670162;23.440236069732478086 1.3171326695686231822 1.1455776236600259121 -7.1260692121481747918 -8.9582981711885576459 -10.888243585098585697 -4.1470913801021414002 25.489362878183406735 -16.724972490149148996 -7.642506586491034426 -24.787623170561154495 6.558212612979091638 -54.11205021325205422 11.000562967366116496 -25.202281333530127228 11.466878392822275501 23.392430048935043629 19.701963274609827437 -7.9406442861304098457 21.925447982997781082 4.5856470260777522086 0.66436080132576491586 0.78645908005938680585 -7.6982761870259963644 -3.4112253657313353905;-9.0244582242589803656 65.994115890866410723 104.16862218028674647 0.078470615312385597062 25.043084632075540696 7.1928200020716319685 -44.764098856115893454 -84.566367722451971645 78.011074513856399903 -99.371634691075243495 -39.653351796719981337 0.5266464844640766918 375.24306084063704247 72.386540152268338488 75.457493589277973456 -11.445281566934204065 34.962372124526147843 141.99894798638231919 -140.02359041093640712 -77.047797117740813633 87.806303174434702896 -37.175415546977532699 28.845170584127551194 7.7874319172462609373 26.083922039483372401];

% Layer 3
b3 = 0.6201255347963537945;
LW3_2 = [0.017039129095416327697 -0.022818949974412621601 -0.012002786730813775479 0.013483666658580112374 -0.013255639250304639479 0.025974516658514870682 0.024532771615670016796 -0.0039033713378068693803 -0.021254619631270184332 -1.5884675025504548795 -0.024904386345592673113 0.021414184755000412186 0.017271762472606609629 0.04251884430306977769 0.011378580920380538394];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.0470499670650231;
y1_step1.xoffset = -0.28;

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
