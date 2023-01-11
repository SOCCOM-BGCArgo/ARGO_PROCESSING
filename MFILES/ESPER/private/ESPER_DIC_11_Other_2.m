function [Y,Xf,Af] = ESPER_DIC_11_Other_2(X,~,~)
%ESPER_DIC_11_OTHER_2 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:14.
% 
% [Y] = ESPER_DIC_11_Other_2(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-1;-0.999998793562952;-77.961;0;29.715;-133.803853032615;-0.1386];
x1_step1.gain = [1.00000000274156;1.00000063368067;0.0138606862425759;0.0002963841138115;0.248478071810163;0.00432440760873659;0.00826143816767911];
x1_step1.ymin = -1;

% Layer 1
b1 = [-1.769377401533372085;-0.85497299759066569536;-0.87549012162530270853;-1.4784892004506591245;0.65420601499672370593;-0.12792619528165108389;-0.93570280245712833533;0.8850789578207042263;0.39190534409987826114;-0.2350241518579044564;0.29860289471246542758;-0.54528168388892050977;0.20331690402734373113;0.26029607357885198304;-0.076204973051021723807;-3.2870940284323144631;-3.4581878781727177774;-3.4973938428183597438;0.46517927442830131479;3.1128017481415000134];
IW1_1 = [-0.11939514632811001893 0.1691252099493028016 -0.31176534464467470231 0.083358980710645222412 -1.2753742239946601611 0.53383763854722776632 -0.026945965687880202993;-0.28369466271510840327 0.38350306238430093009 0.34583542943543288395 -0.80900279681409537069 2.2639439093538471504 -1.0687045127401806344 -0.31560414792959495811;-0.11383056965697435636 0.0097640216096448809169 -1.189521111041507595 -0.73444902755754237766 2.0045914520367076328 -0.74035304542993674026 -0.27339864742600678005;-0.01419327419461110032 0.05198005147525932429 2.4698449844054373159 -1.2396330901654935008 1.9068899369743721461 -1.3891490103532717626 0.27103491913863314577;-0.16429153138048410709 -0.025261018251408810104 -1.0249784944067577008 0.57815198187501470883 0.0040918282393231555222 -0.32576899718559898522 0.22109853594663575471;0.21538515385093340293 -0.27195323684351518922 0.35202685238385950095 -0.022891775155237378131 -0.44112107494784802419 -0.39448573231345557266 -0.39394093461974799553;-0.22799476404148100173 0.065356022941556082828 -2.3263026639685175923 -0.86703544776025598395 0.33788389829232823303 -1.2185381627594595866 0.078922030133660087126;-0.085273737912476429979 0.17363166361435894425 -1.6728081540175687714 0.74546978867290336623 -0.20388859339625597933 0.81545102408233505198 -0.50847723647391840363;-1.0553100832165840384 0.07994696388870670134 0.82410729296828955093 -1.0543864438386023608 -0.099833274013363781929 0.048363767286794492883 3.03542950120911037;-0.15651728575673509725 0.30816832352937845174 0.25440479208255173393 -0.46713516461771775745 -0.86122453885511340843 -0.52075607763143083595 -0.36749393850412759299;0.088639123326121116442 -0.32674142789152360589 0.32738983126558041103 0.11546462968682123795 -0.19494169433945690573 0.36210085954173204525 0.65677999848219725809;0.98205844381067963766 0.054443328462594892114 -1.1705117004632286104 1.2086606137323903098 -0.13506858069956709234 0.43322577623952679415 -3.0003473780665412818;0.45697587920243343262 -0.38283219700566806232 -0.76158315310324953273 0.74985277064274691305 0.97413665431167950626 0.19134858867912807323 0.04563749456484164635;-0.087802826763277067035 0.29796111236790173704 -0.75956311669439813272 -0.065675721471746828506 -0.88839956090961880175 0.98138496932794205829 0.96746016104738619124;-0.023988150146755853748 -0.12111953341113808336 -0.53706441041211139265 0.39552386164929842494 1.5366887656462415723 -0.59970839727613112213 0.29460787081632061613;-0.074966595613551512134 -0.72310914151586902499 0.13770185699642484289 -2.0848097685495288189 2.6289951351496259058 -0.60102944059832752455 -0.54894675259845016058;0.16234957657852094881 0.31919335117190689166 0.39271535063062867899 -0.57183533224204907697 -0.10076112625015908042 -0.16604642517294507242 -0.90859078106679647124;0.11734023862668775673 0.13421036512800219032 0.059845165445635610524 -0.55733378190106785954 0.57942419710312687009 -0.45246719372668298131 -2.2188896911236017928;-0.063205533278394926411 -0.20044867457387818832 0.5074966480406140068 -0.15329807110380938573 1.7396541227269890495 -1.210216830226285456 0.11115089107149979386;1.0639013056238635357 0.036698487324782418095 -1.2626333811640497728 0.99627000525933973041 0.406046588384154028 1.1682286466184437224 0.98934337420429507848];

% Layer 2
b2 = [-38.851625271510634718;44.598595924117688583;61.349385302214621163;-15.326905788249195695;-2.0270133800245786482;25.775562958297498994;7.1223481202357943332;-2.5489821674364905313;0.45901456611666380248;12.097970599818387072;-18.87378883909731897;11.288026486150723571;4.8723052306271910084;-21.571787254747491147;30.961628682286505665;-8.7194887951006503357;-8.5893180501498918034;19.762919495604798215;31.683551023147469294;-17.100690436774740988];
LW2_1 = [-38.544816124905445065 -0.10323193941110353489 -9.2022019595319370922 -7.5290761722138936918 5.1439488365290264582 -1.4389448986791144591 2.6858027815513860759 -13.407082954034880373 3.0731707182134599954 1.2474450489917427998 -5.6252978456104871796 1.4599231566163555218 3.1775101231889353315 -2.9180336366335697384 -6.8340989000582528234 -2.3091257691335647095 14.623210702415045503 -3.7275208870510159187 2.6419282728400470539 14.123426255606206681;20.092770227099325808 10.616226539842212873 -10.267117494550173973 -7.3756251306011213487 -9.267041150454840448 11.882341338662152808 6.1233305136616724695 -5.7643287223719203993 -4.305130001200017098 -3.3554227860707435127 3.1997268118739450316 -2.3556029788450558726 -6.7476148549708740987 2.2536663946366730649 10.971732974632107727 -2.882359187457491867 28.700642010902178924 4.8945152675633005757 0.3577643918588388483 6.5204829816147267252;5.3611473486379983555 -6.8935349774338030571 -1.2733429109905263221 1.4915027022606106577 1.6168722427627328919 -12.101226365461380752 1.243973985579978514 12.254582772693469295 2.2839675887046633207 3.2980962589933224827 4.8617548908722740819 3.6536105785004773239 -1.9311939992743833994 -0.79066266910331495232 -4.793186940447171196 6.7970935169013388943 50.275954487536331783 13.255959031712707841 6.4485373344391678785 -5.6387994118911519692;4.7333213645713270168 3.665012373637747789 13.438425369955799837 8.4810206454507213181 -16.676855877611018286 -9.0329033121579129784 6.2498126833590017526 -7.9968956311616494759 -2.3816629749828570617 -6.733024152854638622 0.87202433583429217467 1.3165354685173871196 3.2986908679911461917 2.0486157756357621373 12.8040439303404181 -3.414455490154719719 -25.225057862355910743 -4.2288352993535092139 -11.206689164721025165 1.9009987541306883063;1.4945794982126743289 -0.039390281884426676939 -0.27539641155989363908 -0.30859326968524652957 -0.44594385020602156366 0.42873175428615289961 0.24835203831631119664 -0.25017942504546186733 0.019462755026620641441 0.21170146707647366502 -0.21797117441097843549 0.012864692074038093852 -0.013268191135864616792 0.30953228357817713645 0.52423567670091220094 -0.033851964145386347071 -3.0682098203082239785 0.12668306898364387059 0.21280897755253427861 0.031771180953206829001;63.133139077368028325 2.7331334618291749372 -2.610785237632384348 4.2916073769809708338 -8.6703938949786731882 4.7712854060903122644 0.062391439020872969945 4.0508538747473812336 1.579148091146941324 10.073194023003882336 8.6160839338063777149 2.3028555836413460156 2.2528899240722672914 -1.2594901541867928874 15.436269666253469524 0.45738488981037817016 -34.029142628588083141 -6.2115783113393643333 -5.807607765520953258 6.3058386776356449843;-30.187221050322627036 2.2608138349558175584 -4.0381953055687329623 0.74830796226751661315 2.8464659720636986684 0.20439339568838463257 2.7277933482775917895 1.8828079131764565712 -0.14511587406658946198 -3.4125858273941105381 0.44213629140386234706 0.43220227730613969275 -1.3193086011519965517 0.71949666449720295347 -2.359035060289435215 0.92854795937722078847 34.861804625868074936 1.5668773329273952832 -1.67499230237847474 -0.030307178014460513238;-0.58150812833632670174 -2.0882329425311842996 1.4481736183652598182 -0.81294390963083074109 -2.0782450423077256652 1.9147276733716469543 -1.2350775029902520963 0.095523852841893963639 2.7888548484430426555 4.4534744043023986038 0.19839898417945112108 3.1885618709790293224 0.13749907231690544629 3.2837807321737919608 1.5043538112270145035 0.51843819883681463434 1.0353581569479914926 -3.2044378449160850764 3.5519125756941520855 -0.41916139134099422714;-5.6325997301128225558 -0.53957328732957599016 3.4376767014036802017 1.6718590319886428297 1.6293777709449741575 -3.4398534361938581938 -2.9187298053418948562 3.8640778530588892892 1.3788340930452094302 2.3643615714149146534 2.6614178829822847661 1.5593826226139571034 1.2850565949111349973 -0.44421321319642442793 -2.5825348338888454158 -0.55190093119969052537 5.6616835793533075716 2.2544514559647237384 1.0432921540674149341 -1.6440740524704788861;8.8316664638729438508 -2.6751838660670248515 -0.93570939926846130597 0.66761120182830957503 -7.847132520737333472 7.8914535554584341881 1.6579190105966106206 6.0415532622763832293 3.2933481412591274484 -1.7596835009718125864 -8.2308759647357767619 4.4645333357064096447 -3.242101581091341167 8.8609456675060638275 -0.87112017420743748097 1.9572330867314948755 13.06514145741267896 -1.814964698257707898 10.956282458318160167 -2.3893618046251652309;-11.882133715214260761 2.402466710940148964 0.90447472600839196044 -3.7673210153266061262 -3.5186632576818053231 -3.3974751176146842369 3.1353824863014905233 -4.4380130667637676112 5.9038147870351913582 -4.6722785126437189618 -2.0999462834529345479 4.5122631883837636124 -0.41512261298186103886 -4.4247087624553618213 4.4396916325924049573 -1.1646752837613372211 -9.2069578132496374678 -1.5223691777905137279 -4.9709924446946516241 3.6687656806808082166;-6.8977941783515985819 -0.49113170947133638622 1.7707254144275119678 1.4250941057571080428 0.98971731314251809586 -1.2387595710014753791 -0.9903940610949811818 2.2431409940614921439 3.2000083340297140744 2.9825164595202249274 2.9568628464493227348 3.8475466274837435421 0.98368515288603763924 0.85778224280962878101 -1.8432598257013013132 0.60822990081147287889 17.758479312278815598 1.108857634561006833 1.5268983041941290857 -0.53799924415830702173;-28.837814546658357528 1.5564239162745561185 -9.0670019923617779511 -0.47307840793825661763 4.4157019068226164293 -2.7308291779018984968 -3.2480947076363086978 -3.621489912752672069 1.6517624491924067787 1.2642935646554540607 -5.1190951433403446558 1.572340268117024209 7.1721357969964225632 0.93091335728060542465 -5.0800787876047115432 1.1797163651044206834 22.124410006647913463 -1.8515474826505955441 -0.91608089717704033195 -1.0239379516701381245;-23.581823750615477309 0.47177250640226758849 1.2546940287132206659 0.026539788457460876392 13.653673072968178559 -1.9513939849470625099 -0.55008891191685582722 -6.9825998671141586982 0.40431587349114789243 1.9894108839835922264 2.0410063895966650271 -0.3699855740168493945 7.5925186423620258935 -1.106062702497789596 -12.930684051514202082 2.1480214705184117463 -8.8137943889106225015 1.8428818745910775512 -9.5647053778347110864 1.347059579161669074;-8.1868315490042800064 -1.9882312226641778441 1.2802389066705857079 -0.67171171100926740927 -0.70781881283792513759 -4.948412239430406423 -1.5535825964713130265 -2.4304939172221367016 -2.1144038273144012408 1.3429354561740909357 -1.0704418644184001685 -2.3808341587423735852 0.54983129592921697881 -1.8762331368428180056 1.6552787698353090828 0.21074447452614056275 37.95256867482261498 1.5518821661809696355 0.39680098413388764644 -0.2287880342116984933;-7.1148760891096918257 -2.9380900054680552103 0.014285793125700357864 -2.5811942847690416691 -6.5612866301501107102 0.26245607376698987245 2.3659867001400614051 1.3198449212798168961 1.7825923201911204607 3.4123048176949675181 -1.9537318749883612057 1.9187087084481950772 -0.91079909773354972291 1.869037331342159769 4.7651913644934591474 0.45577509466784277281 1.6271833186740776167 -1.8045705759146672964 2.9629703914456504421 -2.3221554177671452912;-16.587060284154571121 -0.47996117900735807416 -0.10811707575407984716 0.68256015164326255817 0.98999255045253897745 0.0062090229569600561416 -0.36150932490659876128 1.1205246943379354452 -0.86480784651355302461 -1.0122089252494350653 -1.1240803205018128352 -0.86017351923570117833 -1.2943105561622214594 -0.48423966358145531519 -0.42995886442839714725 0.7422338154837520019 8.1037192151704307008 -1.5246866363860083204 -0.98519774715255303832 -0.17582616238899406724;16.017527906459147147 -2.847666613404812086 10.65169098363718625 7.9752823962947676506 3.5127589620926569225 0.76562646760090502518 -4.3721428606797294592 17.024615510494097492 11.960727533070290463 -0.84791655577693503254 7.7860214381482544255 10.888689035994925547 -8.1982524120208406515 0.98570118255149785647 4.3527995928866731745 -2.0692438162911419575 0.23840264392759766277 -5.1584946537525260268 -1.3343687725608788242 -11.864840983260791063;53.379003496972899256 9.6889805608933734504 2.1866010461817602106 1.6843849149064213933 18.887667805301386892 13.334673498683526915 9.1481547416010648988 -15.556718645468297524 -6.951924768788424025 -24.657241652697830858 -19.355096329147031042 -9.8632112853558009391 5.9236882839071709483 15.367258211280484304 -14.056142607872789796 22.651291671381354575 -16.470314623205013049 -14.797598236484928336 -39.731888131372841144 31.420033214444440972;11.554973522804766972 1.0958640925247271625 0.30436061020726262738 -1.3367985148507708359 -3.1133216664258545059 6.3805618365466871822 -0.92052717938140893317 -3.3896705997074403349 -4.2721743538953118602 -5.1346785615552610338 -3.094781847152969334 -5.1171010352864634285 -3.1538861644619009539 0.56745792724264343398 4.7842047744375788554 -3.4560438554290411339 -31.645044285767063741 5.4844130473896566258 -1.6382243137545615674 1.6501846744449690796];

% Layer 3
b3 = -1.9131972059285053955;
LW3_2 = [0.0078630115618604588118 0.034331301166740768593 0.016388348782599863884 0.0045823298565639642918 -3.7627001427390647237 0.014149515097521500953 0.12880614513663268816 0.19176620380929182419 -0.13311397630830196226 -0.358259113728156664 0.049147645012481137872 -0.28203103386556566123 0.18541045230584818349 -0.039624824371713683602 -0.1212160158074818278 0.28233494435141681667 -0.43250298882103244136 -0.037706701840953163463 -0.00304528804987876274 -0.11106873713668673209];

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
