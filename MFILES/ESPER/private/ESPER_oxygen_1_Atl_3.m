function [Y,Xf,Af] = ESPER_oxygen_1_Atl_3(X,~,~)
%ESPER_OXYGEN_1_ATL_3 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:44.
% 
% [Y] = ESPER_oxygen_1_Atl_3(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timesteps
%   Each X{1,ts} = 9xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 1xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-0.999999996490807;-0.999999994516886;-40;0;9.1619;-1.99149636107278;-0.03;-0.28;-0.52];
x1_step1.gain = [1.0000000017546;1.00003626070655;0.0153846153846154;0.000309885342423303;0.0656325902039552;0.0604107502911417;0.604229607250755;0.0470499670650231;0.0145232735458572];
x1_step1.ymin = -1;

% Layer 1
b1 = [-0.17107895174416085093;2.5175131498405090014;3.568560111301072979;0.32934786267227000867;-1.74325928246735673;0.51050247059573461428;2.6959208446368485568;0.12531642541644194555;-1.5901551691093043228;-1.44813392830202603;-0.18640974853219569041;-0.11817276434500602333;0.95903368040395897687;-0.56275654505722927379;5.4676550382259181404;0.82576745024528541705;-0.38304763131141283949;-1.335765759107687467;4.6011253148788560097;2.4115560177910833062;-3.5265275867928607134;-1.5823897865730416346;-5.9726985281457833921;-2.3455210960413186427;-2.6405617729831334373];
IW1_1 = [-0.1998328566308275045 0.047821231081920413664 0.45876456843539997754 -0.12830147024086752117 -0.37706714121599377298 0.79773243172132635159 -0.10031592787842828762 -0.074772749859136469763 -0.25557344335008241876;-0.12805510362424882276 -0.16912745665778697313 0.55246007496613191989 0.60584287634841882397 0.36106862714733228836 -0.4586046130715781266 -2.9364803678883291127 1.8919224980127176128 0.64265883197013529138;-0.0027733480027028463259 0.37758540716848959518 -0.19971236969597108701 1.4150262187365918454 -1.4068645391138836498 1.0568816557576901083 0.29340917011698242067 0.55096071083173392413 -0.19780117800794339011;-0.034136221462465324294 0.0058461898669341710519 -0.17370832058297203115 -0.32930620911253383198 0.091102584139869061119 0.74661479106351202883 -0.38185505268799035949 0.085007709195047131967 0.28369723531394352989;0.027889549916275981761 -0.046045668434335308705 -0.3579034603715013918 0.2386603814201891105 1.8052983966188618226 -0.80898566662763016222 -0.024216007225377398498 -0.049466381193949686412 0.18242665599988891478;-0.36049475197125963044 0.13592444658060023888 0.024234528454082689825 0.51970741336350079198 -1.4070379857512975708 -0.26641533428492775171 0.3594429804470348766 0.3681473401571992432 -0.62668853189902329337;0.30534889361681105679 0.34818172455371243501 -0.92198227211894268063 0.43534215282147475934 -0.13713304015022362292 0.50265160290620680694 0.10931752963414673174 -0.10819179370081127844 0.2245394521306277591;-0.088580617927734472961 0.33244395645737506717 -0.36524791024386016325 1.5001815736706853865 1.2862234283083038289 -0.0037138022386168541289 0.28643101578131130447 1.1235253326949870001 -1.4683966077840677755;-0.14627703899869234916 -0.21413540223268065499 0.61223469648097694762 -0.045908050567453867474 0.063441982314088984918 -1.1070583363212505024 0.86394894243916542376 -1.7541525589594337831 0.32070761278474585465;0.058343450906898984554 -0.12426627841409242092 0.50200793419206557466 -0.060169955665418217716 3.527455856791288813 0.024434026049328164132 0.013721861570438845593 -0.86828479911708200678 0.34769655416916950275;0.4040660345097434325 0.91594760843296418873 0.61564777817163041096 0.11094471326808116873 -1.2736900179953472012 1.3149690240928244567 -0.67826989562957762381 -0.38301254543976709366 -1.6142178756448670729;-0.24241881630223233168 0.28140608263662497679 -0.015130858210304888589 -0.064337061887577109975 0.39424863677119115213 1.0704619173548053368 -0.61857058267671127805 0.26132013517603874808 -0.43794085944428101653;-0.31051766444985900995 -0.0023184853355330149044 -0.20974408956605072873 0.60509769881150510251 -0.089814105569436475385 -0.91044323273189464718 0.20594360222290045148 0.2558681754909436501 -0.34134021331670288957;0.042772276408928985858 0.042536202208586032514 -0.16477619681850308453 0.055281289000718283888 0.39852802338552278938 -0.004377111976378987579 0.13407060128904407792 0.20493306555572135541 -0.11245966736264209618;-0.5898794834356464456 0.69652127995798185367 -1.0258260516405348639 -0.1054260337997612873 -2.7557385116313390938 1.3879816818689370006 0.70497515587857662034 0.79728921833712274037 1.4845571034086326723;-0.23835142420388436002 0.092601539949328279255 0.21483386568812459783 -0.17488659553659480217 -1.069247267360893261 0.73955046147666425593 0.18219126399966795171 -0.084817020851988417229 0.023094024855328056101;0.20360419672388452295 -0.00013515965062415864588 0.035261907751093723995 -0.3184684438674262319 0.44898856316208785033 0.47492703627601018379 -0.14105466588714074549 -0.077856048541158850984 0.16881269590250636004;-0.29341502852475820351 0.50112908121602539602 0.56993148137804494624 -0.35863678954542710375 2.2610109430104810002 1.4171643720793651866 0.36287887065997515013 -1.4326849739697529174 1.4436114274767708476;0.60889260960615920748 -0.027076612456581085342 -0.76388179716355264492 -0.38517829089009236743 -0.00052929020876042189546 -0.14898591909315955784 -1.285249629970765195 1.4551948987410061598 3.4154315277463482836;-0.36303667537456807413 -0.050867354513651437398 0.285020202899038122 -0.15055608273396445251 -2.5124575188277327165 1.2463270928514564861 0.10101783830701525912 -0.92493522581349374434 1.0186280556105451467;0.088334918107230084305 0.58855206839330376223 -0.18436568689723817349 0.66961287891108722814 5.1938709612804343507 -0.39425138315426144775 -0.17275583145871667434 0.93819075877521407314 -0.46940273686231731576;0.18482807652668525877 0.016155849019720974757 -0.30613229423305754473 0.61805903173223919467 1.9499505848599321123 -0.90828230565255752005 -0.15506194894065405965 0.94397060613796746065 -1.1896543278684741374;0.60304929787754313342 0.58141236920518402087 1.2897497082683204717 -3.2953311471445401182 2.4617734601948715323 -0.78672802589576396048 -0.87953493894445522461 -3.3380459465631453497 0.34484051160230511091;-1.4039203420309285253 -0.52245598324858222306 -1.2573604894199301718 -0.49393494003287186089 -1.7352798066507224561 -1.47980376991429452 0.1744701081237526874 -1.8922918449686383369 1.0117444205404593216;-0.22963389332167663537 0.055715728561385610695 -0.42169215012072669424 0.17455269575182569275 3.3340479483378744341 -0.023718586542060597633 0.65046249876094630604 0.050668119723287968759 -0.0012604040604536591308];

% Layer 2
b2 = [3.5549723769461798462;0.17533835398584093901;-2.9043625284461600522;9.3277899339631620279;1.0962477095314819842;-4.8732496447656235006;21.623647222311532801;-12.687881064661548081;2.532229807997292248;6.6063226017482667984;15.752468833049412922;11.365198218830615318;-28.667642512496247065;6.8322786903343670772;13.464257310699821346];
LW2_1 = [-0.33164212288524208549 4.0722144262189274855 -0.57312470396762404956 -3.7069147674111442292 2.1637421133228142978 1.8750217256547383471 -1.1280841266667323008 -3.6536849521548213637 -2.9439840910578127797 0.53443833553282116267 0.85663013022696432586 -1.4652612357766516826 7.9278327763335356693 12.457445369229613874 -2.9427004571905661656 9.3245266238014217208 7.451493473990844052 -0.45398846885713967092 2.3370070542713974326 6.3664137724071467872 -0.41501381182217539356 7.956166836896585437 -0.04762572001061721183 4.8967923007474345454 -0.14333709759049420573;-17.656465491911063737 -0.18178920735435866174 -0.26269190290482663519 -1.056823617613803723 4.5040451608296194763 -4.5704785595673937948 -0.4012102312042701846 -2.4192290878164928181 -1.678571286989889888 0.1315479488617979642 -1.6217626345862350412 1.5903013246252788448 -7.1283746798082043838 -0.39105491318029161896 -3.2209790754161646831 24.31504294820314982 -5.6082122178211788466 0.71382066480360351068 2.041195749512315416 6.1304832122557462171 2.9512375534735770799 11.887785074702376775 0.17175976991990937748 3.0781668174460534715 -3.8243454595660293371;9.9118395366083724696 -4.0041507041312387472 4.9673149182093689902 14.743857656905914766 -1.7179227773014440483 8.8967551667478286248 -4.4778420909511655879 -6.7945812103719687158 1.0214608949548771299 3.8351682722940236836 1.0639293904635260812 -5.9517935882767414313 14.025843512997127505 4.14503111603901786 -2.4580715661108683001 -6.912992122689680663 8.6697761841084943768 -0.82767955610797461219 7.091096648827523552 5.0379103170850703108 -0.30373409242352206183 11.638767328073425844 1.6753659037331534165 0.36134995008594200661 4.9601535707170736345;40.278786530136819977 19.180592430022343819 -5.073743339178560241 1.6545224901944162177 -4.5375481114817288741 -12.750236454824776189 -6.6346686032412263145 4.3345204356085904607 12.613350663608178692 -3.8610033677845234656 -8.5360769087215402351 -11.241564154906356521 -3.8685917979803057243 31.318578740558137952 -11.417263121597246212 19.446552401224419526 1.1128187575646957264 -31.702596972543940268 -22.129209669714342112 -27.384565740126124211 3.5206392256843184185 -27.655762917027416847 7.8346010063967348458 -25.83193389558848807 8.9202159603500010832;-1.4431068238193094366 5.3930390312575937628 7.6908043198447257893 2.6259195438450242399 -5.0375182890296059668 0.71511058412153338804 -5.1400977831453502631 -5.5384298684962809389 0.19269327904431313159 -2.2095252657877337832 -0.35700164194465416889 -7.1728434818613511226 -12.611686383975596826 -8.4612204185257784417 -0.65259731916808605412 -1.443186018135836024 4.6267129953780781193 -5.7376053938658584386 -5.443192760694248733 -10.557038288360809375 -0.6887442997572340575 -7.0398635988998741198 -2.7818570454842177497 1.5724972418753970782 3.7192441851226800686;12.982789446614631856 -0.42127036569606474137 -7.9390867215795974943 -7.5658632733982837948 -3.4680088188392166337 -2.3519214038171791792 6.7533300966741069971 4.700440825012596413 -2.6028430530391579723 -4.562468291484115035 -0.22228728606791733902 -0.15679748718012387343 12.300264737806365645 8.7133801846688445636 -0.90815364293945788354 -10.494507081115472147 2.1393574179802881119 4.9153544870326655314 1.5882968472049812103 -2.8806192965886796564 -2.5429179042393919019 -9.3921325006675733249 -1.3681725354983118592 0.23616019113207634139 -2.1709322888185926814;-7.7609946649103980221 -12.590531102825021748 4.1480138652908129515 2.2387346440748290277 5.5421573467520568457 11.582616306139357221 -16.803069998460177459 -5.2546587807514129054 -6.5482959082705614051 -3.6216281983785085785 3.2744072548771057107 -3.8838265451488704905 -7.9384469277309248625 -11.571782453646560995 6.2961069457974820907 -23.30756640732498397 9.6297345804710907657 -2.5843225704477226223 -2.277114546069510137 19.685395562222556265 0.18230391403537854478 20.135522412271022574 -0.53676725642248368509 -6.5477060059931000424 -5.9347786882828907906;-2.3389718176287606077 5.8571766977361967577 1.7833113167375844732 2.0493973872925317359 -5.5776767754532805554 -2.460890664179122389 -0.86801380597738708733 -2.239450267046295906 -0.083845204505899675884 3.6509199925733066827 0.87440279229391926208 -1.1869522945325583496 1.7234348120475873056 9.3361370727973689299 -0.63400335468254331861 0.57283338605925993026 -5.6211514584665005501 1.2378679986014287095 2.444012414948195655 1.6879366924417504858 0.33974400495017365031 2.2298347503516922075 -0.026156960565302159183 -2.3995006627860209214 2.5455002002759883339;3.6796642159508969883 3.0880161739039224678 1.5416488972998640961 -8.1547926054038750721 -8.8414393980049705846 -9.5007416362195176163 -8.1278205605774793696 0.048071252488048540896 -3.3480341921960254936 4.2588031646790494023 -3.5654288247442087467 -3.5952905918560764675 -8.2476155066777909752 10.970221289983058455 -0.48497387110611067662 1.4203764922604982246 0.24016102844477099643 11.336292998502482732 5.0416009425728836035 0.72994339084025461339 2.6246054380211556456 3.2298475182422579088 1.5742593092305381397 -5.847971417563572416 2.8935274242149531254;12.362106487881694861 -6.5053221391246003336 -2.4574589859910398459 6.762389652228773862 0.39308333912618054207 4.2881847284354082106 2.8582633839052142832 0.032159419181134504817 2.8080807949174459992 2.6928509482766123107 -0.11513049412698089791 -2.5805024323595588953 -1.517123758449521187 19.968105194953139403 1.3265278213019831011 -8.4555306288239613366 -3.8986114144108738877 -0.9629996524679164116 2.0820019782082193061 -0.045480977982651479841 -0.22162241976468954885 -0.74172934654942690447 -0.40670614970109403297 -2.5402929562620037984 -5.6241507286973053681;6.7020337631230866293 -11.68580509722289662 17.795422032222198538 3.9609238301655360637 6.3093038997112911304 8.8046602959687643164 9.0254189973150076298 6.6028668836725321611 10.907608314616986434 5.0522003161359378964 7.0714047065601421949 -1.7793369824317364358 6.3637992761863770852 21.543670498955293624 -8.6726585738019519312 4.5071057784154806569 -26.879708309313674164 0.8462971756173834903 14.235645790684252177 10.698213988630788762 -9.2731431942321336237 4.0984271794889526319 -1.2566585297146966038 1.3671804475656246236 7.8004347467263359661;-22.404497216903170198 -12.389428290939250132 3.7664808971337930821 4.516495098747265402 -2.4145090260476225197 9.3891707489124947017 -8.6258746327303015278 -0.67857019425606912577 -1.6272125040489646697 2.8382812665607493052 -1.3449449040360053065 -3.6753257000381456265 -0.82154220697109758831 -9.5659013582381202667 -9.7587809736132218319 8.6218757268847685538 21.464397878500786732 -1.3894052545511454966 1.235528426532838342 6.0152757995702614835 -1.0116753336959465326 5.5653138946870983617 -0.21583424313361437852 5.292398310809053541 5.7911364738362918203;0.57990046867824696086 21.214214101661077905 -4.4416136754524755403 -1.8635930585553626138 5.4870987873258236789 11.257410619205423075 7.0239666147967358967 8.5919593371130797976 -4.236707979865065532 2.4796862594252546685 -0.38794701732208436429 -4.7170676355915892586 10.611752993563067449 -16.889140613398872404 -2.3100968413138684987 1.503013147357289192 22.758893567337647568 -0.5505281185442526537 -0.59836110091181760051 -7.8532218081499651774 -1.0237411053283478335 -17.943427753151830473 1.1248398143930191662 0.89227540204693078785 2.9997501876062204751;1.6365401390206659737 6.089488591316708721 -0.75127051432638225492 0.25400755768906069232 -0.80277172155207832116 1.2815155363096979446 -0.72999755698338886223 -1.127219885114193465 1.5605644936476865681 0.7860077328573867872 0.013056167481785955348 -6.8258448526926791544 -0.89410617011404547316 0.98245789718183573758 0.51448331536757452298 4.127297915971717579 -7.8476969561013962817 -2.5049927824216928762 -3.1330783275565203816 2.386457827206319493 0.6099516125268762412 3.4658678706212251441 2.7274135901124219039 -0.20882649475737341538 -1.1548874825838784108;17.151831878124976782 -7.0725812516008543795 -0.2355590172307599206 18.924986850240973268 2.1432521856308257568 5.112630992689293663 1.7497008381393752696 4.8700437556414195583 12.058561046279448803 -1.0920283747871406543 -9.2384973251914317416 -4.692616389941250965 -13.479577866577185929 0.8565830930350945005 -1.650721449606240121 -3.4334550890304504911 -10.473276248275739775 -4.3846084516081100801 -0.79787766864441733716 -4.6583800976667806282 2.4053906947889065471 -0.44523308725732996249 2.283476159278708284 -7.8140534917769608469 -4.1577609269835944517];

% Layer 3
b3 = -0.51651237655418091865;
LW3_2 = [0.07901979365713110437 -0.14091000580059109826 -0.051991716437348470259 -0.0078257061550425534302 -0.030610963375156171468 0.1471188863019716464 0.044169945191853257627 -0.3228119398818075636 -0.013328371379121045637 -0.1860467425120115037 -0.013905900107912415667 -0.32808482235678337879 0.15265479529948172699 -0.078691244268226240299 -0.020532284643649444111];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.00368391969055075;
y1_step1.xoffset = 13.5;

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
