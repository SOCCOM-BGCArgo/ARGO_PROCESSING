function [Y,Xf,Af] = ESPER_DIC_16_Other_4(X,~,~)
%ESPER_DIC_16_OTHER_4 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:16.
% 
% [Y] = ESPER_DIC_16_Other_4(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-1;-0.999998793562952;-77.961;0;29.715];
x1_step1.gain = [1.00000000274156;1.00000063368067;0.0138606862425759;0.0002963841138115;0.248478071810163];
x1_step1.ymin = -1;

% Layer 1
b1 = [1.7646915599501020733;-1.2939312995515730442;13.09799548596159191;-1.8482261629837386518;2.1827138833520596251;-25.418439777465000162;-2.4662851813736885731;-1.6363960676919055981;3.4528988126881126064;1.0345270261409751367;0.25977859950432830027;2.5147965180916300731;-0.65582725182481604875;-1.4641613365291168414;1.2714968999482532119;1.0542507840381940554;-0.38274432790049817132;-1.8989880397309757409;1.2677547214847320056;10.617662586260280477;-0.12535262298180313567;-1.4572712534243550397;-0.1164830803807599402;1.9574365795502040033;-0.26868531923817207119;-0.33686254662596260001;-1.0017560048925731042;8.9305000951456303682;0.68152979749096309625;-21.505820676017513904];
IW1_1 = [-0.4562556715214899361 -0.12059812105833697671 2.7441712955962391796 -0.018286032092939424193 0.61161328709992068831;0.45150343286043953528 -0.064026037502600358087 -1.1672777153811220163 0.88928768473033970032 1.5129134081566273373;-1.2919642839659597477 0.95613410211505212022 3.3781434578485756148 20.113587104571376329 9.9387622424096200291;0.62429979131361890143 -0.16454782619587074444 -2.8537623242794674105 -0.86493849714562498487 0.22918719892211450451;-0.15462078085909897185 0.13896574870632993903 -3.6441913855614647133 0.39392371401080283899 -1.2961063198145263975;1.9146911351437239102 -0.38794558026662573358 -10.945356673736991837 -31.667095721635089234 -10.392873929799234745;-0.58763697986900165215 0.77412444747613673446 1.4343777073655010756 -1.6789959195352179844 -0.71043351148469746636;-1.0604193615641677173 0.62030013336946521196 1.6284409520841922614 1.3162685957765751787 5.9170533979893757959;0.41862672873809547713 -0.16687579039937300296 2.9062264797528825078 0.27450565142693822418 1.6230919689366170022;-0.69002832435639316344 0.45460024969056689326 -4.6945391298135099589 -1.1970099010701096276 0.49005314669094890423;0.36599867996392038139 -1.5096581586727380042 -3.4409208041453074323 -1.0783365954775647744 0.60337320596859289612;-0.62035425602344029805 -0.298662953203657211 -1.5148947400055075008 3.2039714133422494413 6.3923650775066134244;-2.195729504837809376 -0.83666126456292000402 2.963687657531335784 -0.2720071793472447208 1.3155148061776875679;-0.55421030267546733672 -0.057308268702639081349 0.041420217414872652351 -1.4685232481831047568 1.0054514017457887309;-0.84087561648319530949 0.64291664337068887836 0.42230164229358052586 -0.92341437154631877959 -4.3817259323201032117;-0.060811484314660352313 0.22653179026834324716 -2.2540787952863117205 2.4263678961693138092 3.6620965230331976592;-0.086994215215722855983 -0.18189557231229847023 0.74376620776053459227 0.58637082131052586931 0.73241249720151624381;-0.54064009757910513887 -0.56542377869257398615 -0.49884209896661985129 -1.7029656039043659099 0.70231553820299452706;2.7723268891739594189 -3.2771641611395034843 3.6747889418170585074 1.1775602615975830911 -4.1806689744414642362;1.7262050325076332413 0.31518666773465020547 3.9055371102491815982 9.8034994297126676344 1.6799586466967146858;-1.8013958753800940471 -1.2527009791490386625 2.824361325348434093 -0.10896170850785295847 1.748187391316784467;-0.40742546587069322639 -0.093136284477611797072 -0.30739235054725655427 -1.4981906888349925122 0.31625999675376631481;-0.42856250572688897593 0.034142957043457079769 3.2461458672937295589 -0.2383958936215265223 0.97804951243769511215;0.15119386627149394853 0.51963696182776697441 0.25655537078246909433 -1.8820650466569819326 -7.1652919268361339178;0.54110985305879755192 1.487996968036165768 -1.0646283643375975547 -1.4794251260138051585 1.3644201063644605387;-0.12144259739053578362 0.97571898018989999812 0.036807196287610546848 -1.1589363808512960219 -0.30248515366928641068;5.8519857428009984801 -1.0934060223280999136 5.0677505134269020814 2.8793821019588179411 -3.469061444007747852;1.5796365300686550093 -0.86750872452484617536 -1.8293111065615283994 4.6332006792612050461 -9.4562044369023272594;0.56033654812600663941 -1.350469628183012416 -2.8009478476476181719 -0.92469092034652200507 -0.78777490561356777743;0.30605676570723516772 -0.15029538761643457256 0.85671817062795729836 -21.017458618380302937 1.4239747073416739731];

% Layer 2
b2 = [-1.8161078947239039927;-0.77280858535003715826;-9.1549495227619015481;-9.4343657609687490151;-2.102898798989231377;5.6799128526207542933;0.25934284449931005456;-6.2830123429870550211;3.1723142176630250155;11.386550302017599989];
LW2_1 = [3.8913609391182117925 6.296456618324556942 1.0374471750214315779 0.12586691511263214949 -12.911725902069790095 3.4811571050484242384 -2.0970649821291034165 1.298836189610653058 9.7480595147860835681 13.649158909162578013 -13.558835795338154639 -0.80885980024703008073 0.21497163307378694408 1.7615194889635872322 -2.0006384100050760289 1.4354089812982386842 1.0852549086810596712 -1.969530778043547814 1.8672173732968282422 0.076044729098240476661 -3.7261725701831474034 6.0773665165342620753 1.8573625188602957969 2.5428154421585609235 -1.3027232596960558464 -0.88045150253456061495 -0.45256207875229598203 -3.5610192097989243898 12.367925756837136575 0.55638638749355906477;-4.5142918294295961701 -6.1894571750213600581 -0.84440854472209625481 -1.6124132242192437126 0.061287082693144520618 -0.25099115089401380363 1.1600541189888187166 3.9157519721692186643 -3.5618729185172690599 0.75421762027417182939 -0.00087563308537308426338 -0.54089488140747421507 0.91242963541296528973 8.6801890985696399383 -1.1624726896781514274 1.4701862579982203361 1.1305742636154141323 -1.0388801392157103187 -0.1372177376662584436 1.3880736132420385953 1.2933446413022784416 -8.0707106644965964648 -3.1456427705087897095 0.6763852237510208365 -0.39084146216845561117 -2.3865880670112513684 -0.42881927552523357861 2.8041165689542491002 1.9913333770162247038 -0.20398835930771727631;-6.027178886023682125 -3.7856450723227652588 0.99641978136455033876 -8.2166213314312468441 9.3583825446656998537 -0.64732670679347614762 -4.5191428332370229271 -1.5432450395599472959 3.3712451743608387034 -12.289709252001269846 2.4128171061598839842 1.1650331675676697252 -1.2375292522112879556 2.1860126732359534962 -1.8281841340907758475 -3.7722751516258163917 -11.494410835817737038 3.2743299580494009682 -0.00079027401183806057772 0.58053519696034294739 -0.54892676728801748087 -3.6048478347409012024 2.4657750381496303937 -1.4192806124602372897 0.10032319777813228368 -3.6036338766185487259 0.76805767519631040763 -0.68189341594489993881 -6.9195441186258062061 -1.0184915196802193371;0.71305835333338718485 -3.9677809078448946778 -1.8037758553599081424 -1.518094615824697069 8.224569979754100757 -2.3175798928584612213 -8.5643576708257604935 -1.7748268996269991327 10.524844288314907814 -23.047775383236587032 3.8730713970568917404 -0.86354802741560299228 -2.8783057079004237266 0.59621606104522084646 -2.1112599449077604774 -4.9932893846761565726 -5.6150079511932720067 1.8233014228767188136 0.41893410597156277575 0.75983456823166106719 3.1526788077889049688 -4.4501088861328073065 -1.9578111679139127421 -2.4180507763915710662 3.3252018936269029936 0.44407753316099340957 -3.2865908218975432575 -0.40613306903007406756 -6.6089937202431894292 -3.656371117668786308;-7.2267850165617355529 -2.4777396163154405961 -0.89793585921473162426 0.16029173860615900105 -1.4380601101314147705 -0.74582301109532567196 -6.066576517919383349 -0.89947039257453798022 -3.5389202312208203338 -0.64769128915563478621 3.3961028818023231146 0.41012696708019941161 -3.978765993236462073 1.9532048755665356587 -1.4888626022699591456 -2.0265096531070230235 -9.3754744479370693 -2.5883592707585840742 -0.75988041432425512678 2.5975004395370615562 5.4672425693634245647 -0.19546252519777140844 2.2519001424284796897 -4.3253504149094865738 1.0830240345653598943 -0.74885459058811887356 0.11369410333079292552 0.68449223747546472651 -4.042837104713828289 -3.2612332031409185795;5.7334181990004493557 1.414055383887591999 -0.34222480359286988616 1.4896148444362642405 -20.000996562123823708 1.4551176190481414086 2.1421193876978499304 2.0121483723396389465 -10.984044716311128553 20.414174638533800987 -3.3637987584739947877 -2.9264461130379011422 6.4501660878498281448 2.9024714451052706643 1.6713447404843959632 2.9582242412910768614 -21.398167177930446314 14.199642809803286525 1.1452483557397012959 -0.37397348628525828307 -4.7548392012905447146 -18.007414864787481434 -4.4918940338250061828 -7.6377258445706743828 -3.9080833336389639854 10.569943882225798504 3.1626739603436573312 1.9852246209373565211 1.7773980872860946079 0.81091949987307565539;2.2187888771535786425 1.624821141750324216 0.077995784941568710447 0.8621683270308577729 -0.27351504106884050316 -0.11331273588339703384 -0.48916890985111460832 0.21382529991886528586 -0.99489739392974829713 -0.002583998060524841326 0.26308469917065591126 0.44557732135515754068 0.051760069639356441618 -5.0992527860783170013 0.14791965464904191507 -0.42438742693894943026 -2.5724637059317565502 -1.409273743565739867 -0.019164835792945339454 -0.058893267204966237571 0.093953696622156190887 7.3668124917683144304 1.1488775038313536747 -0.51034846556047441801 -0.3682951893215493655 -0.15454404956518350001 -0.011968118260144142151 0.1023737586878128436 -0.49348892027428348239 0.16718402671727131126;0.97050170255023371357 -0.94974407128886917384 0.65071722819457600284 3.5801430377046008502 11.866009979913714645 1.3369689265106305953 3.6639845005558973057 4.9237887262296986535 8.7880268850778193723 -5.8896702019194213662 5.0009645551368278049 -3.7052678569199541947 -3.8609142219670555818 2.372794945052679072 0.90462329926861573259 -5.6970460644163143371 0.70167630517142820512 5.8931035956223398742 -0.3677554151138026084 -0.58042872312012694369 2.8890088969362355265 -9.3608863143117382322 -1.2346584012027952415 -0.666266110761133179 -1.0818199446629439553 4.8302695145418583778 1.4476697919192804687 2.464527926824683135 -1.9632481793016338667 0.087254112100906236482;9.6680760606206650465 7.2929021116316103956 1.1343817049889657333 7.3844821186750895237 -7.6593425772820236119 -0.87182568323445885916 -3.5338769359388968461 -0.89050732110137820197 -1.0615776437007369282 5.0056370897581716761 -3.9460296087236166329 0.95164533564006281718 0.10971097394473708353 -4.5027588728596334988 5.2113890818491466206 1.5469104215130220492 4.916716421418587224 -7.8017643342582960742 1.5674187075048380269 -0.39534984137580536645 0.74060871840636532948 14.757248268917077638 2.507415216381656986 -0.98690942802749637153 2.4431692891596306083 -2.9465208924048025274 -0.66893300191583704439 -1.4964217786898470486 3.1397306131448265987 1.5205372687660001496;2.2969170626731751916 3.2424690292008757631 -0.1575729346206900805 1.4355975088565084974 -1.0557618364786238718 0.11244489391782161447 3.2661992597988125375 3.8950072341331956061 -6.7842115975516863458 -0.13078650065543409586 -2.0542978220494396702 -2.6178683144198835464 0.82295914084337928873 -0.38223553276394578448 1.7259059666509994901 3.1070925961923805403 1.5943545882477196951 -2.3541923037667911167 -0.43250053614425509885 1.0376550576764351774 0.44089655995065962779 2.4201776147701079545 1.3602836281844203903 1.3788158461363813068 1.8742414226177717129 -3.6539871133587511487 -0.74951540662374804835 1.4620177013060695614 7.1994300804183124498 2.1389921868940451333];

% Layer 3
b3 = 0.35857563936581243746;
LW3_2 = [-0.067238308912822883356 0.15373539475609840355 0.10045663355434068797 -0.054887179768284777293 0.13811118130685062777 -0.11955155686273909177 -0.506247193829821307 0.071728218154073403179 0.096661357974590134123 -0.24258040652241377733];

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
