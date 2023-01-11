function [Y,Xf,Af] = ESPER_pH_14_Atl_2(X,~,~)
%ESPER_PH_14_ATL_2 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:25.
% 
% [Y] = ESPER_pH_14_Atl_2(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999025224415;-0.999999994516886;-40;0;11.29;-0.2178];
x1_step1.gain = [1.00000048738803;1.0003169936794;0.0153846153846154;0.000310077519379845;0.0705602618858232;0.0479365703301392];
x1_step1.ymin = -1;

% Layer 1
b1 = [-1.9776544287230239938;3.1182319373516595995;-1.3205905683654577398;0.97309950159889002652;-2.2739620394316997221;-1.955759550172890382;-0.80654001158362909063;0.33559544522261924859;-0.15545678039313318153;-0.51285232213915521449;1.2636639281576789173;-0.59562053048473806438;-0.865393427153428485;1.2820832009768194926;-1.3861845040450357924;3.7980984289253063579;-2.2630636565590074305;1.8773549700433684873;-2.1241206186657608512;2.4518670456229836319];
IW1_1 = [1.4129228264427116102 1.2188373082352947918 2.3561231380262346313 -0.62897234250733446981 -1.5705308241376485512 0.20442918801063067469;0.15146352567390930499 0.023527490159169953121 -0.71241267330403867941 1.29473604961258415 -3.2087913826523344341 -0.28599090254345371775;1.3299367047795598307 1.4594205537846176579 0.74730901088036905122 0.57505105734527106609 -1.1227520698056188131 1.1695280435206536485;-1.4486130377435835204 0.90116139877639334088 0.40479444272784492753 1.5973448774877143652 1.7953045805734120233 -1.0807861056291974133;1.4695733714941487946 2.1529834000929497506 0.40434146008422888752 -0.15832712075892604875 -0.025225270335745956368 0.86039015183635780826;1.5555916037499002424 0.72468158950721195399 0.57456385516065711094 -2.0392209736606323744 -1.4801578103230748074 -1.1931661008527632006;1.1813972309039053599 0.22847965120589300292 0.79841579666202577847 0.9603026366661595592 0.62348238744128314792 0.81944314936436202235;-0.80922751394761027743 0.49370684067526343641 -2.1300040996104097069 0.32252010581087947605 -1.5311217438315773975 0.9864254361005515559;-0.82799688788135927808 -0.18903964065359460389 0.36363987037041711758 -1.0383011498025600527 -0.79990501397638447134 0.36397925485935711531;-0.43475864126205376881 -1.2631625513196922395 2.0546291159371032897 1.1019265266626891986 -0.019007451865921735218 -0.12292707516521636513;1.1107350181060169358 1.4256633770034927711 -3.5117460883317996867 0.053354900106870273291 3.426873921920957855 0.2194108911035639875;-1.2226516845196062633 0.97364028352633835706 -1.4810292567533316355 -0.98738015377963950225 -0.74395391901452179173 -1.1449015531942061674;-0.75355375012545677826 0.13020047651059807325 3.597641727432082881 -0.60246728688065609436 -2.043928061793700568 0.027048036902848082014;1.1320833915079993037 1.9718930640944183352 -0.2440417734400142924 -0.57245901619714989472 -0.1849559535699550894 0.24386195436964921979;-1.072327583002607776 -1.3010897675969670395 -1.1174880706311633993 -1.2071843722606681037 -0.7150216557200702594 0.61897228241763890821;1.2805683214239089729 0.73467718218751387749 -0.26970031333833793807 0.11224912829609015597 -2.5854595684195178507 2.106098844278991411;-1.2493182042127657372 0.40709145745553532159 -1.3510799597077591816 -1.4776728688801141764 1.0843084828521591945 1.1654738434420868209;1.7342543321312733351 0.59681212468060795828 0.50864303551954681559 1.392354173833850739 1.6619151438209813776 -0.28340511440914961661;-1.105292684940154091 1.5708630083877048733 -2.7786264393621724622 -1.772001870310310423 -1.4560809131876388811 0.054543605964259424257;0.36072118075608006205 -1.8958603908279796357 -1.1891336627125446146 1.098120364661258952 -0.0023239254698157948975 -1.5486014809676786896];

% Layer 2
b2 = [2.2910305366200227617;-1.6819752598905544083;-1.9065546964005328778;-1.3831304866609588355;-1.1954493579300791151;-0.99874635814731105299;1.1224328533689542642;-0.15597116568113156276;-0.69432147022348100762;-0.70804493029297344009;-0.0035112173817289514706;0.37368070490306959375;-1.3431824680174104802;-1.103840733965026466;0.19199505517321679982;-1.1608384107472971003;-0.97884581061690412174;-1.910559370340894958;1.4703183177996281827;0.80749983890007015486];
LW2_1 = [-0.080647915169630254861 0.61476477988686040188 0.091411243763614091762 -0.0048477065891416332849 -0.097845773550139333863 -0.2770774252799015902 -0.45623856578175270293 0.65294628486305006376 1.1388453761805623632 0.073205402646558670465 -0.60406179714542862502 2.0317574451832256521 -1.768561554316141704 -0.29085481764285631145 0.29212932212473763371 -0.36606808640279830946 -1.4570787702773519978 -0.68543544352797136643 1.4616121353312514497 0.39426221282720907979;-0.75441483318332158703 0.7000924497224008114 0.86749548159702982719 -0.94015419558097845787 0.73919626410906513758 0.89473038208504063196 0.50989529783004139407 -0.013291712417079160069 0.54500445991429069537 -0.61404335881223182092 1.6265196805105031075 -0.66296868261786157106 0.02106792945822914126 0.84641364188449319528 -0.17258934329702632637 -1.0214372538360849685 0.28619013830448813973 -0.3149196687178852927 0.42330394843919055825 0.087580932822845347641;-0.066156401144296703154 -1.0884847567205715446 0.48708211258477390748 0.20625786094225923017 -0.27971473573581251415 -0.12979102693624819387 0.22422985200923056803 0.2858086353683419456 -0.59686108889554267964 0.63107386791416708594 0.35073363268389773362 -1.1775717556359608373 -0.11684985937341089979 -0.12482219617596183869 0.89096323150078982067 -0.20842704499206782254 -1.4752383796432733121 -0.40787930675936012737 0.062017656146379294968 -0.18810458128591042715;-0.066509804641357514465 -0.49479139550740314224 0.40394121805296656635 -0.62852193540905454405 0.37957552891045948096 0.35309279905086354834 0.29388640109579250437 0.27379081665785698352 -0.19417199392042394646 -0.30536316210742397992 0.063646091547864047655 0.082998142032110061583 -1.4359820145942046032 0.47337526591695705536 -0.37942529549715797721 -0.32456614150464740964 0.27889056630995723296 -0.17754104832862263597 0.17926439762197876604 -0.077785817632882198192;0.56448023969628580154 -0.7055107955052287716 -0.23808986811675267314 0.12040944022960108073 0.055238401414128522615 -0.85656700317001943645 0.38122940957751744184 0.08788269724451255871 0.044765972126957515043 -0.66920055471678596692 -0.52257783955877579185 0.84042249680159764047 -0.39284506527098544959 0.063430621631700309648 1.3064743234857125742 0.2466768141035840245 -0.40522504724085722794 -1.5906055781409269567 -0.53799886423234721811 -0.38289221812353363417;0.6521941964752419052 -0.042066595436300896482 0.24023664182312848281 0.12409249318192884515 0.24065002306431687584 0.35923832100504204945 0.18567277676345703141 -0.47079212413229498679 0.1152970639878860698 -0.90763572196597941399 -1.2616380619993496914 1.0094596001197482238 -1.7777143579782290761 -0.35703699148030532129 -0.76131329513905698914 -0.023504126979929956809 1.2776951830988798609 -0.52269230730837379362 -1.4891650410530468118 0.45773891477478328982;-0.093439518300125715133 0.4318721210855125392 -1.4118205356871720291 0.72885941153622657573 0.86327397742795619617 -0.054308060412227097957 1.2671186512961345461 0.38163800285581511718 1.0259830491345398595 1.759226737583365674 -0.1431217490300615014 0.34267830027374396318 -1.9577561908766880627 0.4093106562714990404 -0.63291704441628193756 -0.041622814180080743018 1.2399807793094037578 1.7155200443041367375 -1.0041971478345179669 0.2751620621420781676;2.1192708281459053943 -1.0870945255774531191 -0.29283887324005142272 -0.54464915471702912697 0.51676515213619134137 0.014286537797489652635 0.20494605896262482747 -0.45640334327804049641 -0.21547622873436542146 0.27801890818756086476 -0.47976742121220206094 -1.2031100497479436662 -1.1215791978022908726 -0.18032943439680798603 -0.31274906379114297827 0.20031425952015810554 0.79866618336025463343 2.2635274567843191029 0.83580272995623960952 1.2168242306645671302;-0.24172235401057984294 -0.48293382487092750832 0.51356480413035365817 -0.10903013755787133165 -0.71421675043556143159 0.070887031300009534229 0.63921202771379181939 -0.69030126174591999177 0.59327842551409704086 0.032751058628783386029 0.49076465268665586672 -1.6080015382744268226 0.80613901706752733034 -0.72975346791199646734 -0.95739300754119405212 -1.0279614010497926113 -0.058815475418243369499 -0.79254061933872110846 -0.94467621340732099711 -0.5706832091918199179;0.38707144225975870233 0.27773686721715240555 0.16551741771204245279 -1.1139654369033498149 1.2666479453139511957 0.037542751638515961532 0.040166385646963252698 -0.39273990113098172783 1.3423188296450629942 -0.11629427659279185925 -0.66177227715106767825 0.00015454094431452403136 -0.7154605393778763478 -1.3018353290100674879 -1.5607485243057563906 1.1950234776318240382 0.20281517382918795578 -0.69213397929903985872 0.19192712763798883135 -0.1454217265277705784;-0.036135140247592983931 -1.2004460957350995276 0.71062299821821672818 0.29894491242478471449 -0.76054069902189780361 -0.41744259835186792174 0.14411178659227993903 0.90867745281867240426 -1.5827083553334082566 0.79428732372938259942 0.40154357189360800895 -0.54068285451447894996 -0.51816291220481625057 0.79878267835119265428 0.88472752844680924156 0.43404306678842363709 0.078684804342440281544 0.35604601922619300014 0.2398000051874693006 -0.044532550309530669064;0.74541540791807747723 0.72566442628783922597 -0.38229851903832651194 0.57591597606529187914 -0.20772119271764785675 -0.02659733982725816237 0.10529320422715729433 0.64205719783182557681 -0.10456498793520019319 0.22030295220048051741 1.0752646760059434428 0.13491853205826975781 1.1988435174987992493 0.67583083605402172189 0.34438728721281264278 -0.65219476264974340918 0.81196141616267414332 1.0314128244269789025 -0.92279222145375006114 0.58717263551579912395;-1.4957145269380289765 1.2932962532863394589 1.5727232408008777753 0.93439944277751241231 0.68083259814989449676 1.2302798530310179448 -0.4363246069894722301 -0.12912679943611424238 -0.089255683315494965813 -0.6302065377406427249 -2.1984737646221823582 -1.7008895836263420165 -0.63270453309675711573 0.57306514616854220545 -1.2422950927299929447 0.0057279353695316934103 -0.50685228919559544014 -0.6785901777177202332 -0.086563533027273958931 -0.33949681831688466316;-0.63834543780013397463 -1.0069604291230618287 -0.012761728638741972469 -0.011822736026990066788 0.063616227902564548291 0.28803997171094514185 0.42833753166707916327 -0.75133054630499118698 -1.4063157625957412034 -0.2566703058696890305 0.8983216476694642072 -0.27374214524317924679 1.069314130203811386 -0.18972263940876055166 0.25200083625001185039 0.27058613728472696858 0.83820342767348388069 -0.099616883629507976816 -0.77509029447351873721 0.043414040131177111803;0.011186458524353187424 -1.0887597812697724642 -0.27685002671749198155 0.25689785251409746891 0.94238624975858842348 1.1372167103243400188 -0.15904250452833204599 0.0031075527111778036837 -0.24034411774346442492 -0.030949114623117876532 0.69602549332603913523 1.4716083695355770544 -0.55206490314768363881 0.16301997255279054855 0.38072734744859704215 1.2897892300100024521 -0.85680369732093797364 0.27668157705729823359 0.099236445003855838776 1.1300405260870167723;-1.5508667544648644387 0.43275029189348324099 0.12307278454180659022 -0.53879181083580529599 0.071876932560762676894 -0.20349874202532561296 0.13947590980613347456 -0.78858066247435876228 -0.073919271978417702695 -1.2855515676433044714 -0.12627116521426726137 0.40997167578084608985 0.19762770534870569006 -2.3794791541142239488 0.86675955123531256419 1.2197828525126244159 -0.47065061113242384616 -1.8102024072451277092 -0.17565877840636681673 2.120869712614729341;-0.018285144577207358002 -0.08489360499073966837 -0.084036316453605902632 -1.2373861599777224995 -0.27109344958426073724 -0.14034926189006310948 -0.6346598198508192068 -0.90163787554540320013 -0.12658395493213978056 0.74186426770129065833 1.2351983007597178155 -0.89301991897232391171 1.6172110042028715604 -0.48974973568589486028 0.11547948171565576569 -1.2163002740618034103 1.0314172795105016611 -0.55268863960803626956 0.24380216500530543788 -1.1194991791598052355;0.33752863974289265547 -0.5541698075064773743 -0.49071826615657709869 -0.28091464819884121029 -0.19843457144891607502 -0.60793224046724980703 -0.53714452216799513096 0.66907458854063750486 1.9493683586611805225 1.353173710951863562 0.20257867700864112814 -1.3954682590382456731 -1.0533384256300881709 0.4696861555935140653 -1.6276238460347416392 -2.1297875792497706904 0.43160840062807614981 1.7024914301225460012 0.32010416138227576477 -0.58803128859578612797;1.1083783246198404271 -1.7472765795839899639 0.10707952463105092911 -0.026716966801087399519 0.077733154945289179016 0.13210267353815355329 0.69710681460214063421 -0.53309539711732045753 -0.93687487432924088093 0.42322723005878110092 -1.0049360611139344091 -1.4627874393596764779 -1.7897933516850952351 -0.52585855557186400588 -0.2246795259239938447 -0.11466402192651223291 -0.9113121624782557495 -2.2910830152783878688 -0.97697030726386002719 -1.1723565983026174031;1.4117243690223160169 0.32107402955778630016 0.18639832422243618471 -0.61620219979194335025 -0.46478158315472428708 -0.27480148713207008804 0.034508079334188423992 0.72488150816110374919 1.1331339980329417916 0.95588410893745512364 0.25990286629502640237 -0.66287639653196339662 -0.69531057883155877875 0.37839259082350679142 -0.0445253475108519578 0.46226582109473990378 -0.85257325627589952699 0.212791166070616794 0.72636885305449050421 -0.63210245571280765375];

% Layer 3
b3 = -0.47471661736674625187;
LW3_2 = [-0.38367102041677914048 1.4122808332961092059 -0.64008009007579713678 -2.1860913181968331109 -1.7074229886380916099 1.5340666307639261312 -3.4280515213271862685 0.51216955959878918136 0.97908791043996767911 -3.9525605284280893592 0.79720072388726637147 -1.0412035993957098334 -0.19056225830288997525 -1.2537671641566401348 0.55664183437371006491 -0.20153580836426038259 1.7309405734239107222 -0.2637084365729834845 -0.47892091403025593976 -0.43541985181577758102];

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