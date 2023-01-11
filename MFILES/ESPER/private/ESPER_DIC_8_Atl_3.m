function [Y,Xf,Af] = ESPER_DIC_8_Atl_3(X,~,~)
%ESPER_DIC_8_ATL_3 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:12.
% 
% [Y] = ESPER_DIC_8_Atl_3(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999999390765;-40;0;9.1619;-1.99149636107278];
x1_step1.gain = [1.0000000017546;1.00004398834998;0.0153846153846154;0.000309885342423303;0.0656325902039552;0.061098595529566];
x1_step1.ymin = -1;

% Layer 1
b1 = [-2.9321543968805587888;-0.63149040087636665319;-0.49229998076110326188;1.9363553394532937446;-2.2154715227311836401;-1.8564699849438801582;-1.1537444800731477823;-1.3991130138968590479;-2.968081120946526319;2.5725134103826139409;2.9741498596502244922;1.0327018774165279691;0.20893249422394755133;3.0145961434480428309;-1.7380752291297778545;-1.0169551809742085613;0.40948769631556852788;2.4439251517176185935;2.4791526537150980758;-0.19218296907927062689;1.2056910603745325528;-1.8729857668277252181;-2.1020747256040586137;0.68757131076309385698;1.7699473360622199891];
IW1_1 = [0.062837885185706765512 -0.3456410529279514976 -0.10224896230354652005 0.66000929950340569974 4.153650984442419869 -0.58235447641521886908;-0.04910243100355833612 0.5652119602978521451 -0.10977485867304789813 0.41747007057193286084 0.73540696508476799398 0.90929437077659736044;-1.3618495736032640853 -1.0866457775733506796 0.089744068260635231904 -0.98771762642464011606 2.3637568394801715677 -1.5111078531684971615;-0.29025519059691790025 -0.7952835735287224761 -1.1044328001157970309 -0.46262325432917228429 -2.3059137825230764385 0.26505645360386920206;-0.0012594394300569442674 -0.095009071615342902128 -0.46812903777314740505 -0.7993993278484463616 1.2555317967163182047 -0.2413177599303316101;-0.086851122656947665757 0.45633710689188705656 0.93572589175377474113 0.055467461899768144562 2.9454383569261448628 0.52624181887652521361;0.059086797092816542576 -0.014127264618528208306 0.6602942632865398398 -0.29514111842294793453 -0.52144924794782554489 -0.07783166022610402135;-0.1602857840043661608 -1.9787619630749602173 -1.5935798403436323145 0.61851915745672192415 1.1506539257077965654 -0.41845091582953897191;0.65946332224800185973 1.7861213381504914555 -0.71166737886561948212 0.60158352850590290917 4.7412221549339914617 -1.1626614466415650107;-0.3522963332499903788 -0.069705626534801168326 -0.081001645188255796959 0.023167016182011203329 -2.7683620499374774759 0.37220644436617222883;-0.22048204716830907546 1.3654474352973848816 2.0201752662518419612 0.058205103796711586683 0.54373372204224335036 -1.0807940519547076796;-1.1929792352651553689 -0.025828443060784467278 -0.45161119922030495966 0.4068122682002616175 -0.34014607580741740289 1.2474248482885290645;0.27765633527206773357 0.36785060722021811763 0.0013377715165680354945 -0.17792743689628245929 -0.18481464059376750964 0.073642794857776000539;-0.15564559591290463847 0.59466523717157504514 1.0117667198141420126 0.52683533206925814429 -3.4878424964238146089 1.943569176866612791;-0.40498460347854819963 0.072424510021631707812 -0.14683594101340835536 -1.2090674847081979149 3.6841346662301814874 0.71569249257385381124;0.7288693083942930917 -0.017014930393454693464 -1.2608139217764098294 -0.62850934419250881202 2.0467411872793617178 -0.10712117964647151214;-1.39012721744050749 -0.41731243163322734313 0.2785143107116482275 0.3075083145211741642 -2.5216932983731745743 -0.85627771253133833085;-0.25193659442985200103 -0.7506083516142256773 -2.4721999880032248242 1.2917112422412779793 -1.9802738454379007482 -0.83557293010820288259;-1.4286383726815736761 -0.30924166103071215606 -0.35583649995428229262 0.45064539923392932241 -2.7118141489544154155 1.5880171799430879975;-2.9165547683698025416 1.5444435759933727859 0.060524602743721449272 0.95194716310820315375 -1.4359537103990334117 2.9086680880080555411;-0.044090137074970323083 0.068299625590693746613 0.41080755804681701182 0.92741020644885596713 0.36637156024426970058 0.12689131272481235513;1.02234395001631162 0.37643730319964807673 -0.83436870222877357328 -0.67813034940323302369 6.6757777738203030182 2.6030262907638266334;0.23913565194798175328 0.56518280910446050402 1.8339148610457485233 -0.5663060245654715219 2.1069102183483652269 0.37507358711512917004;0.77886908146529909569 0.85994927491013339882 -0.91567575687042124244 0.081940061156293167444 0.43887642595361148201 0.11051755778212095371;0.13888742514578056308 0.20129557128619138262 0.0029231257558047125944 0.83632031179637678076 -2.4856782342816887699 0.23217179760977896263];

% Layer 2
b2 = [6.2537972788542344205;9.7594466807676454323;-1.9305577290031661342;-3.0591529912120805079;-19.387262256806391036;19.483936508385969688;11.771124464822792177;0.4212965665719259345;-5.9312256122778777012;-1.5067788921500007149;-3.40366722833705726;2.7720967770430919863;1.643829526050181622;-8.6508127684191986617;-15.429332744934226795];
LW2_1 = [5.4090088189393021878 1.8719777631226035552 -6.5840242698072737326 3.275698502053101091 0.42046211804460564831 7.5868714843416933391 4.6177645654015906374 0.94143995094720633876 1.1566953828425023953 11.605287591293411253 -1.6409608366789951717 1.3119725844831731099 -2.0219547312521894966 -3.3225948434172738466 -2.4457709366451636512 -0.44977828733533448125 -3.9041237478634478464 -1.4379344162717304645 0.88720794635136068784 -2.0046593209120620038 -3.3954286238537907394 -5.0353101091888969876 0.016050003825309470257 6.6080953089201281969 -8.5320168573005439328;-8.7514656422154004645 3.8412136668594007993 1.7237702275626212156 -2.0557285884661573583 5.6956616317439170416 -1.0873665987054297322 -13.454609635818270519 0.082463526798630559789 0.9810801541715877816 3.3959791671014203374 -6.2501345564030170721 1.6101277861563429994 -2.5792205112708614934 -2.4042004362673998763 -11.114596086266438135 -0.92632528373689493328 -1.2591201934951783414 -5.2978179043362265688 -2.187855181414394945 9.4569564989180943115 -0.27579858015353364564 1.3418987978646621517 -0.83372258187659720186 -0.83314851591614202952 -6.1721658091589146622;-1.4014321536742235796 -2.5774952411450731127 -2.9508627962045781423 5.295719298198818592 -9.2309898318926411065 6.1721348065723242726 3.6489142774038456629 3.1099397592566928594 2.2714949682057348213 2.7020965612750775797 -7.2767728572107861496 1.679007388854931726 4.6260620217346781757 0.05498366041209021976 5.1716342928818388458 6.2155398208787566716 3.4052586889270832771 -5.9472013821871785311 -4.5906373723233073747 -2.0276093626179845231 0.17248413238867676966 -3.3692226554084174772 -6.9405540666714244935 -2.6437658506141250569 3.4419569980438518009;2.914701339067284902 -2.9135573734326265516 -0.2274987584940692531 3.6821277358928727885 1.4237594590048896492 1.894167666591249688 0.5397191929060188853 -1.2700836847018923681 -0.61868307457882987244 -1.1191662132127013773 -0.70957846796638146536 -2.4165862752568290084 4.1939921537475850499 -0.029766834891820491416 1.2120166637069182958 -1.8409623233051848334 -2.9733114561138997978 5.8134802047500970801 3.8353088894935996045 -0.70381550805611625687 4.4487777300567943328 -1.5120222138912382537 5.9791995750850182034 0.7926074422355370297 0.062257017696514732608;-4.8966940405643804013 -18.445219611726511033 5.3445970896628605118 -1.965765835554741825 3.1643124563646938263 -2.1197995557024804647 -21.470483399337208397 2.6681704484041799397 -3.0475348678547011794 -5.5208964246477805204 -10.363679587773665602 -6.8735794924015856111 -0.5753440722871726587 4.8199446110422377032 0.80340609076740776562 -6.2372411977897277424 -4.222841092178460265 5.7383075142954238501 10.945858607734917811 0.35114395362548300739 -3.1098486624511592069 -1.068929207663402936 12.851862132069562961 5.3068791982251575234 -1.2357603546988442655;-2.1875614953270892116 6.5124281881198315958 2.1122663566525727852 0.87060210304138530013 9.0610026917060242369 -0.23158938090113040387 13.425039474852306753 -0.23060890848339196668 4.8110600382528332375 13.742360527713753271 4.127558265396375603 0.90874596941504559577 4.8512983253953221308 2.423494454777933349 -8.1566883287905955768 5.0711392634071223995 -4.3951712806649316789 -1.9778277854604409125 -2.4157939351216293211 1.7237366846943065202 -2.6162893259283830361 -0.21998045746371602815 -0.75555111276660147368 -7.3928599740824472164 1.0058158843565820906;0.079334870982640331993 13.522008664939248135 -2.0194725375424238401 0.88249035202119452403 8.6673062582906243989 -7.3803122604556303799 9.8786499880511353666 3.6607951251790229463 3.5747973398094021746 1.7880059945069810556 2.1724790824819502078 -7.2296250750153427944 -10.29506694998404015 -4.1216529776474919444 8.8836754718526069752 -8.6807163154912245773 -4.9552241475603393539 -5.1178194810414421312 5.3896534462933081144 -3.2764289038987275404 12.829640858315043417 0.10849487719053110824 0.11629160889554965297 2.3744545630376889811 11.15903712386994151;0.43628015926735458363 0.48870918098355331516 -1.9526099050687835135 -8.6235433170140058934 7.5076212933239689917 -0.90210169652020388842 0.588854029218443209 -2.7629970618635644186 -6.2723253913672216697 3.1976285379441922174 -4.1664911654482903458 -12.042486247883179118 -9.4726344722497355377 -2.6887252324336530229 8.69553339939361436 -3.5543828206146486082 -4.3835929975881766651 -8.0338658324381988507 9.8615761521120575139 -2.9848836183590705673 8.5401145703006289267 -5.9720329330840673521 -13.942702184056665615 9.6849804743325993428 4.3590685012819134769;1.0334911278193696926 -2.4097185664647331471 1.4648547195089349771 -0.95624588815559685706 0.81835055134443301927 -4.2073816526161369467 -0.61752101720071206348 -3.3038137883663623029 -1.1927328639299770163 -2.3834291198632442388 -4.7975655580725629079 -9.5780858189730277985 1.1392326840506274976 3.6425579084737753632 7.5297390306399565318 -0.5478899341525363953 -4.8345507094726221098 2.0805775392021117476 6.3616021406171983088 2.9876301267587765942 1.2617471708894687499 0.21711453876683350916 1.9568297796630056418 -1.5572261396898947616 -6.6208459678108324553;5.3091423035575022737 -1.5032350196019659272 0.17575063129585910549 3.8543599079572832267 -4.4844818468211835594 2.0309637678109591086 3.5593931877318705226 -1.2614811080348493633 -0.62402358637152299803 0.12348901451398880402 -0.14669137424276637582 -2.1806624130057676858 5.8103631406456006658 -0.73261883942910754897 1.4546659924303271261 -2.45536654400474319 -3.4798260921703856319 5.5385364346240510969 4.2095143949289681728 -0.54379877625426897847 -1.7182034887174417737 -1.4277193092929498786 5.5321322492301545637 1.1057626661578487948 -0.98495753465161117735;9.9650527212323911641 4.9917906497350035266 1.5030282795810285812 -9.7098023997379794281 7.7640256652259767023 3.8851220744707339882 3.1983061825591332372 2.3527461786741614702 -0.43069363669530996708 3.1013189408089418819 8.8629729019344551944 -8.5681895161313228471 7.1401444247985779157 -2.0279469300065517956 -0.22736396794030050428 1.3339398416821990345 3.9139462105267006642 1.2725573285414044644 8.7543501855763370401 -4.3327548428710880302 -0.7513851583928637945 2.6352875660354784237 -5.6633840482193802401 -3.3044638092127409656 -0.083965808013000256538;3.6361765885202377291 -10.754528281952241642 4.0320886083037628111 8.1821388854975634075 -1.5652634311615001295 8.5735216729663932256 4.8123693568954868383 5.4511391586312631929 8.6333154129046825176 -4.4457133620669377905 -0.16243650431330564077 -10.717287756673437471 2.0354509469425061674 1.6839959535164683579 -2.7006804526536529742 3.4244034762554038664 -8.0373473312838576987 8.5316393122672256055 0.6499118093580343869 0.40563264499547885267 1.7618328801287908014 1.777942704744046809 -3.2395540318637188193 -5.9176383779982044686 -1.305977550430382772;7.6136401101146766734 0.53729095552157091209 -7.2654764690658284465 4.5557438786170090239 -0.96482286506339609566 -0.58376675943210132314 6.8636651149885503642 -6.1400685052041197309 -1.948082386090561835 0.44863402860359707924 0.025818014253709772182 0.34752653638778358003 -0.21285575957611088671 -8.7996933717497167748 -4.1885228527594025039 1.0404752760938649203 3.9903676700840318148 -4.0879627963548266223 2.4284321709196521333 -3.1521664749930677907 -5.5418286567691525946 8.2817964779951758203 -1.428231758714034294 0.66594390141786730108 -10.770145395053983606;4.5086681329447486632 -20.632882215461442854 -5.3406327589309086079 4.2223550091729729061 -8.5559335173442647715 0.3424448328647597406 -2.9043715312517113247 -1.0346114375737485158 -0.32499986736792119135 -9.6017399495843402946 -7.6934430356662586448 5.4159808582864430448 5.5331627160117546893 7.9663699048175633521 2.1948578829909783927 0.17540519853459654587 -3.3050533067019811639 0.15806694447456520747 -3.3325270609466688398 -1.880392395748524903 -1.7155525536970608336 -3.0084461527052877017 -2.1790209017848587436 -2.322887983035729409 -1.9570473453729915203;-2.098696596258144087 13.824049201140621079 1.6662256831668740276 -4.9564207064509657741 -8.9376823454536182822 2.3028919291301361838 -13.631284952258306831 3.845228248847674557 5.1034708574084595156 1.015469633770130109 4.9335244960920752177 0.35588126364897554321 6.6329030890366711404 5.5871534255386876566 8.6083487131049931662 -3.7177935104585304593 -2.7283735023213138149 0.43332860838012787763 -1.2944878038648701679 -3.8053616048435578989 4.2494715892833783144 0.76623207670812321091 -1.4561982501950450075 -1.0010959255938107759 2.6236858581777391386];

% Layer 3
b3 = 0.62309666940002816915;
LW3_2 = [0.034520674426012996905 -0.030099398031063779452 -0.021533158229197803257 0.69509348066586917891 -0.017062737514626825225 -0.0068254545096289778033 0.027482797181459968361 -0.013568657958957787771 0.018959928681061440775 -0.65723988426804902385 0.021791189822290255723 0.0070199378262166024778 0.0042887435693114532576 -0.039039713361749543707 -0.0046104907041229752779];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.00115315130997576;
y1_step1.xoffset = 678.286531932822;

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
