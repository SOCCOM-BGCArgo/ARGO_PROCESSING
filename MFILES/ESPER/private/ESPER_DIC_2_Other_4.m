function [Y,Xf,Af] = ESPER_DIC_2_Other_4(X,~,~)
%ESPER_DIC_2_OTHER_4 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:10.
% 
% [Y] = ESPER_DIC_2_Other_4(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timesteps
%   Each X{1,ts} = 8xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 1xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-1;-0.999998793562952;-77.961;0;29.715;-2.08184650291172;-0.9;-0.1386];
x1_step1.gain = [1.00000000274156;1.00000063368067;0.0138606862425759;0.0002963841138115;0.248478071810163;0.0588359962254402;0.0419287211740042;0.00826143816767911];
x1_step1.ymin = -1;

% Layer 1
b1 = [-9.5798777381296158495;-2.1316433970577066681;8.8957846985216342262;1.2106658691894573998;-0.89193952722289193691;-1.1285605426755227487;-1.338281508759216365;1.1470108946242607573;0.49654651174681602166;-0.54456613118148722297;0.9068182019304508179;2.4203946625428196526;0.35591148269392597348;0.49315969661560304926;-1.9310067111537831952;-1.1005309983277893338;-1.2863994544472188153;1.0980250289509254635;0.77164237951126524173;1.6465758879123857383;-0.89774352890866970789;0.56125836086399494285;-2.2475334148373629084;2.4319793781105141228;2.8275355995045519109;-0.64583267706812053088;-0.29729215887601528889;-2.5394200637499513284;2.6450559652249241438;-7.5835312822336469196];
IW1_1 = [1.6567242341334089684 2.1137936393871976648 -2.8599392276820481129 0.35919198106453226194 2.019245590641978616 -1.2976740230480585669 0.71752053614487598665 -6.5371465825778001246;-0.14410683917123776054 -0.026772522110635389725 -0.015793231545261649656 -0.23540282910724061316 -0.149281323782987696 0.024035977678523751111 -1.0778489372619179942 -0.082399220623054256207;-0.037195133899022583823 0.063176215549877789979 0.63168642481781578368 -0.23802656042506439849 -8.9098597242510884087 7.1781284231434030474 1.2529271931364731607 -0.83384517528953161492;-0.57272838406671378042 0.48291829534307556537 -2.1772843020980903184 0.91633156010863936469 3.3977286643826074553 0.99033428777950849309 -0.10579192089188481363 0.36393542164225001212;0.27389207939391885249 -0.20829974209515164651 -1.5489961192760193232 0.46696108714480355495 -0.44354689591927942693 1.5476266076764559898 -0.56019458554433543274 -0.4380926722714394983;-0.20334163602444738084 -0.26652657143438129861 0.90169400506979346943 -0.1766606134094594327 0.77239827028837626521 0.38713011262762619946 1.3894887251570509079 -1.8608988411910500904;0.078663693351989305569 -1.0974077839017921665 0.68064118453955613219 -0.66970697927207312183 2.2157239596280020599 1.7890836847775599683 0.29033915874021770653 0.026617749868558554516;0.066218623969332521262 0.29659430493175553689 -1.239911351896605396 0.25732524067079814989 1.2001728061310532603 1.4412012288789708414 0.31695475780630411045 0.63895756178101814626;-0.08416546284035630876 0.09906535572284662361 -3.5544405690616946991 -0.57653372153895932062 -0.049199329093566099291 -0.045962707099328967142 -0.080408954860731338132 -0.010093088190735483023;0.21823896889933233623 0.32391437789067084596 -2.2973493226780790089 -0.29595598219050228694 -0.11221113915627711766 -0.79787994313685695058 0.25827109304983814075 -0.83738769051085915063;0.086774175427596203503 0.501266389642544552 0.41249898986354299746 -0.17573675229975013745 -2.3386045127467149918 0.55326875373581052564 -0.46956563927624783705 -0.10391990331031784878;0.074247123510850235317 -0.057971564309960574302 0.06054786808203742321 -0.010545298076434951942 0.30979916112583055776 1.2636155541003224467 -0.57695048018735051176 -0.58921799865850466027;0.06150130146598079689 0.50086739337816832407 -0.081219287999406294953 -0.10107190853549598286 -1.5664414492347229491 0.24227867281951739886 -0.44414573172267413081 -0.025579560579186887137;-0.046059454970905348847 0.32745504207684139608 -0.31625800213054361398 -0.044504929053322712929 -1.2066631876853088912 0.13477299951611654283 -0.72366252106512929387 0.82693722062239871651;0.81167880238219114464 0.22798916996500151466 1.1377715733339965887 0.1142379967941524882 0.72554658854535514223 -0.55688363770402271591 -0.15278496807683231151 -0.68435130491677742537;-0.075295578342983698894 0.42129386761707593445 0.06275801328586612382 -0.60518511321902412892 -0.98714945256801278628 -0.95770075768815576822 0.29653741445414560918 -0.54158294761322123279;-1.7095322712496023865 -1.7033529389579151747 1.1623285173819191129 0.62589121954756377075 0.57894084377676280351 -0.25175007311796049336 0.28557896492602258265 -1.5117961674701219099;-0.012849259286403923175 -0.35107296899701995185 -2.3941310724015769651 -0.061854572180806047244 1.8419333268911086687 4.4979206341261717128 3.1457830260253238919 -0.026889994369264071106;-0.39888835234391722961 0.018113798317044776248 -2.7412456301358663069 -0.49983091466761320198 0.52826913503598016408 0.77565281824344389783 0.051493212254615332302 0.17709143956590234215;0.30090208192814560384 0.6213068056426026553 0.16769720210677752092 0.86446948171033721753 2.6317777498939869396 0.085653253374962518008 -0.084103128545647637404 -0.30451952386530112271;-0.12310782095994705909 -0.40490818395031796806 1.5976053315646385755 -0.45837049526924728005 -2.7156107138967384707 -0.87679955912279083918 -0.091346094420629589394 -0.56931137173587054168;0.0060558641142446712272 0.16939863949805855636 -0.10458774886229182943 -0.039850554506218070772 0.58435279824379138258 -0.10265073122393866822 0.27206200421650666987 -0.14650044230588007976;-0.84030736325722277247 -0.0036806878686923734378 -0.481606715433700594 -0.18882797662984579512 -0.22919263637890069374 0.40604303993753643187 -0.94312347799147500105 0.015280339023554929262;-0.73995969570787056835 1.7177896063093129886 -1.8410086361691999635 -0.45722669306614255325 -0.67133761006849412833 1.3143048682220552781 4.100210141341047887 -1.502089170672902485;-0.21579723992506813501 -0.094169011876680946971 -1.0982541897843602907 0.21433451756918811548 0.42372399933620313872 -0.32057360409950313995 0.70762742313772131197 1.2148734497817892386;1.5506402076088288133 -1.1355589621082871066 -1.2464606171053140837 -1.0919238110011475484 4.4916083455772879418 -5.9871062336620619604 -1.4136669256749962997 2.4268790605790666604;0.081749484518501253083 -0.25371752627455318452 -1.5384419779527216399 0.47018568684005523917 4.2601353170664051362 -0.76142832743741484247 -0.55992290553512102935 1.2816836068361743273;0.099707069138773599115 -0.32489165589608076656 1.4073424780987013882 0.23563968896352271987 0.47320855032175163091 0.32760264439765524047 -0.0074274305085632354173 -2.1913994485981986848;0.27641291066794049236 -1.613834850806963761 -0.71164418992646594209 0.50899439556572512178 -0.50909007609889644907 0.7869143958862258037 -0.012656922989082823594 0.70102602448771267341;-2.7440623148846952084 -3.6272591845368991059 0.65100422337204366663 0.71389721140368611607 -2.5259780901154198496 -2.837744682501507576 -0.58413100640554249665 -0.27299409544734465349];

% Layer 2
b2 = [6.4483201303352375433;-3.2861494337061838422;-12.651326086147216543;-7.6083930977555347397;4.0358044083745276964;-1.9541346313198211249;4.2158485761065751518;0.55239705932550176026;0.56189708244792446745;-3.6264080240672540434];
LW2_1 = [1.0374556371734127147 -1.0613574529639109301 -3.9362999153318494194 2.7187955941026626938 2.5059459606942251497 -0.24188390669397880073 1.2866191959589030613 -9.5415234950072029818 2.8308337952064874798 -3.6825099967351118124 -0.23888310665555553824 0.06081279387263172076 3.0528285998978978988 -0.11746620327195256284 1.8930532353943718693 -2.0670002301960539448 3.1683626810946154961 1.7266531909559543756 -0.18858089749134954172 -6.3848599936269101462 -6.1937785471170245799 13.315332049336577214 1.5418957235428021324 2.9049777127512461661 -0.52653760039702302009 1.5309389652887377409 1.634572778472962673 -6.6557369016542233808 -4.944472755554588872 6.958747522167631594;0.45131753345394692856 -15.582619560816263871 -4.4704619330874271199 11.162942361599037966 5.809231198104289895 0.47075476381911846024 -2.9046268865047233376 -7.1523704307360969068 6.7873918095054790456 -2.8355264930742447582 -1.630570878814209701 -5.1165596154051442568 7.1141854435612721375 1.8581083995630220951 5.6408105446950864703 -5.221972504805619586 5.4374788939861495862 -0.83414990628054985766 -3.7035833022119275526 -3.258964608423799536 1.371372993909679483 7.3794270446805443342 -1.6268643262613922129 1.193009623369284844 0.35576323052339942077 -2.1875426510949704806 1.3465524803829642675 -4.5890316401658415302 -4.4633808832969563696 5.9521040843100445628;-0.22635242403857691595 -1.1918221830154254182 0.5433840167016465994 0.36234134923187244492 -0.5996735428400540302 -0.13148610305879368521 0.22666177080039279224 2.8728629468193829233 3.5435699017828552471 -2.344408793080925868 -0.95124001467373953123 4.996157107030780864 0.039200076463527772086 1.2128387862386080798 -1.335328183145269243 -0.047077096665889113514 -0.0070642545525008131752 0.10081771991368534969 -2.2540661996913016196 -0.019699390847026106155 2.3305225347527667168 4.8854179867948488436 -2.0655184896784208881 -0.69505048687402548602 1.6014949796371236168 -0.050723057310377135665 -0.54707849633443916826 1.1662603539668223718 0.34647635906488172264 0.45339855799218869992;0.29922404104598659957 4.1136458455730755546 0.35673175411676122382 0.58762443947839793079 -4.5097439109193340201 -1.950039861974628197 0.92187648799316701531 6.0198660923595017636 0.51468769936290736933 -2.8226372553998775317 4.4592956085945099076 4.7363966216718340618 -2.399401498400485 -0.99278875244499642516 -0.83571251053566708578 4.0673555205338143281 -0.054846534361559944049 0.096999295939393187505 -1.2236228080364948134 -0.32831454546865296296 4.078025749464663896 -4.1237497802155891335 2.345817269443449149 -0.12193714997446765402 6.1610933263812279748 -0.17507118649929603871 2.5518507767048834722 2.6524100944234847432 1.6604638848842647736 -1.0339999412943865487;-1.1586147415364687063 12.711904364153006952 -0.092859611410634257389 -1.504735612052494842 0.60493830278732974026 -1.4101402034642989136 -0.59004083089551018215 3.0955810550894864441 1.9153390354588188949 -1.7521482983859404392 4.906052863698332267 -0.86998116353538679224 -6.8677927649883727668 0.096176185413193657214 -0.45687604615832300148 3.1504510213095469062 1.1174200268265561675 -0.85008272062289536031 0.20652144508866349115 -0.53896381297386475318 -2.0959086186494331194 0.52028570886495417724 -9.7902722795012078194 3.8243632549306849633 4.4465875484883889968 -0.49493606793235556118 -0.48130574573114598413 5.0709745481715629012 -1.2556183718311266784 5.1285420942795463262;0.17605510946692490282 1.206538720242089413 -0.89042408370361436631 -0.34358149571914098352 0.6574212447203680032 -0.33719437937118779125 -0.70714044358721672712 -2.445951897591068569 0.22312614145244008546 -0.53251192235593758362 -0.56670361004646785386 0.6234844754491324359 2.1156811682085847615 0.34737034434991742593 -1.6683607740416688259 -1.8962151617753428834 -0.1020316274881479246 0.22428831894301975991 0.31229048872184966834 -0.046514841905892598262 -0.81007484149428399256 9.0550035400160684418 -0.89366831712841288127 -0.096935083188500942297 -1.5628726332663149812 -0.025722087867924514626 0.15274179457034683938 0.53995118348029436106 0.3904708909588853194 2.0554682590470254588;-0.26661027314570973079 15.454545384793663132 0.88539658295369216834 -1.060598397545596816 0.35205595427323038704 0.19761697717874979197 -0.67931716792133212124 1.4570759519510740887 1.5444950005768343093 -0.66616129390295575341 -2.6077695820056301379 0.65847648382488077612 2.6062119090038549096 -0.5657928539480925334 -0.84266913555084499432 0.47745048685884100559 0.35645057793932805046 -0.53508020051110649629 -1.4811980915617197763 -1.8590566171318405608 -1.7262708328739224672 0.18624675085165684485 -5.9518608314110421631 0.83114545621946367415 6.9992496768138989793 -0.33514057584505080323 -1.1461754083263755444 3.0686332488130547702 0.15669191238540580713 0.72697355039411415412;-0.25266596981192979143 3.2666625633667489836 6.3553031535949608966 -0.25761087627063100625 -1.3919248682161318165 1.0394335475538944991 -0.10084289810396561415 3.4008212368668799819 -2.5316421281914247565 1.9560701119689685168 4.5766824598289019832 -1.3763507730069481916 -6.5957613026852301985 3.6825041444512383038 -1.1759639809773145824 -1.312004473989575315 0.017750404566172734572 -0.32340319340720224961 0.90914023356449713553 2.4084715102133418263 1.8543128291203205826 -0.76248057227545895653 -2.1894610721813401177 0.080368064063750729864 -5.6714396358465863912 -0.33808180865030079687 -0.12529992819093391243 -2.1565351638189991235 -0.0226258051029897371 4.9190313509889502797;-2.7390104450145535964 -2.4865868986724399825 3.8682681433870707011 4.379386970948881519 -1.7510756301079324349 -0.6555122008098898112 0.69563815407822982007 9.206452682074603544 3.188772364597172615 4.1483342021190416915 0.57661513848951262418 -9.0101096948851875368 2.8066792623574556842 -5.7376995370671242114 -3.9960898328700560356 5.4500492953979398081 1.2891943480239811404 -0.76186483914126290617 -4.6609068404185052259 2.4559661392633644894 6.6224552957911848594 -6.2968846884849565271 -3.8920580274922675912 0.56759095155648042752 -4.8801503492645883497 0.63300196548246800976 1.1021955140399264383 1.6075723490478446287 -1.0900091758110574958 0.55845065277000782711;-0.096051050129512538511 1.3680035540532264893 0.9994076200729692383 0.15040124155906417491 -0.024111666772946370452 0.52614284840958081002 0.106268527006641883 0.30884719697010731831 0.17475337323050141314 -0.019384838286373093363 0.51588213775214786239 1.955487331862967082 -1.2877680559697182172 1.4447884293130264233 0.26533158300984849287 0.16219346356583255742 -0.078409392804952834921 0.06163999696206477058 -0.15717548147724402785 -0.095389071177890416564 0.36740236387040020594 0.22480495212487239254 -0.61858921606690164818 0.072766506496893340605 0.2414497307806184323 0.001814972492204005184 -0.064483859458793105213 0.10825238341260640107 0.12413996744880434475 0.10787109153077557855];

% Layer 3
b3 = -3.6768565475058294645;
LW3_2 = [-0.015610184492970052839 -0.0060825908320972675233 0.20440152825519561364 0.062028204641624846982 0.039446608553522125029 -0.34995500317868011297 -0.21255903467191850087 0.11408404985804931897 -0.02305193274038285986 -4.7129909808761576429];

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
