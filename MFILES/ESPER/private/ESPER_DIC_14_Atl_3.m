function [Y,Xf,Af] = ESPER_DIC_14_Atl_3(X,~,~)
%ESPER_DIC_14_ATL_3 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:15.
% 
% [Y] = ESPER_DIC_14_Atl_3(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999999390765;-40;0;9.1619;-0.2178];
x1_step1.gain = [1.0000000017546;1.00004398834998;0.0153846153846154;0.000309885342423303;0.0656325902039552;0.0461352500991908];
x1_step1.ymin = -1;

% Layer 1
b1 = [-2.6864217575510167535;3.4137088466502891748;-3.2395237835648842406;-2.1047510230590846803;3.3776657454695024896;-0.90291318880860216289;-2.4741443174674895289;-1.6588798540694114081;-0.66223615233449018014;0.23768714389385089958;0.06373230809979585354;-0.83106956743461279924;1.2846818889445512646;-1.1410664854502177157;0.60134414400064650774;0.37863981002192870218;1.9956790694212831472;-1.3272483561443213507;1.1438095657444318221;1.3947622660368010905;2.208215683883514302;-2.0721337512578617002;-2.3715498767614762166;-2.4733298129046712432;-2.783033274242204147];
IW1_1 = [-0.10967130662113350192 -0.85413816141454645958 2.072464056729542925 -0.29255303102377389068 0.89492486582811858753 1.0347107510862925128;-0.32759722244442335271 -0.77389450812046800543 2.0369853966539310797 0.24064001152137462491 -3.7184453755781339268 0.9646292889634280332;0.42358815566566715427 -0.52540377404471683143 -0.23066886864487290421 -1.7079771030828934109 -0.9607642415281438808 -1.217212270456833112;2.5687156053489714758 0.41708767797615142747 1.6957747511023593123 -0.024060429785722217688 -1.8479355175969944813 -0.55277183910181015847;-0.94140750359212888032 -0.74752156242091849503 -2.2328302670798958474 -0.08230607200329370654 -0.83205564955699184981 0.27216050556295745544;0.18885399870062899441 1.1993102052809376801 1.3037430866427663911 0.5933510422837889875 0.34169711517501089038 -0.4564424871928262073;0.70908779004219002484 0.08551691547735457144 2.4613919512574433668 -0.032711948117371297817 1.4515369024957105459 0.0017210382485952592166;2.2114371088407187571 0.065192805983068521924 1.1751311325859312706 0.7584802536562821329 -0.88689139233061187184 0.87115247119137084741;1.4805704274819317945 -1.3517999416153365644 0.51542386356692992955 -0.79914881209271970697 -1.4459276625629837731 -0.0086534192031984067534;0.35353307323415095187 0.93564557911061951856 -1.6558459602068524852 -0.32256306503120102391 -0.7211346394203820287 -1.1118127447898555893;-0.17332330501759624086 0.89995873639183732173 -0.03902455013122842209 0.86764616226834634372 0.0014787794852444881982 0.85246088339002246137;-1.0334950029639284885 1.2012927441270326856 1.3203009416749906801 -1.1311243966693180774 0.44831710846157496464 -0.47695148082333221895;1.0962674712666093946 -0.75515101811528884124 -2.5385327994959707709 -0.5608342250648237215 -0.13581062048928668751 -0.11685497995724893061;-0.88933901841083895512 -0.48609624784510357953 -0.84492220822421471027 0.12468708674884382281 2.3167658548825991183 0.36954442185220037764;-0.50258890003283607228 1.2533902230567692815 -0.84251162007569058154 0.76769514176660513805 2.0759734680285433228 -0.016951340902991024473;0.21411085842950999814 1.0518954928548902572 -1.3794477256429618084 1.1604269671087013016 3.131547741409737462 -0.31249605619076881169;0.23020782787086158105 0.47325822995373284519 0.97795076527951751988 0.47765469518312675401 -2.9632998104504726378 -0.17959817919754264559;-1.139694262510773326 -1.1678454333742178672 -0.13987029776521919811 0.14702298295449250931 -1.2510816306508238949 -0.52940383878760766123;1.4281479302258088815 -1.6184502196214267666 0.2409945737338056726 -0.19241188880533130479 0.88506121482649446897 -1.5264826163010309923;0.92565690957028667896 -1.5833859625452457731 -0.2513972585283151262 0.0091929044421481687782 3.9954632758073014287 -0.42287635521375643854;0.39116820765076737443 -1.3605862347789481337 -1.1183408171484383598 -0.73376434110142862988 -2.3129545575809040692 -0.052875147001880375142;-1.2930356543156791282 1.2553430514399324114 -0.30793564938027390143 1.1838215981860076287 0.57356782497016234323 0.50218078108303154661;-0.79716646964982684764 0.15035692116216670122 0.96612139550235476282 1.8240865707498801207 -0.64524439472326611256 0.43822170742752258876;-2.2695674492393207444 -1.8143770048911649351 -1.7969290436991369742 0.59677819353011740677 -0.12537678089797504866 -0.41464013165604196587;0.94582827640115141321 0.7891642758929326984 0.16764858610707467079 -1.3711941264783418237 -0.86884171556885181342 0.78921446130725703139];

% Layer 2
b2 = [0.81590194893481371441;1.5716264770554215779;0.59491388855920046552;1.0937455361110468566;1.5435339007612436379;-0.099210596355600630392;0.8546666634380509775;1.3385707045359080958;-0.12343275386605442878;0.71627647216821666021;0.91660083639522238474;-1.5466371973723145317;1.2551495152543814537;1.6187954609632679848;-1.8254641039433723293];
LW2_1 = [0.60021857003347234283 1.1316558190177357091 -1.9522560538502427008 -2.0388904095662909732 0.92331106208073854269 2.2724446093633341448 0.91529190393531856262 -0.58627147574891402826 -0.8260430667562252971 -1.2813311608965409771 -0.38577969385295213245 -2.0699118355236527655 0.868691373879891926 -0.27761089655387999331 1.0110872874150960143 0.85700782665321517673 1.1378983498253361084 -1.5187474472915674095 -0.10599256843170327957 -0.031535424952826685008 1.1180924194396220095 -0.18757681316639229774 1.1759418877978706774 0.67759473056791474921 2.1065411140300276394;0.70655350337698241692 1.0215593689531403676 0.92761278391933055243 0.17835468963108061513 1.2282514862490714336 -0.62323707901749636484 -0.88272387335349156157 1.2758611704974429646 -0.050696398283190710621 1.120689197006452309 -0.92853537712400735504 0.21960330135588351785 0.801410136730166478 0.1192180274409549573 -1.0901971711546341037 -0.39352098961892334472 0.80312303152666664197 -0.068660553484844080718 0.65041331986867323423 0.18113446193569371356 -0.61323266118979602446 -0.37087715791564690448 0.32062234712600079511 0.50932873075728757328 0.5653757597326027895;0.11495458825619968291 -0.42371658757780478188 -0.31522500482439091973 0.78296432852463215823 0.14448467233472434001 0.71081117628847401146 -0.97065001313036525232 -0.55576207923603226924 0.15369322164025919775 -0.77330648910477661762 0.27931724131306334691 0.50526555655863447658 -0.042494192559055440905 -0.026387373369545299301 -0.77244525378820183548 0.17633558063346049805 0.87874405528884358407 0.027103763142751809839 0.055484573163664938078 -0.23821643106590922878 0.073617081299216693968 -0.1893338240079189827 -0.14347588708773117983 0.15178983026821823987 -0.80668531246857877459;1.6449179951540060518 -1.0341719596551195437 -1.0762271863612260514 -0.070588235599495466333 0.74421535828995999395 0.15549709140710113453 -0.40397955844235394762 0.031890504609202251218 -1.2160852131396537867 -0.24678492287043399256 0.28767621665770493067 0.35489420225281459453 2.0917214858854440429 1.192472726978925035 0.22012447132441656517 0.70203121203447282905 1.156334485432277237 1.635058089408398585 -0.54617274787288372373 1.3575196795284452111 0.36833602242771085189 0.057725091439495432311 -0.0026863911328903894438 -0.18530673784265647153 -0.24839468996044641869;0.73058628181729678985 -0.71044940647774412756 -0.37847742588281862997 0.37941044949900593064 -1.6537472058374607098 0.034492794976067199519 -0.22375141536553452837 -0.29213091156179427088 -0.074799910138805714155 -1.4405280302127971837 0.13341297136731533612 0.21889203148964772594 -1.2460092069823616168 -1.9572744169096745193 0.58067160283481400462 0.22079233977875925565 0.31193814400889613436 0.66305586074607236746 0.24009990815954371235 -0.0012025061606017444643 0.898839329501951112 -0.055336534662074464419 -0.4255009080105434105 -0.5450875052235727436 -0.43643863893543360932;-0.3397387358921115097 -0.23800176404421108378 -0.62890805287528028789 0.27613150632166783005 -0.062760476013552335406 -0.10636630622054996098 -0.48814047660474546575 -0.20537222063653176596 -0.061272329015248823414 -0.57880344288058360736 0.071625303756439051561 0.48401153129459462132 0.46607536557494372298 0.13757890537945152709 -0.91060271740824938469 -0.07039255289640358737 1.054393941854768979 -0.7711433395145931291 -0.017940341792577647884 -0.39494343810671173056 -0.31528944005414188956 -0.5098851861577639033 0.72428033413196835077 0.11161552410896646848 -0.51750436350329642732;0.55073908443303376536 0.55248530330716738579 -0.19863252832528061065 -0.38869174839997472359 -0.94649574860916496633 0.48877846052509044128 0.058558148177180843752 0.025314417570804229052 -0.056610254481264422888 0.82658963488506664863 -0.3446921758318691742 -0.63968254472806107547 -0.62961933420304194708 0.048105950776887389242 0.56509178925536529547 -0.31702225548483864426 -0.92000526831967299746 0.28460178950985237645 0.17924374926108596617 -0.93654959011951710579 0.055596617077031680598 -0.00707714344596823558 -0.61038092907464092107 0.52770070287146986221 -0.2800898424092791994;-0.17896378927550263271 -0.52936517428468865543 0.74501390598385441866 1.1747206044500573263 -0.30708121832015039532 0.5641904638163688368 0.38220444806381426162 -0.043066160840240202967 -0.48355500076285795608 0.70710795450276220375 -0.051142414327123263129 -0.6281593701023090448 0.049873906542155121979 -0.79953435049091114983 -1.0782682723225853749 -0.38941707397689279135 0.10295476785741813908 -0.49087997548375200685 0.59978149317884998926 -1.1525403709482697767 -0.56938661790126232631 -0.8207236076136966707 -0.5529695606918543227 -0.44733459913717077416 0.26875657274188707868;0.33130221222402822523 0.63433610976605425691 -0.3954563951114855036 -0.82646897945818853515 -0.33249650383209999882 -1.0822775214600151017 0.92882588590306081056 0.25925795102377846568 -0.44630585488896923341 0.69920007551393603684 -0.45696486943314379436 -0.33564916200337824304 -0.051477946014693909405 0.073969192578542980465 -0.097931446994670462436 0.48148228811179710762 0.57278047161685352595 -1.2429500729551012927 0.044149273573496142231 -1.9080968956814261528 -0.53032903725649538362 -0.53593865344172841869 -0.48910741254707984504 -0.083679978855723924358 0.001056936482047174386;0.58604937253444822076 -0.74415384470242673665 -0.3695190129584958183 -0.18032340569432844801 1.9741285053772148483 -1.0115618995885133913 0.77461333957459965749 0.01232608922338475875 0.036294148427509352606 -2.8414446682762140028 -0.079783742826504994161 0.30244248993150479876 -0.72924148261500942869 0.84199674374095412333 -0.52363658018564296359 0.45830661852390525146 1.6516529134804738277 0.079274841041592325475 0.52112054603729207614 2.0339938490054150222 -0.81984212705986236891 -0.5297561631104353852 -0.66600869832203291043 -0.026481549336421505297 0.54391613316422926516;0.11247942243589995059 0.30708758543026692944 0.88247989874440047053 -0.5099688118166095796 0.9599951737860753509 -0.808141543654534944 -0.45813344619427398552 0.74220025184621840619 0.15891999782708621236 -0.43686478559546032541 0.2918259450764091878 0.034488398184701243854 -1.2849212768325952716 -0.97701688782601681105 -0.41274222224465900721 1.1304972735794796002 0.021868129402478536649 0.77264764586745560138 -0.17740832108018486646 -0.57195349056663113974 1.2012610228228615128 1.495596680162723402 -0.32607265186036388238 -0.97219923771595528894 -0.23829365080932279897;-0.5619650902356432276 0.21809174266352512883 0.85881347371088523346 -0.09251213097078439529 0.56982427268402957576 -0.78392502283648202166 0.31011534439111759776 0.17229655746703989583 -1.2031826024888521154 0.97730639248987882794 0.51351027793488257522 -0.29891784732972709326 -1.3193930599856580343 0.63371535249162103121 0.9006560876112970071 -0.84715353423910388919 1.1558539010407056224 0.63402746245588648133 -0.94741989951976568474 0.63555888557169426534 0.33315220331918937724 0.65704358261009110365 0.12483159967189137729 0.73851851551021430176 0.31454020641158075833;-0.22649319635097897585 -0.43720983239475985638 -0.053159414198535877294 0.11802845117713975376 -0.029643655274386837561 0.27134174045936432229 -0.0029726742319425951411 0.54942097413434187736 -0.1114370924444582911 1.0770435145029511048 -0.066605029739906396191 -0.38316829955780434824 0.50858517426234228775 -0.95441145563091711601 0.30574015132558002295 -0.38031021190998881565 -0.40888357042128004037 -0.86638931551464004599 0.89670308597898795711 0.61060595751042390233 -1.3035432166439155921 -0.22405706183274817778 -0.19955670936440550456 -0.15436449641334132332 0.5269443936633272374;0.32201875194018736437 -0.91370524384237150617 0.33203775696216591973 -0.12081443583277025333 0.45859424630417999769 1.377544440889298949 0.33240011334731739412 0.361531667617648933 0.28547948004004330969 -0.40627356441359641082 -0.81195062861867306481 -0.53518120971754157278 0.17147073548837321599 -0.29267334246228476902 0.6653316690307178316 -0.21766672685916821361 1.0114972990831427868 -2.586019161278748868 -0.81455021233427005178 0.062231561834613736972 1.0555404438741007311 0.59186850160992021674 0.36682124104860930336 0.80663038284252119858 -0.14487944895214471819;-1.3243220778186568509 -0.77955759501678711931 1.321080071108503784 -0.13267511873347662132 0.55581914366162898933 1.657688924339028036 -0.47707526746399597783 -0.59934089545081803418 0.040452316018037504097 -1.0573377035224178488 0.51721927487057584027 0.1193789350835771923 -0.28749746591255709305 1.5631286783487068348 -1.3279644252880871136 -1.0436915496256575775 0.13356366737616853779 -0.23014543057687131267 0.17700896025120280153 1.0166674295673436834 0.59563634954379895969 -0.10818130888968720271 0.59991745512327543377 0.4165439461058684234 -0.02426195400341144004];

% Layer 3
b3 = -2.2696919562570738726;
LW3_2 = [0.23795612250443115565 -0.21316085989293362402 -1.1648910942216696984 0.29700978956142903842 0.12721596452988817583 1.3383425267976605255 0.4638209339052918323 0.16745714061565442066 -1.148895001102421709 0.83952193268691710504 0.21335029098418936688 -1.4032028556561650934 -0.224768802155895564 0.28542149366503921648 -1.328391249615863634];

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
