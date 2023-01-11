function [Y,Xf,Af] = ESPER_nitrate_1_Atl_2(X,~,~)
%ESPER_NITRATE_1_ATL_2 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:32.
% 
% [Y] = ESPER_nitrate_1_Atl_2(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999994516886;-40;0;9.1619;-1.99149636107278;-0.03;-177.232103709389;-0.52];
x1_step1.gain = [1.0000000017546;1.00003626070655;0.0153846153846154;0.000309885342423303;0.0656325902039552;0.0604107502911417;0.604229607250755;0.00425400504498441;0.0145232735458572];
x1_step1.ymin = -1;

% Layer 1
b1 = [1.0034533336957607297;-3.3981278662410532299;-1.9812698041419141148;1.224493208621600715;-2.0029733322451477306;-3.9227826209925318324;-0.84942524131265872978;4.8092545863568538778;2.4099763596454524262;-2.1730286549846926647;0.36659970460097901501;0.96067144167454221471;1.5307334492290303185;1.6430171334436871611;-2.9935548817626584928;2.1030283925690702596;0.58151347776406292489;0.049657917450925584402;0.88434161676263400853;3.2004479378293573255];
IW1_1 = [-0.25199136478160272246 0.29267438127123507341 0.66613744178422418241 -0.17058490713137500205 2.0785454958111806434 1.4025662598780768153 -1.1640379516713532659 0.66847286583920439984 0.71949650230047945687;0.62653430666987819375 -1.1570728120419448981 -1.2163774090676806594 0.29610723373815550641 2.4231136554417136431 -1.1927118110797541028 -0.5566815824808839519 1.2590959428568559986 -0.41630083954782998168;1.6795129636604522183 0.82039800961617970199 0.12194252470039984582 -0.62219660306986157572 -0.43297981227605503163 0.7522332290167708102 -0.4761956846395419829 0.97623868072588348976 -0.088417621416817776536;-1.0249880503700634637 -0.2456979279642775571 -1.2990442498510959002 0.39092385976754162602 1.7976636893170119258 1.5099361736380438881 -1.1902215800132009971 1.4364301748860717201 -1.3719115833022652584;0.33177542581888153084 0.0042769600725806192557 -0.36615711309629472847 0.98433916714170810458 4.0950039876843158382 -0.09568225804226283826 1.6593979610087119525 -1.0216369797382276463 -0.65976712029416428873;1.5718557432225708226 -0.99480310256257065138 0.7133827332147630873 0.61349613958502402511 -0.081151547314133848365 -1.2155650568246454402 -1.1425134831929535384 0.7910703651741988196 -2.0498632523311202114;0.50764860394219457618 -0.070950429044935134359 -0.2744463867834255888 0.21648894570531457271 0.179323155157761982 1.1593057071740775488 0.54686193853144338117 -0.16463611530736030342 -1.8704279075765199991;-0.44469714732871257201 0.14287035576226247979 0.55761236967667604514 -0.93618255317837117957 -6.2490560047910257424 0.3201426035389927649 -1.1313213739063172625 0.67322729292730520978 0.62521717188263781217;-2.1483707667429388088 1.7890193008098593808 1.2861346669307651691 -0.74011109033743938124 -1.2241648639038746893 0.37944947979549159545 -2.8867572191416916993 -1.4671222450558465411 2.6745074672970421403;1.0992744372371134887 -1.0300299735935056766 -1.0417941875123915008 1.6349794732296833999 3.7661015872463501353 -0.49473891721880319805 1.9591455280122311855 -1.1655379111655139734 0.24146811302627071982;0.11285099331776959897 -0.21762010952773908068 0.40770995596196868416 -0.40928711677191692564 -1.3774762925522467949 -0.49031210697222865624 -0.095528572002282069464 -1.0372864859325316989 -0.43799494395101723487;-0.02987152126003637595 -0.56276681478529733127 -0.96748179542372791939 0.23569729717192594576 0.49724883679668452352 1.7592948533990528936 1.0288453957463026267 -1.8779495562989740876 -0.95382011076750927714;1.0357210228223376891 -0.70630469414635610459 -0.68377682108796189464 0.8575403739168144801 -1.3856386487127059226 0.33462817127939220585 1.2361311176885083629 -1.1260300558215039501 0.19564490851587762466;0.29114637911241958435 -0.19376472770616096075 -0.43245536228441755489 -0.086154437640434206047 -1.5637489579308572196 -0.33106209791204438408 -0.0060838605520329297471 -0.98251904103883058728 0.82196552233004349652;0.46040097580808891919 0.3721965511477279831 0.18919366384051874985 -0.052219208360675811731 -2.4578996723656749701 2.1392639997668334573 -2.4359681756563511534 -0.83164795832785531893 -1.6700048992449740837;0.65536041076158269014 0.59891674685090501296 -1.0057703694978938902 0.16567862953536954818 -0.014470884781486368503 1.477580127469565463 -0.015543795777195532765 0.26844502127681280079 0.79671502487300138728;0.22193670397167564179 -0.04776333453821910574 -0.16068814189263930237 -0.73905935951390056005 0.96336088553802012857 -0.81783381389654141724 -2.3899341511276674943 2.0164089081806886306 2.0991434686767558304;0.4079655274518702579 -0.57954663747544310493 -0.46591692787664329423 1.1363739162058423826 1.6033717352856995397 -0.079902333267529893068 0.41661174088644803426 0.75499899515850876774 -1.1350940897232906224;-0.82213207743587979337 0.75768162498027280538 0.22691250896951731253 -0.3036534224086525624 0.74410632936815690464 -0.49863262490469834054 -1.550544257727691333 2.0648240057128952607 -0.38352356048158720014;-0.034830042548600705088 0.16325338373384326585 -0.046550087424973893291 -0.16356366336461472533 0.4789508672558047131 1.1393860756027853753 0.67841403268581523811 0.13637351219076820907 3.0692251809473418689];

% Layer 2
b2 = [1.1140520839872494196;1.5867903752254635119;3.9685145552001377389;0.68094252390262677288;0.52779839166516762106;0.024659182682506559925;-4.5146429121171962606;1.0628925208524293566;-0.49884574561681183091;0.77681077151172861495;2.5520762631887636829;-0.35194971027157029297;0.9905723247016464672;-0.59654118625894891892;0.66486667452523506849;2.0621671830806960202;-6.4792483962405809095;0.91465977983139634677;1.3340879539817762645;0.8553242581138095213];
LW2_1 = [-1.1714737898705211983 0.006517198249681695546 0.028651795110118852272 2.089958486836586804 0.99301091583250078454 0.30869916607552866372 1.4080815933669210693 -0.88280454494838767054 -0.21164215162391369907 0.9536370228127291071 -1.4265439192028783744 0.51729587286284828274 -2.2905360152173419586 3.9319594765267473768 0.7136330319256661836 -0.20733441399509489056 1.524205961863740022 -1.9013446995656513305 0.19268583931095090156 -0.86297335202139802135;1.4988906661495182693 -0.32279911758941987054 -0.34534462232883611499 -0.40151871112802667296 -0.88570031518429470641 -0.44663011039616173914 0.65310992553714497522 -1.3385002188343113616 -0.31650050505012644608 0.40507539654568036624 -0.48524586262725472041 -0.58589396061427356788 -1.5298931529119861317 2.0047485676774527974 0.88592150658034862776 -0.78931083626141018694 -0.1584714866130962907 1.0522010251203195441 -1.5252961447729225242 -0.08710455170640937772;-1.4342999946809598111 2.6574602637977684516 -1.6494621525251249583 -2.4284609638233454199 5.8017841260772060608 0.87290903141488929151 2.1772801409429813901 2.2929450667486048765 2.0805112748245471188 -2.6352775132769079569 -0.75720706886913036193 -0.093194857238805356436 0.441513686205322331 1.1831935198657683728 2.9979086204935341264 -2.7847743755433471691 0.49040029299475973446 -0.35853467886855455715 3.132634055374948101 -1.0523242229791756053;-0.56180432657621437453 -1.0781092163500503034 0.44571033219847305817 0.4832163804119577466 0.2213267713743921894 -0.29357791087405804564 1.2591032599688469773 1.7293493805723290091 -0.34664937405939039161 2.0913016305503644432 1.0839133014866166871 -1.0971043467528525106 -0.16011244551716244389 -0.14161649867471914899 -0.42141362299431572858 1.8424176873098512264 -0.78196376857527072346 -0.26480793216604986418 0.54210050065889214821 0.79028568424664713543;-0.64162396228235929474 0.2861811686933887211 -0.45596290215648671573 0.70547257347704628305 2.6567361438522976158 0.92573258378966682969 6.2273384186749822788 1.3999532599664796795 1.3978809188988419177 -1.8781916128221090023 -5.900204515976306574 0.56950752372386415079 -0.76733693459109453272 3.4001476528884198025 0.17847830408380355172 -4.0859408313182941086 0.94999191242295311088 -3.5852434507227437521 1.0634430282129498924 -1.3498140324064569207;0.33481611527478144508 -0.21584792477204886585 0.10080304834321035523 0.089690677546018837685 0.12815664681436148009 0.30903927670669506522 -0.61444849519445698949 -0.19384052169340254546 0.44322387549076402991 -0.30284113091463982226 -0.11968531311645411175 0.136298547958960653 -0.17224224335715648815 1.2161850017791124845 0.56120317913274242461 -0.052185268219822854296 -0.1292095386283476921 0.0031924825553983652779 -0.013607450072860023671 -0.13010789800163355401;3.2732884456387147942 -1.1715161601768042665 -0.37943083064985289043 -0.10737473728589246036 1.1145978953709707593 -0.005900141700205316711 -0.45061967148286008467 0.016955299256326876456 0.17045468462761051431 -1.1752048853252672345 -0.74841685964815107823 0.69540678370364727989 0.96277680910761886945 0.74623303682836172879 -0.19449664295234136979 0.44289735285081099914 1.2688978869550640027 -1.1435651963431407996 -0.42844109068019481379 0.44941175523415832771;0.35414403350951167981 1.3676905177830211979 -0.22333028202055521527 1.8499545166267257201 1.4220546800414777255 -0.5437004520088437598 1.74686279970191638 1.5533864226000064779 0.64964043827449446944 -2.1337743099065455965 -0.39077137986860394081 -0.75090235530618243409 2.5699272472401979428 1.7536018208245054417 2.0756968146413319332 -0.67389234746255755049 -1.6091054052081099623 -1.0990824586103189464 0.27983503412232785257 -0.47316318701701304272;0.51931999563596653591 -1.0001847138531074766 -0.33568335602170129572 0.96220395025051586035 0.033206134437973225393 0.1049834917880359958 0.70572068586013814162 -0.59722871537693733579 -0.13503621644305732752 -0.16681888773349645705 -1.2582008807353006752 -0.37580102457695280238 -0.24452556458362142977 2.2954695853723219301 0.5936447192620584179 -0.090704760763850214578 -0.31542325857586456239 -0.81358064510707273076 0.25018387311501361081 0.19638366910608665017;-0.70903551695733746385 -0.21982180843725038222 0.32378919464263278272 -0.21464354976893601501 0.21722897937720148098 0.0076829431765088832679 -0.21946620646650452646 0.075768979682404188902 0.32585115589163721195 -0.26414677520656870691 -0.44719062550782967547 -0.10038213485815701553 0.32069306171588690679 0.50140449756214156007 0.32619518503418387656 0.30874296868190137522 0.075865925622966234876 -0.22389092712408797947 1.163783148355779673 -0.24070870007343694863;-1.3740305359781230798 -0.16947622382353286219 0.35388397204578453392 0.41950553926128858473 -0.75385529561634534446 0.87510705041250247405 0.71086338819435723746 0.48188456128149959046 -0.076323544945769788983 0.09739244830123236496 0.12956023141541594956 -1.9116611329819559462 -0.7499306948472281853 1.1869815691964342719 -0.74848957697596119054 0.090572868655771840429 -3.681146276837270026 1.7584596316360365797 0.60084894026965018909 0.25228523330493712873;-0.43494833494967122256 0.60158386301016220088 -0.3157225883300771585 -1.5043637467496426829 -0.37863996469172300774 -0.12920768536927706993 1.0622216971350804116 -0.046761020253564475335 -0.55921554727027400489 0.29506668957884429894 -0.22856641522827006452 0.23794084275496821257 0.28419467684522065065 -0.62401847437968349652 -0.31961059648494094532 -0.26288727940999062627 1.1891438574787589477 -0.21468461861786888867 -0.3298924192559373747 -0.18936635708041402171;-0.18322410301283287293 -0.17855728180130917004 1.2276287906556322405 -1.7527444487117012084 -1.0439015734397061319 -0.060647901436763240746 -1.5048562136666585598 -1.2853052670919611344 -0.28855766848547220738 1.0670396340897765342 0.74540454165037151046 -0.26716310597001846627 -0.57804334174532123747 -0.48849572610139363471 -1.7721250069316321074 0.34454475781157800629 0.12366570150397998562 1.7978001712859041916 -0.46433236707317726566 0.13846419356394074063;-0.98105849640607301065 0.47797337883493651312 0.01153121046762072173 -1.113230908147733178 -0.45483922380343300151 -0.075142395544453358114 0.49342516409893555362 0.040311660950831926242 -0.41371624372948434534 0.34582198054956508892 0.020930880963305576598 -0.40609701465770775552 0.015053691733308203343 -0.29946345562723442058 0.12229384482920073463 0.2579532799074260585 -0.12616762847827073246 1.3412542488722376532 0.26917436175645692575 -0.2193864195193737654;0.16111795854380822623 -0.77059295465594346286 -0.55876003749637537243 -1.9415320557346140706 0.41562307430961586974 0.21623818819148699077 0.98337299830508539955 -0.18935211476972488298 -0.3975755551170946922 0.40871309366827851539 0.94733399192649814236 -1.1965428977848409087 -0.07246287807301973749 -0.19524963878654924554 -0.89763567586464887427 0.28259572726219517325 0.37913505500026717376 0.42540182211784111077 -0.78186919321557990514 0.91859712960720374397;-0.86598302198178189748 -0.76929000152460713657 1.91815182875963397 2.4181238278181798051 1.7715572061710511687 -1.4354056165627633312 2.140527039123516051 -2.0033101648702023034 0.41384238861396754894 0.54867456776184930778 -0.82619250114225351034 0.15193775169241099943 -2.8162373871879498033 3.6522168934267975615 2.2916455594792877193 -2.3066972895018293954 4.1501153391816547611 0.89098187805604001088 -1.5221563833743281258 0.35664998999137792479;0.47408774651756813157 1.5158322487368476317 -0.83640892802577038179 -0.61740234114008230026 4.1682332274915978587 1.4490434152783315813 1.5454329437850218287 6.412888716292608926 -0.51413299746312179916 1.0319262191291438402 1.9674028433173862407 0.83684271202805471379 -5.0145295247316941101 1.7852687036872880721 -0.4835534769894962448 0.21442644489953049947 0.51025076181488493443 -1.7844654541159568062 -2.9941210163154892143 1.7996066271986523866;0.74356413221127681812 -0.20410091055934714843 -0.22856096113049609309 -2.0584067707621871435 0.34883049363787016484 -0.18618483918547931721 0.89553502060802114837 0.10414043118863650084 -0.27408227160274911682 -0.19617430093203980235 -0.32249115863073385402 -0.38979715525933605891 0.5046510623660060979 -1.0680710983602503816 -1.0433327697120322597 -0.40499606663881060342 0.79072082514939667774 0.19936759495866346326 -0.94293946402731088074 0.44556098067603250312;-0.99466380316506897152 0.77544762875085138099 0.2381119582514574784 0.38883544149833337844 -0.18043714258288609531 0.76870605633407018686 0.40223783578551824958 0.029723350467656284302 0.52603988040581139973 -0.0338132641430715572 0.34090676169681377594 0.5752749858197635513 0.17164052624522541102 0.48209268424239076767 1.2269655336595162165 -1.0901121001757523299 -0.5959829497639945739 -0.077589235207845289044 2.6863986418483807306 -0.35992656756774932258;2.1513995906654352908 -0.68386628028714724348 -1.4109143461291886457 -1.4905857448395369858 4.2871871773391880822 0.025330900764707842132 -0.9905542162746130419 1.8091623993657264258 0.37042316196133484585 0.78798312878096465894 -0.70331667898850724896 -1.035088613229812049 -2.9519554940016639044 4.6908423311059666005 2.3517025064564069936 -1.435231784467741134 3.436192106551775538 -1.5657891582646354944 0.97961434807722314044 3.3202559220119294281];

% Layer 3
b3 = 3.221763917488563056;
LW3_2 = [-0.14495730030570444513 -0.78374797890799663413 0.16940216371393593264 0.11270498685089472279 -0.10498480692655387525 0.86037925841132756322 -0.29155517397795321921 0.22176465934887631959 0.47177453393207602339 -1.0782483217343192816 -0.22562697986527124883 -0.74981439794851456604 0.34009015960320138072 1.5800020730636330502 -0.38520202833474909143 0.21674913696260900142 2.5116243111345144534 0.65357813038644962944 -0.43030930751096302389 0.26409363832380416959];

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
