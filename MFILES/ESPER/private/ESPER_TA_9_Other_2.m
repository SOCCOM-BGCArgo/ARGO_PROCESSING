function [Y,Xf,Af] = ESPER_TA_9_Other_2(X,~,~)
%ESPER_TA_9_OTHER_2 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:02.
% 
% [Y] = ESPER_TA_9_Other_2(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999540494077341;-75.751;0;29.715;-0.9;-133.803853032615;-0.1386];
x1_step1.gain = [1.00000001035699;1.00022983623563;0.0140762793578401;0.000302617642608564;0.251952632905014;0.0422119037568594;0.00432440760873659;0.00826143816767911];
x1_step1.ymin = -1;

% Layer 1
b1 = [-2.0965834311956559333;1.0553015445261137994;2.5490703127156599805;-0.83474280126284361447;-0.3783929124174461589;2.8193774000666578239;-0.59804303990269414015;-0.90163491661118788656;1.2538333606424569844;-0.50210917195408122371;0.30189280062961770801;0.25583988332333124083;-0.78661509508892801801;-1.1117637837841782833;-0.25367614652210612736;1.5478759283271856795;0.60976972288588038307;-1.1271807778150615409;0.20954008312713537365;-0.59500613943899083313];
IW1_1 = [0.76199905657109956447 0.60393893601305781438 -0.63846970731254670817 -0.015337655329754068187 0.39141620261274207015 -0.46822713087484235661 0.69974928150948279004 -0.45851999334146986698;0.60067546405207938864 0.0055169576480667411808 2.0404606825011679661 -0.36965899940292656911 -0.078879767895587094007 0.56661125030455827822 -0.97825843799448208848 -0.13961429514735509971;0.22099701185006123683 1.0560513367578574506 -1.910903345151636179 0.22739237816444543561 -2.0959774771119614734 -0.20560722354599919592 0.22260522781920616397 1.0739331963484641808;0.49376890639486320511 -0.038666378916898280771 -2.0531630188621523558 0.72496373688374249333 0.64845340906748272225 -1.4536499699148870413 1.0378129794605590597 -1.2634069814949639987;0.086141658689967551266 -0.039377111584634527652 0.27823765892517637877 -0.040054110535563049178 0.44290398068577546686 -0.077061520477956682584 -0.0047444561115159234257 0.03654981583249655025;-1.6231892490590065226 -0.40973521333072604467 2.7639056783881401635 2.435001652183665044 -3.0448358701098734436 -2.4394242888721326956 1.1224307494988692291 -1.858653959405880407;0.055097076221360584414 -0.075507459763015469245 1.0068907495987939882 0.1784504961992712524 1.141108954437271894 -0.14779142934404124188 0.062123968924967520877 -0.45815860479018388762;-0.39053779703784641253 0.28448454193626826481 -0.015376313311893185717 -0.23965868702512541089 2.2321071961195273481 -0.59981088771535751025 2.938043693955543656 -1.5705332334785908621;-1.0025977249556683457 -0.28218213369964534332 -1.5210651002240354046 -0.029050598439849081162 -1.2285544484028316781 -1.1511707980398229179 1.5231629825710040205 -0.28581989628119591762;-0.04117431946748335847 -0.013690656678950542322 0.15801805170734722905 -0.36692572345550272805 2.7080368976780544443 -0.027327385233577902723 0.61911643482243206016 -1.3689829216950850643;-0.12829760903330147959 0.011194480138823178356 -0.79258425916536645595 0.16300008819619793421 -1.4556364487175088573 0.63983317070993950271 -1.842047182880371281 0.73330446580259789346;0.31656414916550518202 -0.5009571445227962716 1.3466721194039708021 -0.47104945334575892835 -0.084067582478338107044 0.44579418118989638797 -0.99565468515205879374 0.45718455421322395926;-0.14995919885586456166 0.11580987756020837132 -0.10450741234023414039 0.61384327449620623707 1.2222828912305787608 0.18707339461954089477 0.45752454980532442086 1.0937671824542893884;-0.098204751356861935685 0.2332046659157262225 0.2610889038260520989 -0.51173056000453931436 -1.6149005886814107313 -0.21449996898427514669 0.41786859155126898635 -0.27136143403493362403;-0.46825061822726782301 0.81575358760259653046 0.19799996819165247008 0.026465640563434617016 -2.0327240023771753208 -0.15891412794043885404 0.38867448441948831972 0.17296719942514218737;-0.36014397991235092 0.30320092806499499005 -1.6748025523194076669 0.52094593681970957633 -0.78551286233727701802 0.3731500353775733525 -0.19466368729060265497 0.67272671043246812239;0.1191572779351355027 0.056583316516095188975 -0.42308581155487112735 0.50710075637500151569 -1.6935814105890119041 0.45570617632021070964 -1.5692887947993361397 0.83271699537720511586;-0.63425111148263524363 0.17935113933864807456 -0.52126248853210099821 -0.95960287376332753517 -2.8210894203141729797 0.85859723018342537504 -0.47030187695609215748 0.25321808627507513823;-0.046057378484095350124 -0.042132419294048352454 -0.056793127563552044113 -0.32846466691638093183 -0.9530657995179978581 0.32488097324110992714 -0.067024815204231938393 -0.57681863601056915236;-0.37494630907837062095 -0.5780544064165226148 -0.18966679523397875218 -0.1379336714481080306 0.70694534603879044532 0.12413957738780055651 -0.34737285100343540112 0.10406481844879576781];

% Layer 2
b2 = [3.2006326020374902264;6.4373639597350535979;3.0372361743570728265;-0.43541538045035044213;3.4664584220310028506;-0.084145267545312776858;0.98144017372227032592;10.153221932256547788;3.5305276968145382632;1.1439566858326188914;-9.1239379278518608629;-11.664314031763929336;2.4047335409232197811;-2.6179724340414538553;7.8272808335470696051;-0.73975833243482325141;3.0739653323614231084;-1.1122152249619037256;0.46654166285745229947;-2.8886726916613145733];
LW2_1 = [-3.7122999172442949778 -3.2469185227545405148 -2.1791258008516951605 -0.44676941991492330475 8.4275963793750552355 -0.15603892180055095595 1.9442289092803328732 0.56782201609028337685 -0.40857859673241930798 2.1554654154617511885 2.2964848043267980238 3.6249919749910004541 6.5724468351818048362 1.0713252407350442041 -6.0131740098032464559 0.87183827088725507348 -1.4428485380456250731 2.1890967210060248149 9.2791438427415489087 -8.5136243948349008548;2.2648192381694625119 3.9718710927357587792 1.0857153750423023109 0.99858562989974097857 -5.6263034842746284525 0.34275339242133628925 -6.0307014015808153573 -0.95104428790384809567 -5.9024303001362206089 0.75661331897444872308 -1.8474882057540786828 -3.3620333166381417556 0.65137674316914784267 -2.8938570886014924888 4.6130665482869375182 -7.5096254859925766212 0.29672503747723599776 0.12985375731605386385 -0.6036783672271738288 7.3578303991499192094;3.2379188102022355977 0.46518653872130277804 -0.14131249674805612671 -0.51443268933661490117 -11.302001891693281621 0.95936538184995889811 -2.4330254087306082411 0.34026354289501498096 -4.7242981369048910878 -0.18961711699916139473 0.071724642113410216893 -1.1886939892194952062 -2.1849937518154161076 0.84495241850626245927 0.24444007622247340694 -2.7344658181442000888 -2.3873944625439413869 -0.93841729813186247622 -5.6471701894286745471 3.0220122698099918246;2.5512440912153571126 -0.23239116290947064991 -0.27652226643945310158 0.35461176615522571609 -18.256335940755910485 0.49696501022084010879 11.624677857296086714 -0.52343088713172014614 -0.099090139099791835209 -8.1288027894172749654 0.47913681942752994747 -3.7010994256619000886 -1.8204306665664828913 -2.0181112676818862006 -1.9183667816857217225 5.8828256761479185144 -8.6665645104266584298 6.6869528399570583588 -2.2285922189105926705 -0.39613984642067556763;0.71291194512253008853 -0.10510811433580651009 -0.0049158587391908165759 -0.17412373642770959115 19.480997065705260241 0.13779417038602559398 -3.2446050979758855881 0.85193085018349434812 0.80141811645939098341 -1.3769636480481250995 -3.0510907351705927049 -1.2876775682711487647 2.656979062313010953 -2.2508731130450008173 1.5232146486557078813 0.23815964363449435015 3.5670647088102995781 0.018188136615711211763 8.3768226010877366861 1.7153559536788649531;0.21149785607830726142 -0.035164240484148900912 0.0073670168662938722431 -0.082518218855851172 -0.38726263626889101399 -0.044133111311418934952 -0.0092452546858226829912 -0.035545856327265180907 0.099929774836882923061 0.18095855336399907509 -0.10370349779258546119 0.070049664945905754099 -0.0088956123653141858099 0.16159974227679790859 0.058487639501939220343 0.034163562499626771995 0.16490741868002134485 -0.078266311948632802387 -0.12195391693254321874 0.12660347536895555365;0.89516043379683807313 -0.04701598548596732513 0.074225162243481884738 -0.29829685851787168804 -1.0172933092351708062 -0.18930491955021599981 0.35982561787263522834 -0.027402560033267193973 0.18738877093581773892 0.32168962498365760716 -0.60820126090431425592 0.21495429023747636421 -0.022345318031467265724 0.41945791888357625288 0.32735243919473988594 -0.049745821216395279263 0.76270243869085430077 -0.25012394649420066628 -0.011504609776847092573 0.64344579543047708547;6.0750206242108788501 -0.86343522575685116571 -3.0753672492055481591 -1.8562061961621110218 6.6636610607113153648 -0.011122089945668045541 0.62837339486694065993 1.7701390941406587309 -1.8833235240135779431 -1.5712675911312323773 4.5970742973501463879 -0.42441328229924030202 -1.1284661109540787827 4.0253444438590912924 3.8367094951091624822 3.433181980725679594 -1.4255509979825817801 -3.7368819040086949101 -0.030823883409461013544 3.1789544495435508864;-0.14373435644919221876 0.6779776060848280439 0.26730429764740193388 1.5505663715527571167 -17.770810242277004676 -0.75783845654401393865 3.6559248067192435805 4.4127636830671113088 -4.722172266918189365 -1.8791568486282825834 1.3007312662837644535 -0.17623463831503505261 -0.50804983416662485851 9.3548317450850930754 -3.2808635562639212324 -0.23233346011823549637 -1.427799391141130414 -0.47088610729485225193 -6.8685167216252800415 1.4677558954317044826;2.7798872618809618018 14.816056102534625794 4.3009813938893710628 -1.9892126789488004945 -21.468918134175694945 1.989532898381839443 3.3953785109433418654 -10.25689865812290158 -7.3167841615774120712 7.6951798274590519355 -5.2478592253445386362 -14.837384865612079921 1.6362303413160563803 -12.37988061590273503 7.0295804159791792642 -9.7816959884919985058 2.2180152401128863104 5.1826812470549077361 -5.6848080036750303634 13.863916813393904803;-8.5912623388082316467 2.8398294165733641847 1.6213133723862245272 1.5340802567715436933 5.1505175197451968927 0.21000568330072219569 -5.1848594586187415345 -2.5804524365153440613 -6.4611677522854451539 -3.3810572906093541157 3.6028246219951647333 -7.1748688778975635572 -4.4563222783127685034 -7.7929800658741230635 5.0375029519725265104 -4.9186762011785933879 -6.4466849867398892115 -0.69753847024146486966 -5.5416133772267830437 5.429477159171738343;-6.2004542099635147068 4.0790918071787629628 -0.36863237250414376822 0.37399636923230522578 -1.5470339745489181116 0.00034208708848775535097 -3.1422960390767262773 -0.41036079858324586711 4.0829114298076616407 1.9159055464066003793 -0.071345062335804496079 -2.6373951957516319489 2.2597183633743043174 2.7498023863835352465 -1.6835547325196842383 -1.6184199496063575996 0.073610577230442841845 -1.0254549773887875475 -1.6264191894195978438 0.07191386098391987225;-1.3610713372806333421 0.23989084505968985561 -0.31119361531919087493 0.10680853026179720766 11.54415071015512062 -0.52860506540167029321 -1.4859770866246613252 1.7606835753817398071 1.2942087954668437533 4.0575379158097621968 -3.6698485518507868441 -0.95662082763789346718 3.21545185546673018 1.5845846806674384055 -0.66695299336327684703 1.3151417145313146762 8.2674238191805287101 1.332467417174112434 3.3675321302545806468 0.69834394477917360877;9.9459992381026154362 3.0982165348319723464 0.64113441458667297379 -2.7281453486951194876 -31.768804132042461674 -2.1115573360358275146 1.2388376371473917015 0.92513649292871591445 0.12197290701227435361 0.66444854730709890234 1.0728651749788777892 2.9404083078520644179 -2.1309078910896133152 -1.0896074139980003181 4.7445422108348553891 -4.550201199763887594 1.3949102104196560159 -5.6779906556174282173 -4.5016613062296890746 5.9214723310908050635;1.7410525299686707701 0.76651632366428590348 -0.98191750124101584429 1.9570218920665243267 6.9180785844832000109 -0.63339041782916727374 0.65857770671629767278 0.5450548231613637773 -5.4945312593014241642 -1.89371699890505929 -0.5190666507202291724 -2.6736953455254308665 0.85079564790081296799 1.6987097819124921116 1.5722464522385126706 3.7279260045577795601 0.017177435084002167648 -1.0911447183586957532 4.3439959745684317838 4.7159407389970757762;-3.2455664679931919281 0.62114365451749631308 0.87503004893045788215 0.25443888563051059171 9.0485527247354227143 -0.82777913199497243912 -1.0294479606316511155 0.45070694763557261897 1.5807693034206877769 -0.43238127620485977687 -1.2008880573529052693 -0.58536362060012525621 0.42415708869693924399 0.10281130967179749613 -1.1139476564348949061 -0.53979023349977561796 3.1603116106531663476 1.1614551544293523211 2.331220214743738417 -1.7825890996647777431;-1.9746423420763326018 0.85513410483931862061 0.038444941980777289081 -0.14076746319165486798 12.71765660034114731 -0.81491681108641567111 -0.56315621895976764044 1.5642363712043638735 1.0092509496057973095 3.0955673251326221163 -2.2768895508896020985 -2.0913600688330462596 2.9617773156842761573 1.1783126443790077342 -1.4326995089132075201 0.147894680546887064 6.721159245152229289 1.74120215886747709 3.1490411571367693533 -0.010923391707928198011;-2.0125928131188390857 0.057612710759083592127 1.4168596395234440699 -0.16715081047243518508 -2.6403484541516952788 0.16135954212738515268 3.5704617551815349685 -0.94789908898885177901 0.17217787156835340223 4.5693766721030426581 3.6881398019222264573 -0.62212690817863669501 0.23657465405026475991 5.4942306026807807839 -1.6750728946992150536 2.1411668788512105976 0.13197388379431695959 0.051780236206770122176 -5.0833420386951830849 -1.1474056757743416046;-0.1298727973975743466 0.078813331235130429375 0.054498655922122656026 0.036917671786574179915 1.8983950938118532292 0.038726046889863974254 0.23113836287672018255 0.15410136041541835916 -0.31239894709572291198 -0.31460852915193615598 -0.15475420534355532887 -0.12542183483102561925 0.0468913362546380455 -0.45704092150094088876 0.1239942279889604676 -0.13133915823903249964 0.22951509369181319342 0.20193010756698334673 0.42909851332673765167 0.044244936506669101106;-2.1005558538256563494 2.9180143798781115194 0.59906832602997828197 2.1913197448749057727 -17.255059385520681303 -0.44391634731385898327 -1.4731904201073351413 -1.1253661260087517437 -2.9169726363971104632 -3.5925894389502452952 2.6065869468051672264 -0.43701295380224258746 -2.7525919072322300174 -0.15640902882247126326 -0.32398596421119424393 -2.3202718642326884968 -8.1346796289169525096 -0.39294234867202809669 -7.2285295073745539085 1.5566982271026474915];

% Layer 3
b3 = 4.1027280589771004671;
LW3_2 = [0.067232103422231057066 -0.068661427377255557225 -0.1354905184104393534 -0.013538670689092676194 0.22650984153661238296 15.230913974886851392 -3.2980766292885701496 -0.1773798349846543998 0.077191451411763395418 -0.010818083359824042997 0.15693756814970760693 0.24332322618605434217 -0.29258612502413322565 -0.03163186506025739142 0.064510464499289960072 -0.25841576294065354835 0.27748520780540619768 -0.058677416745623565963 3.1030481276533650181 0.089453209143604939246];

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
