function [Y,Xf,Af] = ESPER_TA_3_Other_3(X,~,~)
%ESPER_TA_3_OTHER_3 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:00.
% 
% [Y] = ESPER_TA_3_Other_3(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999540494077341;-75.751;0;29.715;-2.05132667284751;-133.803853032615;-0.1386];
x1_step1.gain = [1.00000001035699;1.00022983623563;0.0140762793578401;0.000302617642608564;0.251952632905014;0.0584431850247228;0.00432440760873659;0.00826143816767911];
x1_step1.ymin = -1;

% Layer 1
b1 = [-0.079407541840380679865;-1.4430160881479157098;-0.68783862770941206488;-1.5381802101584602926;-1.0621986672204573043;-0.4806892613228141431;0.41539518141725745304;0.87396026804632487295;0.36438488494951842833;-3.9062410564093137388;1.6791701988642140364;-6.4884229132624504643;-3.6332962303398406156;0.28364350925815245974;0.45830333126229672036;0.84009300247950990936;0.41609842064014523899;-0.62723346041095817682;0.95314538862908859684;0.16763258132063110772;-0.57848232960556655868;0.49289785791528739622;-2.4964073753240532128;3.486739679154268412;-2.9752243028808691427];
IW1_1 = [0.018172891684385836242 -0.36688427414126728587 0.27270898053228298519 0.54858724148587556524 -1.8306108274318204643 1.0263761968486295739 0.29027472679405713629 0.52121013581429820327;0.023119216137888653345 0.13732540507053395129 -0.14995704229484138614 0.095114427660147410348 2.1950620703927326005 -0.73810513632391749894 -1.3345522012947417245 0.72609338268642742431;0.42653350049648863784 0.5041102100471164249 -0.66449177507823187216 0.28524166250748789198 1.7018021724174261511 -0.75502913434975926776 -0.49845843299027625672 -0.9718836447954447566;-0.068266442159014992241 0.13814296329511452366 -0.086340782917491823434 0.10599529131960538919 -0.030525699494940151252 -0.66719014723219083329 -0.45549752224227219033 0.64227526240226062715;-0.054905718991557529718 0.12715147631377057835 -0.12734886592528868454 0.14055894228092002951 1.2727756343365250746 -0.43669463686577997263 -0.94206803535452876552 0.81916052024308649848;0.89009567983293524929 -0.75161444812559330408 0.85830006810525927285 0.49216476237144163308 -0.36672867100999401169 1.1135049447282716617 0.22287690783737482136 0.1501296307504731975;-0.83187635474067789243 -0.52594857768071356929 0.40201570915325601785 -0.016169910641813958763 0.22396482517498672582 0.8546252583963886984 0.33819759230817425522 -0.21080931160857299678;-0.57825047342512891912 0.25339870287798016424 0.060951111956351521071 -0.41944153412917140056 -1.1297454540178117899 -0.43101389512927090752 0.25888567721443428438 -0.10970291767878119493;-0.60650819593405258523 0.31524527128511337093 0.045568052472866965585 -0.49023766560350701127 -2.1912105802639825214 -0.57771690374041928528 -0.040016419687466601995 -0.84358259500636845196;-0.020717702526806491825 0.15274567182501297413 1.2172152725972871234 0.65093157995042949526 2.7477677359378822608 -2.4413171285901618113 1.242450027155130865 -0.67952057381254049595;0.15374499261708723941 0.33761201757911973864 -1.2784849337309691641 0.40806716529978892538 -2.583411043994385814 0.82236798823028667282 0.2920024658744742041 0.16697721910673729573;-0.28203760581842207689 0.16461188502076326556 1.141191575083672749 0.70638378577838778494 7.9361196377451328487 -5.229483579328556786 1.0524475336572971518 -1.4070979590875349885;-0.10791310559745594111 0.14647451860871318896 2.3325312368691006881 0.11527799905685304216 6.096783118977964655 -4.2934389307820772785 -0.86756562714964358118 -0.39983245245972642845;0.61424787697549398935 0.43225748881075254371 -1.2303412885670050869 0.15678064467022465034 0.83306953840409592171 -0.4576585950123476243 -0.53323646815920333619 1.1689650114395975233;-0.084762127167801570859 -0.030233239735430513656 0.022287309865571439116 -0.21852800321374821202 -1.2180738488248530338 0.11536569419435664663 0.11144792577270143386 -0.1507551140525977551;-0.31994451823664699086 0.41067046851723398859 0.79329843238317498955 0.080355480398415210508 -0.35040254908125545441 1.5363662593664186851 1.2585429611644303183 1.2662708400039999201;-0.25714690415228214082 0.36504986635113434446 -0.91993832583536605529 -1.2291654831282616822 -2.795434247832121244 2.9759581773910155356 0.60158841310452981777 -0.22264477069558949451;0.11289275549080196592 -0.0051356639529380083975 -1.4072434526703891944 -0.45902530342121056206 -0.81237563315327443458 0.44990248596289444949 -0.55952509247505055612 0.1620845043973081212;0.44024795561729207494 0.043711653203688097258 -0.82914533442356064263 -0.51927299074179245064 -2.0778795318771452472 -0.033854399193498807552 -0.21591013707715464576 1.1848012043504292379;0.0225189970954608569 -0.082286720595379725385 -0.70211950896174157855 0.17800834888128486133 1.0106784763295997287 -0.1809338449049937847 -0.45263094510982493368 -0.10912895208720405082;-0.70570183654393581207 -0.1023284944009286157 0.10699361222270037386 0.33010159597777899432 1.7617339881376377075 0.48902810556708614653 0.019576966375066066212 -0.52402678797838719849;-0.091788729100026805741 -0.49766841374880277371 -2.2220814461734232914 -1.0003477093217363336 0.36798340082312008814 0.44699073811017453473 0.42464930002248757468 -0.25223662581751210565;0.17809684743111617311 -0.016590024019120661741 -0.36807834531046396132 -0.27653523004427721377 1.6828938902436738623 -0.30907108908403302383 -0.8912657432002513147 -0.13264198321145034831;2.2266474314025606063 1.1670520950359801748 -1.2440016855479369973 0.053459746259334035845 0.41728413932726293734 -1.1608373896389969637 0.44973782415578200666 -0.96007152096459857926;-2.0892733473543780853 -2.1578955568997719716 1.3144130574096699604 -1.0039521608425936439 -0.13466052656513557984 -0.94792919824250188032 -2.6888875443076245553 2.8872338612216128517];

% Layer 2
b2 = [-2.5602846633208624638;-15.694668127924849443;4.4378834484552163175;-19.405345707811044065;14.588504363649709816;16.606004012602081588;-0.85821233626768711478;3.249616974147564008;-1.167922660461063078;-23.963893999404060509;-3.6782766588819479914;20.289989381805444424;2.0978361111129419925;1.5187337940066238939;1.8815305594687548041];
LW2_1 = [0.076949094678738366726 -1.3083391986647221472 0.46519298692149035901 -6.5014887704198978113 3.7756994601606717055 0.87225540996320671816 -0.5081360728290603701 3.3791313363681418913 -2.092947400339858266 -0.45906685588567369694 -0.49231775169355046495 0.64985174276701340723 -0.35039413915762118723 0.094721276658527503445 5.7744673863500945998 0.37159549524104246165 -0.057074853158325195013 -0.94236036240749898596 0.32398433197045578558 0.94159761848973233533 0.24178803441936425056 -0.50154381978623863159 6.1408237044057916165 1.3859431385709002349 0.09987352513243603902;2.7603313578308492637 1.100131338172385842 0.27997944147205333643 -7.6214611381184287353 -1.0601376223520124764 0.11824834535465743313 -4.2084208887782592612 -2.5425443305587007714 2.5406424863958423366 -2.3936746025527173565 -0.91746930679890037474 0.3604830958936312002 -0.26732909026694023913 -1.4470749808209970588 -2.2948654003379616029 1.1224961174059600033 -0.91970051229468896548 1.1131717369801796824 1.5083328703449672137 -1.4853280825019963274 3.6923666164798785161 0.40222114709650402808 -8.5315080980038988656 5.1157895678883642532 0.35361363393193512294;1.3419051758724325207 -0.91139393859829542777 2.4078291107332243648 0.85948199619330434285 0.0073815253888410807537 0.13113543032188870741 2.3589503023863920284 1.4539555817218072775 -2.1298962185357472876 -0.78508918987402376288 0.68467683526643163816 1.2158947638117068024 -0.73367504367683178401 2.7597850071903251568 2.8636121773477491992 0.04016594872611863043 -1.0568519830329994313 1.5038467356208067294 -3.7672854615599256078 -3.2659282925434114375 -3.0301706457562214148 1.1963519046570234572 3.3002444491094835577 -1.080568006237619505 1.0948602980875508095;7.9496036998604493462 -12.922275484995308759 -2.6983209991895016344 -37.928994274913378604 9.8916301609372450088 -3.2390527975335561273 -7.4905111915243969278 -4.4244647161922951994 -6.8974009057521064037 -0.22741483918657445784 2.9168463271660951541 6.3609479277783878004 -1.2589145747472911818 4.8039512437182274596 14.563776586355153597 2.3825662850214679445 -4.4160152315302765302 2.9038996706825357386 -4.7304773264259338816 -1.4975238068075102937 8.1683657346614904071 4.8744863576242609682 -5.8281120074501391315 -13.100397221972714945 0.48947802544150409432;-17.833679832895306561 4.7044699040525799205 -2.73930781210782337 10.456791487895484849 13.882134144586350644 2.5236945709027880369 8.3457690639907422536 10.268031544184880843 -12.0388659190275078 14.406655570851263803 -14.166817293617151918 -0.060625312895099078836 0.78795627431603054713 8.8669560816004899806 5.8507006197454325402 1.6891742097834137049 -0.92114992791151872886 10.960306496395887166 -20.742013410604634771 -21.309256833744967707 -11.677945017094655711 18.196054623278691054 -2.4569635079765532737 15.375619087975962884 -5.0606086041450666357;-2.348346703275223657 -3.6337673570386224142 0.49362635210234973027 9.3979386195814509364 0.16335294305454553521 2.5485055043457824908 0.92657149623427781382 0.26077919417787748246 0.06577262314011252109 1.4817574943922167119 1.0819862872190206815 -0.44451790954474568185 -0.56087933083395546952 1.4358947136828938973 0.88053341208482638347 -0.60408094041557269804 -0.21210438019788360742 -1.4063122670012879745 -0.098646631336238604915 -3.1524512679081073685 0.84268315581236163947 1.1597418716402152494 7.7823508468430624418 -3.6213607404515988364 0.46355499304184760723;-1.2305492233329506924 -6.8124202445671384609 -0.39308661737361694222 -5.1805420982035821353 8.9045567815579804716 0.30503449908279917091 -0.42051797963314846429 -2.118596934864301673 0.53822896056030233858 -1.5973918367217541103 -0.21461926656658197876 0.15069646264825756199 -0.014764259985099172437 0.48697051208286440405 -1.7508130995032529231 -0.2005306765610838371 -0.61225948845751010374 1.3963546731989020522 -1.5106675233509565537 -2.3024399730524223173 -0.84118161302179850214 0.21857099177951685998 2.2072520623449722876 1.5053857059226960491 0.070067356180959716738;6.6703100808459021565 2.6747091699078775484 -3.5684849075123574735 0.56151093340379454322 4.2456630894699802781 4.8078633181107832684 -8.4590850909052104356 5.7144253711463433021 -4.9084556476287453108 0.12258169659710002608 4.4569767709123597754 0.36426714361414941346 5.2568658632045730528 -4.2185297985667533638 3.0109628245436179661 -9.5713253829840336806 14.818380172393419869 -5.9044457417300382218 2.7670022261140529984 -3.4081118401318377664 -7.601599979105178484 -2.0023197054602035472 0.18173167620964056357 -1.9108266132271654048 -8.1475959837494205118;0.86929444937986710684 1.9467857678643321062 -0.21767470703593405634 0.15548004574538729083 -4.1641519760591387822 0.95161119334981658646 -0.32324713298767004943 -1.2861227889148845183 -0.021241124474107354186 3.6771277623051439321 -0.67961660204391793894 -1.6238166529875133204 0.34091409087555812896 0.35445371846309325914 7.0691073880891179826 2.3981984557455029439 0.077696650328836255883 -2.9631448766690038354 0.38660782021103995421 7.1553800830516198772 -0.73413295183790083787 -0.2755730366937145237 3.253843187116704172 2.9438954507183519027 0.51594285141142126605;-1.0356935461930401043 -2.8137524663107522649 -1.436431787831868645 5.75310748562416574 3.4429307106702471408 -2.3218625549619753556 -2.8017871662985505132 1.6086906989568179327 -1.9840007952655485113 0.52197190764788370387 4.4270942441311156301 3.1534415254821519525 -0.22428630025343793131 -0.91853840056416624638 4.6838107025203061795 0.39738090840219125743 -0.053723321872831630852 -0.089607273718578600574 -2.7587811233638657171 -2.2180619558750733766 -5.1825385406213024808 1.1639109859411065351 5.6593372819678053887 25.379860905144070671 -1.1323582100416016427;-1.6357061007674391639 2.8163540595654565912 0.61571631007383631484 4.5851679873253248232 -0.86440503901130871789 0.09548134754526489687 0.086244562663811960568 0.37600753322525093392 -0.28970181855766580226 -1.1788030319185887063 -0.23079114107011564072 -0.84621850868148240643 -0.054329255025733193263 -0.047857336012165144212 0.086853055728198055863 -0.0082137988925982431676 -0.2151113049528282084 1.2200192579744009169 -1.2387640221769196103 -1.2751641779434861323 -1.1533374673521175069 -0.35211397761763507663 -6.1063742024332299252 2.605372792119734271 0.25230844693925075228;-11.660235141352632837 -1.3932569800558802786 10.835972443851764169 11.913888563691353184 -5.7548218756133815432 -14.491945158262842241 1.8935762925609334051 -0.24719474245863407025 -6.5828912505311336645 -0.42132787505265834049 2.2231480296649412942 7.2039678322899645124 -2.9511136152254442422 -5.0424153406871035088 4.004217321948769559 -2.088580093743626076 -3.2159653643035519544 -3.3364049105061766021 1.609792444402249334 -15.378920690952732286 -4.1837416342373332512 2.2762244680170269362 33.549504406027551795 -5.0444809290717866901 1.1566434249773918097;2.9476769182687685422 -2.0674986875409784304 0.54192022576464549122 7.7208997107024917028 -3.8548069761131311672 -4.5411695039214974301 0.61423135924642213723 1.5788626835511121449 0.85420776807421439081 -0.0068472695462666521501 0.10547664785551354172 0.1009759534240215334 -0.23153676310564266783 -0.43306122142605313341 -1.5628365813357281144 0.42063662436236332054 0.060954946733840498219 -1.2491794842280885547 -0.28438127564438930861 2.4533525921063046304 -2.0055065044134821051 -0.81614717481725718162 9.8779163812910777409 4.3143535081327923919 0.69679197204125553622;0.10878470307475852918 1.60653539971159387 -0.54071981082761699788 5.5412351408403823783 -3.7039772980292466897 -1.0210738977438749142 0.46205983342545375869 -3.794627727793216998 2.0695266553785085506 0.33328384881112155247 0.43092647023321295485 -0.5929874656701363822 0.36560494491670858519 -0.13501284427926030229 -5.0495183672360770544 -0.32567550897757363559 0.05695448819893887682 0.99172918736887316626 -0.47055751570411097351 -0.67772348898782031057 -0.43576457805312052152 0.41863800571391895922 -6.5844849481173577388 -1.2528039018034791141 -0.13167679895313633143;0.035608642075444028874 -0.92189352561749315917 -1.1781977820095075327 5.4531338321997608887 -1.7492742648107844072 0.31183183619473786585 -1.2823784711356340793 -4.0364027058532041536 1.9440414304497497522 -3.0820941119363132188 -2.3045747261286386731 -2.3147825792477543239 0.69422485045138138116 -0.37132415263043161424 -1.0536356412243246261 0.61526883375266028242 0.78299452185739271126 -0.47176852470731689415 -2.2588522811779583321 1.9341571903127012799 0.3420214962790308344 -1.3062257168187749734 -3.8317092778365355876 0.62309445811133457038 -0.4139675055496123135];

% Layer 3
b3 = 0.12723438321306654308;
LW3_2 = [-2.4043557519179579174 -0.079152756627465972739 -0.12683128750143010666 0.00752347428878032045 -0.0064296897429182145109 -0.11419597422884217308 0.26885696112441698213 0.003582334777299332057 0.14899876569364101098 -0.06237521549350957123 -0.16583381810784139643 0.0065014712769472281576 -0.24945499004124258025 -2.4969525316470981835 -0.17803651912606513474];

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
