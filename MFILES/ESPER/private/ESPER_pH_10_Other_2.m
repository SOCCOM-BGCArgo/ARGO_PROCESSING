function [Y,Xf,Af] = ESPER_pH_10_Other_2(X,~,~)
%ESPER_PH_10_OTHER_2 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:23.
% 
% [Y] = ESPER_pH_10_Other_2(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999540494077341;-75.751;0;29.715;-0.9;-0.1386];
x1_step1.gain = [1.0000770825255;1.00022983623563;0.0140762793578401;0.000302617642608564;0.251952632905014;0.0422119037568594;0.00826143816767911];
x1_step1.ymin = -1;

% Layer 1
b1 = [1.6543301549792324145;2.0014559324772229765;2.3681103243729020313;-4.3864915898698919605;2.062118420858130019;0.96638843776185334722;0.27675888997893993704;-0.55810780415693694856;0.61041850127310726126;0.57600243045872923098;0.66939803836929223646;-0.85138139644991228661;-1.4731959878218201521;0.060838817947323896307;0.16105311607039013433;-0.89808110829049625501;1.4432034978287895388;3.5057454513270815966;2.7004579364257614849;4.1395206126753301135];
IW1_1 = [-1.8045680847324476748 -0.48343901924560772621 0.33273926065597431156 -0.45595262040154599381 -0.50281537974822509973 0.7864834140156663489 -2.1194413444849748984;-0.77069555866054695326 0.20081537254334508313 -1.0989784183905717541 -0.31901546191283214338 -1.887836844729234187 -0.58549999581166167761 2.0710795930122616326;-1.429173652196081834 -1.6282384454812561803 -0.75541704595949354939 0.51588617158203353696 0.33031844377785413158 -0.47995779856218268034 0.15308067988364512657;0.45661676141318813249 -0.51656842617941800544 1.4697121508533534939 -3.3838443580404220512 2.2588905105848966492 -0.32681041503910734924 -1.1716183463956011401;0.21374166004289815723 -0.016216168093676316259 0.32914162253980921991 0.46308122319248939958 -0.52918114274510974493 0.5610113723211667125 1.6621483813145576658;-0.30002959451931576762 0.70486676599625852635 0.69709818019328173566 -0.29412807567370130579 0.6104287304237767664 -0.24195018101664117482 0.020877962464395459496;1.570644570330438361 -0.88590100839969720958 0.91171630057771435585 -0.69966722267820924408 0.83446544298839431253 0.74684811184159993758 1.1648878081306988097;0.2087738824066595611 -0.23508785866120687769 -0.055356464468896859898 0.36037392044501825872 0.86456106299998658482 0.99687724564609114797 0.26663176341340794062;-0.43784062983707172645 -0.37989157453761812633 1.5885872536268867972 -0.69282478558133142599 -0.79000090592541960532 1.8773394260439850711 -0.63868327439213634644;0.27206443863450974208 -0.46870761072334216069 -0.16626888582574428366 0.63365570457882569588 -0.35628351086274229198 -0.37528537239835879635 1.159645603772767597;-0.40975968291852360448 0.020926919178765081286 -0.67392960058064088624 0.53221625803665528753 -0.76901265319445843982 1.4417558733102173729 0.48970240345370102952;-0.16376316101104052003 0.27457627125691275172 -1.3659973183090681204 -0.48238568175117962911 -0.75014457541017109943 -0.041803095781166726974 -1.4103870027016820465;-1.7876571210745377805 -0.06299183964182616946 -1.0893941747381905572 1.0111597029561125716 0.72820100232969653398 2.1772201763387171169 -0.52329143727037918943;0.29042055257392029866 -0.11003388758951149207 0.39617539234240412016 -0.7296790066621671178 -1.2099401907361824193 1.2986266090593658884 -1.8169411027415818261;-0.24281428865971174669 -0.25391647923305143975 1.7771680251301107401 -0.66559605842131674969 -0.65994888557534869822 0.18502709062455077738 0.64033709947944006391;0.25177773146277698801 -0.41689309976771221633 -0.13245867193237670367 0.0024330314466655322478 3.5904077544005899014 0.18136594631537955591 -0.094995971227502779;0.92760081217842249579 -0.80335814586651388414 -1.4809250578999735293 0.33927455487547342194 -0.10350485195303189834 -0.63422400219741259519 0.022786230196177252594;0.46507261010520328792 -0.22873527927670384008 2.4861707643981887372 0.42152077233538864931 -0.16887906495484197866 2.4495749011368403991 0.90071205294459644186;0.23153498171240799453 0.13505229375740984699 -0.26716950716749238159 0.96452626460649426576 -0.049209627940098477072 0.17560425961712350373 1.6861294667738944142;0.045395088909512541775 -2.0203345913996790806 -2.981883967751368214 0.55915409259515724827 -0.983023071928111003 -0.79401833176810254766 0.61578227372164195952];

% Layer 2
b2 = [2.4272718763792502372;-2.4197505340106917338;-1.9755237870023838465;1.9117109693662663883;1.4081443086135190246;0.33916038450906677015;0.3430035895296407733;-1.403253816304146806;-0.62130770822986547408;0.43596546176501288494;-1.1635253147416819708;-0.58130939256460434361;0.40461863058772351565;0.71542656428884654485;0.084071687956332893932;1.0251031257448341982;0.28904424622223262586;1.5847266765626339957;-2.3282088765886910231;1.3506267021591404909];
LW2_1 = [-0.20667411210602137928 -0.08083018457845107807 0.14092745329296071932 0.5388829155941776694 0.20696849237835540647 -0.37446062529555346199 0.61820832548429982545 1.619584277094701763 -0.29846305849769083629 -0.43920032028049488249 -0.43090369144171675186 1.3078453561921392634 0.067382727451795909612 -0.9220099868835827106 1.5599632268754279885 -1.2995319762575154865 0.65661645990509032167 -0.039916215299119296722 -0.064036278580207550704 -0.38777282566159876342;0.16924077509005142073 -0.7049323543373383627 0.57354690153331844282 1.7299699849820489028 -1.2046628941866341922 -0.16927719674351643375 -0.93014270962335154902 1.3681419103245553437 0.62398084117857477437 0.94035410777362393731 1.0378735453849943227 0.37474057699331370097 -1.1830477672452752724 0.17746964047379060192 0.50490281754216226329 -1.458263070933133676 -0.22998495627462486879 -0.87310770128121761857 4.3468931009198623627 0.21972829391374859576;-0.11266596404771071771 0.66409413379885895434 0.46169545512495252515 1.2325602005432740871 -1.4903435342293200883 0.26169986899825609106 -0.079804119721501462492 0.72105503442008822557 -0.53518414523702007024 -0.095170900418815485011 1.3894101250939057124 1.2110695000879794225 -0.29560350863358086926 -0.34083049397704023864 1.7710130389071363499 0.70255410632366555035 0.010962512468145380745 1.3129484322498723436 1.3683322732721396875 -0.39550894364282168647;-0.29105439210858413679 -0.26081959086324252306 -0.39231947982419046461 -0.032847466005517501897 -1.2906055181939746213 -1.1286846073045277983 0.18172016997774625979 0.83137723153478393368 0.23298076224868030648 0.089306996503271424137 -0.67967498212403543167 0.64106217463142423085 0.30139547074638356472 0.12806929580364759258 0.48536854771245258044 -1.6420514929713121433 0.55711068324406620711 -0.35154636541893019563 1.1545684822389783442 0.3868445861463361557;-0.042577002500650495964 0.61175335413811038343 0.6070826969126547068 -0.069973995213787651415 1.572461033294602295 -0.57621115810951706315 0.65758162396214403511 -2.0224009025720892652 -0.42471055218727205016 -0.014944423327447761696 0.14697976784098160108 -1.4039811266941022438 -0.10686517018391705758 -0.52435559806558773577 -1.7334182901340688954 0.25069016278020878064 0.54698131308548081364 -0.51151124341528764194 -0.38125988261363519261 -0.056927235719188687413;-0.45505893418079729207 0.173119871878038456 0.063011306930255889247 0.20796218303711058728 2.8696151577097106511 1.5235537602677433622 0.29634674323760590298 0.51571450502024163143 -0.25748802544195142916 0.80406636920242857691 -0.62531333711569392353 0.49963320002016192722 -0.16225853780759363909 0.093055984004850628155 0.6439263694205924482 -0.11681952261465053389 -0.84447940873333415812 0.89364842948674727552 -3.3868358246235019671 -0.046863998577234819432;0.38085889410921280884 0.044282915551799480192 -0.45666796973651763158 0.59683726607099396588 -0.13101667746198686348 -0.28278987368236863098 -0.025723223708596486964 0.3552095891476440559 0.14451054331886362125 0.75448228094762803675 -0.72145598652950893026 -0.082179900542706282862 0.10148033966810729445 0.58019376167332592509 -0.23788707012437959598 -1.3994304834212032418 0.13515638350990219485 -0.23897549350109467681 0.50714186581593934111 0.31455807813592118194;0.45767363803788446708 0.318163463791348744 0.22786386177338102677 0.43082118737233515704 0.42612553070229630681 0.35662225001476149933 0.060286144722061146828 0.52839864019516868954 0.15324667568006783425 0.34693521244153058269 0.05463417113743553738 0.095477387060743101888 -0.23592723624673755944 -0.098521315551387891873 -0.15139052242648165514 0.33169664724332753591 -0.3602684538270680692 0.23418726253299776419 0.0016612844927788752818 -0.027287927397410982611;0.093194508260717859316 0.13388590297017258224 0.265083813186287931 -0.57580043453087426109 0.22883431712967280958 0.57456261526862839606 0.070870083402044092846 -0.22594240517443020533 0.10125087920320112456 -0.42839127773496482599 0.29112874676632183313 -0.14230847411555161419 0.062452361728087560366 -0.5497749844248474016 -0.17117227379228117101 1.3220547675262430687 -0.32137552183301332231 0.22243854466561266414 -0.56165819590212595624 -0.13288670024058429897;-0.58730745306103071535 0.41903397956152077519 -0.0085275033421911838882 0.97866579926207231122 1.0061293650479197215 -1.5608069798630481717 0.346856556491132062 -1.5335366140608770458 -0.84349870056953624164 -0.69459765511189530685 0.61240033928208936764 -0.74672599473846679174 0.40064062105803605052 -1.4182143244431326679 -0.27494207053887043113 -1.0866650556181878251 0.79996627406793241466 0.34286976643871780457 1.4900896450733849274 0.02773681464054312748;0.10643055294626405094 0.96489850929471110419 -0.32813745786478148991 0.64902889123507401248 0.33842522470172836657 -0.11516552427188293806 -0.33308579679722133404 0.93697238860745424294 -0.64001985813525141822 -0.14214620656777834928 -0.0608783984835977246 1.3853348618184671004 -0.19864442195137818614 -0.52969416466440066849 1.5103111807365454666 1.0593875399143743987 -0.71846619571265490922 1.1875984014529288313 -0.065008628251266797893 -0.22476016430918746569;-1.3019090039702090333 1.2213829630365860002 -0.06738381133672856349 0.90675490099216726758 0.55219557474568936595 -0.8578580881572245298 0.19808492328018123541 -0.40727113966188638328 0.30403367964332778106 0.025208086048366482035 -0.95754363780275708606 0.34423044939295693867 0.67534517451471609562 0.063533999920616743284 -0.19767485750321184268 -0.17624251366885027026 1.1772399334889789291 0.50421090933753098806 0.83240386859034332101 0.11743953672787460496;-0.67008794094057655411 -0.11595289882022383199 0.20451703263393516963 -0.13580831429843379099 -0.72632807977734648031 0.79510234081953679031 -0.24105989349391535037 -0.83679899291716153265 0.2625236922360194014 0.46574410679578381878 0.45555671095690764005 -0.46895830154911449794 0.19521477028453945723 0.78677770984375272612 -0.8459737785803048693 0.19970084857830261438 0.24796368691165179832 0.006477655409742613668 0.22406968581928393691 -0.03725532778012021623;1.1877513335919454018 0.072840147391951021172 -0.88331326992021841349 0.36348360635161186538 1.5729945407707901506 -2.5715218310387042422 -0.015787114423237387506 1.7189303422362076645 -1.0158759699999122539 -2.3809833029603675314 -1.2340052665241056928 -0.83995108141543217695 -0.65525665707487623646 -1.5460903686885820729 1.3755804787848804338 -1.3578526721542174549 -0.57026529493470567989 0.069448015623244374606 -1.1246307562582504858 0.62479870304291429406;-0.08345067075066681106 -0.45468912262937599689 0.085112587802456812769 -0.68738154783190341668 -0.64616381878646766967 -0.095595007202323903361 -0.16480656725895861947 -0.53992230708157273877 -0.20694858097782284356 -0.50787866752535626436 -0.058678186102261442603 -0.42399112081587791989 0.32129245251331156608 0.22515939551761082749 -0.36961516755257789457 -0.24196392547727907374 0.98243338768373389414 -0.28437400356630276876 -0.14277730624598192199 -0.098552483256808737222;-0.77165912503042466408 -0.015842100035927880047 0.55409969138225956264 -0.36267256542018538878 -0.75355967874098506165 1.4276980312024289965 -0.26189124457129281387 -0.1298465904937566584 0.33936117848372487904 0.38127590121748711205 -0.4698013409952980779 -0.82267805262264648292 0.18583757390943750343 0.84508554333249974455 -1.3368875187723900932 0.20278096041770130031 -0.27276052479747020518 1.0590741195954787024 -0.68853298008362451288 -0.31018107852582216255;0.19583347352700647415 -0.24991101342205901403 0.63570817278987989951 -0.27572817640413954665 1.5177080103357356577 1.8995766132445766061 -1.0064720723732663021 -1.2662325878448124072 -0.53470728255790045402 0.42057084138014672048 1.7126912153268063488 -0.35296411447354064661 0.74363320986906611409 0.68179559167711101253 1.2005177410791960657 0.16568376076601501468 1.5381515876142732857 -0.23185513591028641955 -3.241493501183378978 0.10753700262798134013;0.230889274052110266 0.45140224705118281445 -0.39250578824102699826 0.36173443391378490297 0.452134776278513284 -0.23336952013712050635 0.2027837224982507458 -0.29620297054245720014 0.11935845133226477555 -0.12747264098357852236 -0.75377433646023050962 -0.13328039364120530808 0.35031628848643375784 -0.27082250343979380469 0.38330001404587676905 -0.67838440954725398591 0.17531757882542919891 0.020356578566053978169 -0.41928403351281573475 -0.3190284242407750992;-1.7517090921295148664 0.29079862703009912162 1.3155824631229970922 1.7866676348085150927 -0.053343124426741099497 -2.1211472762509537837 -1.6697828093526010385 -0.47916566680765959596 -0.37073213771666896355 0.35654044275769164463 -0.51743521722427454623 0.26257798049104136284 -0.24932981160639253182 0.8379750420436872238 -0.27082816012057409072 -0.43774830671509851188 -0.071427557210164027701 0.44099552581477990953 -1.1981461785832185907 -1.0013639764189974368;-0.20626942333998746215 -0.51802854509262974414 -1.2137142951181627559 -1.603010753405938793 0.41817738179852387281 2.0281308018301498208 1.03242007168271277 -2.0215124760423184291 0.10893489554460071556 -0.41559642183305733898 0.39728128785199962048 -0.23194832713669616409 0.98941064121900790251 0.28958808594681267889 0.50483092434778675806 0.71391969351379980235 0.039548386341861543891 -0.64256299036177177975 1.1043942492154656954 1.0172217605148636821];

% Layer 3
b3 = 1.602871133549290672;
LW3_2 = [-0.48792222042796740133 -0.46573442319426028302 -1.3824615937641582253 1.1834662871602781298 0.62924887550886721144 0.86417476053976838291 -3.3069856090664191406 4.0035452844354582069 -2.8554273466172288565 0.30006684473484040998 2.0550765777786637045 -0.49868138842272335953 3.2664043632040886678 0.53695906887111055639 1.8401550766106220536 -2.5179893557792181191 -0.30990496877708839651 1.6714566412355154057 -2.7056145242020135377 -3.1132082870130850516];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 1.90239167587468;
y1_step1.xoffset = 7.36300175998845;

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
