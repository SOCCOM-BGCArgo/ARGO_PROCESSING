function [Y,Xf,Af] = ESPER_pH_2_Other_3(X,~,~)
%ESPER_PH_2_OTHER_3 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:21.
% 
% [Y] = ESPER_pH_2_Other_3(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999540494077341;-75.751;0;29.715;-2.05132667284751;-0.9;-0.1386];
x1_step1.gain = [1.0000770825255;1.00022983623563;0.0140762793578401;0.000302617642608564;0.251952632905014;0.0588888685539285;0.0422119037568594;0.00826143816767911];
x1_step1.ymin = -1;

% Layer 1
b1 = [-2.3676205261211098474;-3.4645812033587275458;-2.9411638310243599648;0.89504235196530479435;0.034014226234523207981;1.8302699546323368818;-1.5819266120819248655;0.73267597803644624044;-0.68796684629393078936;1.5522107946844783122;-0.017885812408463826195;0.34639051324082459349;-0.47310110834579000016;-0.023412974213867027845;0.70438824084924822611;0.16411472659211354119;0.52459685307388426168;-0.56127158109206332082;-2.1062508997342277439;2.2591831629862633157;1.6049910867003629455;-1.1155637196781589271;1.064583095867367657;-1.3722219030492903702;-2.8947669286678343603];
IW1_1 = [0.32194493551296099465 -0.36236319126000310353 -1.272584524652735638 -0.3357688709820488282 0.90329548808413795413 -0.10171463587972898179 -1.209251430519521664 -0.0002108131022470118715;1.0015083870423151513 2.2271317118383300304 -1.0059523839044022253 -0.14401134611131566343 1.4984649057538379502 -0.043910793180527048218 0.86733712537769747097 -0.13344713551518633254;0.7637905448075359427 1.2363496113849559155 1.5759053107685343509 -0.013941508157841430923 -0.47637451857548124279 0.21341233682299753927 0.58660324374524053237 -1.1065230564800316859;-0.80431508379674088172 -0.29735158761927626969 0.98159641100737071717 -0.47558891745120213113 -0.12778537183596011539 0.38470929512830465136 0.14342546711014550476 0.06651947914302630549;-0.128249367365672301 0.055433294270881700028 0.35861477246189921164 0.42760903020890950188 1.7531819465191991192 -1.2565020939134277089 0.61004387204340937423 -0.34863850816599667448;0.3369146097149722241 0.057042886046117013965 -0.15144257616677361811 -1.0660035759842461101 0.77618174221429869775 -0.10822610886908443084 1.7121004850852952206 0.54760799211554411769;0.77361254230174802959 0.61093144571208379823 -0.27303295530162752236 -0.23993446006693072703 -0.39111425052241782785 -1.2188721404558566785 -0.86648984826563724138 -0.99655322263045409414;-0.36194137324175307269 -0.19962036220805448594 3.0888897415464486862 0.27215201461665256 0.48923464799635218991 -0.66399803611412200333 0.18484077755805075616 1.0217649216853288419;-0.34865468790925835396 -0.4547199216883155426 1.7695306062381668255 -0.26139623544133688826 0.75055058464250379302 0.30863630604890252007 -0.17712493241294052781 -0.52971947977838740318;-2.1953527949080862136 0.6281828430395406393 1.6543581650732772648 -0.68731190026284072481 0.16867850522505176136 0.90077263583883326881 -0.082303086698899907692 0.21383817595110135579;0.47298583422208773097 0.046294386765509037973 -2.0619964388577312064 -0.5438650221552778774 -0.68310026710962690188 1.3197304405388148307 -0.20777236246945107734 -0.79597481495111999106;0.38923986942197125183 0.91807921276312076664 -0.88566704327577505662 -1.0084193512849544128 -0.34830946768951004389 -1.2888229836644571336 -0.52502696945366100767 0.87026279427863983962;-0.10375699966024053145 0.20057773390318700546 1.4170925555353743874 -0.26855012478919387808 0.42338767771636781623 0.49450458776095623792 -0.45250877298213326716 0.40871943539490096509;0.27472774116255899157 0.33843133470245279071 -0.53180932827052140865 -0.56868590980962019366 0.96789202316249967062 -0.66367947741142097051 0.3225060800653536397 0.34249109410400746523;-0.033553699740202475033 0.87391678327552968497 1.3924794416215802162 -0.22631588367522878502 -0.82515610832716335121 1.3068014467387727429 0.25510999244287224652 0.20706194777414710373;-0.27897920261568481815 -1.147404549716762201 -1.4038380787904050084 0.022504859680688188428 1.5765665252595812973 0.22628012229282687406 -0.00014628257753197620138 0.83593791744557077728;-0.12760155119417149194 0.034404442812769421567 0.52091495496303341639 -0.29548793070951140383 -1.188169724110299974 1.2122203469147203592 0.059122706879572224647 0.28250977075884564638;-0.29450781416997118267 0.043860007456078403676 0.35846478833915901285 -0.24134979471295428533 0.62674127466640106032 -0.77186068584815903382 0.18676854199753542551 -0.32687705072310524468;-0.22998977279439514199 0.039747716293656150388 0.87457743086813422018 0.35169944614405301841 -1.3581000226627502059 -1.7250472560607639227 0.03010954440522206968 0.098608705109234279229;1.0900820549559546713 -0.94556848875631116513 0.082733081707627401102 1.0559529877450093416 0.2989133287445291387 -1.5719082349191040304 0.085941198402255530864 1.9291198210853346406;0.46796100992600092594 -0.37423822570554415856 -0.035286419795268558763 -0.57444041568943959852 -1.0012024373624461848 1.6651873841594426739 0.14044345552037185154 -0.2155130000403411572;-0.8396106016567665975 -1.2329451927036694237 1.4231452956513164221 0.88892952577012374071 -0.31066439384202232388 0.39313359366332772948 -0.20131410515792039218 -0.12554973377075021035;-0.048968144660956156033 0.28617169799619807113 0.89718750332805830006 -0.089013792133200472789 -0.69371838873846125573 1.5169756872903716882 0.3444191729125281598 0.78030826133708164338;0.076849519286572645638 -0.23617133936189471655 -0.11521733718296076066 -0.45559742841009914782 1.0351829823577423095 0.036191680292997886748 0.9993234186767826488 -1.1818440127710048149;-0.49552150309843845388 0.79503409060186680524 1.852459384122849162 0.26537835480830240797 -0.72456823474908149674 -1.7937942300028828502 -0.48148339621191799109 1.2384465025606574695];

% Layer 2
b2 = [-1.157077921998051373;-0.92025678354970674722;-1.4209847528549977014;1.2238796555291142631;-0.75512959525219880419;0.66757472789416494319;0.9572126890899795848;-0.57979311962643875766;1.222595250815096346;-0.059992603674453193197;1.2240998858802700777;-1.5736179875303215603;-0.61583905351037582232;1.1413700686118934957;0.97224224446453377801];
LW2_1 = [0.32419467891847414398 -0.15073139414363051491 -0.66254051859378881861 -0.36680592650033705393 -0.75130796882234340828 0.15665887616642784619 -0.25399311384504846822 -0.52472414721920046521 0.1544952286700946309 0.033894576676547910621 -0.62848449641258086196 1.8797453136048605238 1.0631913862961583916 -1.908731761833848406 -1.8856136250203499038 -0.77048687180823483001 -2.7025309405113446815 0.55655226207240371394 0.86084737883294537575 0.17371520203122986503 0.20071928376568018426 0.39459625626790650887 1.9238148377139652201 0.2219923510148702761 -1.2531417877165911268;-0.70387778119701993607 0.11812174723931859188 -0.81458703546701782816 -0.99792262974854428759 -0.33539402255704126032 -1.9699654017976024889 -0.94506495484270136576 -0.084510506541080998577 -0.00038612938408479676594 -0.27296568001246579449 -0.13118986694076681387 -0.66019227749256970128 0.87714064130898772298 1.8776442154036696586 -0.095429987966310220182 0.24415162607167867348 -0.5532672588817646453 0.82290828647595126721 1.6013166476552991302 -0.022352083274268536722 -0.12086259060527457243 -1.0167285640596428475 0.051701938575819637867 0.28061025449219229033 -0.8012065758220282552;0.78188909709134779558 -0.53867329046539713655 -0.64380074551394039784 -1.1089632609254209683 0.690209680428884087 0.27681620543185136984 0.17012072081270016022 0.80157435743933169725 0.57277490439070422035 0.39167042767208293252 0.77665689955079497686 -0.52254031979323045665 -0.30645212372080393726 -0.68125707273590774271 -0.06472784472815255874 -1.157513103821595335 1.2207268444701109278 -0.81862319809453154562 -1.0538732118965876694 -0.22379338733921835813 -0.83348581868801885353 0.075151072339409194512 -0.99411853498102176196 0.37700466882226807419 0.52817501982019210427;-0.74070430924022545582 -0.50378148681256140762 -0.16958274825519439877 -0.14275378771906760411 0.47317072759374168012 0.50374006270830684606 0.64704300530148517456 0.018683150603843519388 -0.088119919373322766476 0.69614353239035198229 -1.043295403447354941 -1.2694533783949750738 -1.8261729306279674834 1.1657473629925070302 -0.75107467547005257735 0.11622382338597556095 0.84392589812940588256 0.93848180330236308233 2.2436135235397958176 -0.58090246527902400508 -0.50289565185899665067 -0.67869671844754286827 1.7358269745944814311 0.020489793870987718649 -1.0575194538540659828;-0.69879883597222614 0.50489568539243157019 -0.085207220748257234133 0.97350294343132259822 -0.37073031849174231667 0.44583838216328242066 0.1628155553427285962 0.1144123698039409831 0.71825455791278314255 -0.49283407036187643069 0.072350047618926757731 0.3240609150130447369 -0.25627310627260191378 -0.29000736232934282155 1.6371335232258870107 0.87725816881307316919 -0.81120948715762020154 -1.4397507846710020374 0.14280223379315562515 0.45014729405378040017 -0.074062560165173693671 -0.31567530824713174731 -1.2407701836534070772 -0.24247969894727974816 0.40446937209261657475;0.64837350612367727898 0.14172796964196932024 1.285569919408278583 0.10133538844867939765 -0.25041682444361362814 -0.23938358134767864782 0.11506509582687275695 -0.12620666183285636497 -0.061784356613237660649 0.10132633299702241414 0.26783560993244154291 -0.034653915396293519247 -0.42390189695150354288 -0.51092600845827285116 -0.38621538070339811455 -0.44482376587751193409 -0.41636815393842735444 0.59924416611668962496 -1.1864420590720592408 -0.00044647747717661823413 -0.067797520684584416939 0.20198744299385898437 -0.50848807264432904773 -0.028312667674109632393 -0.0084747180531074238485;0.36734201888467171715 -0.095248713507070842388 0.25402857070727880018 0.023186126258483829438 0.07946715967160675187 0.13558084088260213829 0.30649900960708392761 -0.036919908045983981792 -0.11612353139520797951 0.14062978179872509532 -0.088694519243499636563 -0.1653485073641051406 -0.26183733439727624459 -0.24508793188323207968 0.12702127343813784943 -0.17989851650485488843 0.61523023061418502522 -0.2525396515324906388 -0.2866437543202796645 0.0013749471498817529037 0.26745146155238685992 0.20408724141743159186 -0.059124736026483201701 -0.041920264452954811296 -0.022763578406195410231;0.081752578934966166346 0.29531625899342789365 0.37074250740594200115 -0.19709469690610409431 -0.67656456197153591603 -0.25719663411823201482 -0.74909378333664755623 0.56413527084788928079 0.61526691973779812361 -0.083146601099517741806 0.83686971557177181058 -0.23019989676648192067 0.7399520629473298472 0.7981792687292985633 0.31589647091717770744 -0.36863614735216199136 -2.159442357749521868 -0.70372762198599214667 -0.17827432394514267089 0.081828547328552561635 0.35453186990185203697 -0.48193727838274041941 0.3771753236348859839 -0.31183317793457143718 -0.9381702580717331541;0.046125653021122585418 0.30813612805090323521 -0.19458356973996684958 -0.37091790457372886269 -0.77778932489570395958 -0.55804375028167407802 0.21034335418285965424 0.28768211689913569806 0.6144919042042729096 -0.15830617299260832387 0.18173821063233278061 0.4003532225355663221 -0.15129467889522199564 0.18462007938726587852 1.0457097393201357782 -0.087488253546673847016 -0.71089509985721033303 0.19328885794813074894 0.065519356322388852543 -0.18900468131357237023 0.072513527437958982813 -0.11374161437664762064 -0.95745007301992224313 0.31426700210546909542 -0.40647408403184226433;-0.1054601452832546471 -0.37803837904626369504 -0.32103997470213746501 -0.29946784458020719866 -0.37035790981803412558 0.40236462701698427757 0.56252590401802926579 -1.4150407950097829701 -0.92654170033984128985 -0.17405411251292696173 -1.4795265647348132543 -0.35119159660285442781 -0.16380294386981988919 -1.379965092008866856 0.074777979697217333843 0.68924528540719487424 0.32867297416853069603 1.9951331158124676612 0.97865391161504666773 0.24405123571154285766 0.62301032709200954152 -0.10204499869154208536 0.53960205717550413507 -0.46545394265746792284 -0.039798683778948851231;1.090665708034338266 0.25447644137388852892 0.15114786628472826058 0.011802825251586979843 0.31973875689803038513 0.93680970806705643472 0.73952774700940104857 -0.27352135060749843687 -1.0625478648855766206 -0.056111879756464165536 -0.25546457604612565762 -0.87054141120144445942 0.064106162332088656441 -0.85675373056655756976 0.74316673843884395456 -0.3732051520674435352 2.9994730084076071108 -0.029119517678830156293 -0.041917677374792218503 -0.01071550913283178183 0.019013561065173687631 -0.57173806653394143584 -2.1263198844189754588 -0.069410359882032632295 0.19088154423463291054;-0.85330638395409508234 -0.082554523344142163865 1.0314928056180443505 0.6953377346634348477 -0.68722064002595151244 -0.40389492996378117962 0.88587146597380250768 0.21616363453982601062 0.5207647187389787824 0.26002354855911596898 0.20632319709968655053 -0.52939608097241952489 -1.5155642900441017584 0.3613407086260725376 1.4941251643559561035 0.53110754244440905314 -0.20770204006166290145 -0.57978005109555008634 0.97029110329080525688 0.39275886872354376367 0.69869961878942854483 -0.17832883328186016625 -1.2275206494010453273 0.16284305821672237502 -0.83095217562437906356;-0.63612329048545512844 -0.24253284030859845188 -0.21121442761481210026 0.51166127804731731388 0.42659983986493910102 -0.55155170366215100497 0.79080768934770173662 -0.17941209642040592542 -0.21535239363371772026 0.48565084447639578213 -0.1292764245359609987 0.013380296055815085049 -0.69329124762346616695 0.25354432893533845883 -0.72600457548589247825 -0.26433532938617043673 0.27950190848393274834 0.025568283778812317397 0.50637380878006987839 -0.53444038618620870462 -1.0338630537160848899 -0.27127077133307603285 1.2870965328762249413 0.75472757327007033989 -0.14639358487006837173;-0.66024768333997974068 0.16538916962213717765 -0.53045264865246810615 -0.57634569985062833464 -0.098843099019527785876 -0.91608250168288274384 -0.35393219660604441357 0.50332012730994513117 -0.45181115701431406873 -0.37003372223085440051 0.33248272652029731145 -0.63865406319006445823 0.63710452879710044449 1.0883100640081435628 0.63038958093688912676 0.61203406941753535353 0.94598549374602791673 0.23727528775316522602 1.4062207445001084771 0.15147638961303136473 -0.27666200312928385108 -0.84615120378464503403 -0.027751480587828298435 0.30861896175361669403 -0.69293588537066674871;0.042901228123191130492 -0.14230537201148424353 1.7230112514913151767 -0.56788968416890472923 -0.32201234725126787195 0.13864735467341329778 -1.2202092458911881767 -1.4824803293947796945 -0.45616208377652101014 -0.048385798531364002106 -0.75515997754891028837 0.45672045778777081715 0.10629491782890407303 0.97581498498992802482 0.93004526067516390953 0.91661190656524460518 -0.88062876480285856573 2.0490305387262797154 1.525625951417705739 0.58289083382183759507 1.5136281851623387684 0.19084825999166349364 0.77286688909899159405 -1.5713526080579451882 -1.4867254327803958169];

% Layer 3
b3 = 0.8841159202074982959;
LW3_2 = [-0.78193221632124565446 -1.4218971469326762147 -1.9273121109717712418 -0.70965230872526263717 -1.5100977340963614992 -1.4264032959963486302 -4.7696957841047815663 0.82993630576419163969 1.0955953372242706045 0.83653629698041742646 0.81640044417313295533 1.0647534885064100507 -0.71445838802105343657 -1.0688705228332202601 -1.2430801446806671873];

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
