function [Y,Xf,Af] = ESPER_TA_2_Other_3(X,~,~)
%ESPER_TA_2_OTHER_3 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:00.
% 
% [Y] = ESPER_TA_2_Other_3(X,~,~) takes these arguments:
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
x1_step1.gain = [1.00000001035699;1.00022983623563;0.0140762793578401;0.000302617642608564;0.251952632905014;0.0584431850247228;0.0422119037568594;0.00826143816767911];
x1_step1.ymin = -1;

% Layer 1
b1 = [-0.72018417925208877417;-0.64518316949080800615;-2.5551735217573963155;0.10474346708650207616;-0.058085387184249742676;-0.15371290978548235451;-0.99535958251675737962;-2.0670670830595665457;-2.2565190895809292826;-0.2873014287004502787;0.54487925355514199577;0.28810786091257128971;-1.2448370411179658923;0.67621889428142156753;-0.0062893497777060311527;0.14915164447481552012;-2.4699041931457990096;1.7169796251818039146;0.48736372180134091003;0.90357754320131933312;1.9547315043230613085;-0.37700456906910828625;-1.3447681967713176476;1.9697652503095566257;-1.6911659316182447999];
IW1_1 = [-0.12384954917213543413 -0.33321168991005084559 0.34840504949442341331 0.049453593681784611957 0.95176160005050480795 0.8225967462077665715 0.38108237615704443346 -0.19640806242944586213;-0.64347405606588758964 0.31331319507520655199 -0.69467087480245459652 -0.10492778533755378612 2.7782224361556897563 -0.90374598355590562981 0.24444468954288875451 -1.7452735238964409703;1.3837262103508860545 -1.1421048167621965508 -0.80723795784712037982 0.71547565364784038167 -0.80755929185497798972 -0.95760451126114776166 -1.4598747767617714644 -0.22631429919128034922;0.30982990101711066488 -0.11746565494240973626 0.9593269175227048251 -0.0011304823345527413407 0.050207099456604994958 0.61834893035557891849 0.24056476764018763781 -0.07726832839446923229;0.050070135191033295075 0.15271334641843181679 -1.5092878840392389694 -0.21459695059967504038 -0.77669970810801547412 0.51945853623964788337 0.50539377069240809526 -0.66663310091015293146;-0.015970175972350363069 0.44425696606093662 -2.808762082878531352 0.34265734919082435628 0.3974436847228213221 1.7132655007651986789 0.98167079139116519837 1.0116391615732165388;0.85605252314271751413 -0.37214376507471469413 1.0768904923544777574 0.67274935969006577352 -0.14118014766069486576 -1.0066834576748824315 -1.3873298475431332655 0.92872826373318173587;0.24516997932140524585 0.22374358568197244557 -0.52632076951674078291 -0.17702526771411947193 1.8010086140421177614 -1.6424292111283003148 -0.70288123762381182758 -0.69414515717657010452;-0.057712043789354253165 0.91526248826806777537 -0.87309437776826159805 -0.7204318783296086437 1.861981033324347834 -1.0492501998165704169 0.66669187505631943846 -1.4147125230672790153;0.35603746272963177066 0.70978656218976632353 0.50297690007883544272 0.34795650617861773046 -0.20132684178097653982 -0.25266555244099736521 0.047965570627908904466 -0.472949344434769825;-0.21914874923654809047 -0.14850415700533262098 -0.32031839731593342302 -0.25670392168012845513 0.15181698759914319385 1.0256924381763963261 0.40957888882667481401 0.19473155154594914973;0.51706501417169481449 1.1053001065604579622 0.48332137032655314579 1.1639177686601636097 -0.58237987279378378158 -0.46066100473586052422 -0.26858968235594732343 -0.49763990151148557173;0.16857096152116624044 0.10894185493125299602 -0.12149083173029898153 -0.25775394062663697925 1.0308905739755291897 -2.2072497002054833892 -0.45345981839283555637 -1.0101650942149655865;-0.017430145692276698827 -0.14351145313674448789 -0.73067081230904540323 -0.36700177031694020124 0.32259022944054926985 0.439345619551547939 0.27468410335763043806 0.68682636422818121869;-0.11559745633871931769 -0.45488821897455722354 0.10952180630944016726 -0.21092755875923366315 -1.8593450653135452644 -0.50548198768532714986 -0.36441951378725206956 -0.41880701414667820259;-0.22910332894878565613 0.064897683891309468351 0.21217470163820825735 0.1543427911093229743 0.033095011465927574867 -0.83768780729880132352 -0.6900411451411460062 0.16732791502850782872;0.30692596190186033356 -0.06631922034976490965 -0.35092653475236967875 -0.44623152954725298081 1.4611654178540927074 -1.2350300100446938156 -0.48544773408593711039 -1.3213350239695176302;-0.26076734595117190274 -0.64378014014815765353 0.3891707831688596686 -0.15049168904557050364 0.027526006741854638915 -0.17409943756762444322 0.23573343039215247097 1.1715146915038745501;0.8075793270460120965 -0.1807652330032980259 -0.012380794784783077586 0.61835238176083484607 1.3231498406372137477 1.3703767297476603382 1.1197953817488586203 -0.21259908077395495152;-0.033460103215050072956 0.1624154812499611078 0.35335576297324933215 0.049755252243354412822 0.2037172211495736629 0.04890340297656487617 -0.342881258274918721 -0.6250300601746918705;1.7000251701513711122 -0.43057811136753088022 -0.31390006310198576012 -0.51028291937544623114 0.3507045064541647883 1.0631353573521165057 0.74222528351782846023 -0.22072761547119709724;-0.19295145940622976388 0.064029004977986136993 -0.8744675298329602775 -0.044112604416047473954 -0.01050156941241467759 -0.20978062910169356647 -0.04402131619945426344 -0.29294710087665032372;0.91813439917517980504 -0.50939517738814132919 -0.088577257547499274604 0.65239299398941208441 -0.49185631306500116056 -1.0589315195301578498 -1.4063239175479764764 0.34919389344154166954;1.1977668504233864866 -0.56113658315702941159 0.1561427954557071951 -0.31071190055049446865 -1.921565322662000197 0.050600719917468041376 -0.56972463330515765723 0.32527402790469189808;-0.16203569340657750231 0.21548815997249490262 -0.652736412649094766 0.33388148708324211222 -3.4204911251161189689 -1.9738493296424395318 0.11675902622666059494 1.0273151162274831716];

% Layer 2
b2 = [-4.3546441566464553219;-1.7859878286459696106;-3.0451300262302440913;2.8037031551991349332;5.6359686536420685599;1.3261703184360040808;-18.389392400082400059;-4.6142009762403599638;-3.6083910910601799671;-2.3672303683496838111;2.9414417364657690079;-0.26416592206810862731;2.4439447807274627955;3.2962326362205662811;23.681990096141078084];
LW2_1 = [-8.4504544518097546302 4.3674173739098671732 1.5681858509684032388 3.2135589542817255548 4.1978306614056375423 -0.43482968858552067948 2.1605452968249960755 -5.4687380138364867221 -2.5192967493527929079 3.3903455542770828401 -4.6220206780535093571 -3.2286732637150734604 -5.5591025886800284894 2.7056358876255588441 0.54592900187737170015 3.0987940996369411906 5.0618391287026973657 2.9926634391446116368 -1.7024648596112843535 -7.3929157819656401429 4.069011216204502901 -1.0978266556563804901 -2.5513365875756308654 -3.429990298080210831 1.2141383199931816694;0.56315318181327089775 0.44589779122257672261 1.6482642700104614963 -3.3022222421902434064 1.8769735052695224464 0.6795917436968036407 0.60380270253863976837 -1.0335167390377406615 0.1419199607711613409 -0.43695755900823340934 0.42136708996971977692 0.22357399067164257511 -0.10335407211280049633 0.25743990837803976124 0.82732801768748820415 2.1118688625652448643 1.7122806424568737249 -0.37798487127087238013 0.67754880975496090123 0.83769405514054140482 0.5127542748758280311 -7.400454576955993069 -2.7785766810889773737 -0.65491998898425141107 0.18465221848856022691;-1.7827082005506964535 0.24111104114228792961 -2.1313378238192037806 6.8450076758894811135 -0.42826785028283642642 0.43752174469161364545 0.28844627238379427769 1.546837487032456826 -0.29674868811197518959 0.31485798771978090738 -0.14109206496908408002 0.83808188135838401767 -1.6441045620430971841 0.76460173758742444594 0.84532602247150789854 1.0679144176083128492 -0.70075014954280212542 2.7949789588296565235 0.071642650401925092596 0.16030277448939037055 0.31571875438587770901 7.0805709161658718642 1.7157255470738475012 -0.42579972242193042664 -0.62196174319840424882;-11.16129295457474413 -16.754304812352245335 3.2851302520990945943 -2.4553914590826892983 8.260100746507307079 -5.7072256761128672409 6.8546693701026368828 -17.447207669356984638 13.95205372513101949 2.0243076138993849611 3.1945264568626727808 20.009012942436292803 2.0225919211029705913 3.2476304648551406196 23.195022755257298286 8.330526946272680533 -31.691269819959995147 -15.23183796103897869 -8.8158038781412866314 -13.55507652213804981 16.267093173951515439 6.5417135055123614151 5.7880371459618764618 13.55529950742047518 9.9878548019124817614;0.53040900723595130195 -1.222203578496097176 -3.1798238696054754726 -3.8888243632667549576 -2.3657983305613932146 0.49242784714155873882 -1.2696053626037937523 11.996561693794308923 -1.9807465026109301043 -0.18592747628709607732 -0.44015061380765840848 -1.9604955733502780202 -3.2112717943156834366 -0.33607729278110987092 1.7988166165235741634 -2.4262509366917761611 -2.2681380053988902112 -1.0189678155960235628 -1.7055863856605479878 -0.05375699461154328368 -4.241328349703606726 -3.584644637683988222 4.7523148084391095836 1.9612870900373435834 -0.21674592106640677702;2.9039542785401026848 -5.8382224403185096762 -3.0180876508921148371 -13.953453202090667418 -2.5807057226176306486 3.6656131500457416195 -2.6392720663220297084 2.3846590032887955424 -0.28899037217266809741 4.3061937893946327094 11.665965848388676562 2.525769910266093099 2.4602045214285670482 0.34451714285822876427 3.0040236097118722114 -0.53670790311988336718 -4.1350962449684276478 -0.35082200955886483351 1.6348027994644689365 1.4315728443807944625 -0.41184881603569806208 -4.2864919390823787992 3.2707099493905089638 0.93392752298431358504 -2.3232581690248279571;11.959282082265360714 -6.3564589582982131688 -20.310590180056259157 4.8254328344579660026 5.2420326246652528823 4.7303377411099836181 -18.176495626953727225 4.257850164592045239 -19.616915526766391054 -6.2763524570727380336 4.2306764881460940941 15.596176686977392123 15.273688218225631275 -5.6067280466640001535 3.6376332629664211993 -2.9901658028519411126 5.1889191805026024795 17.580230799463212321 -7.8137935019711193263 -11.267671255101616268 19.61176499868209433 3.398294556028079505 5.7108790648250673527 -17.140987067853455272 3.0184803039627485788;0.27052250187420229244 0.2275001393594053567 -2.9058154476471949046 3.226575109733648361 -0.25572474531636835149 0.48329091728566436759 -0.81454981099218493057 -0.82341978684437200364 -0.29974014371036328264 0.65758709811599602713 1.9603209213024543889 0.88721984405596454604 1.6971780699086560507 2.4796065520661998427 1.2034249329303465803 2.2357680769882173344 0.99349488877578651103 -0.095491470563573294106 0.90265214406196536689 3.2714391045319093365 0.048097278657767521259 3.7431295234269916783 3.4296647599322023225 -0.37558379672424824225 0.38155380882380723895;-10.668535169497818771 0.012134250445534644838 0.76578253623903547176 2.7913934667703124148 -4.4811069746203777697 2.3610785917971606196 -2.8052918533383768285 -0.38794336708238214007 6.1443686316446184392 1.4036809356720834074 -3.7174499926937540906 -3.572177157239979195 -1.6591161304821571143 0.98769674725774292412 -1.7453590707836883489 0.72451555127409383861 0.22300908509156877613 2.0519226012590960906 -1.2281534568610736891 -4.5313315244479985822 4.1034581509098941865 -8.5606119732451677606 -0.2480699174052820033 1.2677477198111231083 2.5101960098847571068;-0.12722381738826446518 0.29278420285187378713 1.6096339668717747351 -4.1512208243037438393 2.276782147271532164 0.74159534249829484054 0.213918254384471207 -0.45124481393893883796 0.17984994898449008227 -1.1698274411253046789 0.82329905430448779846 0.42622102533722194817 -0.58188737024412262588 -1.2426695806460175131 1.101547233224884792 1.1962153064320313334 1.4428947610336855245 -0.061184778978325540799 0.61242741479624107637 1.9897108016534454933 0.76224775528507382472 -8.9386448924990791198 -2.8598318775926347968 -0.71523236655993194244 0.71729500111029387277;0.7132606600742568137 4.7881428641116761113 0.076663084119154970009 3.412193819169341058 -2.6486889622156537527 4.8600643644554413214 -9.7299096741474464523 8.2501779930355176873 -5.126772129105986231 8.5015954530207071826 4.0285303274285233499 1.3443894539921341647 6.7113248193003274267 4.9591146023205530113 5.7263356972003682799 3.1828718380582041725 -4.1087966415258208031 -4.7512127272801070177 -3.0282331496093179979 -7.6573378722726266332 -2.6304255581191022806 -17.015136809770648796 -0.83835129127187779652 1.0642441066410175399 1.6933074390999400816;-5.5163601315639008149 7.5216994118259066227 6.7697890791803869703 -1.2545741054138546122 -5.9142644046661105861 -1.4883730370389762765 1.0901032162980988094 -5.4939127858408616234 0.61933157060982391151 -3.2069742180430327316 1.6227309180198474792 0.77928657780898291652 1.7465486433319836745 3.7477905325840392692 -0.46508877520216024282 -0.12604470076682583457 3.1504198479217664364 0.030670364523209852486 -0.29724424519626230223 -4.3121927118328464701 0.1794076429775574022 -3.8681781715010568057 0.35808380202495659894 -0.17813906135011001175 -2.019584135929460178;-2.9012033864953630236 2.9376620100209378528 2.5153586628806530179 7.5213502873393576209 -0.17705574609801655317 0.49148008912385543479 8.4815762443204665288 2.5072416620705189061 0.78702357917404563992 5.1917286612534683243 -0.81874796943062788657 -0.81377137664533649719 0.47283719504955129542 9.8158891170441542329 8.0280915850027163572 3.9830621511599746398 -1.670584381697796994 -0.8369780937110996577 1.848969257105665065 -1.8159854441877512077 -0.33792409043361831822 -1.2454776664370013073 -5.199591683389589214 -2.9272638286501408977 -1.6008104273144692087;5.8943330735750834037 0.53380830853531902047 -0.64983773127013422233 1.8690014708816815592 1.8412662365437111323 -1.7244802696995238378 -1.8324278990034530423 4.9119812682093932921 -4.3517351336746781243 -6.7024396922536704224 -12.197334563831333298 3.3767271815465718099 1.8264075414268126085 2.868676088300710969 2.1373577978220197338 -11.392998464791235236 -3.5242266835523010649 -1.1076127369843744574 -5.3967523988031045334 -1.2160675642074052494 2.1252078234693327907 4.1827732671061754743 -1.3065845573377685707 -3.2670207425765065423 1.2725731355801042088;16.055644272413811535 2.6200152262619380572 1.2848461036424521087 9.3650642648681063918 -1.4626866046705584434 3.0588813606432703196 0.73221794798357064948 -0.70583625626054047331 4.3818682845692977779 -3.8677828818227246721 -7.0061562683882359082 3.0197565851150547367 -1.5311428828902464705 2.6682705640409141168 3.0724930064487065984 2.5610991119999169108 -1.6835631528286600922 -0.16819055076600911147 -2.9060451448941124752 -4.9843839348794913846 -3.2057574315566723655 1.2496749825562392111 1.0250733364583155538 1.8333410107805496736 -0.52898204681138538685];

% Layer 3
b3 = 0.34055777657035185557;
LW3_2 = [-0.01341208431597792558 1.2418440316997554351 0.19686547442342769432 0.0020030697064918343356 0.012337962221210141411 0.031714622323396834958 0.0075448064677181821483 -0.13228472412456313001 0.0074968409244020907389 -0.95556958229803390381 0.0093276351794912819998 0.01557879813451783918 -0.032732674272462743137 -0.021039638789430877702 0.20653561534808698164];

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
