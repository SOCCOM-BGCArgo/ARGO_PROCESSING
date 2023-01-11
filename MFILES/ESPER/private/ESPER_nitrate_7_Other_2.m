function [Y,Xf,Af] = ESPER_nitrate_7_Other_2(X,~,~)
%ESPER_NITRATE_7_OTHER_2 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:34.
% 
% [Y] = ESPER_nitrate_7_Other_2(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999998793562952;-77.961;0;29.715;-2.18999082656054;-138.904784718693];
x1_step1.gain = [1.00003298043486;1.00000060321889;0.0138606862425759;0.000266418009857466;0.193068829037552;0.0582073296585987;0.0042772329208737];
x1_step1.ymin = -1;

% Layer 1
b1 = [0.93613431349429976347;4.1994367298538355371;2.8004995625692146177;2.9590130051157568225;1.5097573955679590707;-0.076500137366742770872;2.1274388303044875492;-0.13481846853733153724;-0.9270541150094303573;-1.2560624013003562904;-0.19046861349601953983;0.28483123198695969869;-0.61134313279132357355;-0.4747960360303606353;0.47459830074065850747;-0.68560039525561888851;1.0035718130074269272;-8.3830545488985706015;-0.50435256043464216091;-4.1666728904410899048];
IW1_1 = [-0.08520163925729988541 -0.046897960307565760008 1.3130608065552487052 -0.11264396455916343698 -1.6683931722048124247 -1.6620914646384463698 -0.59369051385884674676;-1.3609099963681086098 1.0598360806410553892 -2.1744401740426106961 1.1605444771735387022 4.938390312269932636 2.8358245500421377372 -0.99174518991776527965;0.0051081626868001066624 -0.29117560680119131122 0.65866256762331887309 0.69175413458687917689 0.066405752101199061466 -0.058751869635110136325 2.422153530956050016;-0.36461405518645056967 0.68375190768174187728 -0.37980634917708144815 0.49785026213334721001 -3.8237655626585609348 -0.65340926626267059607 -2.3217606671348143266;-0.91543755319751474175 -0.31260646234893696427 1.775244858233075762 -0.55775950638986881991 6.5184414824362510998 1.1219103098648162131 -1.1549145247228083111;0.46503562031259754184 0.40926630633154131855 -0.35211045086883907551 0.17639571765895711852 -1.08408260703833359 -0.5372726001827057285 0.088706023795110680585;-0.31394187163615877578 0.22011966258237214422 1.3782946442925381181 0.064614571556799671459 6.2370365490452801538 1.1690144005479925493 -0.89207327604316488667;0.36857884077504970977 0.26073086160334207895 0.85775696035861281885 -0.25354110324234474527 -0.021960635358722301047 0.76114054736024283621 -0.067381087346279996964;0.068253600238607839423 -0.61102342646874763243 -0.24384555681493430557 -0.014763236263184294564 -0.068129241122526967844 1.1029716283890051987 -0.53497200438802305911;0.52206059598685494727 0.034393772848673141718 -0.9143850467857806974 0.25706825734891108892 0.47201739218787541752 -1.0193785385887137451 1.2197604156584824775;0.48448987246208102819 0.27135582508864936058 1.3161473263791199084 -0.24903781987875478476 -1.2508047420327352839 0.42719448412084082278 -0.65800018158151252212;-0.22305631718059393265 -0.30619783195688349275 -0.5808942116302630021 0.21691568535742611346 -0.60632552170608744468 -0.73006302486204410762 -0.39043669677597003798;-0.702907540179849466 -0.23891001022488841943 0.73564991695942072436 -0.29655044369841115914 -0.2585299247630719055 -1.0808121355136253161 -0.38704754785624484459;0.11140718306443493646 0.18020013138502014027 0.30634855071975125584 -0.011642047956067940948 -0.30307966586659290664 0.26872661610830067547 0.80624771955736052043;0.40408425727737729005 -0.2114348414794432196 0.10718881597407410589 -0.51216867138556720906 2.3565359742262179665 0.94452637731625177508 -1.2066623833983818326;-0.47038760357155151581 0.057542127115014010674 -0.42695903364463694274 -0.0345561448686896866 0.24201957862273720812 0.37583876685339934198 -0.48319898167372127284;0.45856469680024675162 0.12284248835655929566 -0.15189031427649885009 0.021654835496013012663 -0.91924323816431530165 -0.69037979326209863995 -0.20486842524575915592;-1.3415730642748990853 -4.8052086656860026892 4.9466948045662926248 -2.0931335454021366971 2.6788073557291456339 -1.4184906700595520945 0.25884486588184502098;-0.30448278164652592404 0.046865904443108577504 -0.25901659848339009251 0.010318320739673418618 0.58384284011223419597 -1.1987601787103612239 -0.10756353542925572997;-0.16427623475542996179 -0.90707816401491658809 -4.027311500096423913 -1.1209258771760488571 4.2811091242651304611 1.3461279299131150289 -3.1751043704365429221];

% Layer 2
b2 = [-11.737735018703503442;-3.9865265822634889936;-2.8141823368894183233;0.98759602309283478494;3.0021863769253132048;1.5572157790674530098;1.910385647439519996;0.92315313302788759575;-1.3606913610748780918;-12.04111033695099664;0.45797872880155365838;-0.040984309421195165568;4.4459009368748843016;0.84472233493587189912;1.6901411120388221132;-5.1907720187964478953;-3.0606215549072928184;-6.3751829807529718508;-3.6627631502972488597;-7.4190807647136773895];
LW2_1 = [-0.081266875407358024619 -0.095774403786125672511 1.29587931020278746 -0.38765452870020533638 0.42677761474385628171 -3.6309360628958291528 -0.46876542601634269003 -3.6011752556002045722 -1.8554219499930351134 0.33816711538587040931 1.285054052057077234 -3.7355742452556919275 -1.1275879367036756307 -1.525037937706784863 -1.0839398263926567356 -0.88669011208245662736 10.342430932959585377 0.32636481822976726175 -0.86292720197389516379 0.32324392114274763665;0.76702822617155841378 0.39308744476439499982 0.54916385258505395761 -0.073259635429480371438 -0.32302945907437741191 1.3911319632017267178 -0.068065682618419554029 1.8331975763355208286 0.76552509293175186045 -0.23629060787562491952 -2.3023252100152320487 -1.2133566790693683934 0.51764641115270848548 -0.29616328405646408228 0.56879035267622501948 1.7840274114374901071 1.2111850905783900778 -1.7058644633098245524 -0.59837590807713070262 -0.24516151613573106127;-1.0423423762475751975 -0.83897007555312352878 -0.19964051857185638461 -4.1336533792026841638 -0.25690913467381748525 -1.9391594174519652949 0.02348206940281319377 -8.051273970320485418 -1.0497582851761781964 -0.70905682482011900269 3.555012957387821082 -5.582461622323598327 -1.866076223439284032 -3.2982296084386679702 -2.7516331002371581249 0.035714251267872705597 8.9334860631882140325 0.21806896023394098383 -0.052546245902991228771 1.3848969603630465475;-2.1040326851724571178 0.34358573054036795069 0.14471795711054494249 0.052387444492321667833 -0.13622370675461573097 1.2972122448782910098 0.75916307868988763818 -0.47258844992183152156 1.0386686826064157607 -0.072141094990945453524 -1.1939462303597065951 -1.7094414550399701369 1.1636366156217534673 -2.4187381246726884498 0.10024451056146538308 1.0899900547098282377 2.1267325408679313981 -0.10437339269649256734 -3.1112927654329052274 -0.25955425880602295274;-2.5254540757701420262 0.86080727325165951047 0.53881053090534225536 0.43579501598428727061 0.43108823545098151619 -0.44582438210007369062 0.39902655953885135531 -3.2447349523835682383 2.2133929288377074407 -0.24472435971269865251 1.5177844864763665456 -0.034641385147224021401 -0.48928625428246197426 2.4771054083037125082 0.63722218518769147 0.31594105369225222724 0.8440726266852145443 0.13525089999660067086 1.0273447435029234853 -0.88415951036726003753;-0.023476747856658888858 -0.88782674649208714879 0.82787980566814745043 0.20669930918780021756 -0.88443572747701870718 0.21859090232434202905 0.85731653234626736815 4.5299896305970559496 1.4551397420241163339 -0.32873715990194646297 -1.1144491455032652372 4.8256890200562363447 -0.86816880513577476819 2.4264861062253784141 2.0438770745240608662 -1.331849122780719874 -2.8805046184890961136 0.024369082695987295817 3.3677075931863140745 -0.80517558082985840606;0.49312783588887215691 -0.61257358634116132468 -0.70316841961382758885 0.033713365838648950978 0.85251742858924928203 -1.2422177222554602327 -0.37003541113063670442 -2.774387328035782474 4.0878164816114122004 0.86747510645043002242 1.0327299912886276356 -1.5227892408063985474 -1.6535923625014290472 -1.4852345731725609124 0.41473063205281895938 -0.73621924669878191505 0.23273408979496834625 -1.2815984650658380772 0.70533686325753142121 -1.4810668162407079063;-1.3956405851647386562 0.52987741996258630461 1.0058309491005605896 0.37120198880957100407 -0.1918537406133185852 -1.5735157993788988495 0.4010470944389403769 2.6534216577855715968 -0.28191722035013272452 -0.51103915037845504266 -0.39850099684416340207 3.0182634925425562145 -0.52173570560802118301 2.1103144559555171433 0.078034453586674906722 -0.86353012640219062312 0.13984640738946099714 1.2844107651913503343 1.5798626189044346901 -0.16175509019100006713;0.3154197844604659795 0.12962405628669482205 -0.63372747524407280206 -0.49263906877638719006 -0.61425504784498408117 -0.83078759088970877134 0.34550752381151111914 -2.0697395385981849358 -0.59605008821209792824 -0.69513377407575616118 0.38734479695340323646 -2.4701933078559146217 0.026930164268779539782 -1.6517658571875715179 -0.63160518418984556988 0.7102241015600528895 0.56862244838414710557 -0.56467129937778559601 -2.1996812153730349593 -0.2500316236426488592;-1.2205368697357039665 0.43255285832744932994 2.242856741712992541 -1.1970893300460290209 -0.14769600933635051998 -1.1592193532618650931 -0.051925677569137046463 -8.8946749692304667434 0.066885176261177667345 -0.096330822542348273907 2.7987667412912213827 -4.8846298499245328628 -1.6247042321816163035 1.3861303092523835634 1.4317329809339565116 1.5610205733150572804 2.168949238752284181 -12.132821817328183656 1.6396935133726799894 1.0115811106764285299;0.17593392812668723346 0.50716582822831346444 0.64882659535547382301 -1.4746903361406760169 -0.43477532735274160292 0.73491926359764203625 0.088009670777262627439 -7.9644969034277055542 -3.0220835828717835803 0.3278860438435185487 3.2891293348005108399 -5.8346863796816910153 1.2865444412326174639 -1.4475774071454936553 -0.040763156005918427249 -2.0167729524666140328 -3.3708392550061851978 0.56878822529193018287 -1.0438401025948822465 1.358180898908267098;-0.30632563003568935578 0.58498785690591192932 -0.049961684070497384969 -0.98016211044592682633 -0.073037104253758353134 0.15129846671360847044 -0.095328694116269818704 -1.881696549025096088 -0.47482308520076527136 -0.14903104262357130816 -0.02196050582953550101 -2.6462908656738726165 1.0959057640223826002 -0.75657339789916966399 -0.7962260536752507889 1.3183887796588993346 1.8638059357857528209 0.11599163868542311517 -1.9011676919277584741 0.15222959635979613213;0.11837721782181803709 -0.8519753175454587435 -1.3148473588886937513 -1.914206149774397181 1.1141238772763184528 0.95687071956858049671 0.12349132729747618054 -4.7556563636587450006 -0.73411731118787726391 2.2067608138246908034 2.9393517050383177747 -2.6107189338900811393 1.7270972528299437521 -1.8253445577700417868 -0.15112078588030672388 -0.30241086886076651297 -0.30153951167746767537 0.55007435263891069077 -1.5175274786923509396 -0.1180216960619546035;-0.75250581166499508345 -1.0785960021375131657 0.44755121371546985332 -3.7781371776821273123 1.4089242199000979738 -2.2557745658000665223 0.40961366992483194682 2.578456893492936608 0.89262802346100800222 1.8774667405479741333 0.71959629636417521414 1.5764045758324432001 -3.0152024361284301968 -2.1804304724069321075 0.19760715723294319801 0.028936451065791089776 5.4466042888092944807 0.72053057558903355062 3.8098903147855942564 1.01843733454157781;1.4004191651443109734 -0.43654295478401472508 2.6143917727919689931 -2.2562934906784586353 -0.096902928657861348016 3.7343884302741199299 -0.87783160752021582862 5.933652203572107986 -1.8607075536249157999 -2.1679594154074912637 -1.3084321672743322029 7.2132652434643622641 -0.23290218075831611855 3.1788621246875847781 0.33966220365097032907 1.6419610305913145964 -6.0760260400641605827 -1.6374767763950712762 1.2049575473267297721 1.1194447243262528602;-2.073742361020257885 -0.39581772678707893842 2.1561081813804099383 -2.1446023760637671352 -2.6606817295957361225 6.5800666029395680212 4.0773390909703133289 1.8832854629087545373 2.4575859464750151595 -5.5290654068148912259 -4.1985455379040699242 0.51680401295095124148 6.8838709402136606386 1.7048113613765607521 3.3898386531812532141 -1.5884750090124466482 0.74386523062749321245 0.43052225712008318492 -1.8658986589677668544 0.4273887263994879393;-1.8095094114700569499 -0.072386631512141849654 1.0622662690603736912 2.9699659872385222315 0.28701536008571926217 -0.72459787911192752308 1.1895836485269297711 1.6070605977801792896 1.4383078159782856087 0.75425746491516398873 -1.4106702073306329748 0.77374764827888764529 -0.0361118920429658527 2.7658466760116953154 0.78773326732682102946 -1.7798052990118831129 1.8823443450238146468 -1.1480273505774685194 2.5813451530109121634 0.10118245841933352913;-2.0346896347527771987 1.1472487549834962106 2.5113063649465345861 5.0785891445346251771 0.60909927801478080944 0.30356278380016571417 -1.8862246336609589381 -10.487236936041224311 -0.05623784118017430278 1.4160691158698899272 6.3139777471822799981 -2.7709478660143562756 0.69644488323088149517 1.3719183158490899643 2.9142302830830102423 3.2470070556751351454 -1.1845378309490934843 -0.20581093887571924594 0.23792860855026101707 -0.15891983246519042483;2.4816076033476783813 -0.18365119668992146984 -0.6952303932158953037 -0.61752462659856499627 0.13458329556689074824 -1.608491324945916956 -0.70199818749096787318 0.22267596375104806161 -2.2583046460269597944 0.10682070647096129912 0.95336837654057726699 1.120244456144827927 -0.80057307016900913776 2.3068391299483552714 0.017165685403215794058 -1.0309630385060377655 -0.020399586479035225484 0.21806200826854102437 2.6761725411678138897 -0.36116758640583418449;-0.10799812788328844171 -0.86753445119604311042 0.9016590943251115986 0.75376587004249040458 -0.62112901103967987915 -0.80091353344462945518 -0.02695361285492695394 9.1228016384961367891 -3.876517957328827535 -1.1326864388264528127 -3.4744683073970890064 6.2835292074600639367 -0.42512015139323433566 1.6598084289766144028 0.019755717130807898263 1.0403563463698703995 2.1996503041245238919 -0.64510912607183978906 0.62097970200174723576 -0.059815487951782882192];

% Layer 3
b3 = -0.81981175897331393099;
LW3_2 = [1.0565914949503496523 -1.1219043685761254103 -0.18442986221004631364 -1.7717073515172667708 1.1807371096892975704 0.49830757810909553918 0.3562003660497549129 -1.1175639503787013762 -0.9931400735660471657 -0.28125836919317376283 0.41830987594660179463 1.3889247185881217206 -0.41662736091960866913 -0.15494009514743317668 -0.15875242536084288592 -0.4792296344807081554 0.22154505200433358136 -0.65543966741295567324 -1.8335870559766731525 1.1143540636260671484];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.0415843642790311;
y1_step1.xoffset = -0.9;

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