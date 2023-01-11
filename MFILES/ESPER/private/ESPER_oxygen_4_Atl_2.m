function [Y,Xf,Af] = ESPER_oxygen_4_Atl_2(X,~,~)
%ESPER_OXYGEN_4_ATL_2 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:45.
% 
% [Y] = ESPER_oxygen_4_Atl_2(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999994516886;-40;0;9.1619;-1.99149636107278;-0.52];
x1_step1.gain = [1.0000000017546;1.00003626070655;0.0153846153846154;0.000309885342423303;0.0656325902039552;0.0604107502911417;0.0145232735458572];
x1_step1.ymin = -1;

% Layer 1
b1 = [-5.507391307529474922;1.0267743716371693363;4.5964026700609696974;6.1724159310175155113;-2.8835582942066628753;3.0594619977865868954;-7.0977813295767706592;-0.32203869212825986246;-2.9759874751943002913;-1.9931869574042055948;3.5131313866333129248;-0.7401176148353394435;-3.3010424407206704345;-6.2989929022314123941;-2.7890643024860182209;3.8089051841487702177;0.86055803152869247974;-0.42143263914129375225;-7.2055641874850193673;1.1918582050764510161];
IW1_1 = [0.53919996099656786193 1.5377175234154594197 -0.31319497475333146763 -1.5373730847167614844 -0.22544349534153287906 0.2666740634548552813 -1.5914810922503239699;-0.26390903877113708331 -2.1220547972922907398 2.8503946812864686677 -0.46125894725672067498 -4.3945254013240431235 1.7639594877671094242 -1.0408634638095277669;-1.456368233410640034 0.59733799496399753171 -0.61123734991905287472 1.1818378407309040323 -0.71675466583931901621 0.29940310581846518723 0.51609506693735052174;-0.28689862490004580664 0.57873268662321764033 -0.30158715557881121683 1.1590938071881384808 -0.36262985328979174815 0.74757876911676135467 3.479777904986892878;0.33225719361561023568 -0.22058758004746775727 0.62693118120926882852 0.84844150629681236442 4.3509692959714234561 -1.9974728717236498454 0.61264190721045752586;1.2789908186281788716 -0.30318228056315377206 -6.0926283972862904648 0.53067319473886109815 1.7409344187143396265 -3.4809249091628440276 3.6883084018445599916;1.7261088924275329592 -0.46490644263657232571 -0.1554528012090293787 0.67281306615714286945 3.525901749171149735 0.089724458087170860798 -2.8706701909963334707;1.4374260718117979341 -0.92216983668985519351 -0.61363222581147358614 0.4862599770509242747 -4.6221698895541960539 1.9505670747410432941 -3.8620776805328622672;0.18432091989020296663 0.037638005919414713185 -0.083970664508471279386 0.0015027543499848683688 3.60293133337894389 -1.1698795688543202065 -0.34638429045250507254;-0.03093841265935677598 -1.8059055538688544296 -0.010515849396015160624 -0.029273331797099705343 -0.069425459334337399908 -1.6278244436367019787 -1.1936620527265198177;-0.0075934486032317874762 -0.3195282580781683146 -0.47286595658993241242 0.35981234445704080027 -6.8903198498686544227 -1.0710207756034080262 0.28381300173326740044;2.6681694642093343184 1.0033734781230316813 2.7815643538124721346 0.59761987422390583369 0.096277648201040358367 0.24853928106559519673 -1.7987965105094452234;-0.14723141042495499486 0.18889271697143728002 0.47762215945130903627 1.7000122700723323543 6.1595486935538694695 -1.7412257292244113227 0.99407437033910883706;-1.6833122046968977958 0.84274818830210285192 -3.3629519412266843936 -2.6928900832632760576 4.0943946691160908458 -0.71502693105338510637 -0.43594676476187654846;-1.8832890603948742037 1.3680193626200412016 -3.2672807862466979323 -2.0531485873759738325 0.30535851971321870124 -0.2480133468622205406 -0.47314043430242269661;2.071605433160580656 -3.4964566593919244042 2.2149275340917853327 3.6814068156901722517 -3.5845012761956200897 1.7821814920270817595 1.009946784657479002;0.010759219881924858217 -0.8483458774694382365 -0.84210342524444625756 0.34648627137495685568 -0.20823785455992108173 1.6408655906902607224 0.43007510998767911614;0.13985409157044792483 0.60000289814741247785 0.25601247360547624288 4.4942766095482804545 3.0781074225437756731 -0.99918648307138147047 -3.2924340626240149454;-0.32241298044446714544 0.4857780785069473195 1.8016190954460342688 -0.23522780550339575378 -0.24168692570696709288 -3.5676099302951418402 -1.9173134072100892489;-0.089888890463622556082 0.39027005834482980662 0.52426594334007703413 -0.47696270341509827695 -0.37602801521961237707 1.404295163632772292 -0.65687102372312899945];

% Layer 2
b2 = [-3.9241251901799514279;3.1941624258147789384;-1.9710254343523450071;3.540932168621367726;0.73163819734638113612;-6.1519003922842738774;0.10090882776546886146;0.5035575341540086125;-2.6810708643005840557;0.056430195568663146788;0.73048294665233881595;2.1075703306927899838;-2.5910100210039170854;-1.0816623125759987634;-0.41522221713557266876;0.55098948172236072995;-3.0423032051139506926;-1.8584753239113775969;0.33553232175379305025;-0.1502440574252178751];
LW2_1 = [4.0512060344923739308 -1.1695072966770021505 0.14518531444193735758 -3.8950081675137999682 -1.3753103885698401143 0.91722627046959903652 -6.1993593436864165014 -0.86643899460757534126 2.0989432615599006304 -1.8605196441782874128 -5.4790303685122658806 0.71589958383765961436 -4.5444859646771318396 2.3611731726782969609 -4.171473154063778388 0.43801600694742964803 1.2665934856442100376 0.31704987929329325747 2.0554820041548680365 -5.9677924876094481377;-2.1644951067409499146 0.10432739032260857082 -0.81557479729782533795 0.7451335422916245399 1.5638025491810856327 -0.61652464777804327767 0.41324614362705019976 -0.12933790348873094866 -0.77603075272838184517 -0.40593500087486755756 0.48882669787681398699 -0.50868978974537726323 0.52124240262249832423 -2.3574518428985458129 2.7031795222713288496 0.53063253312834524067 2.611732743167987536 0.32905502604992592497 3.9207061734382739893 1.5803111983252609019;1.3217871963672980229 0.20958805990365000915 -1.6271626002715982029 -2.8745470427988064621 -3.5770798263089638169 -0.52904394934125376082 -3.3136960062240290448 -1.0232401804744555296 4.4991162580365902102 -0.17970755565810550558 -1.1098809041399164688 0.97376190083250946383 -2.7245527624931336064 -0.69197461695501116807 0.69190088594291498936 -0.54658563833335260806 -2.4541826512872795263 -4.2147924482012264491 2.3946630940484014971 0.82386582210506487112;-2.1737893352971071614 1.4443349996709389238 1.9026417748854949075 3.4593422026486582155 2.2314546444176848183 -0.33521264307112308423 1.4738159085311961505 1.7331107636409313066 -2.2785390632098976305 -0.13991884064255291631 -0.80120638655157394581 -0.64179575895904394933 0.8565268581554900118 2.4667960239849318071 -1.3800229394578007458 0.28993768931373342657 -1.0026930610692426082 0.68879455976621772884 3.8649920890129676465 -1.4650146361342752943;0.42170708470976064897 0.10031858969212573141 -0.32970809589869631839 0.85350039870492111493 1.4158856097292153997 -0.18671949275109930855 -0.057834213057832692395 0.66526653772233990125 -2.5920033950444798521 -0.59882128614218443907 -0.845806674751739207 -0.0018886099545622914847 0.67667849185575568161 1.4613359125696130381 -1.3903856260026612457 0.22139556670681190265 0.9901611948931869911 -0.094910154217326597736 -0.51662411937366947701 -0.10416436169429717229;5.8268540151489398582 1.5955778463150731561 0.54221983381729021989 -0.81796695634707738964 4.0246886516926743127 0.32005423540253447046 2.6710302924725826124 -3.801639002922757804 -3.7518684939019055768 -2.0860678401770913482 -0.34962708669241648884 0.13146169533661083983 -0.44526638667020229656 -2.2627934806625980002 1.6181260080161066206 0.78972657384320077423 -5.4698044254137334619 0.10231926294607361638 -4.5264409921190962294 -0.83035365286142415187;0.053146982359088303094 -0.44743161442125622962 0.95327881878645448488 -1.0660610075759204829 0.040609455184775350434 -0.31580601104663874601 0.57726732944966585581 -0.39156776221187084008 -0.56931282155377094245 0.26806182627639402183 -1.8180742270332190813 -0.22255668853717275324 -1.6734381287503206082 1.9096494234129248735 -1.2543897115244007079 0.31604835491049126928 -0.51308198418102712868 -0.06690797185963670568 -0.16490492498663494869 -1.0258425646527185027;-1.227279609134399152 1.6520895509869188977 -3.6306890881150875039 -2.9095371836913175301 -10.65145340371904048 -0.89195932907518538091 7.4778646794228826877 -1.2494153299135435375 5.7172510170270225416 -0.18403021292127874631 -7.4767886264078411429 -0.04669914640846658388 4.9617203653433898936 -0.55811902881107700036 1.7638778454182857391 -0.52290986545102646765 -3.5740616645394021056 -6.1892574705968348425 -2.2007775052055245091 0.58558932811915587191;1.1173405238929168259 -1.4136341442133597557 4.9038518867714611105 -0.78968040799642291727 0.28047887784388014776 0.13337456624874477473 -2.0364287361668074183 1.1534134018711938996 -1.6608570747087609831 0.38142187725468690296 0.18431632383219903604 0.10464748542120662356 -0.44404367962575119755 -1.0407364545052961713 1.3722653837721241477 -1.040340191425851657 0.28574264236755936297 -0.12588726206188044499 1.0548747837143441775 -1.4467633155913366405;0.46834190404387499029 -4.5008233942416433493 -4.8707400666875262374 -2.4031967985234148344 0.62543596778670751402 -2.2897372076036122124 -1.5628417077465213136 -3.0854921015532084105 -2.4634471610405230813 2.0008672990022149385 0.4309728108917911249 1.8290134213678301833 1.733208690148970943 -2.998399079426121272 0.42851556402843438143 -0.84752145670815437484 5.8823380213002494088 -2.3407372814306945941 0.37697735017619349307 4.2436955438860319845;0.55012574605354369606 -0.4508595222593639873 0.18621607807740092611 -2.3733376622695074865 -0.82129372893835539404 0.098177145926531136921 -0.99695022501611318955 0.19533817608172446278 0.34792699733725829248 0.91398517042073224026 -1.8899995791895307651 0.18400019380391452462 0.8024503806428773256 -2.8991478911097581239 2.5750311410913604426 -0.16350504636325050023 0.55366542753607572713 -0.43290432030650627393 0.022394925015382712019 0.37532886206419230346;-3.076266887211661416 0.33826904767187138123 2.4865540276124451502 -3.9909216700885785656 -3.2707399798999303719 -1.0672045594105368238 1.4549982330440085487 0.37876739991173968836 0.77797504757487478955 0.40219013335358422268 7.1616733525038407393 0.57360274670174105083 -0.30986465029602272558 3.7121584664195639292 -3.482182943980358214 0.10334471250339398452 -2.6905220535145830674 0.25464078973686071761 -1.7395598856237144414 1.7432061104044560551;1.5177559262012114694 -0.28298435157738405232 -2.8363340805819903245 -0.88734651015880627867 -1.0758160005874919385 -0.37175538083816095369 -0.35846153187925161587 0.71714898588541942104 -1.2824400973960230488 0.45096336102194228834 -0.70233402217956830693 -0.1792703957346366006 1.8829300793343861642 -2.2067546551878822747 2.3994186488958706427 -0.44618081998305897384 0.48350791483010774341 0.57504489327077223937 -4.0699196912420871897 2.5722996097028403284;0.12641006559208672888 -0.88260583579369034091 4.237479616955410755 -1.9829888940010154119 -2.6902714627196506747 -0.32634300632944845955 1.4048089340112832613 -1.1209858249834430133 -0.068879411653295113949 1.019674793068979568 -1.5423710922771856069 -0.30805174212925356958 0.21922282817778696029 -1.4550957331948166651 2.0056344136703616599 0.08762307907649723171 -0.49086058022086942465 -0.54103801841978638176 0.56578570106834324971 -1.2138929413652197464;0.11519036093947443433 -0.91625210141716273426 -0.71462886837029249154 2.7076736837396131463 1.3313102276393202761 -0.2787556654437104986 0.082095953797755014114 1.054057471000923174 -2.2408374733511311483 -0.90558883946450852154 -1.0166370515515021911 -0.073855399226542745028 -1.247926928417196768 -1.7126702866448488827 0.85387068955255074343 -0.2327875878986274405 0.84110040192878088305 -0.66335441399755867842 0.14420351831411148957 -1.8933390025785550304;-1.277051872044292713 -0.66145744114903781696 0.91023361763435095551 -2.2850689817706122398 0.43196303326393120647 -0.85940178441593650938 1.3667120459343746752 1.3403810243101381694 0.22661872378382122428 0.89669927553770456541 -0.10945503368612079231 0.18883417810140076809 -1.4761210353387756999 1.185924968433537785 -1.3757184510731019778 -0.48265443166560700883 -1.6554252260136033126 -0.51962007918117725236 -0.040282714239367550801 -1.049579525256908008;-1.6743090796753765659 0.29578749924145580907 -0.14577998825005944461 0.20458608122574373223 -0.70911263211017061536 0.18274920126938246079 -1.3380280603362642022 0.42765678627922315558 1.254308341830876028 -0.6062600618110819406 0.75969677943133029352 0.083020420605245115464 1.3118802775571019037 0.25854002273388743127 -0.92546209680529278874 0.028236884580366154729 0.22496637007605366221 -0.069325407968392593894 -0.25878109749619659929 0.57506399601386115616;1.3020471136030731252 0.442959443323871771 0.81441182135594625269 0.36534428561249748801 -0.1015008144997888162 -0.59257625305663463866 -1.2953295368606860816 0.50511550242118186116 2.7196632430313680828 0.052962396932090115287 -0.076993912244589207683 -0.25264936032023804779 0.50248555197456690191 -0.79623071653845811646 0.20472816061309581648 -0.15370025986386026684 -0.98640422999048094077 -0.54152885202353773231 -0.92637340579580573774 -1.1429309858238823061;0.42823475963377488629 0.16084899878345110258 0.42057769075514311874 -0.34729336597473253745 -1.1538231411115504876 0.57297948113366814482 0.052892796577973646976 0.14444712768649667356 0.72402294023390967315 0.10453966476947706976 0.056606811458989984842 -0.034719593012369258378 -0.45258220877455618281 0.57655965972513067985 -1.177708563040867018 -0.69089314715564098535 0.03668477236605023617 -0.1343552502564641371 0.28135362498042615353 -1.728127476314958999;1.1611800009204871209 0.37213817251125352881 0.030516161057882081054 -0.07727866615385717286 -0.8996980320914872431 -0.19493648996490942071 -1.3539869975458371787 0.28986385434505135272 1.4090689065840962257 -0.64056362233263963724 0.28404879314776687593 -0.021789479746526074083 1.1835769723524018904 0.35567169395107128826 -0.92710137010112658729 0.091704424109150112221 -0.10369665405248527923 -0.18827599046872323107 -0.42921811437373613707 -0.091908121536972697818];

% Layer 3
b3 = 0.06058169563095349569;
LW3_2 = [-0.53630523869113855273 0.17581253002477997827 0.23425384991297812376 -0.50656318133208277921 -0.25450829654856738093 -0.25268316183011196863 0.42364496055730122004 0.065283646490014959896 0.1681810016610310321 -0.029885955698675307352 -0.1472137672851540624 0.10677936345899861903 0.29511230398515092554 0.20825510789198428685 -0.15079653530748843293 0.14370948350360468471 1.5210968093316374894 0.38740451497711042572 -0.26871641065606327592 -0.952320159434316027];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.00368391969055075;
y1_step1.xoffset = 13.5;

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
