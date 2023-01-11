function [Y,Xf,Af] = ESPER_silicate_4_Atl_3(X,~,~)
%ESPER_SILICATE_4_ATL_3 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:39.
% 
% [Y] = ESPER_silicate_4_Atl_3(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999994516886;-40;0;9.1619;-1.99149636107278;-0.28];
x1_step1.gain = [1.0000000017546;1.00003626070655;0.0153846153846154;0.000309885342423303;0.0656325902039552;0.0604107502911417;0.0470499670650231];
x1_step1.ymin = -1;

% Layer 1
b1 = [-0.59296927018529210773;1.1343839011106084946;2.0480624827074995764;4.9830973501297357586;1.5650622493676171754;4.4617281940741628787;-0.46350962159910347182;2.2866009314732580471;-6.8450729889913946025;0.67432087325806622236;-4.0984876788873432218;-3.8199724633216884939;-1.0699227452232626945;-3.8091640060330176354;4.3062167034380420816;-0.84764718430002494731;4.8876401784140250939;-0.74549828062758405345;7.6902638762254760252;-0.36809750961381237921;1.6891236542982772662;-1.8799943231327413962;-14.038876295605746236;2.2876357113413541278;6.3282612099459276678];
IW1_1 = [-0.83412528725700962884 0.11488768039936965437 0.66671033591211115166 -1.412581434569289307 1.4213560523663455459 0.01705813323650300084 0.46889794206296198098;0.84676132152218597238 -0.30755969761674673002 -4.4047638595572431441 -0.33860376148683352104 -8.5662707740719135074 0.44395256569345081799 0.34577534282797828258;-0.19054857141338785165 -0.67753714535404374608 -0.66880340736820309555 -1.0199730060534966114 -2.6009149745941466669 -1.7931055214618039084 -0.72058879226449279454;2.696517584402526424 -2.7551058269497112896 3.2004411063526863046 3.50734625533344202 -1.6966056027149445384 0.72605434108535349402 -0.090040587539221533797;-0.88350689843209684327 -0.29164733215316379633 0.68695677562320245357 -0.013397933296544966425 -0.69106103656168427207 0.35091900913179785437 -0.61771418852751203143;-0.63689599132433094564 -0.38055726494068359278 -0.65091065455133190198 -0.76004911922279139169 -9.1621574135767396285 -0.48766177631213747379 0.96418130066986595939;0.10272755640394332088 -0.0076275706682991583837 0.39364648085740250316 -0.73849615816387237377 -0.66895104229490287562 -0.41282476588916633409 -0.77633302694377415154;0.22874081443603017094 1.528776897252993594 0.052212722903814873243 1.3509447656490700673 0.90124239228874769569 1.1229050943133236373 0.041224134552130801967;0.63646172782833299042 0.14550216219303946374 -0.12132498347534459437 0.93329082022599851864 15.126931612907995373 4.6195384665893630682 -0.18765936599389379524;-0.036747458599910604038 0.42539013871414155243 -0.78445275371239908679 -0.15603880296775773928 0.52376902271531411426 1.5726268723357932355 1.6426135684775371448;1.7849948144629836921 -0.56391053099862842313 -0.084188161212294346236 0.20604660559487600335 3.6726732959144721136 -0.73177533523509208369 -0.74791269420866901552;0.17431326766537963424 0.46612010311445589705 -0.29423174022315962128 0.96768756505684350167 6.4856508170243678535 -1.5671423779454412006 1.032466116349378682;-0.4732561433318046018 -0.08664907831235417468 0.40294260260571640453 0.37174557373900646295 0.58353678846363676858 -0.75567307199163191722 -1.7042865878828405979;0.2373516542827518172 0.48905593025018717634 0.60597825024502605107 0.9458822161214928137 3.405805141519345991 -0.81543772704745354218 -0.92336771324097965952;0.7728296162624400667 -0.92054968460072705838 0.084433830402500673062 -0.23263005263285413404 -1.9800243944108584238 4.6332005974552243899 -0.81590603872257116969;-0.187077362930225366 0.35801683342697543377 0.46347841878192924669 0.89012116243664862747 0.65821022623598079981 -1.4309998167888888787 0.2089526602056504534;-0.33616583706466202086 -0.60932619204787841216 0.48458342432586121173 0.3546764184934187436 -6.0992774045944981864 1.6540098036641857071 -1.109927214703909204;1.3150923660489315381 -0.26527150171118196731 0.27073975777243269869 0.49085548753164620628 -0.22810402851311550876 -1.3180084816589412444 -1.2088664615042521433;-0.83093925652729305753 1.1330742709833765414 7.371619664109057446 -1.2302513343482883723 -8.4485817346352760637 2.6765798754030392104 1.6257800400415407527;-0.50465015926529965817 -0.63668016746563593777 0.32598667748448417747 0.82718487312391819621 -6.0805005367797662075 -3.5488615869271349368 -1.1442029165175424676;4.4372018699988542423 -1.4873827540860364849 -8.5892166902633046988 -0.65854387849774731478 -0.61532592195234547816 -6.311063507990924748 -3.5062213686317877936;-2.4004389215984947015 -0.84652011379401403968 0.0090706700294087404951 -0.02968254401659780814 1.0063395315493686155 -1.6687866019292465491 0.2837144010658979032;-10.788803823290651351 -5.5253631467811992195 6.7153864376665293179 0.77983970296596427652 4.9628575684043028104 5.1062343013604829522 1.7481601653862426371;-1.4082394365344499398 -0.94702230441305124131 0.46875578695327896339 -1.4776358328608669712 -0.17322597565807668074 -0.294120752023305565 1.4413027887692386741;-1.6824297032085668846 -0.81322213810336418938 0.82394716080260699265 -1.6015505096785871419 -5.9827331466013431438 -1.0380399192151021914 1.207715102801373952];

% Layer 2
b2 = [-0.58668870147947571603;-3.0048739787457763306;3.553514189182378491;0.33033074600594425352;-0.69637248898675618847;-1.9077825575289590887;2.4988600364191846381;2.5114014684655101561;0.72063906096831487069;-2.0850147392093356657;1.0118097288854572824;-3.4183651010383386826;2.4092543868535236662;12.243828184031887218;-8.395564862088999547];
LW2_1 = [1.3652601648827566994 -1.4517913014317289644 2.2354605590280915983 0.21671668293580925724 -0.43283511920826112096 0.8200013876746832242 -4.3507273615429795655 -0.18153243871688029554 1.484938070734364457 0.80338685156189215686 1.3336356637338029341 -3.6247992896431755483 3.851369853149770428 -6.0593339210457344279 0.37404464907518342187 2.6066257417478149172 -3.2572566382065595292 -0.093267758653429719273 -0.031540637124256767909 -1.180406212949620004 -0.23219859985138258107 0.147327394899013836 -0.48975353803162458011 -4.7611808379211204567 -1.6225881066679059828;-0.57729227053795095781 1.8003802484085862368 -2.9058375282952941809 -0.3241405727693783545 0.84882023114621651949 -0.71604849691891936381 3.1895961517841628918 0.3823597309348804596 -1.7500518212508000548 -1.6006886032827156718 -0.63217772232772528884 3.5578513771524269416 -3.4783701226750225466 4.239740832725222397 2.5755863865600518992 -2.7999672566642876959 3.1398897361693411945 -0.80312546035115761622 0.14002315938020781783 0.8290085387221923785 0.34220021682588297463 -0.35105140611543222962 -1.3710074318565905749 1.4509356389691978606 1.8583886373131406344;-2.8310197305587547767 0.16072253998437924238 -5.871493092732767316 -1.4133630532321581619 -3.0603380986435837485 -2.4968325926477197285 5.5773855035497961197 -2.2240787990012740316 -1.9334475729464977167 -1.660411271111635223 3.7359935807877513447 9.3263037677968654293 -6.4477268588093403778 -7.1342781391218696996 0.19197775892319726698 1.5136740556954209058 7.3944868553481830631 -5.7257710335255875833 -0.85007727286154566304 6.39348310145502019 0.64221084961874197461 0.61657873213101177523 -0.21909236453096264663 -8.3139740407008044798 6.4024415524171365277;1.0616027285015010584 1.4508709342352024141 3.4090195354847825371 -0.18762190278995696802 2.5957162045212989199 0.9649263551144620088 0.012023437173209548801 0.17305596262616790293 -0.36036094286518272067 1.8225837620235845193 0.64920592506310836889 1.064086889596070673 -1.2135800718256306574 1.2203558497728446142 0.23123615401588229479 2.1398610436199372309 0.84931411447471527598 1.6267027105730658132 0.53690906984366748578 1.0523677372546424991 -0.57875892445544319109 -0.12748437908691490072 -0.42398722977016639613 -7.1581224375567602536 3.8473822584045120898;3.4443941891151435009 -0.40078135064154235545 -1.6835755443414848198 -0.035769233969347002944 -0.99602821429146848331 0.72811156797380327532 1.9354529956940274005 0.24537038723886445135 0.010620604603759716336 1.5338263276839900762 0.38792883476209427673 2.2021891698993600883 1.5668577691841110155 -0.97640755987063010402 -1.6163133466257249005 -3.4803288090202579852 2.5439848085315972526 0.16682299717688914287 0.077829616298030201982 -1.1258157757494948736 0.7796009877779627173 -0.063308698484048173127 -0.10333940753304801663 -4.5826289801579651595 1.2316023472080246393;2.7627115484103952348 -0.0060067717052446345516 5.2896672642348558924 4.3370976838313124446 -2.960963574260027098 0.22163307354259553783 8.1657941950177956869 5.4070022853931742901 3.7971348038443069051 5.9707128825532977956 0.10491027740603645069 0.15213795419377432183 1.002980268249938467 3.4103758473699667952 -2.5514255413990283294 -4.1713576833894041584 1.5528006539177301448 0.71358034094620714782 -0.98759836502807285541 -1.4489555469725523551 -0.8506650824489085938 0.23599877555476608815 -1.6270648489339094578 3.5525004083730724602 -0.98098855544066199119;1.6272520843133744251 0.18584059648707029755 -0.027540533255536084334 -0.019075542075915177809 -1.1282782946955780545 0.47441826470903469826 0.11760738183916402688 -2.0888565825400693399 -0.045248572000224862588 2.1962694085449192549 -1.5362144574402152308 3.5697494341349096381 -1.2818145863188561417 6.9130940874741471092 -2.1148394432564341372 -2.2009434241920819986 2.5988983476969043096 4.2040370389774750493 -0.61704517291772653653 0.524808637785730836 -4.9466232853777434286 0.068598572389627410217 -2.4123727818735147821 4.2329239040553101603 1.9021992364847246648;2.3625456598060683433 1.3064529107558322529 -0.1121883471490714157 -0.34146005727053097001 0.4773791237287137168 -1.078420551098196345 -2.425404787285087238 0.71557955395891670314 -1.1766294951685831727 -0.39549755533172714639 -0.45296875246481854704 -0.078052528614897073722 -1.5392369152054854453 2.9716198290719413322 -0.92815047614899759409 -1.3173777578685936795 1.6017654218705079661 0.79151476060154635217 -0.2610266922679015944 0.39335799107369912031 0.90918524570804015728 -0.5282618677846091515 -0.34967500644348831207 -1.3205360789015261336 0.25707160284884905677;-3.3631176173627514991 0.37023446373793927888 1.6957634923968802898 0.0035799466417004129418 1.0071605177325924085 -0.78074510370486860378 -1.9966621407130862309 -0.19748900669589877777 -0.024558218017674202094 -1.6073374867143825906 -0.43646802848461896396 -2.2614878232728732321 -1.6040253410411910995 0.99319766170309409237 1.6375384658844485131 3.3669690238421328132 -2.5833948342068726589 -0.12660574066842958385 -0.08659373618651834692 1.1310657118171583146 -0.70223124163587491431 0.099005615769998320252 0.30878358263286875829 4.5125333138209660078 -1.2347953456306559605;-5.5107137189618047302 -0.65610902266421922135 -0.013896232411930914408 -3.8300883012365138924 -2.4744183426932497838 3.2608051503076902478 4.7307506456720771482 -2.9644454776572328747 0.47454736599002245168 -2.5257981574778525236 -3.1963985338128306779 10.968973934696778727 -6.4744979462322405084 -4.8316364442995451967 -3.2052444093250049839 -4.024056573374959811 4.7019146764036516828 3.4092610246504624349 2.6027807316211717392 1.4970360115534548928 -0.12510204995354928803 -1.0943132107101951611 2.5182552015505046406 0.87075675975182675881 -0.49259158202322561593;0.39630088084788311908 -0.38068195695216827668 1.1584318545795457833 0.048117499390989006136 -0.88329129170187015152 0.70676839597415430649 1.4354711152046466527 -0.39731359362290680215 -0.017143044177879061818 0.36149604133441437526 0.53336075299763563429 1.5245180950416690191 -0.90792508613510514159 1.8396120908634823188 0.13613608301257329458 1.1397205499269973394 1.179149117456150142 -0.50473901323744452263 0.19683019960856878994 0.56265332503979348377 -0.094211035744596841046 -0.13086237268995165595 0.00058320116745456951346 -0.5066079630539958023 0.33157174704697572887;2.7964814768194519168 1.9976628561313334309 0.27427735981468870374 4.0707229997094147578 -0.34534382063687069087 -0.78797052941441847462 -1.8215096783661557378 0.85966089652353949724 -1.0209177311463784044 0.1739079711455230326 -0.15644112245409619311 0.83904294509038923611 -1.8340714856711994862 2.266884788114801097 1.8207369497748320253 -1.0916397943780487001 2.7490315567920498019 0.78575205609834997134 -0.36365229228882550716 -0.028222585608941027496 1.2686238906956783801 -0.71645212515133416264 -0.28873638393909450972 -3.5912508081914213953 0.51698321049857143272;0.6318949993959810385 -0.39776340214302580156 0.40549320311041664455 0.34785117399372872393 0.28317856005598240898 1.1070549146611823321 -2.6171134888302787225 0.11513963804411832259 0.5778986618922207219 1.4857244992944078632 0.58006109242995396169 -1.84676116318075767 1.471122490149931572 -1.509670670467071707 -1.1703777949402323078 0.56557838996122922204 -1.2307832577833417442 1.2976063437361302633 0.02980985107597680997 -0.63937951417793059683 -0.12051687366639037002 0.47269456634768552972 0.032068798666303341416 -0.41296152871238578363 0.23205318127019874952;-4.5501164919662677733 -4.2692494208290909796 -6.131617739384737753 1.0672793615700777448 -2.8678501827661655454 7.5616695188962719243 4.2736006980957279566 2.3955212195182244983 -2.1685719780831638204 2.4988218315889350052 -6.8200003807794420396 1.2287681631573099672 2.7019451414426320923 -0.85588137993147828109 -3.5447207144076378071 -7.3335432911603692929 6.0629448477796437089 3.114695001647496575 0.67259940509056626023 1.5902412642928016595 0.56733745307109773037 -3.9180581027242915582 -1.1115233222277796088 1.8627751050368703822 -0.88966159467000716088;-1.0812431466246594791 0.12107699013987655023 6.0200221923909378319 -0.2336744477590971536 2.8657130249204851857 -0.29932305385182228497 0.20494301895699151839 -0.41146911685739678832 -0.032307159526264513094 -0.087261798576139915085 -0.071843609139507297323 -2.6958842885081621965 -0.97070381385661497209 -1.6895013054679628439 1.3748356742215512849 4.3357342027941676577 -1.8468695301386592877 -0.94029159811518725487 0.32347761791525070363 -0.77900154327939474808 -0.032417288636980116923 -0.76189702854799667975 0.0055714946399771959895 -3.27154192136936528 0.92998142450290632866];

% Layer 3
b3 = -0.45099798929236506462;
LW3_2 = [0.2039294553617477801 0.20375279086607625545 1.1326328404199885735 -0.14608487697717045806 1.6196047254199834597 -0.29710454349745984981 0.036229349679998371081 0.15002121289003755855 1.6542350567687076346 -0.26939315517704631731 0.29477363153205649038 -0.07958140794552377717 1.8131823466224139541 -0.08378256733734983086 0.95255538061866040778];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.0145232735458572;
y1_step1.xoffset = -0.52;

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