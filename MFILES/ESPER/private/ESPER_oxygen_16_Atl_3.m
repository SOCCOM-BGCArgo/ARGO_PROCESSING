function [Y,Xf,Af] = ESPER_oxygen_16_Atl_3(X,~,~)
%ESPER_OXYGEN_16_ATL_3 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:49.
% 
% [Y] = ESPER_oxygen_16_Atl_3(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timesteps
%   Each X{1,ts} = 5xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 1xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-0.999999996490807;-0.999999994516886;-40;0;9.1619];
x1_step1.gain = [1.0000000017546;1.00003626070655;0.0153846153846154;0.000309885342423303;0.0656325902039552];
x1_step1.ymin = -1;

% Layer 1
b1 = [-3.1710939554588351008;2.7826786471616045127;1.5634755618358411322;1.8738425848846849142;-4.5114970936838636462;-1.7225764839607311529;0.81470167614081134833;-5.1592365822876464776;0.79353890123796422262;0.95443529016043204827;-1.0514471286394453386;-2.6627419198542123802;0.80510952974076310795;0.16365951589438382241;1.076710818864076602;0.39349940979174435896;1.9926112231640134009;-3.7269906208536025538;2.8784176733415698379;-1.2070572927597238433;1.4434110914017612703;-3.0039245677151709124;-3.1995551726691902061;-2.249753623126164559;-2.7517692540740408269];
IW1_1 = [0.40354083357079106875 -0.41552353812335851568 2.7745138459729759006 -0.38683530647471564556 0.84671075668813644555;-0.8474272839582538408 -0.16585362117144433758 0.8721372055515713706 3.8613806612647452177 1.9354495947652370358;0.32178009011390235061 -0.34457873755556217565 1.7511722898841144413 -1.1948534363823575255 -2.0526308118478433506;-1.1121119719473737675 -0.79240747825691659756 -1.4980830589440599088 -0.90221013223419677107 -0.91105122946864591782;2.5462141455839941706 0.78390582022476640311 2.9259372346839751522 -0.86350953232256966174 1.0211337613027173798;2.5289161680801384158 -0.42769693390533009003 1.6511749462163800573 -0.6917718605677556365 0.80967905057422573378;-1.2266800152004855295 0.95013475973311689593 -1.8604688767367125735 -0.36209519261791411893 0.36435725340960900542;1.2759971995084555818 0.097225448189547064581 -2.265423550009062037 0.93204156164509777138 9.0131513018815763161;0.57966168828046393813 -2.8739299709453036691 -1.1309957286125711384 -0.59264715714630589183 -2.8036100967026476027;0.85724695880285117777 -2.0040930036516040325 0.17865426944017770428 3.5758442599663116113 -0.43932664752649164264;0.88536605516381172709 3.1282394769182881866 -0.30523987071354641287 -0.53121429050785085035 -1.0424775854754479987;0.8741465811657663787 0.59218926478419964265 -0.019239168407994649418 -1.1057931468200439262 4.3587749294468887129;-1.1806157874138165642 -1.8537404446963297033 2.0689876671815983222 1.7563413094751849641 -2.6775009143506549592;0.50870268638975524311 0.9710075733285677968 -2.0736103236309642028 -2.4548187914418053701 -1.2593310634760728206;0.15480884489292001671 -3.2281629553303377023 -4.6713500329795039079 1.6592339305407874939 -1.2813505520147945305;-0.29039230119363512728 0.16319324232899923288 -1.8389217050898953865 0.94657482800243830656 -0.90519785496215143805;0.43365788099752061724 0.053899511693946973445 0.54269475022438884437 3.0312817732883172717 1.4444563993916608169;-0.48568430298402032186 0.70682753602403225468 2.6399633208285897545 -2.1971298245980279695 4.3774947912155894869;0.51800717393172057701 -0.96526250436823890499 -1.273023257766423777 2.2515162304081028921 -2.130004839445371001;-0.9298970948348566079 0.93894755697704679331 -2.1308772061924687868 -1.1842305790508491636 0.025154389001671710602;2.3646968789771563024 0.1588917583660324484 2.9011588785082382813 2.0003380411750599777 1.6668405214064969666;-1.3268291935116023339 -0.17832758692906780507 -1.9520633581961144287 0.5385128607697002856 4.7931950480922829172;-0.84246124101926445071 -0.063870803627201905295 -3.9373814577774526313 -2.9230590967010683379 0.094623971054076272402;-1.5734409195656000602 -0.1350035693255877689 -0.65496059990102617476 0.22464317624822788266 2.3982705211227117736;-1.4036256438839580518 0.59348366227141424734 -1.6615213451116144228 1.7858227626904554786 -0.04756522813514525716];

% Layer 2
b2 = [-1.2316403810623834936;-2.8524876914278065598;2.4005134024712533325;-0.50908394734930240144;1.8668417713993250562;0.12000655266034737645;0.64664635994294872301;0.13444954845362350193;-0.45997626408310260393;0.35101198526207588735;-0.2178012791758044131;-1.3924441552360404817;1.5820356019069814479;-1.4374523023382137499;-3.1126873368219389882];
LW2_1 = [-0.00089108534072966996686 -1.4659559352571474022 -0.53164448094197958028 -1.455983112095676546 0.28034355301454771947 0.22977109429218023351 -0.53177883240136869514 0.33836995238842709766 1.3063658074427568767 1.3270303522186011058 -0.84404580700057574294 0.97713927045909498847 0.0552952697189335457 -0.80982528119896612928 0.8851761502307897711 -1.7453770617487995853 -1.3375695908676261681 -0.26971350298064911133 1.8824615239424100643 1.5225824738067688102 0.12740224136164246538 1.2936809924539254268 -0.41921842377246543165 1.9374196766734652631 0.082949513077657771598;-0.54478494354246798093 1.1167346698375673064 -0.969827511692098021 -2.6148545677189214409 -0.17619260180752621681 0.1852097208983392862 0.15011163071070973696 0.031332331558985682118 0.44814835077942727359 -0.11132954586426822841 -0.69076183162342519317 -0.27724833000216048129 -0.23558253420272426415 -0.12462414442657232527 -0.04605795736927143641 0.14709366793593345313 -0.0014324492925971016311 1.7520514854963602502 -0.7516654938403197761 -0.9887105817787475992 0.1836356646635543266 -0.88796360156011844289 -3.1379654450976905977 0.80889112391693007709 1.794284597936570913;1.0273232298011882868 -2.5454540944456671703 0.30150831521241172606 0.27388119973713170863 -0.55182446532806461725 0.72812379738277321906 2.3148804723330997923 3.0656784420009288716 -0.11108766911884053463 1.5740269802244086783 1.2201731497707724028 0.89270178723170312551 0.84482269589226133277 -0.95385742737291123738 -0.21619972997138650617 -0.54731741436546921253 0.32923215616572354403 0.023175965087951722055 -3.7902715924461185537 -2.5053796846949496491 -1.3985748840294127682 0.6649011984041224288 1.2570121301295738014 0.098536583470218874758 -0.78733237596515648349;0.33963179400036475064 -0.38840792089000170106 -1.1062340360950737583 0.2617710864357110645 -0.928841127606012118 0.69078916456158045989 -1.1611300725021640456 -0.049025789914932151048 -0.91851050412665125044 -0.23744324948243747397 -0.043676394528117806759 -1.3298578100733073892 -1.605341978623329835 0.73909014238110370965 -1.0837074136526845258 -1.7310608917029459697 -1.6192750797650188321 -1.873239543713637989 -0.2752771813694289138 0.74855850494409492679 0.17522157477382627189 0.11758788117297704068 -0.50039363015167315218 0.82707306158881643832 -0.81147990194430374977;1.0564114141727154461 -1.7151919986649264427 -0.7707707313987329778 1.1192097270249790064 -0.050195850980554665433 -0.16482465426077835824 0.68496628875343545939 0.24281284932248739517 0.52391002890657811797 -0.48698770060024348938 0.59022685205605551761 -0.1240605728959191939 0.23390515273791401274 -0.18544960272064189843 -0.56608290748686251881 -1.477817498367494542 0.94646770214932296472 -0.15300189832394431311 -0.048267103877711631243 -1.568579366361576799 0.28061769768418876048 -0.0018809187813179084614 0.96762532670939949142 1.4202644023261514761 -0.44836280150265656053;1.2948826020797836733 -0.40977265854174355075 1.9069894093910570998 -2.0538105139054119164 -1.1214446796653922966 -0.10842947412876496882 -1.3731551559010208852 -1.1575690849132160842 -0.3113150518430441438 -0.094146391052222277995 -1.070381670042896971 -1.612242929401227709 -0.47837076570694275768 -0.11104983377066665606 1.4233696991630919726 0.79811123058207511516 1.3095416985115413411 0.96388250222369398212 -1.89404478117284647 -1.6498169802103306303 -0.37488075588124475379 1.556374700422003432 -0.015825343817284200781 -1.2126897643341700839 0.18718365505065112653;1.9318558877421565523 -0.30821335591394150022 0.52876476957952389668 -0.056061544521825212861 -0.88770021788368203008 -0.32999077417934341572 1.011537824926212803 0.4933120218076237351 -0.18930769928411883551 -0.29363151884477406428 0.22845254846810175708 0.97317666559986781838 0.3508402112079256896 -1.0949667044903910185 -0.19030476600510762286 -0.20227445639175370418 0.029221339322611995098 -0.79818380813678235786 -0.82120529258495478242 0.13046997684480371449 -0.13833477776317312791 -0.00046777742939231592434 0.27920009293289238528 0.021070355070736446895 -0.82587678224866922783;-1.4607733152237463159 -0.54802162106273388797 -1.5400601889236695818 0.037592437586488489099 0.51649608554502302038 0.14981470550007700004 -0.83948265568152702354 -0.57463029166308576823 -0.057933459101519750589 0.15293661009680253882 -0.13147537560845218296 0.92417619657712879455 -0.13377288784492211149 0.67337617675030092546 0.41684225101055893958 -0.72243080333530695913 -0.17885439319697304938 0.068651848233234238905 0.5118574206043984276 -0.63079960597664508359 -0.76845331453156429102 -0.088020341340761004689 -0.36079065724189723907 0.25422579874881612305 0.64778693890492466956;-0.44016276310340640698 0.0016074421117488041626 -0.63587153168520660351 -0.011748770816356533114 0.3541028211669365322 -0.51617738369217847083 -0.55620012194617285406 -0.21520891025912380901 0.030688574675396772945 -0.010228698940489805436 -0.23223265719650801331 0.98595108072480830419 -0.085451572516792681511 0.3767844348270715793 -0.78200284108999063015 -0.23481059898379999518 -1.3343405497060516485 0.21836122165593868649 0.99630365645339358238 0.069847058746359921844 1.1474599007943164786 -1.1074390303953396408 1.0413649237592430552 0.28622160914428651513 0.82282932601014580065;1.1343165284576475216 -2.4205018467140120286 0.7916501655904387702 0.23743895773878817557 -0.73464651976350980878 -0.21260425414881703454 1.9723215616073257461 2.4736519078354364254 -0.16659364731307027552 0.59518498497092642285 0.79548054789418809474 0.85988623938037889172 0.66567762984209888355 -0.81753672909470331831 -0.11648106079331096985 -0.13586097808261937692 1.6212126795662695766 0.76355049702407473688 -1.7342379512760313087 -2.5761642549502914434 -1.2086578550520648001 0.42711453662212151006 1.2053058909968836776 -1.1555179897976513015 0.21157271311299236105;0.52449837069654126065 -1.633880338379337438 -1.1858661086731134837 -0.12262665412328199732 0.7840754818557164052 0.764102102887484258 -0.11194788188881354962 -0.2983192379653524795 1.4309814967583482925 0.021421536326066664691 -0.31078211871987199144 -0.41530537750865997859 -0.029077862483327329007 0.90498675898001867957 -0.75350899807238802541 -3.7204710922135810769 2.0692214241482989046 0.1584013305869019228 -0.2236350670605433133 0.49615663140182503232 -0.74814218156710110552 2.4429974977115342938 0.30466684633470009214 -0.054273951810363887938 0.0034258293725563833398;0.32418707491047554781 0.69707252053411161707 -1.5666708602256589344 0.78664060169342475692 -0.27411599841332212613 1.117852932574084468 1.6063786418734016248 -0.20228659750693953057 1.7568969115195531838 -0.90800732311368237415 0.23683819430810032491 -0.355170414825671521 -0.29248926329752716891 0.80744048934414580643 0.47555158722185869591 1.0188458983838970529 -0.72071965759625100656 0.25531224526875312408 -2.4336142962421445368 -1.4041266220498576089 0.67820346272463538639 -1.7467070336930492491 1.2545385485768125733 1.6836289504752797175 0.38678237176083463122;0.51135157221690474838 -0.42833456532735342259 1.4136105949064312615 0.16409944709821047693 -0.33405488920863096736 0.49219169036405119311 -0.76665125002411915744 -0.070221482635776100167 0.083595655847706812347 0.20702916235754514895 -0.52440088120413319572 -1.2735801203260213033 -0.11547576767883141469 0.88897965834540371155 0.16970293810295153381 0.97642932927512515739 -0.96351573615556684871 0.16828123768012570594 -0.49289414012900040207 2.263791607218689883 0.40112258732835143205 1.2126620540524593039 -0.72890638056148948642 -0.039163754414531737336 -0.30201969978692455587;-1.1909013527258638376 0.13845391691341260576 0.47557293988769822057 1.0253554774893403057 1.1404834221010755613 -0.71716089450588271781 0.0021930748105311714119 -0.73294194437912951745 -0.1695071485756210472 -0.64321048732208174759 0.47956975399279944128 1.6155938980361388424 0.13926759561918414243 -1.0612507062149476589 0.19483115925998428786 1.5647744456835801508 0.68155721986800421597 0.23254157726067825607 1.5357136907571067308 -1.4129463027664228569 0.17472323800901101709 -1.4860901352639304118 0.8113082232461126253 0.74361264796092385598 -0.39935133545885670436;0.21848355906416880301 1.1059023397101777064 2.2899904691795542888 0.46966775185558096162 -0.11405365032161768335 0.44763481104847763214 -0.013051506931304399095 -0.6764542102397178569 -0.70119890999060263681 -0.057825003100178393178 -1.647211930513160949 -1.1236208217320866609 -0.31631319098546600843 -0.30984146979241394826 1.9817485378261823215 1.5653055004072284895 -1.0461377239953699725 0.35703375289550642879 3.225071357418346274 3.4163726842420314433 1.3089519405942102726 1.6109937003300667246 -0.33279393966230119695 0.10794361641019428744 1.456355063216933754];

% Layer 3
b3 = -0.78211550214679503856;
LW3_2 = [-0.73658745230190336706 -0.51741692759918433975 0.67742617979447750098 -0.68280416685001976695 0.88614838506571813603 0.55024129869783644509 1.0568065273906155355 1.219251053694494713 1.5066053936998855534 -1.022407808269363283 -0.64294413770692127041 -0.45045498528012345441 -1.2869995656093069503 -1.3809960361190121159 0.58638158792940553443];

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
