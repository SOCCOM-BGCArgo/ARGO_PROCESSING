function [Y,Xf,Af] = ESPER_talk_7_Atl_4(X,~,~)
%ESPER_TALK_7_ATL_4 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 12-Jan-2021 17:03:47.
% 
% [Y] = ESPER_talk_7_Atl_4(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999999390765;-40;0;11.29;-1.99149636107278;-143.241523403244];
x1_step1.gain = [1.0000000017546;1.00031699124092;0.0153846153846154;0.000310077519379845;0.0705602618858232;0.0604107502911417;0.00458552991639557];
x1_step1.ymin = -1;

% Layer 1
b1 = [-8.6847590802606174;1.7284916047405616;-10.535389873343387;-6.1370834295335612;-9.9132979114493196;4.772747693367104;13.285956712258244;0.7895135771357199;-0.67780050624037202;1.3363456055204512;-6.3840315127699343;1.3648411713246107;-1.1570079433583924;-11.132138075528005;4.0385743012636786;-6.0516349033933086;7.3179671792691963;0.84863040815654667;2.1051474160006918;-3.1017643182418482;10.053574778334726;-1.4192359191921389;-1.1094207230570521;-5.9185720930858308;1.3857992981599354;-3.0473733829241003;4.0104112752533903;0.73605369002055887;-2.2023262111839386;0.49867861574505989];
IW1_1 = [1.3096370311908851 4.2658918017071832 2.4083474882475917 -0.56945495972364002 10.946965453352886 10.456319164602911 1.812885435786459;-0.72699206789074322 -0.4393151720411882 -0.39317317197583473 0.21253244018296538 -2.1524886669970336 0.68135537516696254 0.16957852075625923;-7.464149599104708 3.5298016971177337 -5.1898176132831368 -0.83311496921743278 8.6042929917920148 -9.3475568997896339 -4.0020504100029424;4.5958635602102502 -6.2717302174322196 4.5040357600451975 1.5695505393668969 1.0575409628978027 -0.69346731111665194 1.1512542023130106;-0.39846710432663518 -4.7171188863577322 -3.2910052908044163 -1.892026716174757 17.647664491803884 1.2666516048197176 -6.4735624072256837;-0.19476725307402393 -0.070904523964993305 -0.11386818630672037 -1.4671001487270761 -4.812941055656232 1.686108623346924 -1.3901808598902645;-14.866025952154528 0.62827518855645115 -0.69364610797210746 -0.26383596692293704 4.1923255686180552 -1.5311983109834684 -0.48495315463854705;-0.17835679231953697 0.22752227844722822 1.1939762013457915 0.15410377129927544 -2.3223123322077739 -0.38099648700023309 -0.33392408807251517;0.50089702980453044 0.0013691104859576272 0.2076172873869559 -1.7126580392997286 4.6697573425899401 3.2050425621077134 0.0082863992937019671;-0.25797952987859318 -0.059182887240243454 -0.40539768380355451 0.015911526739867492 -1.3366862547266822 0.28716980029847816 0.067248365424254872;-5.6828319793260995 4.8032625063298608 4.2427281627226154 -2.0952926526690834 14.062445858017849 -6.0481809937656834 0.46140192333468744;-0.18099846511262987 -0.37449147850815467 2.3630475372065378 0.51266293289195775 -4.4825767877507712 -0.44474760044355288 -0.62156722564886391;1.410985146596645 -0.87036704093029615 0.86133904827382002 1.0119369035323407 -1.2717102920679819 0.3720557757975263 -1.1449750061438453;4.1040884325421336 -1.4515870398520963 -1.0848828316011883 1.340574586813976 13.618831101893225 1.08603794730661 0.69558435225672388;0.89503307729387827 0.81738066572011725 0.26250855171642623 0.85920862420395827 -3.1848141780829864 4.3278483880570855 -0.25101682756576194;-1.8173388812943645 0.18043946880460995 -0.13667281648816912 -1.1221123603887764 8.7752030857702383 -2.2341441950380636 0.077854601936783091;-1.1525638282874404 0.9563833317648458 -2.7227630452515705 -1.6850414456197924 -10.748065823401991 1.5726050174764643 -1.1305945654549832;-0.4730789757647893 5.4026814555611749 3.0358196420440535 -0.94684890914287123 3.1376833406319453 -2.3446517879957649 -3.3735175759333704;-0.76547736613822204 -0.47986020513979527 -0.39866698155131586 0.26971352785358527 -2.6011642243275181 0.79819473109122674 0.19189186788156568;-0.068010931854567289 0.35573019655236321 0.81782711011708675 -0.36786886044670969 5.1902704342907686 -0.98267226676461261 -0.67424188733994606;-1.0893364139529595 2.434600723239011 -4.1196214591161482 -3.0694030894156299 -19.621805282106607 3.7016930249776467 2.1979696257221688;1.2838885212153859 0.078734232962842904 0.1674263395147233 -0.017369639613361451 0.91927151172566934 -0.29962478389183933 0.19109931165240227;-0.016585651022152519 -0.5954746540455631 -0.6506212696732373 0.027384791960837544 0.59696096801781884 0.39516840373920209 0.25139727274926771;-3.4503331426056851 -4.1802960892143695 5.2187001136229609 0.10047368628531356 4.6890188011847398 3.4471935686610138 -0.91358874770115306;-0.2699130413928002 0.66566696549387794 0.50602385930747362 0.058195988946992475 -1.7410990270243321 -0.11091405881568712 -0.20041885966709672;0.31199527154946105 -0.60770774634226288 0.43686261747026445 -0.1228539318650487 2.8385726989554443 -0.55045677765459455 -0.042632359198254303;-0.17391411694707415 0.10133611081587196 -0.079723410764948821 0.98669824406118201 -8.2916785068747245 -1.8394720626168086 1.0542443091632703;-1.0192529390664151 0.38497816298416632 0.71453855320918558 -0.15827969066517925 -3.9106442410924496 -1.3502336229971981 0.81461707452172993;-1.6491286387882509 2.3454907895373553 0.14227492384240512 -1.3214126264469093 7.5310191015681553 4.4605855542362054 3.05218519445866;0.70133470301829481 0.21327870364803783 -0.15022832203447545 -0.10428830759497729 -0.89909506626694946 -0.103413878754868 -0.011057850824234666];

% Layer 2
b2 = [107.38747565594336;152.97048529223531;-6.1047463383540368;1.0235116348634923;15.533140166868085;-10.736753161543724;-10.187707531117436;-59.530215640740927;-1.6046070200233595;-2.5207468899532071];
LW2_1 = [-8.1578643945044789 21.527393838426292 52.847724383926881 -1.5864452876769366 -32.115855648628205 -26.631039049671543 37.007281606036457 -5.499254709008973 1.8886883182744163 46.163203818939458 3.4463538694134597 6.5767083926821783 10.279122226005436 8.0602732750475141 -4.2234424283083056 3.4477143059410262 6.6682998780457723 -6.4720744574993523 -0.83412586305234326 -11.952124903017562 -3.7308477267482059 57.31556746523674 7.3723985326900303 8.7572696445646443 -7.9324705489312706 0.55598939094790512 -15.990165149452425 42.555484987400419 -25.614935765578029 -40.842719872817923;38.804206802685563 -7.8308769071955231 0.234671017540139 -45.813132136496769 -30.798619234132058 183.87411662368737 68.230238394833265 -96.323915294611169 -4.9861829256451307 119.12069659385575 41.042191310013258 29.991209173130883 98.143971449851065 32.831443874562801 -2.6431385595306525 -19.284419686539895 -48.035677057900422 -75.262728504459389 -52.187752095608872 -17.773313447611844 19.601603326248089 49.705234837534391 97.806867391503232 -2.2555622919397509 -52.115279307909908 -59.933168871383096 71.039358564492005 -97.80135421164681 90.627895622042914 -44.03729286551183;-35.183460689661388 -28.746333819786436 -8.62077331675445 -0.015512547508496174 -10.708939151900129 48.539347106338376 -8.7072086388119452 27.021641819966256 16.10972364398463 12.513912757440757 -14.87512005066969 17.000619656007444 -21.992582383025152 -5.7010755018043842 0.73865724516990894 -0.95831814266091608 30.806812968804039 5.3531494271294067 5.5495432758540151 2.5394958717140015 -14.793095782475504 -65.221941626963954 105.1213619369925 20.802514836744702 -20.690434370527068 -39.701851053264896 0.63362071030397149 -24.091708013249487 3.2796220490485846 15.227240087436099;-1.1401353016896016 -1.4284177366277908 -6.6813202945568584 -6.7018059706400779 9.7982932746903035 8.9645918176605051 -7.2081938453071288 3.3478823741860988 -13.892652086174666 -2.403735520782416 -9.9993777455318842 -4.5903212568782843 2.7452056614739511 -11.755458934817447 -5.520630770758463 -2.3193568995560727 -9.2204476953461079 6.9282823688495974 -8.0548121186885133 -11.257216819540638 -5.5274463044740276 -2.7366570309483249 8.5657223284647284 5.0446583365800706 1.0038818494208905 -14.270305479522079 10.829683327156578 3.4863325929484059 -2.7118999239740673 19.556341120114833;-18.181289740788266 9.8714404915757914 18.736256965599829 1.1642379587416549 -1.9112809118705993 15.112731785282751 -3.8127738945278549 9.3692112920423618 -6.2829141637927863 -5.7743034189245037 -1.4719495720116866 -20.836573424994221 10.091490476902896 18.942382259868978 3.9410106714751629 -20.718727394805011 6.7098178034913882 11.238001158476468 19.537158127334781 3.568406157180906 13.63013310854191 -23.770526819918079 -23.220214884757475 -5.2394384489407013 -11.174022411955701 11.865090481461493 -3.1900759026715848 -11.000230428205271 3.7022128213189678 -34.151802845687818;-40.825756899053495 9.5028353287831724 -4.2617456482724423 -4.5011667654705834 -13.847203415999653 -7.3341668976070729 23.923119056717859 -7.3361196188422264 -17.726121611898439 48.51683192280349 -29.324885603252071 -9.2406080672737527 4.8184654692004223 9.036174737822753 -4.4973555923449773 -4.9064608998907513 -19.344259856733231 30.405024753099426 23.95690763606127 -13.256758717451024 -5.106989834089056 12.428112691433579 3.3440254924584805 -16.726228831686452 22.312702610648611 -9.1092324582228041 -2.8269236693108928 0.41764835250566079 -3.427593737545318 3.7227022718595664;0.80426448613381474 -9.8669115901909219 -1.776822909650172 2.5632059315881146 5.6682264545587868 -4.1465669893174661 13.061876627636083 31.402913410138567 18.590854473255039 -10.923993961696581 8.2762527851924208 8.909446077161741 -9.8183002940076758 12.272708044362622 1.6794778662098788 -0.4911125762385426 -6.2438038716752651 5.3500129964909959 -0.13209614840249093 10.271763243669225 -1.5505664148158151 -2.6525695859772589 58.262503893616497 -7.2147448876217215 0.73078030548982098 -29.507808457128196 -2.8615536617065547 5.3547048391549126 12.79369607012443 -0.26733407982364282;-37.228787369801502 -8.2263916548702376 -25.948443311917867 1.3454457936716235 3.1313705674522967 3.2715461068140099 0.84339655476694564 6.7678226270678303 -3.8152010720951082 8.1974481970763531 5.2721187387630426 7.9039772861792956 -9.8812139128666505 -3.1554727371069724 7.5228348079415372 -8.0025094382459461 -1.1212071066398492 0.30217116544671874 10.60751779468869 6.3318974815659255 -0.0079073207493095925 0.038538856711089309 -0.65009071816238817 -1.5601493189472051 -19.249834600244618 -19.845682734583242 -4.869810123805479 7.8615321465394512 1.2675498115088319 -7.552385823992787;0.0064657489325140611 -3.7339718520142879 0.04368408360928925 -0.02093351338009895 -0.075225313920275344 -0.22046813361114051 1.1223027093663769 -0.77155012418664215 0.32273735995856168 -3.1535596127776477 0.035557344543938718 -0.38025427572275033 -0.45726468327238351 0.073348270382514111 0.0078534892787051787 0.040145409998686786 -0.010356088829517083 0.0061456202885908709 2.1424211013895813 0.10631338548713069 0.041409223383143658 -0.99402084795320322 2.9662363680543566 -0.023033576046075546 3.5007499747271043 -0.92782051420154377 0.60070601820617497 -0.94723902870491361 0.023532699518333249 -0.055212507171909259;0.037785645994402629 11.299220344039746 0.10758153038460171 0.0029801824653337156 -0.079894005891395503 -0.27339466724277139 0.9433896162643427 1.6299768612747092 0.26718341293868836 -7.2018097517128776 0.013280690240806347 -0.84180069801308754 -0.39551787330813548 0.11132016929320256 0.044464051793567307 0.19828508333128189 -0.088307052831979857 0.056812898440642308 -9.4647362062302118 -0.61054347158237487 0.010766474361191575 -1.6233294935072677 -0.9517235133711176 -0.053384241725893761 -0.73338899450306028 -3.0789418269746256 0.4098377271261669 -0.82315789243351911 0.048714337587625507 2.182969411447885];

% Layer 3
b3 = 0.40460825637540354;
LW3_2 = [-0.0012760847812544863 0.0020341127831345392 -0.00069185861007384989 -0.00049308913770284473 -0.001144572637401688 0.00026823481292034406 0.00094466206124155206 -0.0019277741381943526 1.7347992936099139 -1.7118057348597462];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.00111018595614765;
y1_step1.xoffset = 1025.5;

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