function [Y,Xf,Af] = ESPER_phosphate_7_Other_2(X,~,~)
%ESPER_PHOSPHATE_7_OTHER_2 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:28.
% 
% [Y] = ESPER_phosphate_7_Other_2(X,~,~) takes these arguments:
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
b1 = [3.8304577004838309406;-2.2424265847241899863;-7.7919095030031950699;1.3749227502998080386;-1.5983060332346634702;-1.0680337119310558958;0.84693936235113986832;-0.53384368592313791169;0.26602627808399881282;-0.94685025973900394103;2.736725472698776862;0.71998063831533698487;0.7319820592755144073;-0.94295295723146255185;-1.3698479544937876007;-0.23856969486854459617;-0.034103078340119802325;-5.6173400063938609961;2.7579711104743522831;1.4855174749695370107];
IW1_1 = [-1.8099740557860106716 -0.98561517902151263293 1.4828262073849605152 0.45217009995774926079 0.50568293840043532228 -0.22207416431099818999 0.69652599485431521131;0.075395987882734633034 -0.24859292248746711618 -0.97705568580733670458 -0.13100456284369504933 -3.966347623595163796 -0.15801006224365785058 -1.1124101077019614703;3.0648808858508793662 4.322788666090944254 4.8701761424762404218 -1.4628794086133485308 0.79713779062671119569 -2.5153247875382485077 -2.5049439566546611857;-0.14464376520366345469 0.068705850777166455057 0.19996501712547648522 0.45096063096195410091 -2.096891075514665026 -1.8660358425900416268 -1.0635854795676273898;0.24978229217919561966 1.5358461782828136677 -0.12024400023469571097 -0.051823950160022801137 1.2715762025213082076 1.0874612573432347329 0.0084817079779255719063;0.25746862778292656238 0.83260911197266906925 1.1307396968059277853 -0.58901902602020672628 2.4679167191259567815 -0.2477930509314935914 -0.13101201302666135029;-0.05322599596421687812 0.00094758835011709403005 0.0080075999709589527897 0.033241238484902610462 2.0299418854409845281 0.44978771750485441183 -0.48005473333596759833;0.034656082817492370185 -0.16642111773509821182 0.21640447704923604366 -0.16883488612910299564 -0.42825297426983754967 -1.0554563310159836842 0.63595649089288674993;-0.34198016319701934895 -0.21611949970798732812 1.6260262674148351358 -0.46832414984424852733 0.92317159406211413142 -1.8564423675823917836 -0.36312667926802161267;1.246206860159110752 -0.30013111021994010041 -1.3623536552304893466 0.057875413613789061795 0.16955131128854544764 0.23232567067869866473 0.42725366458374613998;-0.92005984563000819954 0.61306727852396269984 -1.716348136541061109 2.5768777286707456042 0.52983333974495505636 -1.7759559920334624383 1.0855812982301853431;0.23644297921097445769 0.25666758578754922393 -0.71853680712026768429 -0.20080507393062976829 2.0612038656157958272 -1.5311930901314201314 -0.28540932468833435376;-0.5425688826564990519 0.41636171449953179868 -1.1682695807482734729 0.245900778237311185 -0.87577818148477115656 1.0169016816564715455 1.3812765686014456534;-0.10367874328468063005 0.72299192064582629147 1.3867891614681304713 -0.76898465252143788362 -1.7464733872922928448 -0.012785513898317120429 -0.58022237825086264706;0.14903711963187307732 -0.48970936902188072493 0.10964345482537861398 -0.26159025376983324751 0.30884838347367815636 -1.6922475577035360761 -0.74886096173216709637;-0.22932676582027389234 -0.52508217362108844384 -0.14557056516936997537 0.14650820986458013007 -0.79302478574353951757 0.96066612745147594854 -0.49412646539385657318;-0.16875796206416068812 0.40695494895660316859 -0.52895586263124516613 0.22598512597250472456 1.6869499998078869574 0.52531601093398705427 -0.81287274204247728093;-1.3666297754511813611 0.57382203649786389654 3.5435392782612069595 -0.22771937955908938545 -2.8596771419130218028 -0.54350057251460981433 -0.12330210292373582071;0.059714565682376825884 0.040079307939473590716 0.053396206820538270887 0.1344929641623854033 -0.0088086012258257132684 3.0000649980008735618 0.55974635996234023416;-0.17717033405596110951 0.057927140281075061023 -0.28576775017190292338 0.77240420213751714851 -1.8766807255936441834 0.27386034106467932459 1.7739348771809462235];

% Layer 2
b2 = [-2.1562141432816761544;3.3426993398269160274;-3.3357390965209936518;2.0848627278541389529;1.3153130892040050259;-0.60511076609608727139;-1.3144336841884791678;1.5970534998421561657;1.7548093155836514967;-3.6142288982266275177;-0.19285835973339646165;0.20999594000017438034;-0.60157528432370666049;-0.47227539993602746238;-1.0007958061370114677;4.1645570910565874456;1.7238104742316420204;0.17060789857121433366;4.6909074840721611821;-2.3753246054364720585];
LW2_1 = [-1.2035893247829012598 -1.191444974268985213 0.17235319441144783537 -0.40509985812094362334 -0.78753181174100539419 0.093699758670292759843 0.63765647785868029285 1.4617385765524049113 -1.7425944217254458746 -3.531256873048394862 -0.15594676290886375081 0.0016484883438588580906 -0.48934268416547771485 -1.8759391857938818582 3.1731661695541810708 -1.2205109281854520376 1.0529014692484268689 -1.6760092379965401133 1.9267307804456894615 1.2368547495364388134;-0.91332249794048425695 -1.2995256982732314999 0.67123593045345908426 0.95328677425562302616 -0.21271827681051039272 -0.21036824943563534052 0.89525689856328449778 0.40776186366644134873 -0.36642719371925325245 -1.0411549888010716014 -0.22015219489863505076 -0.33977604935227034844 -0.89810987145907961171 0.097812511323752293624 1.6790208930571874912 -0.46976852679656344591 0.94589172261927501051 0.12720797578638512459 1.5381369942350562319 -2.7509671865178453309;1.5357494985562907264 0.37429388855222456289 0.002139379708890127435 0.42584191959229628077 1.0327624558677088729 -0.78398289999040904341 1.7168835299045159815 0.27141032497632283649 0.48519795484873579072 0.99841785736505683779 0.14733445604146322427 0.52572482051552793525 0.62714154084465001393 0.79786786526767816685 -1.0752590094292593381 -1.4314713403486585275 -2.2691243946553050215 0.71652200692633050139 -1.4160931960853970768 -1.3916831896944985569;-1.2960244431527294662 0.79155888993093304062 -0.3269979265550894687 1.2834339670415766399 -1.1608113135491764378 1.0966214176648554357 -0.79336155936253838217 -0.63170337485105165776 -0.30683070540183760011 -0.46792646870226500821 -0.055932383119484736789 -0.37407837803549470745 0.74222406901520399103 -0.67099871397669319517 1.7263014475699869354 0.32326018038944198674 0.90490178188124081871 0.64348469459586987718 0.30141008785004963499 0.56035574138832833224;0.66395956662428468054 -0.88597683651397907134 -1.6363685875427307348 -0.023305354742972954341 -0.89209124611635670643 2.2813065063115995912 -3.1262361244957976858 -2.6806098681898240521 -0.46492230301139914062 -0.25102795929197224689 1.610683877221096294 -0.59185716564891588387 0.85009616714966995676 -0.092786982643822721317 0.1100606436303836061 0.25123140875858679388 -1.0220181846184532848 0.47210559209559377347 -1.5138883119370221131 -1.5933013742322914563;0.075818627358902379587 0.58252902633031689028 0.20107404098251749347 -1.2846595744313027421 -0.2870065107796406223 -1.2311842638634897806 0.003409400335913897323 -0.17753631357293431225 0.33076808043031191708 0.14539623253209899012 0.36013789934425399775 0.8348867955847228961 0.18567571948018943284 0.60653231364329052155 -1.1698278734636371468 0.23440150122341493755 -0.27528463228104410598 0.56822632563520769633 0.44802147845100026569 -0.6123816751393101887;0.47729356949799872645 0.47482384427313001352 0.082184417462666997656 -0.085029083618405609801 0.18617085937398253348 -0.7812468437531960852 1.0987832667835990552 0.85437137248837491388 0.27922531106647791033 0.13604793998407546662 -0.14520247580204415083 -0.20518086584998349031 0.011679292485795281578 0.19475710116709285535 -1.0705622776372343274 -0.58391730517704432124 -0.36125035581822817488 0.8706525058397741601 0.58256849770526020205 -0.79095956735713246477;0.14789905834722560929 -0.31017094087860919149 -0.22974396478875566685 -0.74187049928690407263 0.24297465504495713295 -0.15307790082994621184 -0.13616857382879990479 -0.86091113129904284218 0.18784795148526964303 0.39650563836614832747 0.55534056222653949142 -1.0305348280609836209 0.51492218943709944234 -0.054268989436943707028 -0.19092990381344832063 -0.098014433801307240812 -0.18689560390136772394 0.2940859698857236082 -0.58385240612337541144 -1.1296225699077282467;-0.78118985942734964212 2.759564927654004407 1.2787518467072798156 -1.4936146082756958098 0.31717990798806067199 -1.3269726034477919985 2.1984186619912242655 2.7446187999966156212 1.2146965096490625058 0.28013773011108972222 -0.55536030194573737795 0.62670040169093454985 0.12107603566728251943 0.62229609394468388128 -0.00069212827974546569515 0.70053836581997819266 1.2656027930627322675 -0.547698726319208129 0.89354779273136863971 1.3980975742102308423;-0.68415939442607542897 2.1502791296629260032 -0.49979403526292948134 0.37188896637789203847 -0.0031485747560078027052 -2.0438984421252430757 1.7561283872647883708 1.6432434685950536579 -1.7881042754747928925 -0.59334962140652913742 0.36120331641138891943 3.8965433741947941293 -2.1747926531944226092 0.036886419590614540065 -2.5260454125840339223 1.696942000691441077 -0.33693511739241027714 1.3451460472770968568 -0.53209607780318612935 1.8104324663359010472;-0.06652130379776832747 1.0666258517080751123 0.74840805511878694123 -0.74165780506126965932 -0.73462112155884629061 -1.1203466705782763047 3.3552407889871558311 1.0160963006883689808 0.93907829719190949458 -0.20046956500532661649 -0.61956619358871201086 0.069885431318110050047 0.60664552731155940712 0.72623519402689185043 -1.5825964285793541197 -1.718869549860085888 -1.6084919675267272954 1.9442016416394691802 -1.0793381522788321192 -1.7220514675893865419;-1.4917216627445977206 -0.39915489867954673242 -0.49193357948865595874 -0.27910561254767951578 -0.41203213232280172518 0.71564717486073581298 0.2234977078032280573 0.96754503643207556696 -0.52461719126145267555 -1.1090681221185771932 0.40867171088868131656 0.04543949369336180294 -1.3598679197052583323 0.064205311336328979577 1.2929876322092328511 0.60505625718074218522 2.0226632050679858033 -1.144790718290931375 1.8792278287359984734 2.1794151004750470513;1.4304807112791551038 1.1143144259544863761 0.69955738351833940403 0.45946362562340825608 -0.23700915907709621533 -0.96053198925037219791 2.1887900075578494707 -0.097567839437486078569 -0.22022797613493166891 0.40696084688786404726 -0.42363956003566477504 0.27023462764219874366 0.23257747225437130911 0.7097965774532992711 -0.81489287128864928977 -0.23221966357100071177 -1.6028438820854584712 -0.64717476895345449339 -0.75127816747608788805 -1.0691610079549069834;1.0226885051295397933 -0.10615288600703119692 -0.027796670249043486239 -0.75127126327446047682 1.0255046115029602394 -0.70605639318355817302 1.9121148141338224491 1.2939265304718474692 0.073304262917232673735 0.56694592390552212002 0.035653971322699956614 0.41488068848600367922 -0.030561783270159671028 0.76683553377345681223 -1.0657938616720494807 -0.041476845478452137317 -1.2628797059515985968 -0.34529481089877661626 -1.2094955256768731289 -0.5781594658223736527;0.14087040898446409609 -0.72960377101269258127 0.030657932100095684064 -0.58390036742952100557 0.55171940071605030109 -0.56461891430421828897 2.2139991017562290487 1.4290326061558586979 0.37422445170181589669 0.44456862937714586304 0.22456406224951522277 -0.79988398161658702534 0.64886608908989285016 0.082832954501183705953 0.57693659849353551294 -1.2320434709913150595 0.94918107270608409465 2.0064171002164186142 1.2801293577873142926 -0.35751946346348340011;0.80462853107944598552 -0.60978145773038761046 0.70324559560667931901 1.6294548622777371705 1.0260386453063972123 0.75848371276502268845 -2.8522204521617422479 -3.2102825731942910004 0.84829094686621153709 0.34362743315488253337 -0.60972867118878359971 1.0910342392280647505 2.5048863249909061324 -0.71805692674125587427 1.8955541624015836977 2.4615900437889002283 -1.0530860176563630048 2.4550835350428790171 -1.7204307034368471196 0.69563322659143056459;0.40112182513718469545 -0.6751513068757866165 -0.57917419248781620844 -0.59486796220794624368 0.55036912636435764146 -0.11930585737788075573 -0.70799853161508330412 -0.380988194211736797 0.3091280832561532077 0.64932499574183011326 0.60442460967608924172 -0.68637676122270174783 0.36596429841980859798 0.0620310401708023193 -0.37736070798152959727 0.68095528981398134594 -0.0077319677606869597827 0.57011298337549576321 -0.55249067282158925707 -1.0827498681336924502;-0.41658713216226139764 0.76862020065050196127 -0.049457407402440096977 0.69472940315106068709 0.59593252850175904012 0.35072500900701680804 3.0084557041360771024 1.5311948568786375446 -0.063900348181728988672 -0.24270718755663805988 -0.31022909838664664006 -1.1401220648226826349 -0.056195294928432115711 -0.11294549680693090343 -0.29954326804326264488 -0.90975601920091619945 -0.43556900537884729596 -0.71449647684064421238 -0.8496845528121781399 0.046557524758405688381;-1.6684935193905299222 -1.1167825643354090115 0.12715145554103796099 0.55878688557494449185 1.5530042129137691109 -0.83995090732173771642 3.5688094719312930536 1.3198460732501697823 2.4060236372346213152 0.22003167087196967566 -1.6005880154644314128 -2.3362271944430026416 1.9509216100896638046 -0.20627466066690808288 1.5401043840646819749 -1.8441300205723882844 -0.79524371335760002388 4.7258174826602648011 0.32975247148840619582 0.57212629651573254641;-1.7607425495649398073 -1.2947118738466587562 2.2111577283047978426 -0.047934251219374383879 1.9587566463158632146 -2.4292069929064101785 1.2129577537210751714 2.9809019211989999931 1.3280332306858682045 0.62852556688190763801 -2.2358097422009604038 1.9003883846896709731 -1.1532936458200770957 0.048787596763775337161 -1.7092236677809000689 1.537222138137332772 1.0631035943756781403 -2.0396143696060025086 0.11404549754097076386 2.5384267576033585634];

% Layer 3
b3 = 1.3618782292947428925;
LW3_2 = [0.19922333327505759026 -0.31809837075666286266 1.308237852757084374 0.7988482232770510727 0.47348106925259658695 -1.1079391269204399428 2.4770425030440543779 2.3625655043182525183 0.3625066614485079608 0.36901719820179357257 -0.46373804250236316804 1.2167134240640427034 -0.71279319536964202975 1.2975057216507523972 -0.79817464451539310577 -0.16569981958964197233 -2.3007307173450031179 -1.2516189297969566585 0.21594341816251300115 0.29827270893643820227];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.568038194888224;
y1_step1.xoffset = -0.04979;

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
