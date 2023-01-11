function [Y,Xf,Af] = ESPER_TA_11_Other_4(X,~,~)
%ESPER_TA_11_OTHER_4 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:03.
% 
% [Y] = ESPER_TA_11_Other_4(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999540494077341;-75.751;0;29.715;-133.803853032615;-0.1386];
x1_step1.gain = [1.00000001035699;1.00022983623563;0.0140762793578401;0.000302617642608564;0.251952632905014;0.00432440760873659;0.00826143816767911];
x1_step1.ymin = -1;

% Layer 1
b1 = [-0.67955050623545443234;-2.5145293396052950108;5.4393372499351597682;-1.543754833971539675;-3.4846923589851188652;-0.081500022678825639888;1.3132657336846020879;0.59879732661199702726;-0.197347336080903768;-0.84718288401070018789;-0.78543667312360221366;-0.75411015227028221908;0.067375546965412411171;-1.1188569936412211536;0.065467162529390138248;0.63104061544000422668;-0.16514338559686791985;2.447617654255109354;0.82903254032329243284;-0.66914758344965541426;0.64540768089731659352;0.99298386185458298314;0.2487716343368425731;-0.030232618485881580794;-1.5346742243934519134;0.012447603260406227624;3.2758073745955473832;-1.3835498654987719025;1.1672619871314386941;0.29124188163809222596];
IW1_1 = [0.66953746707700756424 -1.649257216284427674 1.0002951985488512943 0.32100130138512317224 0.042283655936269058662 -0.87537436058206852252 0.31885236438941666481;0.12362807141564194358 -0.918024801465376461 -0.65378629501011609193 0.24471332004104009994 0.30306749321519388385 -0.30793898016448278954 -1.7660055938222909422;-2.5807031836325999841 1.4227213657389179424 3.4598793701850447668 -1.2345753034234083767 -1.7435772875853785546 -3.4527099055732866084 1.6281550099326258518;1.2100512879721994075 -0.9570609339717054409 0.42952833061610623533 0.049338660194115320501 1.6454701887781282288 1.0324496920126391419 0.70216751001155031098;1.2440579830817930507 -1.5789928676366944416 0.055184115440302133893 0.70868573756140640363 3.8681527154573314142 0.2384675567166002208 1.2160886665956154307;-0.17375690936865240133 -0.27812871087699259087 0.55841990557801990391 0.59501576373252185181 -2.827495661182574338 0.15957725383150142551 0.57512642741813846392;-0.71837619981796441149 -0.23963056826663070575 -0.00192772283090390361 0.36434131243120526644 0.2877849539575960125 0.26063609161749662357 0.61682142460019062646;1.583955847106879089 0.53158745863890799477 1.4644853357153606588 -0.06462754206460934181 -4.2886844882147450164 -1.1044044227812992442 0.86199356425454065622;0.26747757397945365332 0.075386777757650749732 0.10338292785090420767 -0.19268801720543937694 0.015908925500584245039 -0.58041202656510437752 -0.38194590296388664097;0.15604693044938303115 -0.042327647194242423723 -0.64832256575580737845 -0.66535621463499672323 -0.8830218688427227125 -0.34795778660000348248 -0.17264686959741973316;0.13478028866246091755 0.028178165917809552127 0.88871049499505694325 0.18684286765389571117 2.2471730355669281565 0.30756189864984984084 -1.2610093699523718946;0.65167745400246956944 1.212490018656384283 -1.4378182915897970151 0.28185131180888978175 -0.35501684024907032011 0.1525623787229412609 0.86601230208232682362;0.5054002204282781463 -0.42218333090536175334 1.8553922938445126523 0.56643120415804237844 -1.8883733689591446492 -0.017407699341702363022 -0.66304617076671268361;-0.19361792621743026399 0.34138520440819897228 -1.0599263723459109166 0.53343907274386692663 0.29572305407170079716 -0.91076931251031012504 0.47394207496822166048;-0.26332496342405958822 0.018441838928555734367 -0.18443422912939833047 -0.12540545703398267574 -1.0476559577995450123 -0.038504995761532566367 -0.64149779773583193165;-1.1744644034555082879 -0.61496030535666657002 -0.1069481945918858562 -0.034934205515465131808 0.31464635512197725031 -0.1340226113137422681 -0.010054856937542985124;-0.3393388013127739411 0.44908346292530404442 -1.6076528963876790534 -0.49146694568442728102 1.6164348012240923325 -0.017795834117890402259 0.18052457694451676251;0.11346779407495656555 -0.32320804745977965711 -0.91464122465013752983 -0.5960649667162533083 -4.6956606214133698174 -0.94824864938299691097 -1.5956160475120217956;0.53443044953499307947 -0.085960382603073884344 0.79106656117441764842 0.44372745829922077831 -0.23917491332482757049 0.08477903468618490368 -0.52658798730386968678;2.655389104506745479 -0.6191466315247462493 -1.3554076543361448515 0.0047020759965595640417 -2.4796707102775874532 -0.91110279543723127205 1.2756245344893455407;-1.8872342600061025042 -0.29051738481388805768 1.1332995810050670027 -0.29523647117653251115 1.6001587666226722018 0.22763184679820719358 -1.9754615136641335837;0.5939147951598088282 -0.14531819961047107159 0.64110273295866471166 0.082559131868463378479 -0.65602439310643489989 0.22107914518097301904 -0.19858152603983128537;-1.8025069950653620499 0.24481138933970258975 -2.492939093520565752 -1.0631068126188958356 0.50238528186738329318 0.84203342858225127543 2.1618861650423601972;-0.34858832927510929878 -0.12629368760794978943 0.10704400411231068668 0.29583921388741357683 -0.28030215261625018863 0.43625195916176356947 0.70224556153150885507;-2.0639928063096566468 0.87063246006308769864 0.50070336699482087539 -1.5006508063186239088 -0.95568754262447330738 -1.9046106660684962719 2.7785008631631571419;1.460135650320788514 -0.56099367576827718906 3.8309590425238515543 0.81753717110763979292 1.6480306234976105362 1.0288109043432183132 -2.3858968810717597897;0.013375073504955784995 0.11289254809249053835 0.55526205853565846482 2.0231629731076163914 -0.25880478769288345342 0.38391806253356652867 1.0097166644553148274;-1.2708414772157325956 1.5000736453466343434 -0.49863016993568176627 -0.3962628649312521234 -2.2536604975948946183 -1.214743678704969243 1.9236975967229530404;-1.0337999244738096394 0.34016859539551386193 -1.9078542646969023355 -1.4922425376946657494 1.6887304815733503727 -0.82965420270354517385 0.1461310229602627575;0.5461924870228755946 -0.40656369396005259675 0.045114563843630897666 0.026450489770058455141 1.8230369095747114283 1.4424070147829926025 -1.641405511982693799];

% Layer 2
b2 = [8.770330182533824015;-9.4673827644262154024;0.45504811498008351789;-2.7363008317501140709;-0.41380242185241744579;-1.6850251581758577402;6.245740711974762327;-0.77464652774529019208;4.0652755672705298551;-7.6348252095443065457];
LW2_1 = [1.3765631394972726032 2.751180468716255767 0.77884085726802010363 -1.5789741006167203885 0.0016373866945813449203 0.59048887177771491874 -0.70919267297162702501 -0.61066145896458012832 -4.573410317281553894 3.1036452657923918608 -0.7940682207098759271 0.84450784296061487399 1.2830579201533152123 1.5235159894281562476 -1.3460419799288811937 0.48033504462128268075 1.5762819463625259431 -1.1678892413048400645 8.1899910626238092703 -0.032670526362911775364 -3.6902928326148396465 -8.9678866550664526613 -0.33274295260525371676 -1.2201920831258932232 -0.15869936222693001304 0.25550530237078189488 1.7848417305744759265 0.60985200791884397464 1.7454148095694876908 0.67116935748665051964;0.92434347005462813573 -7.5203222943169221182 -2.9839515844836195235 -2.706732692623907699 -0.25474760587524380995 2.2904055436441193017 7.5592846572771854596 -1.3188038076788610731 7.7310040591647846142 2.6672711622662355069 2.6861512547070902635 -2.9696360066460782967 1.3281609294727512172 -1.8749532633099721401 -2.2083494402842571702 -5.5990909413920926241 0.28271979806832781001 -1.8947425913462148106 -6.0271070832933490991 -0.041419360655170003871 -0.74758856079947078044 1.8196357000865950315 0.10595495156228071365 4.7461758700216174489 -2.8831467609202778846 -0.35325062205160429052 -1.3180987771801264685 -3.5491959607768048812 5.7917610313699769975 -3.3623188442271190191;-1.1310207665443907654 -0.25086396654721498933 -0.32444779880074503531 -0.79310703351062405986 1.2211066490891071545 0.98215416605940986106 3.396066600687861925 0.43597800807909442566 2.2067808523303540369 1.6819979623083400888 0.32568133191075165911 -1.171028693676980259 -2.4726364672253056121 3.4165386144318032891 -1.6292177285412148269 -0.71058899337347380953 -2.0495804369889865093 0.63767099827028284054 -4.6696889890788240152 0.83795746356021416634 1.9822101578494324414 6.2896880423902183921 0.25806021326809946981 -2.6907466190359872193 -0.038566664902409113214 0.29834826899204003325 2.2706246050097100841 -0.13857189693951241227 -1.5073720080067996019 -1.3027357623900230887;2.0181668421231369059 -2.5081222909291689405 0.029622690199159369684 -0.96715394032245838307 -2.0344474717995124635 2.4824723980587677374 -1.8225931212592245689 -0.99814590067656050465 -2.7764674992573588064 -0.40275380326054321944 -1.9090051442201456666 -1.158290619679541944 -0.34468664805384346961 0.65884460607050221448 -1.9852137054764578128 -2.6505181761003155927 1.3641839716954533568 -0.37347708998152284732 -1.0816102839902153043 -0.16485057050503115761 1.5746671461922363555 0.6901425195097335985 0.24718566096493077566 -2.8861514932709022574 1.4256073018755002213 1.2898645871905922711 -1.5134355086011650471 0.47669368949868379248 -0.84586627050680107676 5.4223286521515934311;-5.0579526315465610509 2.4787391993832730286 1.5035975745945002746 5.2072369184935851649 -5.8447743412740305047 2.2721760381643956173 -4.7051271089627890731 -1.690363617081911185 -7.2806573000599943413 -4.914257270936303712 -2.376334860539913052 0.38810841685667701872 6.8861557802593882371 2.4204763179680393037 3.4984879118418534816 -4.3722640594984758522 3.3008299197430601701 -1.0458297930536797349 -0.51597358970274243006 -1.3231355310372001632 -0.02371705579480605855 -2.8749258102682144767 0.19753636167739593854 -4.4239633268696900714 -1.2826993687530028332 1.012189409347837632 -3.6167919660269243209 -0.19249164453332504032 5.2246754727422448639 -3.4894797931638912303;3.2788322616749052685 -3.4579544835823377014 -0.73087913048114017123 -1.7250252668677357182 -4.8259185844618102124 5.5243862021761298564 6.4463665148297542018 -0.68246023193242344451 6.2076051422929490187 3.7713947247209320679 2.4129083322688398283 -1.1879464889839828068 -2.1353617826897264287 3.7156551806871118693 -0.4580256528953366435 -6.4323695467608352772 -3.1602411176284008754 1.7993505746617053376 -1.4156765560480926158 -0.088639333490215743705 1.4063114002435814776 -4.4424217941692276668 -0.108321738802623177 3.9303308849006430847 0.36085845552963258243 0.11168356488109026925 1.6804061663701579388 0.10145998529960265078 2.9119638820402560242 2.215324975565083232;0.44663333230795065898 1.5933437425481993266 1.2467671486532387881 4.3325510438320060658 -1.6341825328561805364 5.6333617138077185515 -1.3464059940843802199 -1.0984789162170109922 12.529946704039899075 -9.6276088491874833863 -1.569633829876282638 -1.4548795510811234344 8.0563452767389112097 1.2842430861170806011 8.523245748852454895 -2.7257220594599540497 11.916627347386262059 -9.6707340224425770714 -11.372788998000803673 -1.307625843182084413 -2.6564990215938735219 11.409174515813386819 2.8481469807204145006 5.3373567552748175657 0.87108161255693072889 1.1791303752192625609 1.060867762286973548 3.213867317596423856 0.90563975410271657651 7.5612123975764431449;0.20515773315341612193 -0.19331395742672077831 0.10051964293854652455 0.24271196767977826214 0.096254931772370447285 0.52961652399070013875 0.52972035244849824398 -0.10091250712275560975 -0.12291750107295580652 0.4807812028885308786 0.16912404114901105623 0.030797130745320901168 1.2949490553581837204 -0.41989065046474061127 1.0932898454409680777 -0.49307436359763401335 1.6504390449718076006 -0.27121517541013878683 0.59369410988525950224 -0.059253791345129036561 -0.15451908750340004328 -0.16073392888919690868 0.02467590282618716685 0.64994554962035222534 -0.13404413081021460119 -0.00048072776648688179783 -0.19736818774224670903 0.28205071881809845591 0.49706883543024948935 0.10510706229595258299;3.0231108941083069119 0.23266188532742493322 -2.6253674751006834676 -1.7742801895426654202 -6.2718876247926056422 4.4120377060675153302 1.4614154541669859277 -3.1581604451242157872 12.393123447239238644 -3.4366035557686309687 -1.552389190565452104 -1.4274896116991304673 4.4277452878832264105 3.5278033887550490633 0.40949103162770844611 -2.1868027719367395534 7.6592345202684786187 -2.0703439702215722384 -6.3925016891256118967 -0.18975892002666730485 -0.48388772701607507543 5.275361497161290103 3.3686840582775321451 1.0439837286003013084 0.29783804279985498509 2.5827464561389619213 1.88457149246080502 -1.4569550377865359003 -6.5771200086794170403 -0.8611927606844874683;-1.1629712497973141438 -3.9636895061332841905 -1.3550809394741500924 -0.83196493780560210052 1.7845615054414603051 -7.2196981431659983386 -6.9538265231481286577 -0.61111225896877319208 -3.9110458345598266838 -0.38456594279276501425 -2.5917855358519821252 0.8261237040917518426 -1.1053773328418070587 -0.93317601198697097686 -4.0202355471075303939 4.8589941658711968131 -2.1356536756580362102 -0.71737055553395912 2.6719899532656632424 -1.6909327342679618678 0.66824208171066801221 2.0855486778498049105 0.66274626187766882968 -2.7494275728060384445 -0.49910483971001617931 0.19640426930204141209 0.78001657781131528147 -2.0752867311353218049 -4.4496118345068698119 -2.9940443816070647642];

% Layer 3
b3 = -0.26696268957668406463;
LW3_2 = [0.1245676654028753938 0.029917029344260531037 -0.23549800600485043534 0.054241100049206278366 0.033441942566686193306 -0.043061782434100581707 0.060063181271746646228 -1.1897354493513307983 -0.02027510549872861903 0.11053027776541540783];

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
