function [Y,Xf,Af] = ESPER_silicate_15_Atl_4(X,~,~)
%ESPER_SILICATE_15_ATL_4 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:43.
% 
% [Y] = ESPER_silicate_15_Atl_4(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timesteps
%   Each X{1,ts} = 6xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 1xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-0.999999996490807;-0.999999994516886;-40;0;9.1619;-177.232103709389];
x1_step1.gain = [1.0000000017546;1.00003626070655;0.0153846153846154;0.000309885342423303;0.0656325902039552;0.00425400504498441];
x1_step1.ymin = -1;

% Layer 1
b1 = [3.9172079974091791321;-1.6397886675322814742;2.9892900945062788054;5.7140822092093159768;6.0162457516258944779;0.84487106624672569932;-0.51540750878667973556;2.0580166133742965329;4.39468582496317417;0.49877492804972678453;-2.132397720001781849;-0.48839891875905305785;0.21041132979911864043;3.4861403045171153181;8.7065192061822109793;1.4313793578561391229;-2.0425888670778733491;-3.3417875484262125063;-1.4852023847155673142;-1.7489374523612581935;-0.6477427071260429825;-3.1190881120989866204;0.86179139292073658751;-1.5402976939356591846;6.6654637904152425065;0.17564078857885143115;-2.0335130877777141123;-1.7885051736003407807;7.1547615982825130132;-9.4408349820298376187];
IW1_1 = [-0.97519509758453160053 0.27605712880085969685 -5.1812116784407100667 0.19388814199321419141 0.4197530832933566991 1.1373056203812490494;0.12610076848091952839 -0.31185535493509164429 0.076481395929633488873 -1.4249344422568297475 1.0720813208187456311 0.79914658730439303014;0.0073252282463045536101 -0.018487073882764729071 0.28345345790248477336 0.34046801478634342031 -4.2479884637164087025 0.50909386432847514303;1.0161157400125335748 -0.57939210155298326299 1.1200177033885971856 4.5082948340216617567 -1.5766908278528355414 -3.737116450237511156;-2.3538891571910984801 0.22688695301608788335 -3.3174045736471020618 2.5000572657755770756 -4.224015615379593136 -3.6115505275041748234;0.26169813464252517088 -0.15323771816346840047 2.5143298095984909502 -1.3169306750003930695 -0.61011629603774719666 -0.83879166616922073985;-0.84910716963493204634 0.13872289556321856785 2.0673180289467247839 -0.32973121226575929921 0.98676489886357054981 0.072253857535450466409;-0.67332762112603228299 1.0256596124280852322 -2.7318001391824404678 1.3578387483082763332 -6.1595303989593039873 3.4957975088015973419;-0.25140276944641520851 0.010118987878859952989 -1.9202915891811114246 -0.64881563913600737159 -4.8892587181582749523 0.88027854361894330992;0.32886813755256316671 -0.15862473779558605713 0.3983016748871103796 -0.90419417088195519483 -2.5564968419222489082 -1.6439219727661222858;-0.13586186831042298961 0.10050917435830157565 1.1634775860580321982 -0.25136383766751985158 2.5532568196648526282 0.23022198793784751092;0.427261743494307622 -0.54911360050378821285 1.370960117230735742 1.1039052989669806415 -2.1897956063613825961 0.54027786911952069282;-1.7576997846117055069 -3.5293533560699392737 -2.4193464606864938382 -2.3403864431129939305 1.4194423512138281129 -2.2256837200399854915;0.13557295834529969891 0.020490439157128167297 0.089501011821470244367 1.087281453013959398 -5.7863722671586899082 1.2283406195346275869;-1.4317796123530008856 0.5030102448488467104 -0.65414940959202527715 -3.712733160441203939 -14.641671445614427682 -1.5067013867834833096;0.33438919338907296952 -0.0098752743152031204066 -0.21802563616029888016 1.6020346904230700069 -0.92329233115420616951 -1.1110604714555207817;0.27280388627109286581 -0.262798357825043849 1.7903672341055769657 0.12738000931863327492 1.8602556589987118141 1.5178907961145087224;0.11669728916415354758 -0.35995324640718434628 -3.0954678528677388449 -3.18308430694878286 -2.0022133450062380433 -1.0858763881608797242;0.018980542086286422004 -0.062930954985262202439 1.2150722748800488748 -2.3482487725279765911 1.6469691433042856765 -1.5723657499589593911;1.1105926791466944614 0.04680060642081738248 -3.0671036263833477697 -0.98211134464055449556 4.0198305419420554685 0.92778510885005682773;-0.90503986664698288234 0.51471518822661332937 0.36386003022452717559 -1.0512618274440308497 -3.1389277308344247785 0.54974728132425465521;0.39461782949565660239 1.0699266974260648322 3.9939005455066727279 1.0480130046953590739 -0.10250673886515450361 0.24337685227251862341;0.62076431886161709262 0.60040242500128726544 -4.1178918944807012537 0.30812423064032395947 3.0713887952915399637 1.0391375685365364934;-0.033709437867215717721 0.35643899615897739475 2.1243576206860019795 0.65949301902489565474 -0.11962785794257707606 0.28097039328294270977;4.0317470853605854586 -0.92218715477163959449 3.0539582737903656628 2.582559478161457811 -2.2093025004131812139 -0.48132249087075640315;-4.821458219201899098 -1.9592487155614040084 1.1974933063445323267 2.8030879738894789455 -6.2431402073549495668 -1.3865998044961751479;-0.21631354092145596568 0.12161846484118286893 -2.0106527042783950421 0.16242015821338184889 0.89277717559137215186 0.1621935973384122931;-0.035556936984487676201 0.055694957291455300663 -0.050460485392976393848 0.84327762507526238256 2.4699085961736630956 0.45308401150791255407;-0.04070320489072973974 0.89471424835658952546 -5.3786604027470978906 1.6039914663471497303 -3.3136350698682313087 1.9764450843704113314;-0.11726392787635021853 1.1245424117057092328 4.7775377858243093243 -4.7104665595566572733 1.9056174473985032503 -0.52733745080861427201];

% Layer 2
b2 = [6.594363380031925459;-2.1356303631278361621;-1.0083693822578165822;2.075441887234568128;-3.0006769822231702527;-1.2266222938746507154;10.278298319685777784;-3.9839581988353240938;-0.072765142405706195738;9.2598290709103387996];
LW2_1 = [-4.866561026751326402 -1.1431742359967900224 -3.0001981658646128537 -1.3998676213269583002 -0.40005682867644654088 -0.22238456756923294377 1.9561791367163636757 -0.31433422233940772061 -4.4807159347793508175 -0.79106408641500514101 -3.1688002743062448907 -0.27580400287782974589 -0.071899119175076742216 1.6515134482563873508 -0.5383179075956171511 -0.21332588776844321754 0.5508824621032155644 0.07536616677755506033 0.69122459710716865899 2.0227668389275916816 0.8361842632792629848 -2.5351173187293563061 -4.6954110669041968507 -0.87800044537116883614 -2.8389449720054762416 -0.28678829687431284245 -0.82368358985626688451 -3.9075241384771821629 3.8122943049160227247 -1.0244776630714722554;-0.49272741092216321146 -7.8105124595254231679 -4.2021639171963824921 0.33184838159692436887 -0.2487614700280979807 0.94639713735886665802 0.31418308994410326651 0.30325917195368939927 -0.69187584253171985349 -1.4992331246879884166 -5.8695081540401217168 -1.4102671125206116542 -0.49218962303122643087 2.5654910477153682891 -0.18713889818908718832 -4.8020342323818736574 2.9509005636727052568 0.52131972855722918947 1.5445438756085168031 0.18066659920929373229 -2.5653777106102255701 3.1496021528992552696 -2.1209636182728872633 1.6775038679529661323 3.1648947958926072488 1.1340523090650656002 -1.8689570689246317592 -4.7498128999970052888 1.1559386710538319676 -1.557753677352681887;-0.69631459814354612625 -2.8252447229115689176 -2.1484476401104042331 -1.1622880603934266919 -0.15167470254157755516 0.83974547091572426805 0.12970336398601248207 1.5607274142367102421 1.3746160719867270839 0.24313840334906836049 -2.6814568865048467039 -5.1074671549235040047 -0.076444748625151520693 -0.16304446873115080452 -0.054471972406518853749 -1.8073986213896777731 1.7739228815777097203 -0.4725368473692906357 1.2192010615010699937 -2.496872534003684585 -3.4007757395897382224 1.2743314257577766924 1.2272937305004876229 1.5141157192232186368 1.6340335716399330135 0.4604902143479243648 1.296295212477345915 2.9783151577166768753 0.37778918066556277378 0.27345200983800282968;-0.21282038749125845034 1.3537592560654720941 2.1077105004685927625 -0.3872604384523757326 0.1268427719871376147 -2.1368678124203341895 -1.0844107155497075823 -0.25674110362978697086 -0.10380189246169505002 -0.74773722754601457119 4.4602012324063800008 -1.9895948184730414265 0.024088334980625531606 1.7987617473802428592 0.51185088478628848474 2.1932402385631650787 -1.3951119730272358854 -0.57845229883519300174 0.55304622045302276323 0.95227673279243352233 -1.6046833202205925684 -0.5229234060763334746 -0.44247114509657148318 4.711213545963196303 -0.22147288105538176506 0.4690264779031593223 -2.0917102874780177935 1.1080832806894078679 0.0047985925876395437578 -0.38126972012095039499;-0.19594273162175362968 1.3134933285431162897 2.9130997647131620099 -0.2763874455898512883 0.29191412157276547257 -0.39392313680309076762 -1.5218532061877703132 -0.2712702334123771597 -2.2924277300905488985 -0.54613006717056999939 2.0210626249721226344 4.6553816772495029142 0.17509797334377669165 -0.95750857692157609335 0.16521513079105040744 0.64299288532613052549 -0.72131351600966753246 -1.7119678773347397627 0.16606084397891818227 1.8694163576004194649 4.3623566482245985654 1.1084271578189308638 -1.4059928742947320757 -6.4340007986236589232 3.1379205760119504021 0.013005605383316392459 -0.54305120463272316123 -1.6664457987036465081 0.61903619277125898623 -0.25142007801651444776;-0.39663385562871333878 2.1687520622299936157 2.8272034005177837379 -0.16566099262516986723 0.13013566077855362457 -0.51175368676480892383 -1.2022237892122664693 -0.77297660552169200621 -1.1598997006271796195 -0.73015927251393264452 3.1708236235528661773 -2.9547756430842824216 -0.24347240523028307324 0.69973651301236428068 0.5823577772537740227 2.0904450652808841404 -0.94331231075927157725 -2.2004188430509028507 0.78136835461532139213 0.23170490087252446076 -0.96386329670953896365 -2.5578745286896644728 0.056368954723407672935 4.4207539314772832384 1.8554302352058742454 0.71440114077006011684 3.1419117103913856148 0.87980682858284497971 0.16671531603679914935 0.04897636521297041301;0.73283577892294859968 3.1863747732538016955 -0.19645286730375277617 0.78557710461896523846 -0.97286443319444215749 1.4094305201166192898 -0.73765741347107904868 -2.3020536805087568055 -1.3572161511225144803 2.7150960763615872295 -5.189810922055801079 -2.0542387715807413073 -0.073690802192999158016 2.273581000711071276 -0.043883060475147638024 -2.0724496034880135475 -2.2511071272282663358 -0.38240054452333388513 3.6688922796951231042 -0.90429885836049572756 2.1398159445948738977 -1.2306145311734335568 1.6846243390250930094 3.7385220044555227403 -0.39395508842576754116 0.20796906075950108428 1.6600719356297499729 11.529881859054578186 -0.63136772369211813238 -0.3727166570123959799;0.46265998557801418789 -1.1352204714503775218 1.4817336513219199201 0.14229149222582940237 0.52148263188648658062 0.0058924895598691743454 0.48523444856190484087 -0.23536997806917903109 0.97940338350207067908 -2.5902489394790406685 -2.5283918419250102083 0.56304010162092910097 -0.06688115156265579675 0.37366215589917456397 -0.26400581154976165976 -0.8396534581431587263 0.27723356300211815517 -0.74468018967975357914 0.066000384410099935395 -0.45406339201296053165 -0.86822238618548863709 0.37729878036731206992 -0.33062498493220376927 -0.5117299965415762486 -2.5714943351359873525 -0.26924583870294072074 0.35608410739154094937 -1.7630696901707485491 -0.47585771258515768833 0.78493583148570500452;0.17592347404879690731 3.2407323912886831785 2.2130746408241486201 -0.32862268862911941936 0.14221511110796586363 0.63595608492701949377 -0.84824136273076677917 0.64656330681999274468 -0.15535726087799067718 -1.2035587151690301244 3.9801218298407454732 1.2888850606805892784 0.30195211600066990121 0.22155870487055032148 0.23332885700526520556 2.9803096591177573593 -1.4925690884497031519 -0.12689706101167905206 -0.80692063233434263569 0.1434931972123915378 -0.27760082566783661484 1.4970863266882898923 0.35592423195150341142 -3.8777571366959424637 -1.1551097031101122958 0.63709327813829175202 0.046942692796950354073 -1.1985940204967882039 -0.17570009725172883175 0.26745180913773103359;-0.092355766501567612692 -6.9674534298134735977 -2.894625392794786034 0.55648358562061828181 -0.11936857228448115309 0.91638998405361671296 0.1944327554832224747 -2.1866623106323772596 -0.74566564749247510324 1.2201958529847818369 -4.4222179614148116045 4.7355900039380065536 -0.14960236447005906268 -0.66521577159182765371 -0.13377170133989405021 -5.3852405254534758683 2.8338517133257448144 -0.24418999844054295467 -1.6583976835359888202 0.41643795511409109622 4.267933297934636272 -0.78013320723240797694 -1.4320894789722997675 2.8099250253214300699 1.7169251984313538273 -1.5597838886031056838 3.3431388145388747368 0.25689561182965370278 1.5053664290452428709 -0.479301380786266773];

% Layer 3
b3 = 1.6768277347327316118;
LW3_2 = [0.74191192854573040361 0.31779939282855140803 -0.40902161102336370702 -0.65444574353812967615 -0.53023595768444109488 0.72478357751789401764 5.6519919239938340993 8.7448597138713068944 0.87436445092993797523 0.20414276503418687936];

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
