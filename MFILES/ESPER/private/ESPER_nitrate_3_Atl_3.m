function [Y,Xf,Af] = ESPER_nitrate_3_Atl_3(X,~,~)
%ESPER_NITRATE_3_ATL_3 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:33.
% 
% [Y] = ESPER_nitrate_3_Atl_3(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999994516886;-40;0;9.1619;-1.99149636107278;-177.232103709389;-0.52];
x1_step1.gain = [1.0000000017546;1.00003626070655;0.0153846153846154;0.000309885342423303;0.0656325902039552;0.0604107502911417;0.00425400504498441;0.0145232735458572];
x1_step1.ymin = -1;

% Layer 1
b1 = [0.66506735922328152899;-1.0341796627883814708;-4.0437147118845331306;4.6353310476089175651;-9.741965298150006447;3.5860172415253321887;-7.2545030977767090974;-6.2618688664193786764;-5.3965796833539112498;0.49058022621092756133;-4.775615927511744907;-5.3028913090800333663;6.4239083966331111952;3.721724231386840831;-1.2793843493197827321;4.2827687530066560484;4.796244448821678219;-1.0626395652103770573;-1.9466376550162969306;2.2025196900758716723;11.424519478805491346;-12.846885862532177569;-8.4982726727129165312;-2.7670787925281583597;0.93996091296856443265];
IW1_1 = [0.085716622324254226895 0.25086389362443023421 -0.33792027974231941823 0.93600913108156391651 -0.6360727535441361713 -0.68950322134849229005 1.4366445749792686737 0.67790074325094706875;0.071694879967637126228 0.2761609218615208583 2.6994576684484616891 -2.1075893151900841893 -2.6772378087399459901 -1.5357568788523732195 -1.3145919522879208063 -1.299808425072251028;-0.46549523634424500029 0.078297219852803495721 0.72937595755650397322 -1.0422946974411475551 1.8363218402148471231 -0.72316723186315878813 -0.5219196864465337482 -0.27151054329512191998;0.62721911995455659206 -0.023833979658295537801 -2.9364850147109704537 1.2477774061735953381 1.6761943657527214047 0.34522706307137773774 2.4245788389489768555 1.6308389827237572689;0.16337704700663818747 0.65405906100449562324 -0.56187334334277860837 0.043035011503239378738 10.726963578324870241 -3.1603986324065713731 -0.85069594874980536403 0.73622934291718900024;-0.08387551941856975779 -0.069739337583621496819 0.25253773242673688282 -0.27402515038556080063 -4.872000494074372412 0.53118302531216710083 0.033108898018704126243 -0.70324467620185071581;0.87184923912316525385 -0.030624489639081306502 1.6690301220485523892 -1.4273365111776143976 10.348226046077531493 1.2906953302100048209 1.5425737432574906371 -2.2864803600189103072;0.28874899831206823908 -0.76896359022791382909 -0.53873277643131656234 -0.54550834027339123455 6.5810084349082647037 -1.0238978087514900572 0.5416796756974675553 0.033405579522860336106;-0.51632254471016281183 1.7670719951094102651 -1.326337246836974515 -2.2307228886521888178 10.76323335715821905 -1.0176016672320553624 3.0334996621306986597 0.79682510952661755876;-0.034851127505560013098 -0.07068870330753546094 0.016550222882165607774 -0.18150592779440322122 -0.20147391750342485972 0.98681308531085232261 -0.19167259892825869172 -0.42945807468902436055;-0.11372347957514079841 0.053327412068054094729 0.0022399945248491866072 -0.7943996882915033142 0.6465714955893080429 0.6573534586796582202 -3.1585399058576437703 -2.5180986543146315348;5.3172853047405341798 5.0340237763350001998 0.51172108626664081221 -3.1814679778340329008 -10.178002651870148298 -6.0805291065594113675 3.9680532555705272379 -6.1378735239197421691;3.6540863550415578764 3.7142911314221040975 -1.5222302935612386232 -0.31595254339274314592 1.5725706062269995122 7.5999149558398704229 -6.0841913625462407111 11.512546269550828271;0.17835465031694422433 -0.053716394966981662051 -0.46266294294105159146 -0.11714877152546408767 -1.919816051141067792 3.2044380888013308528 -0.43516940803382681002 0.92796243686009383111;-3.9012714103598233883 2.0264588324167065458 -8.1206906190412251334 -0.20085340210997651433 1.3392626411558952171 -2.2372291214118011915 -1.2122921494894813854 3.9942449753459201567;-0.21383296419877736994 1.1466239149246015838 0.016687570562216806624 -0.21065796569860731347 -5.9462518537080253367 1.0587767823625906516 0.88038341997395630489 -1.5579060848208399204;0.12391330094337374423 -0.16940122610238395051 0.083857958871380075538 1.1682142490579083827 -1.6495830216388343015 -0.60234181870714464946 3.8844450024725256654 2.0503639660958761226;0.12615307466644659629 -0.20790417993784740536 -1.2254421043308041384 -1.3301134756682904126 0.43076634337199654068 -0.47231568884719005164 0.76008371612483727198 0.65044044815901502776;0.10721996068183257955 0.39071237360385269177 -0.31338175518914851425 0.78042567930165929191 3.7868757544339071153 0.35611823582529544829 1.5302083384217370909 0.625446850246269892;0.34032722236207413324 -0.19980055914569178066 -0.35078096080810788893 1.1502076643901966158 1.2611083740474566817 1.4500139108447216429 0.22533294521874389083 -0.13995046391639515138;-0.115430612196081403 0.84411365033838603456 -1.1117369750023029518 -0.84445531898892289835 -5.1983366339208663121 6.2470388945572885575 -2.1239232929699203289 2.3772705139064527202;3.0151646391824078464 -2.194999546136738644 3.3331382273914815784 -2.443697587505280211 11.454032899630712805 -4.7001358561979600381 -0.52188925312514200527 -3.8253493097089341646;-1.316646823898510954 -0.077696349963944505435 11.260779393666792814 -0.44924527476932823999 6.5564004980440779846 -1.4818918231727760126 -7.4185779498778785523 7.4893630430633031381;0.22689999899334820088 -0.0913057744467001553 -0.29782837814682061595 -0.70547789445865227709 5.2185978030222175406 0.75594172783309598884 -1.8575955976663143243 -0.21857181841230899355;0.44457458085062240283 -0.52884009263615494589 -0.19999591303363159489 0.57696173732174604876 -3.7277830927733677235 2.5887509857139807679 0.67491862205008790188 -2.8854015781206996394];

% Layer 2
b2 = [3.6094803687935006486;-5.9119365750269725623;4.533462360160116944;6.1697977719339904823;2.2637348477582555084;2.6426396822689119936;0.5334669964757178251;-1.3631803493346794642;4.5459968079712513855;-9.0061203676402019624;-1.1375345801886060926;6.6844393398791632777;7.5077055870653746084;-2.9480940594917717235;-10.244720711458718654];
LW2_1 = [-0.29533259055838245732 0.61432411570916622967 8.2146890008865174337 11.399824305770238198 -6.8809240452805706667 2.6532919269749104529 -3.0986422541921645468 4.9435093124462294867 -1.6851557718167409217 11.112209917459418662 -5.7957621208960636849 0.74909481669264577075 1.3051669771821210464 5.4461662859267523373 -0.56400320818684845925 4.0706427660999331408 -3.756886356419486539 -3.3930957509824692941 6.1419336229277918804 -1.0278431174729809605 -6.1992740495206524187 -1.0539576917642783371 3.9818824703360471418 5.218876956062218575 -0.86105148260121100456;-1.2252290966514562953 3.1718392909438128946 1.0849590246170370822 2.5185323271829824776 0.66559084094665355025 1.513735283694084055 -4.4850436303038074826 0.21377178556888357108 6.136396354357547267 6.0248726393479614671 -3.3226014017875513673 0.33567327363081972402 2.9737164524695716139 -0.9071297051847386772 -4.9555604426351997915 -8.0374128751415998551 3.6240048669295781281 5.0763736902976193122 2.0821055276797815203 2.0994348678836449729 -3.9296580458575931516 -1.1646188495317861822 2.1103115656941655587 -6.7915326609199304642 6.2786178861544899021;4.8394685583757244274 -2.4546495906863237124 10.140488005714622943 6.0725766846289568335 1.4368304503942421757 -0.72268633174380092488 11.028546454730333792 -9.2652156687286968406 -5.2691274634807303556 -3.1672349833182420475 -13.679209515877044367 5.4474672775190722263 10.049852314122979635 5.3076308849354836283 1.0326091981410616683 -1.8536474039295081084 -3.886840703549069076 8.7509660398530986214 -1.4697312099728156998 -4.1427406966689632739 6.2444289242877948354 -1.0009261979157897393 -5.0975785968326023934 2.6998303370885508023 -6.081158467499742315;4.361568134868190505 -4.6755021380228445338 4.6164601751154838993 -6.7799511827939946684 4.2340593575999942999 6.1250565150518578861 -7.8228133732711420123 8.3997361885875978516 -0.61848118936906593124 -7.18268788257512103 0.89603745425398051427 1.3641309004777331104 -0.33113615176867183809 2.7685262791101896696 1.5761721566198507727 12.011278365767783072 -2.0642239322319997363 -4.26306170479927804 -9.9033751428952747631 -4.4711437289065685263 1.4142140982181548026 -0.28082907536641432689 -11.731992327724096725 4.2296810213914799448 -0.087435420017703616158;-0.39249845692163437549 -0.067812081070571170183 1.0003213775584514078 0.91652431769484776769 1.3652804258684352323 -4.3993963986442414793 2.7206075747976852242 1.9191257596911994732 2.2873982216138335488 1.1063450087009656375 3.6460355034585489875 0.16520166139285144546 0.40308964084780990467 -2.2594323759855288181 -0.16465320742117239838 1.9283270076607359478 1.4694111294523835021 -0.59527692778202578872 -2.7085106499851190343 0.20550353813847163642 0.70046593834159709413 -3.020258078655866818 0.1267730453835570481 -1.7532827310597480786 2.6428680434494866702;3.1470633727912549205 3.2933033925882053694 2.5802386332840887562 -2.1029120868591797411 0.53826992442457621824 -8.1468702851512198038 0.79712864231391078462 -2.9650174055618743196 -0.18578039586386746462 5.1973426582960220443 -1.0280248273696670491 -0.37287426754346042479 0.74562151654530250955 -3.8841377885899293076 4.7649427824796317665 -0.32870534083014774085 -2.208418955918046489 1.8802209126480327051 -8.793747719786265904 3.3638962604695601399 -0.54949731566669435257 -0.93562265176969561242 -0.41897090878428572269 -4.7005363684429797289 1.3982622503332360342;17.02555274531238183 -12.813784453386128703 10.309902747138153245 6.4958298043906843589 8.1733526158274223405 6.9425993988669914003 3.2599569818042875013 4.0263739465361023662 2.2020968855321014779 12.96000955924980147 3.8184081794285726907 0.98537991846916750394 2.8575308407982218384 1.7475403509137024116 -11.585973573457961194 -4.7885582057295641434 4.5608747447499959904 -4.2412958628411931983 -14.423666622220663314 0.27933986523533810287 1.4177205345793812974 -1.0641549619771097923 0.50059334054067727848 3.9899609715743764582 -3.217815119720246031;3.0897888575642689091 1.0095000147555315273 4.1263723414263573375 6.5457957258379595444 -2.2414953334671512408 3.2643115064320329388 -9.9611104579918965385 6.0961922251267806061 -1.1931042268607281454 0.15935900834132615822 -5.7709315237919165398 -2.7550099310871756231 3.5901430246827050752 2.3842935222358212499 -4.8820519930079804283 3.6637410632685019074 1.2482658454943695858 -4.9725070937292761286 -4.1491473485684711164 -1.9022281468991855125 -5.8541940398517278155 -2.1259997998820980847 -1.6046219843116911186 3.1526390850694094858 4.0191435648153417404;0.23858062389938203274 0.069347780874045347743 2.3929948612509854833 1.0187708156326431563 0.39028487194298955298 4.2644796223528542711 -13.999355165972543702 0.87366852944156281335 -1.1685539503212460932 -0.85910521284057062008 -3.2710246105842193387 -0.70463203965219156633 -0.48405931089931125699 3.6634777486587997153 3.6040527936402684261 0.75020945395523130905 0.15227873948789943936 -0.84314188398035572192 -1.0758758924610209995 2.2736855699736864267 -1.4343367459553848509 8.009682864660451429 -0.61743976534276845403 0.82770388735973188243 -2.9770898175483595516;-0.406407596457714404 -0.095519918863178121526 2.2244250073412135116 0.35398664559619752445 0.33919516520624137801 3.0076072219601486246 -0.063881715260146923474 0.55125674213864406248 -0.17311582844651873803 -2.791025219432034632 -10.230255442183999293 -0.012884522249305278493 -0.019029998595601359257 0.67285242182965598889 0.0086349920831326181653 0.084033552308965434041 -2.6061879729506731707 -0.16419001396194568554 0.65087668629691508571 2.1784343820724076934 0.45642503878681622842 -0.0070892037216854766543 -0.011661514789760218588 -0.30403928508558358823 -0.1596809381701779218;-5.5658985940806235959 0.27933160459653633145 -2.7766007260713463367 -0.37817512323031066135 -0.3163111365153113419 -1.8681512156036361727 2.2334451307087261895 -2.8197486281166939115 1.7630770898110614819 0.87763241037047390325 2.2711202261880822917 -0.51547180299031170048 -0.12390415348089395264 -2.3409547361912936836 -1.6971706097400431901 -1.308064488919936208 -0.39014085472673964983 0.79154746996700453554 1.7385270305763387277 -2.6996958184856030272 1.8474281753690080521 -2.4301182290567560074 0.39896941851744738283 -4.8881761961390211013 1.7963163026603796713;8.252708561712591262 -2.8381133960740640987 -3.0295297894617370105 3.2923031840972867279 -2.1955401943176644508 -3.185849447200378215 1.6915598751924958432 8.6586351487744845912 -0.8465472280888745038 -2.0812018619135366215 -2.7254232089325554078 -0.76689148527907990971 1.0834883996531328165 2.7697601619356380098 4.6300628920566264668 4.3005896771168901083 -2.1942468236268544146 -6.0403457328134875226 -4.190079351976470079 -7.6634667202393487173 0.21421674202913845719 -1.1691914096006528201 0.62183623645210184261 7.2769821953302447071 -2.7707909380039952296;-1.024380569506608607 -1.7969168295373443289 2.2019101394459981158 -4.0827382414890260875 -1.2272074180660976772 -5.2682940530363060461 2.5485155323123014348 -0.93761268293391197837 -6.6133850146235682033 1.0526841262000778965 16.284576735254397306 -2.8602799397903186929 -1.453264486377902065 2.5466282860827389456 -2.7808697396285761982 0.67966205697531067997 6.3544971631558757608 -3.8694341819256594661 3.1427022152397978694 1.4782547937636307811 3.1733971381197463302 1.9914378028330164128 -2.0160701807249061801 0.53562460092862174132 -6.3402182167399052659;4.0666872489006369307 -2.9976633342868281851 2.0141676702213038297 6.7895624976121524696 4.8477915705402052637 3.9350974673696863526 -1.8761167668526517183 3.825298138859347219 -3.7276580521212112629 -1.0158156569933756153 0.32260523826814452297 -0.43024745091382760043 -1.014889293653962854 -0.64061019170720445715 0.13861503009432168954 9.0086682191307012602 -4.6768822168124275507 -1.2566469874326540168 -7.3875655201235970182 2.9390382404970787356 3.8418757803760259506 -3.734362230533223137 -0.99553127812844244282 -2.2689851167451542402 -1.5752852708519788916;-12.904847118609048096 4.2816291369152974511 -6.2715913925239181737 -8.1687653887041768996 0.27163606761463104933 -7.0669400405906328544 3.8925954460281699809 3.2030280973074725104 -0.85331051593924567289 -5.7738631707645646784 -7.3133297694408341272 1.907251620132563108 -1.5049692619899113222 -3.0019828189762045589 -2.2694226670964798842 4.801317345668271841 -6.4463296876791744694 7.6849587817537416967 -1.3356800663686509978 1.8363940195905168462 0.12016513428660882679 -1.5813246823585476175 -1.0560939819084986802 -5.0512575981907792411 5.804901686082047263];

% Layer 3
b3 = -0.01594246687376000074;
LW3_2 = [0.016822780904881901032 0.017345564063333044486 0.010739634400177769033 -0.017377421496521279487 -0.22220753699878381671 0.093897600036949385727 0.038363285603557971404 -0.011145413821161144138 -0.062549914762889319286 0.84247639582743860664 -0.079080698004821647906 0.15594098198336303618 -0.019758969750887742101 0.01538718365870480817 -0.033782082745637279086];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.0470499670650231;
y1_step1.xoffset = -0.28;

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
