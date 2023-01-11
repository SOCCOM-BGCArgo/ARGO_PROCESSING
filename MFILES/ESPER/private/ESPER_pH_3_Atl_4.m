function [Y,Xf,Af] = ESPER_pH_3_Atl_4(X,~,~)
%ESPER_PH_3_ATL_4 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:21.
% 
% [Y] = ESPER_pH_3_Atl_4(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999025224415;-0.999999994516886;-40;0;11.29;-1.99149636107278;-101.945317545701;-0.12];
x1_step1.gain = [1.00000048738803;1.0003169936794;0.0153846153846154;0.000310077519379845;0.0705602618858232;0.0604107502911417;0.0054210015363163;0.0145655815308426];
x1_step1.ymin = -1;

% Layer 1
b1 = [-2.3878550025800575618;2.0069910728733764138;-1.9984123277828824694;1.9956206213682226025;1.2353862455031838774;1.4581143055511871509;-1.6565718930651318619;-0.99035037287498317582;-0.75933708380937303684;0.66198255020899554601;-0.9365720582406753536;0.68864093001661752513;0.91302006442531269403;-0.090078829670572124755;0.067163035087732927964;-0.25089829828844034676;0.2341577816312764293;0.37175233990566375342;-0.77250139079146329379;-0.73142551802526012406;1.056344671427307258;1.6142821388054664489;-1.2129875028143901616;1.3467096881498261673;1.950773945509634455;1.718939388464093998;-1.6092963934334514331;-1.9554322407290791741;2.0789194050600645625;2.2940623888794342555];
IW1_1 = [0.57823409250154333971 -1.0689858976335420149 -0.6493246103199644681 -0.34748184245319801988 -0.0011718381328671642356 0.13875872689728890874 -0.69396614085545338391 -0.72769525030197390603;-0.83286432813477295767 -0.35665934286953837606 0.39270236476078396182 -0.48710478138969220963 -1.3954081607849695423 1.0034266732829952851 0.77284082991449354338 -0.2046154729452358545;0.7601440301929246468 0.36439012521103564257 -0.0084356812014409893186 0.021794192988664019134 0.61440816078686799617 -1.4464068601767055 0.9323944672223080099 0.25023733265753289601;-0.92020260329579273151 -0.40314035899170219679 -0.67959263906276956568 0.45760955445217460502 0.095941486165614667248 0.22085629863785216642 0.98277322974763969832 -0.84874276874627785361;-0.53346715139825351404 -0.6758689461888536254 -0.77824536155201939902 1.0061391998233644163 -1.6821259136007669088 0.24384129858654782219 -0.23989777945203197285 -0.27405074434195220157;-0.57317714043316347894 0.10927064619634506426 1.435192716011730818 0.44738803102005453294 -0.53283548192832697055 -0.5442023954756176618 0.87274359788275757133 -0.052873523621351149748;0.48307030972525760726 0.69944571036074587589 0.53544474646209538893 -0.58275033813980292585 -1.2086371322708819864 1.0894278299736099846 -0.13442096240603890189 -0.90051874623034544598;0.58303574634781829555 -0.061378650029678738509 0.66286907505155889808 -0.13520288165008118031 -1.5718727493189013078 0.27511312050660868822 -0.96038266618731238378 -0.23818280704003302883;0.70591949128421227488 -0.76386826622253112262 1.108611621885962073 -0.89353130529049862307 -0.26768505763111172246 -0.64337889859998620423 -0.46622204800900196586 -0.2465216919224533787;-0.16850474065835463477 -1.1671916502049564102 -0.23172267397837109648 0.70065304492439473005 -0.0068206852583105751289 -1.5460726484918525525 -0.59589458584553678389 0.58076386392963552652;0.57345750672532680436 0.70903373448670281576 -0.55073152123713786654 0.19378911822097352746 -1.3194185243737397428 -0.62590656430288893475 -0.13017759991749777715 0.91067014899879938916;-0.2622757696969885699 0.33413203406410740826 0.040092757158289797637 1.118261595675567488 1.3733558612762584161 1.1323454654556419641 -0.11580758675850436057 0.27894927519684570294;-0.26141294171014678449 -0.26658616026427101042 -0.046129981320551648183 0.58317924418951949672 -0.14600226056894197169 -1.1686394561379047641 -0.46893586628853772824 1.4395865066669966126;-1.3254400132776655852 -0.22747710682747571131 -0.071981577663510556975 0.2203261160397168017 0.53612271829401147993 0.54523842963047264298 0.048662916725241978089 -1.6953252004172985057;0.62466582930345537772 0.41613823783300096792 0.50446371869879247374 -0.14701251435346163188 -0.53377810527803770668 1.0738761171896773483 1.0854651076211578875 -0.22528596211700580021;-0.51217974713757064009 0.85005325999785652513 -0.11921291509947018472 0.075421655617177271735 1.3978018191125618497 -1.3504254035603904249 -0.41601087550345389188 0.13209044662816224469;0.25100126304472630823 0.34035866413322235502 -1.481118594983954484 0.056118846599600494263 0.76042197649537768278 0.24834339178699266215 0.96867036442035203692 -1.0949188930272193421;0.55732511675666185358 1.0074374136634978161 0.96972803703161991606 -0.2705918873547561998 -1.8141127637013008123 0.11580624436498895835 0.36825139040534476953 -0.31981633049715546546;-0.28485767403060296976 1.2157897865396622983 0.5991914806710226804 -0.70548684645695003415 -0.88739342520472197506 -0.78503879808530563356 0.20764550253311050976 0.5937070237684324292;-0.97453174599266878531 -1.1064758200287465595 -1.3429417985345455389 -0.63390716335430430384 0.84029866323869750744 -0.027409515232161129478 0.70278394106124653096 -0.063108536774229326283;1.2108689507273682207 -0.67537174865977056548 -0.61691713024245675712 -0.70701799661952546838 -0.082509454476187724103 0.81042000020963533835 -0.29223260100713216447 0.80522950161281492409;0.15926115304940308137 -0.10521201267772489762 -0.64963694196193699781 0.861030011372630244 1.0831306295774887349 -0.83011866312207960306 2.2746639091999014681 0.08906923036323043108;-1.0903884144370168396 0.52300586736654031128 -0.7922142160984564363 -0.33760288788368214563 -1.1081248503248373183 -0.45111413378895326431 -0.52350185076666189143 1.0213833744380833579;0.20941479418071470842 -0.80989678518773100535 0.33820464534102440135 0.73816724934432742522 -0.34112363773641130038 1.4111354219515668262 -0.52398229110807725828 -0.18194065399397230109;0.28454966750602284087 0.77282712659060759375 -0.77451371349718400872 0.24858141888333099012 -1.4169185679891158625 1.033289146252019286 0.21530705971693475664 -0.31531277122494416609;0.56039071459047662405 -0.41496870511992911545 0.032760122502464128491 -0.6911784557580169519 -1.1051802972657589486 -1.3870815368658497846 0.86836739738957335888 -0.058343294304412290918;-0.80434340557112704584 -0.73891226756419969846 -0.6563361106120974009 -0.72641556100214532687 0.44544986166240396752 1.1462611322934013192 -0.45297642300269658211 1.1561141999138724312;-0.078139781252218154517 0.59604265935893507855 1.4435346107656503722 -0.60330602753329276577 0.29070843683840885685 -0.99697247142248335283 1.2498227718178844814 -0.06577321676979330789;1.0613302759463456582 0.18680035418736914976 0.79238627752559831485 -0.6683555359370750848 0.62790926208544572962 -0.27000890444480829533 0.93183657902824512664 -0.89995291777552255397;-0.00084774203087478912366 0.16011785912773543461 -0.31601103855740697579 0.1794610009430746278 -0.451464175091528086 1.0015771667303419168 -1.5028649655106223193 -0.49533163166996979632];

% Layer 2
b2 = [1.4049727307090300599;-1.0043903507325784297;-0.76421832558341673547;-0.57617663408058483743;0.27852096183024838139;-0.25525544112137782404;-0.46916035678420831001;-0.81951731987097653498;-1.1939890983487053688;1.4078120009872945317];
LW2_1 = [-0.021638728997988347647 -0.38906408073098069922 0.037856236369802918895 -0.13743180472752028476 0.35308563317900959477 0.2687119744942247368 0.27215514859533218583 -0.15630762618786261942 -0.1528670915788069018 -0.074385357655262690502 0.11544970291594440492 0.28735363215692122374 0.092368254503003621725 -0.14149365003067942981 -0.013265901345486268148 0.50036900018993935024 -0.51223442824403653262 0.78279967093076996409 -0.18498949142564305448 0.10012269862273165399 0.30649649731267469832 -0.00098341700224672222441 0.4880761354283377762 0.14486347963948409401 -0.57202110628339575271 0.51561405746094879365 0.031500778493727003837 -0.15559141515495958363 -0.35094201934803787024 -0.055152352186467343276;0.16956295842181229316 0.091878853326030990312 -0.26836601721792613384 0.015401294538599868966 0.31076070713281001012 0.29996407710525291312 -0.62981423468508634933 0.35023048343104662417 0.042345498824904841129 -0.54929429478719937752 0.1902604017621289112 -0.015213650549927577871 0.29250099728942235977 0.29791920392515996996 -0.12948172019887338458 -0.42780442636533311251 -0.74721711126553758131 0.012709411823201909736 0.17207619372682361747 -0.013009154223116317292 -0.073176800588557425198 0.033644437262592766236 -0.45785258744434093048 0.41291892578885464049 0.22442687222596935381 0.23221024087153172011 -0.37066981616889427276 1.0063905245038753389 0.39343098128285353132 -0.52709443141837386548;0.69965164072881003143 -0.23323651397326652934 0.024458789846662463974 -0.17444207662448871221 0.077724953618620945539 0.16099268037390127017 0.6516428153288895686 -0.10816908845097221448 -0.12114902659248172334 0.10953466643429801697 0.12902497404818061288 -0.19818391120791437054 -0.66208072952073249162 1.1754600314855241372 -0.18616871873966267525 0.17245646697410296766 -0.074055363114100297706 0.75667241876471769135 -0.59333404409134848922 0.28668444131881920534 0.18154896186968397243 -0.03089339924644735752 0.47772181582128891986 0.02067887675202489603 -0.26862916315835233538 0.49025313795175806675 0.17997304033704161763 0.96789913934171389887 0.0026049938112028716945 -0.051527459782697675184;-0.052630120627334027461 -0.74711274324740650776 -0.36500710020745852358 0.18304022430681834988 -0.0028634610333896663367 -0.34689099115714910893 -0.4129229381662943843 -0.3555678616117262103 -0.0049452908935512188832 -0.19226146354385922899 -0.075283510336844175481 0.27344013340208339624 -0.51254446270066578162 -0.30231740388010913678 0.029286754670639472686 -0.38512492138558196419 0.8603802125381534438 -0.14800052154428999795 -0.40046888974775896042 0.20629368214324633746 0.13716461342787228395 0.13130157978544051978 0.0031167773344537561854 -0.24428869734562411176 0.044353604727381090722 0.050889804351829097395 -0.61250736791964199845 -0.44311905513216542918 0.081472525834679399015 0.024799463917159959847;-0.41289317034715089516 -0.0023088168395955606035 0.27743960864075362105 0.067252256000142840442 0.67353122072966431499 0.55161246201648161946 0.067922724250712121496 -0.37612375826603139384 0.26186052787965419242 0.40524672946684942021 0.11838840906231724859 -0.16284148791426611425 0.38341220052026253606 -0.046848065788869013015 -0.28954155330594610884 -0.38906876764837272642 0.16487281406927545291 -0.14076213991641969847 0.058506352720121751476 0.43495136143056745137 -0.11965738065136531654 0.73252525071840235427 -0.53723330464644247773 0.07346893919782376825 0.14215335432405534766 0.35359450432223321625 -0.31560572173348167002 -0.030899587593299714106 0.54750897430517986297 -0.10935145195339525392;-0.030388385522029926628 -0.15168436573219098995 -0.20081046626132120769 -0.47601973208202713606 -0.47540740851278545653 -0.48371439524089854256 -0.13453693213721790167 -0.36981359136701796864 -0.074473649250303655678 -0.41726019900104688309 -0.40619803672649751336 -0.62521884533967397868 -0.059052791818606155394 0.21956233717229420299 0.36480511065242632318 0.15480261973585407453 -0.55180954962502826167 0.5559605705940233733 -0.10763455881320434238 0.35600854271904436299 0.21082104480161778515 0.59966991530321400727 0.24523407150693679601 -0.16330660824273962595 -0.52046266144118047414 -0.19517975344071361588 0.60111309978964932998 0.11330180030312775319 -0.087969644086063933375 -0.60856010281020267438;0.012939711195161730101 0.14009667545606063754 0.20464356600769775807 -0.24826724508098063637 0.73747843742285124069 0.2404237735855549174 -0.029397059563572872087 -0.29760768887273347794 -0.18319080365539081034 0.23056697143510676562 0.25719452015788896171 -0.29335307589874642931 -0.015446808693627930506 -0.15191449188174371843 0.56818788794275687515 -0.29022762316687922279 0.35158503362908022361 0.22407962927489347149 0.0022566026978739447441 0.18463828128057677147 0.1486080733341855753 0.56271571885166860039 0.10928410289392330446 -0.33694746886141763165 0.38614463985176727956 0.087816613701073131804 0.13744521593872566223 -0.30840186882914910615 -0.37416008080631307786 -0.37374200480908720268;-0.33814918925477027711 -0.43777851099603093399 -0.14745355981522556554 0.2384314096542996253 0.32575083820852029914 0.23157823125011572918 -0.18651027395470728965 0.66214092614044561813 -0.25549299168251887737 0.65701037917377058672 -0.019813737165276439772 0.15268637880589891465 -0.25500147926033545076 -0.28112764704135401184 0.28055491536349502768 0.90006742762260594226 -0.25115695595453890743 0.52356013529304057386 -0.0030029446381821013612 -0.021743551463443710342 0.28464171842317731453 0.033947838545268400889 -0.34459509483724803935 0.49404479679230556277 0.49634123839461286032 -0.10340870515700963106 -0.31896399324287544319 0.26591754192866184914 -0.41108543693182703072 0.13468892653618244615;-0.29715617229375318464 0.21969717003508090758 -0.41187054762995306545 0.092728289517983994439 -0.33534667519196220686 0.38795648917280095569 0.07741793363920898563 -0.051154534616794387902 -0.10912673483038157207 0.29195675077966498634 0.086083952153991602496 -0.25707979357658611308 0.17525988061648717364 -0.099385325200917681854 -0.0095386618239285611565 0.59627130407817485747 -0.75975680910527898515 0.26572153547682503261 0.05667849827515824418 -0.15533600983321957378 0.050217924364781944824 0.6536760154440389714 0.31533500822726362456 0.58715420455199940353 0.00059366165109881255391 0.13327566454281797315 0.06247791825307995317 -0.37820731352406611325 -0.11632545777832530021 -0.033387610846700110923;0.19229071724990115522 -0.17698747629120581415 0.20400980308847863065 -0.058960184143606521323 0.1798134394615502063 0.22238777854131669409 0.070288550617726694436 -0.24784657753732439178 -0.092619332386501385113 0.19475954335495915259 0.068373959155245730979 0.38802213461751677848 0.099116598092011384824 -0.050859564541468782473 -0.46044690020342571302 0.021151411136650227762 0.048339834404788444011 -0.49934042774623665917 0.33009077407255743575 -0.029807041830800129478 -0.73930833160212083843 -0.12101971602699838904 -0.41440561891210148104 -0.1110225418313714324 -0.9052460299422079304 0.13536996708617182739 0.21366426723948361843 0.47182265068443846046 -0.32807991944288089625 -0.30558218204760956294];

% Layer 3
b3 = -0.5634229115422256795;
LW3_2 = [0.64162906372887928974 -0.93425742943207923386 1.2675463811403353542 -0.86428845591290359707 0.5120100553383130082 -0.79435954474733228015 -0.53699431877158865234 -0.33211932280327249867 -0.36478318540413190441 -0.27958698149033339719];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 2.03830178109291;
y1_step1.xoffset = 7.47181919175367;

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
