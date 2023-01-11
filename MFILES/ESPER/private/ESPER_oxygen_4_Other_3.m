function [Y,Xf,Af] = ESPER_oxygen_4_Other_3(X,~,~)
%ESPER_OXYGEN_4_OTHER_3 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:45.
% 
% [Y] = ESPER_oxygen_4_Other_3(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999998793562952;-77.961;0;29.715;-2.18999082656054;-2.02];
x1_step1.gain = [1.00003298043486;1.00000060321889;0.0138606862425759;0.000266418009857466;0.193068829037552;0.0582073296585987;0.00819772922900357];
x1_step1.ymin = -1;

% Layer 1
b1 = [-0.030959199332437629254;0.33051990323106611447;-1.2423486910836931685;2.2310007345215252528;1.0451697725489859092;-11.482388125996635608;3.7796002518968307982;1.9617762315697748754;0.92717738170727104219;1.8676819552713601613;-2.1312305867344858434;1.1855880058296566215;-21.766280954510268231;-2.8909003615994675052;-0.024149833487026872714;-0.5145672063429224119;-1.925384912193339737;-2.0194546387337721782;-0.35024808059234491253;-2.4426813716243951369;5.2255189817431135779;-1.6846265461752143455;-28.372479699973954581;-3.6553475227587184371;3.6937813329603161883];
IW1_1 = [-0.041459024865292126716 0.53914791037897591419 -1.2348575868305236014 0.30796360398066652264 -2.2794022912789384527 2.0142440530737437676 -0.87498233019028204271;0.33529945795588073088 0.46837873741371421898 -1.1628800068237430398 0.053068715334511608328 2.9099656077432531909 -1.3154239432107286412 0.51680882684937590188;0.025994705057747587856 0.46206429485632233289 0.061906945596731231818 -0.48817950704607271373 1.3901923184311304471 -0.40912746219993068042 -1.8414010796495137168;0.079077610842517284673 0.48449035923385586067 0.82366418523110207506 0.36459559340576391806 0.10980748629388628612 2.1612201216778714752 0.37549050099255892876;-0.57938751222427864462 -0.27416249633363193317 0.33836105507474750764 -0.83215464260020854947 -2.6095588864800505213 0.91483216620998308866 1.3844227177683858798;1.2899381380273575726 -1.4683249098845774228 0.60355877446929351304 -11.487930343930116805 1.6008006864579535922 -0.05080320835729962986 -0.62556729068165572993;-1.8469865116732540056 0.24398864930468150414 -2.1802684233860643026 2.4797698280541022875 -1.3667386183632908025 3.429274527479081236 2.18867284340683943;-0.9913260926156116648 -0.026990551937862750781 -1.1251939486961959336 1.6254793517760899313 -1.6617168024165143159 1.273409226915057868 1.3382556910528164984;-0.028509367645030328497 0.48552626489351863315 0.050444488639992020318 -0.27264694054516497301 0.21998468628015729887 -0.3324077206149457786 1.8210835606988069468;-0.12690755474711590267 -0.021737338189076870626 0.024045287777342987112 -0.076603291823842525887 -3.0311578695080032553 2.1560891407813276288 0.047982376928752464018;-0.39497398209891282361 -0.51229290746002698054 0.034726566471547384851 -1.0441590419739716911 -0.71271466131343663708 1.371612279837044257 -1.003254491279185201;-0.22172080238960523157 -0.11553670620460317064 -0.023549711273322455896 -0.2644107530338472567 1.2563865672700731402 0.73350322951071078581 1.5328522214603614859;0.22829192176311111084 -0.36027025935091805398 -0.32001015036957602966 -14.522044475275544428 -0.0095846125710279082843 2.2535941288287362205 -6.0416935532172209022;-0.33184972123782446873 0.001494882666095828129 -0.14004361650557300978 -0.26883858438376584399 1.2640395920889007098 0.7056605422370401115 -2.6846969932980342399;0.96190381884743825225 0.10849039832334643474 -1.6848261015958130304 0.34420502717517920122 4.6155125092227065053 -0.91075818986963397172 0.60755304285771971617;0.86002815380468400352 -1.4469204874664878258 -1.4548619535232152344 -0.95405804567129903759 -1.2402028325264715658 -1.1580261908773257229 -0.35707266596299858108;-1.0408619752613170739 -0.52432623335708927748 1.4817277230238805608 -0.15305717527924614374 1.6266102274976528275 -1.0089365942123909026 -0.78489975663870437117;0.21635957513179976397 0.24923721588753117162 1.0659171299902894248 -1.2655165527263838587 3.3450613092316969777 -1.8288134051141271552 -2.1503300370355695748;0.016031994379828766462 -0.29031136337452251306 -0.34898093240135069637 0.42366296156399890771 1.5820251242080738407 0.85458399175623922517 -2.0332130052130281683;-0.33158324552806084373 1.1633234296585532874 1.6673202524242687339 0.58129497651066686092 0.28957435406611820028 0.92649556305624514163 -1.1115023233439427219;-0.55317808533989654141 -1.4304199099061047917 -4.5770389678066125683 2.1616697715132899837 -2.3250931177494975621 0.84951393578964196074 -1.9534375026116626017;-0.3461355650409333462 0.45102914107086872209 -2.209868451372192677 0.26868524670314480884 -0.55259802029660864964 -0.61405639140154011368 -2.1175414285176548113;0.10437592931232358984 0.0034095660526980171162 -0.13997613645749948375 -29.642083709004534597 -2.0806529259993240544 1.7118346589451596262 0.8607889425321534338;-1.7324348432285106014 -1.0485713990098843862 -1.9861642710501017373 0.13990561021710354317 2.6471595154798763794 -1.330109448915583803 -3.0844388778542413121;-0.3306158294127641839 -0.21512368800671471947 0.36665441300486317777 -0.0052271846623195786924 -3.4498516049324470245 3.5956527281261752016 -0.086947364320956135875];

% Layer 2
b2 = [15.427848197884742376;4.62321708952191468;0.77179747668350096035;0.34936147902357483908;-2.2699310118323667851;-3.2472154349850681498;1.5153988352564513598;-5.8542499749908252227;13.156197471776435393;-3.851940328704918759;5.6644638129234952117;5.946703770222530494;12.089465532693116856;2.1918438403146511817;-3.6681646586098026752];
LW2_1 = [0.46747390799449567433 1.4963554536387624339 -5.4447409658110244024 0.054423094371931796109 5.1848682499266960377 -0.21184931908436646641 -0.10432881201921687109 0.20310719069650398638 -1.1469936617042519877 -2.1659094209655802565 -2.1513569196296664643 -1.4062840590836358512 1.0824135308013487222 1.4554839663581760068 2.3531186276487598974 -0.97182303173877082614 -1.1963658195080628666 1.778897699972928903 1.4366069027318966178 0.062556188358588449283 -1.0975639591276182827 0.30487177255249192642 -0.38713442734778796828 1.8565028768374984836 -8.7918347747580636309;0.93907834713322047016 1.3737889523322772689 -1.866117856204894121 -2.5381106818898842015 -3.2839385014246920136 0.42411962952690546169 0.27532297695413476823 -1.141726169238208799 0.55943268390889822061 2.2666666818826408658 1.3772107236479667058 -2.3223487295124534135 -1.162531995878939739 3.4752393273972188581 -0.57919102398528143905 -0.36012471281180857829 3.2233151675292859295 -0.13182735121014144131 -3.1585372364460795325 1.1891201954687746767 -1.4129160267271714257 1.2900166499175780288 2.7443373437487879052 0.9053372152898300218 3.1721651438394045108;-6.0651195641418524929 1.9373559121675205663 -3.8355196786864986791 3.6044325225345490438 0.072189981505032418951 0.77651890565717951453 -2.20859825911352603 4.4436559873817538246 3.7600777969166565562 11.563008309483791436 2.1047155394048546739 -2.3633704191520230609 0.71725275128758747556 5.4903130089408289649 0.75406123544591985297 0.18506723527480128944 -1.0167096562956960959 2.7058377438719336538 1.3823761961597382797 4.8231127316144153383 -0.64257951152151626673 1.3408982486378950849 -1.0750900852390794515 1.0762642636697135057 -3.4815388981327552997;0.71463495943024624157 -0.07549754912293248077 1.6826170371280970706 0.16777464262458660627 0.022936466971688260208 0.026893453185689360913 -0.66840593188733510299 0.26776242766274538365 -0.69157920584485954585 -7.5824118567952449865 -0.41312004372277943975 2.4339570579244806936 -2.8829375114084134779 -1.2605461773498920497 -0.089750799322586269358 0.26590776885452388489 -0.16394413952573647086 0.44876409285398943805 -1.6248783041998151599 2.2644456646792447962 0.40192527156568497171 1.8179455200332326203 0.098710479607908463429 0.16968150155077701968 2.9546979180142964871;-1.5165096618063387268 -3.1046000879874013911 -1.5023382018969033691 0.43856372076895433487 1.6037839908149849055 -0.039333306512500797181 -0.24241284376046587368 1.2184641288379862711 3.1638885063406814169 4.679162300885872483 -1.3300550686362835329 -3.3248584959705858033 0.08417469930139347134 3.7724623236133596471 1.573317451180498372 0.37217609769623483107 -2.1619840691876901495 1.6748878825106381374 2.6306713715716982804 -1.5818672169879253175 -0.058493302178180252293 -0.56065024498300797262 -1.3604112992536792071 -0.3659303706354020691 -5.1751831629254407829;0.90219737338143113359 -1.9801473391336144037 -0.0099316480354864778574 -0.29165144117513136512 -2.5528938995167509418 -0.10760550112447388182 0.98505243873377035246 -1.8829006850149130337 -1.0832785474349486687 -2.5096242983448164843 0.20373037648687417067 3.9339152890258630713 1.8238444109346909627 -1.3829881967997541548 -0.83517035991936738881 0.18289418025606724028 -2.1711704507229874217 -0.11556215702563582548 -2.7315915944024053985 -2.1528602110263030944 -0.15663989414941126288 1.9148712958455040845 0.27729532730647465932 -0.045845144122893005667 4.3169105566305026045;-5.941699596962670249 -0.99902474452509593128 -4.7028441873379431826 3.9124350104273881712 1.3780562453873466122 0.45921842735337148333 -2.2395465255822188944 4.8530865180135380754 -0.79502550996209875311 -7.1006103252129042502 2.5574309567441080304 -0.18153567019335156707 6.2274986479390506133 -3.0644871452369661924 2.3086167358639504421 2.4323588632027850487 -0.24224616090342807051 4.4847769560079431272 2.534697206717742457 5.3399946152706565528 -0.65672109772936382299 3.7531577095752521878 0.043150458363720493871 -2.2076921776889815163 6.989379943917331417;1.4557417667923955307 -2.5148358618101984696 1.8900446167300839484 -3.3271092674028195368 -0.89846412637148354019 -5.0628228248605209316 0.87463344585426716637 -0.10417102709158648777 4.9259814385012976601 6.1063780940520384632 0.74900971160834073359 0.18078815002922646316 -0.13113302654335762254 0.62589476465976523478 0.97161302725880982489 -1.5703189363009035784 -0.90568195052839750581 -0.28552122120924861015 -0.40311495068763025795 -1.7709149794128808963 -1.062082307512778323 -2.8591765802450579415 -0.022428470898044162302 0.82928322462072356913 -5.302061111244643854;-0.95566691788180357925 1.6882830632554870487 -2.7861997452177238443 -0.28404710947686928613 1.9829074464749063544 0.29005945436140301297 -1.3144965813278430833 3.6251697659727781797 4.3572187255314505094 5.1951353390920260722 5.7720130527073756355 -3.3492994158944831007 7.0368232319806818253 -1.5074177120982197309 0.75298024125816365171 -0.19319237262774668906 0.075864725909105368284 3.1537032260246373738 4.8974472683891221436 4.2281328728829148389 0.76946925664341181417 -0.62120400896389205148 -0.72177266759757530146 -0.69482560009145888369 -5.6993683220974746817;-0.77648358265212880092 3.6426118518820942427 0.92027753490877306497 -0.13269512994397722472 -4.6340883542846045273 -0.34839728931429081449 0.074721144891098545404 -2.2179062323918494748 -4.3672449777880997246 6.7376438309672037619 3.3606408960337748049 6.2505787368597935938 -1.6340927361790802408 -8.7089884538128323754 -2.8101675882804015849 -0.69102694124612495941 -1.8771060669537746879 -1.2750612440648341206 0.64782235359422091214 1.5874239800843528148 -3.1374189345212672997 -0.59226272534022073035 0.65907077159674387623 -1.4295815153222577454 -1.4751442430258137062;-1.4701359522558989923 3.4178502440369031845 0.92668264411881029474 -0.99432938605972598367 1.3760921608683038819 -0.1988646545919318942 -0.62581819511874003936 0.52632542644730706538 -1.5692030462339028318 -1.0167985339677130252 -0.076334581946920582185 -4.4859612944853060057 1.570235739807270603 3.6059615014134211997 0.52080277426616317094 -0.58341507448849838635 1.6559863593236532342 -1.2742904603815885523 0.055459729363791787637 0.63429505028585442261 -0.67716307877449022357 -1.10589838079709879 0.96020057968289629358 -0.9163772179895327552 0.64417513965093942119;0.46302042357202805034 4.973622086074069415 -0.24856667593295175833 0.39168131484526974573 2.34082389802503954 0.20196324752681910075 -0.75375077625378572588 1.581259405665720319 -5.2982544950661649708 0.96556386800780746871 0.40857081199361644508 -5.4631190173152788248 -1.2533528130226987418 2.9170282438111412304 0.18139255187884267828 -0.37596260765894956579 2.4311435197236446193 0.2023172715077612982 -5.1834117758907645168 4.2798411520042085598 0.27179313861732024016 -1.6844068483279128756 0.37965798612481971297 0.37273759020587521329 -1.1022447293282349179;-2.0448796148007075679 2.1893090261729879487 1.0467324702240090595 0.31130968116361401599 1.8742360747709296831 0.029132650311903967449 -0.56896830718424351581 1.1940905196879609829 5.6764313056957425019 10.94226964391927126 -0.053539687817714218254 -1.7692602509288182233 2.6939586711370782091 -1.5160638491113245241 0.46985848816959896457 0.1831151584392443632 2.1596364362999644904 -0.60733858866850143343 6.1852657015022689535 4.5832404087023705941 -6.558570386423379972 -1.2484623757723583015 -0.26828172455403853647 -0.69572783189208220467 -8.2637302929161009502;1.0008152888758052868 2.3401092938257463949 -3.8480613262810576813 1.7493706942906372959 1.6129308398981141437 -0.29026361413459494942 1.2460264486963350539 -2.0506835129719376276 -3.8364691743093293219 1.6357778844982411304 -1.2902447684548199458 -1.9895565554758984916 2.3771053548173370906 2.1900445979522951134 -0.9928205704355643757 0.85143124314850049394 0.54616392488127651905 0.51776949724857324053 -2.9260926878358652914 1.2537397236978837789 1.2268675343047943382 2.0259732143520130521 0.17054781534831753298 -0.9072678454413599658 -0.81467345981290861801;-0.70124481525350279565 -0.95027270056701618195 1.6688978711823985179 -1.4676139319595793431 -1.3064678023080622538 0.11293511737509734361 -0.21752407961734487873 0.58418991599763736566 2.6911860198279793366 1.7194714302552098761 2.571189655081197234 1.2315039134466621018 -1.6973724144700821359 -2.9164288356848504513 0.59533646901583692213 -0.76408157912913154952 -0.87004004755276209693 -0.011713811373786249537 1.6205057452380262628 -0.048471132202097880204 -0.62926747096779256374 -0.92845909933194581409 -0.40659914968472626873 0.34457015598775919685 0.50527556292052344666];

% Layer 3
b3 = 0.38959368631394486471;
LW3_2 = [-0.36514721176684522552 0.11627223537186477664 -0.066355236665769506965 0.20710954626558500391 0.22446559300841303908 0.2044056164133255804 -0.068048397539299920744 0.22880116248113097077 0.12442941002918651849 -0.047039446012959347521 0.11484314118620324308 0.15952948880642953133 -0.13026604688714132907 -0.18469234422100841431 -0.44415368173986191636];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.00454442172233583;
y1_step1.xoffset = 0.3;

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
