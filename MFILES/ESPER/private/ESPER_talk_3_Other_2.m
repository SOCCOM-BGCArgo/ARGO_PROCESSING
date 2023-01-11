function [Y,Xf,Af] = ESPER_talk_3_Other_2(X,~,~)
%ESPER_TALK_3_OTHER_2 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 12-Jan-2021 17:03:44.
% 
% [Y] = ESPER_talk_3_Other_2(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999540494077341;-75.751;0;29.715;-2.05132667284751;-133.803853032615;-0.1386];
x1_step1.gain = [1.00000001035699;1.00022983623563;0.0140762793578401;0.000302617642608564;0.251952632905014;0.0584431850247228;0.00432440760873659;0.00826143816767911];
x1_step1.ymin = -1;

% Layer 1
b1 = [-1.3299480779813813;-0.69465651205894341;1.7427119577212438;-0.15915495199132182;-0.3897384887196656;2.842927311304456;-0.73916832421375689;0.64680549644631502;-1.2857831516460467;-0.48491844050687;0.60655396733054368;0.13276010337410005;0.71790375556740249;-0.1159181874340663;-0.82925265987906704;1.3114067983787743;-0.079332447937908737;-2.3421896699626958;-0.9086534447582072;-1.9520353080373336];
IW1_1 = [-0.077330322149114686 -0.13602129450620626 -1.1237418877852592 0.12418143561132644 -0.41734804816018878 0.48593205028071113 -0.90874021613039047 0.18452495439940861;-0.089070267326826671 -0.29964234107141263 -0.47865464851964112 -0.37016919961795525 -1.0472911292044846 0.60106172729904017 0.40272282125700204 -0.0061758399500992933;-0.57242214469093544 -0.045679000753777128 -0.05779684604220342 -0.019198089377599127 -0.19716149064105046 -0.005391536346958041 0.21539657789404917 -0.069612473712322945;0.020111381698929346 -0.17052078837565007 -0.17656620974980375 -0.23122234256415181 -1.377436780502737 0.30889615872418896 -0.02026588837755457 0.4192984993629883;0.3034912194327653 -0.34530714935777179 -0.25489097689375317 0.44378810453596484 0.35032403271270735 0.02340081924386625 0.041400156813286286 0.7447658249257193;0.43193802695579703 -0.79744727206234223 -0.76714986101341065 -0.18128598160957818 -1.2468411598069551 2.0028132674701258 1.0718208837337688 0.51992845901343754;0.77636680024831872 -0.14433924859695454 0.55842350351064596 -0.54442285842654259 0.25334463760709552 -0.06117799819690687 0.10475721995338663 -0.42125396228832934;-0.21011212393420323 0.2525216133255937 0.15672486896889626 -0.23652136632805487 0.60321465532318463 -0.27991528170194452 1.1864825031855324 -0.4413471890689098;-0.47603711853835112 0.15423789132864063 0.62230730738128237 0.42424992396850841 2.5575873097970288 0.081777888701935211 1.0384986182476215 0.68529696592764044;0.31714652851381026 0.018356434397381684 -0.11291239383329058 -0.90315193116111248 2.8176530619174391 -0.21709582979955269 0.82507837047197718 0.18245775859986502;-0.076288834679844847 0.091885813044221984 0.23126987167555976 0.037748161600116244 -0.44625999881200362 -0.38898426706702544 -0.37168944389996683 0.64620005039041784;0.027499746144948607 0.0090419501805099559 -0.29157687146911448 0.095304767176040653 -0.92918669094157791 0.031187136408184407 -0.30408527815296621 0.52131176975607496;-0.26573144624542733 0.086371819143220604 0.010048641783124783 0.13962365095382925 -0.52097315539644573 -0.16319029080781641 0.054535989964105532 0.014288833232545999;-0.31395710410967725 0.17555069530823325 -0.54591481754654536 0.28351171416438231 0.78647108063696236 -0.11460018342083013 -0.58673454994901797 0.049306781849888022;-0.40493887086317398 0.078185768176367798 1.2066308229126859 -0.79113042539504996 0.43119929570746357 1.1275293632990613 0.69646702056321486 0.86554417769313463;0.059926182352650181 0.60429470186259582 -0.56828325733287022 -0.52327359109611848 0.19348355319899524 -0.46776909341010658 0.3081977828582419 0.26154519092585599;0.39801806836354281 0.017044643827798377 0.82647639329398837 -0.29574854579475945 2.8081738497247839 -0.47136726005498408 -0.25100772216630779 0.16403681939797279;-1.1038347238885537 0.3649867096163813 -0.89070392681894206 -0.94845404474706907 2.4676231696226254 -0.53852592129182442 0.57235259474778899 0.29967719003166565;-0.23386420524639506 0.50852180653474199 0.20764948903109384 -1.4319532220936453 -1.141278043227153 0.53592763388173281 -1.1566495310214757 0.081216600810583059;-0.72531168066010776 0.72773303582513038 0.48604407402531319 -0.28080274041778258 -1.5273564680169784 -0.13675796338210783 0.56680154767353763 -1.3444601308060606];

% Layer 2
b2 = [17.616113831411049;-4.0530387191908739;0.22658711755385577;5.4477711077573572;14.894852067135451;4.5533976535468588;-6.5397445349410592;-12.27526950114018;-38.203275679646886;8.9561632365651107;8.0437283240085868;0.91140322913141769;-1.0652390596257135;15.606615137000881;3.2443189257285558;2.6970186665664762;-6.5700270946260959;10.166513059434861;43.866424945343844;-13.256545954384901];
LW2_1 = [2.5448435075952736 18.127828670370782 19.785352148637919 2.6114134563384153 0.49253708371147881 -10.910684718534441 -9.9513566438859868 -48.736168606981394 -27.888106686930449 -2.0456523865451408 20.549736830691714 1.1022589932054698 -16.593526104381123 -6.8246940424769083 11.648252167988097 21.817799821658355 19.811713592282171 12.71007554771451 -11.698483875535697 11.916403292825695;-12.544503934992454 5.825069762489437 -6.1853160046449434 -14.542228990219327 2.9327635234043976 -3.6002838385180937 0.31236177442089791 -17.399277930739746 -6.6445103435415014 2.5557464083673347 -0.49187594720695343 5.8859681620825741 4.7854016618562873 -3.1057836096730433 1.8654366543071668 9.1725921872959528 -5.4962773304746353 1.6391759056607464 2.9490647007408173 -4.9251060395967512;13.65680356308258 -5.2497843235718245 12.368559291751373 -9.5266735577015247 0.70940301231871783 -2.1273891274277981 2.0218281834523641 -2.2241723355849339 -0.28905823145409254 -5.4491059442264893 -5.2000279380413836 6.7147040004482932 6.6179145490908038 -8.4081564764166465 -0.11827234309359556 2.2061348892444084 -5.1705354212650034 -0.34734868847981293 1.134913675388409 -3.5471569371222156;6.1419433730848043 -3.9623526823735951 5.7577117902842163 7.3097567501606395 -0.97172467018646325 -1.9823733911469326 2.9554571484354688 7.5123479224585896 -3.0531839502813969 -0.96785391013458244 -2.7395738301804955 1.3073641379607528 5.9673680446233393 4.6443754754715094 0.21175503479382299 -15.900109680282979 3.5923009588530928 0.35467547363219321 -3.4057933541521792 1.1012038597050171;4.8344750103242795 -2.0953377766179977 -13.508844366808018 0.25553027277246537 -1.9465174002240546 0.71916446400563794 9.3025409947089059 7.5378648307480223 -3.9729670327766082 -0.011120641597323094 1.8111844055812976 4.1840980504683198 16.520804078686464 5.2265730137776272 3.4577161748869583 -12.52609868106439 -1.7027195038480412 1.6731147017676149 -3.5420009696813017 0.31107947260026225;0.91518559295170132 -0.79816437996712775 -3.6357148756993718 6.441033525617966 -3.582525102388876 -0.055735826368624272 -0.77330658358006066 0.99042263427910016 2.4195948034757331 3.0328287253470561 4.1875668531703401 -6.9878289699511411 -7.9772022869299137 1.0521136713449044 -1.6455233990478015 0.3215008885363061 -3.9804113403812811 -1.5471178721999308 1.5323374389536353 0.36314546615048826;-3.5226087803965282 1.6946549320192525 4.2072234106844562 -0.81461504205787805 0.97637460000212661 -0.090755036650066248 2.3407179648967458 -2.3746708969535653 -0.50668474641258909 -0.24398598494603019 0.83622518993444506 -0.43750496384657372 3.2411073569599349 1.613466859914553 0.52744365970968132 1.9047331091628998 -0.45687444986745002 0.3168092252987893 0.18840871184654856 -0.36671278944236058;4.7573736308334125 -8.0603823196291096 7.7327429564795773 0.62845973238755015 4.2987174165210833 0.56368940164518555 3.0343770090557372 -2.0492813694593011 -3.4947095497613638 -6.9179151464994399 -2.0416722179101092 -0.6893771546680213 6.8752607925974303 0.11063379044090829 4.8592853159475684 10.770256752857334 0.34384352847164829 1.0103798242546411 0.32321020478994528 -2.2553099622269492;0.56724683933639508 1.2114739511301955 32.804959743431056 2.8405967201861602 -0.19387404117927226 2.8292695388361051 4.325194766856832 4.9497572314502856 2.4229396579965123 -1.6526428085773128 -4.3814580038802449 14.398367451427728 -0.0060839162949900705 -0.034634749290702316 -1.972568927486799 11.421653444827452 5.7758709305510241 -0.4734926757486187 0.038927498496756197 2.513619049981473;2.4742162833245684 -8.9755049588627784 -5.1612090552538064 9.2477043952840816 -7.3940608868321105 -1.1109621591687426 -3.5058201474457462 -6.6522806338291245 2.3581030766232187 -0.11685560761001627 -2.4832463635502524 -4.8733965377012609 -16.889407609422555 -1.4101987145721899 0.01878933413488584 11.29641347513674 -2.4853595986753327 -1.7729838186676175 -5.4146354663285807 3.7960974863016053;2.8350397029069678 -2.7363258169937841 -6.8647649893044642 2.9883905779488775 -0.9494288378933885 -0.11077493266147134 -2.4001359937581395 1.594079897043194 0.16258198547794275 0.13494203287073861 -2.8794486288178671 0.90822136212745419 -0.1947916220352722 -1.7495379835602478 -0.565526971428717 -1.3988650469913753 0.96034841052951714 -0.32123602160885584 -0.28298221133526158 0.50539566381037881;-1.6252391464606077 -3.4539621741080819 3.7585617101267443 2.9104605018460625 -1.0035762262415988 -1.1415653506997259 2.8008820866746293 0.77830659838035254 -1.0252504690284527 -0.099205824680948004 1.4990009547452245 1.0400584561719202 1.9297736508827243 -0.0043088716759383822 2.8215770604802142 -5.8330320160549558 -0.21557104116786963 0.78079558411810845 -1.1630837132886631 -0.04591694248284673;1.3490257105865764 -11.224151180782235 3.9329297155300487 8.5825091742706956 2.8977246357417838 -1.1601152744661303 2.0689464142776646 1.3547792387987008 -0.40836152096151196 1.692016878657302 -9.6692038030120226 4.5200215505233459 13.217596999105732 -0.60399393651134536 -0.46951618430382758 -5.5816012876671932 -1.5097717699346902 0.45003125570062669 0.52344751981724413 -1.1739593170585259;14.23959206488488 30.318196959626754 -25.783185051316821 -17.030259165829964 -42.802428502155003 -15.432590390011699 20.891563582002458 -14.398892287078887 -37.886704189402131 35.934515407092576 -31.40986225175498 14.625343767309412 2.8145178113277107 -9.0897486362026907 -17.065778250651324 -8.9707180957446671 4.5952475081926085 -0.7811551272230588 -5.0299862352950173 31.368449315028677;7.8943623197197343 5.2611232672407526 0.32100499325069765 -1.7109438348586947 1.7194976227960341 3.1738599753681034 5.7770992242998425 7.3158866429582403 -1.4349061721149259 -2.2971211602027468 6.7294559950701798 0.09826770131794553 2.6797140978532545 4.9834437493772903 0.51765772293133583 -1.1784712396972956 -0.012712332021023041 -0.6941240219870658 -0.30103807392383286 1.3270535050487957;1.1201193173269857 3.9022046879273775 1.3213532059868374 -6.9656042693134603 -2.68432954233689 -0.94704226044312712 -4.5933181385616244 -0.89677629200352782 1.8688246566490003 0.87929524114966007 -2.8051821001096138 8.7039259725011302 -4.1305700143061443 -6.3363263712745974 -0.71164968059497247 1.021500786492767 0.039858826705531619 0.52793742706139946 2.4935567463071346 -1.3760665477925937;-1.054154473794958 1.0212289348821948 -1.241601362568719 -2.1819182315523276 -0.92986792316223899 1.9187121173990032 -0.65208202907172919 0.42533016578524463 0.93210397626107233 -0.19801731256130584 -1.2987746205703938 1.6616631276674827 2.5705008706453296 -0.13916715286256159 -1.9183304117171713 1.7409152448134964 -0.42348685029136496 -0.31822351695040862 1.061275032982135 -0.79714529082147523;7.1476676944582795 -24.541012043530394 45.30033361140903 31.573863902106801 19.119595730488165 -19.639992593673981 -9.7096623924514205 -31.981115499027101 -26.353262126638896 -14.431240602820299 -24.232349159088265 -31.784620927466676 -1.1369232800813993 4.4526219371430535 25.061821979122126 15.672694930278007 -6.7170762443395198 8.0252253026743698 -10.362767648026285 -31.99927557529066;1.8144113406929321 11.526366970306727 -12.548194539345497 -4.6549749094934558 -0.96056829075364625 -1.6277908104112702 -0.14915602540171619 13.179357535081847 4.3068083934897361 4.7066473397898561 2.9841540506893498 3.5606961652482285 -8.3924844954044957 3.4708315394801148 0.88605883405142083 -29.060437879501446 -3.8470594904775095 -0.96479418415097085 3.1603620419934715 1.464708048271594;1.476180316054877 -6.689778166896553 -0.99695081165215271 4.9403541711888055 2.7058420058886958 1.2377408075176284 2.2804065198001759 -1.9148444621415737 0.98462442653801663 2.5072718061930241 -4.0010195789456668 0.25347579404873183 13.189807425732212 1.0208970040474281 -1.0185535968905941 8.7002544194026488 -1.9050535458019762 -0.99364312012662104 0.68399641927526766 -0.85452484340734558];

% Layer 3
b3 = 0.15961893026606722;
LW3_2 = [-0.0045426063482539823 -0.044939655086641772 -0.02200880578081563 -0.064245152647864484 -0.029808667634031193 -0.10330416074449326 -0.76090154267596632 -0.022620262245693183 -0.03428615379102222 0.020168748295152025 -0.75128854323114602 -0.060924858059808323 0.12261387365997432 -0.0015963672413407935 0.058842409135116479 0.066579006231266738 -0.15994870728162172 -0.0059165943603555881 -0.026719876650547297 0.095408277815759751];

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
