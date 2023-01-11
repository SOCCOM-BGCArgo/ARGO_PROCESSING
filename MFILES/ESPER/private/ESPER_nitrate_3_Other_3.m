function [Y,Xf,Af] = ESPER_nitrate_3_Other_3(X,~,~)
%ESPER_NITRATE_3_OTHER_3 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:33.
% 
% [Y] = ESPER_nitrate_3_Other_3(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999998793562952;-77.961;0;29.715;-2.18999082656054;-138.904784718693;-2.02];
x1_step1.gain = [1.00003298043486;1.00000060321889;0.0138606862425759;0.000266418009857466;0.193068829037552;0.0582073296585987;0.0042772329208737;0.00819772922900357];
x1_step1.ymin = -1;

% Layer 1
b1 = [-2.742329328597930882;-4.7237953136414683186;-0.92108430262812468481;-1.1702455608775217399;-2.0603400037140380263;3.3932090163746271294;-0.75665543300726145759;0.59817670289784907744;-0.76132377151724128783;2.4061792915181268171;-4.0828672615503496246;0.70660369337751838259;-0.65764674378435605817;-0.37959734075490941807;0.41413928327992538536;-0.4186925315183640528;-0.89947147153301887368;-1.6314208790570285945;0.16469973659695430479;-1.5668189925958893038;0.19596263530115010099;-3.3134117496815029824;-1.8891115176329886793;-1.876370460668040252;-1.9168592527698549155];
IW1_1 = [0.37037446091101711154 -0.1662183111422418369 -2.1172736668038663055 -0.91761132356799146592 -0.90356181290911896653 0.79094140270778579005 0.11792682892803731109 -0.32348518229087569464;-0.38845606249496034135 -1.3642479970332375139 0.24668953750943994896 -3.2315211347669232644 -2.3295572861563820233 0.0088912266283339008721 -0.076350194568209114876 0.16689390729192143947;-0.17326516821473486596 0.16696142318987317044 0.88683943667351106122 0.30338726627327450291 2.0027795132355974417 0.17192808668668516092 0.19036050185860611816 -1.2327538865431395543;0.38250441303494570056 -0.11972233310529160899 -1.8419836309751145276 -0.66063789355577495765 0.80996472852941070464 0.10152529477702204597 -0.65599564414651700961 0.22878171785556541962;0.43251601860027655277 -0.57977645674084288263 -2.9330467722257429486 -1.4509027747808846964 -1.3877210526760825893 1.7810208936671720892 0.6530134383069631232 0.8277363561885092702;-0.19204597299047315717 -0.62561430181770727454 -2.7713135646133628676 -0.58720108872312115622 1.8810634616371770811 0.13346142960747892081 -0.31801183887650941173 1.1246353173668759418;0.23321038744561278788 -0.46174586565739184074 -0.46669427972044402742 0.062954808456961636187 1.1918272369319669046 0.7082799558477055335 0.28121136105233396352 -1.0203434278571208438;0.092309365455624195484 0.13436852073691935017 -0.46051474002792941453 0.81039616396262859066 0.35138295899757138185 -1.3702054808254813256 0.12005439979191825928 1.0061406283302043452;0.10612540682278413795 -0.284771221659301077 0.57683198459491435983 -0.10649408164827024736 2.5232699524147719217 1.6696395019560636985 2.6278606408044664633 0.034711673955315935913;-0.51238443027820301801 0.14274287703547189343 -1.0094044843170926207 -0.085545488584762754969 2.6645428484223345755 1.4233164330965433209 0.26397672855455206475 0.48025379774504289854;-0.38908977823570561894 -0.26973972538143087263 1.1492883458114402373 -0.62866777674618568028 -1.7512329717563954734 1.9739011111803175247 -1.2555808406246438125 -1.4116859565962140977;-0.83659034818324318206 0.59847529813883315608 -0.63346823929779594309 0.13968399112962570641 2.0889583683311045981 -0.47649811657892338079 -0.17491620297517299076 0.059904861152614237441;0.50100812477915057652 -0.29123667560743360383 -5.6224942450903672153 0.7152673242390925612 1.329904735408507177 1.2456672330864249076 0.18047096105076188 -2.1358782295508467186;0.79182804640183535039 -0.13541967388999370803 -0.14917342827223634139 -0.31328241602305023639 1.0783958925751579994 0.75763290100363434032 0.20532868395766643022 -1.056451618878984311;-0.5439692039647646693 0.10000358547207782722 0.99341869357620860015 0.075084011892142199507 -2.0766232301060894017 -0.5347412999939378242 0.15583350874379647832 0.65678375317170301084;0.42102413575728814266 -0.41539460739611389073 -0.17025781542152987336 -0.78501199187115455747 -0.3866540548275814726 -0.36248404739694406462 0.07199339899698641021 0.67642364740917027;0.26171620613285839774 -0.30216193749367248333 -0.72224988314836169323 0.36015900998981392744 -2.4536010543973225673 0.39127442721702299089 3.4562966004884061455 -3.4640548258399519455;-0.34069432104750291268 -0.01730090811522354971 0.66172710752535346401 -1.0466517640405936795 0.423567214637940892 1.0307016513173892225 -0.86284929663826748225 0.17633835166183584042;1.1767265184235218278 0.37459142235505060725 0.30542449282049055626 0.52217163395933996473 -2.9578532220909798944 1.5916111887666604119 0.75845276424730889708 -1.0113619600788765762;-0.23695337363218865745 0.43191113873914688259 1.370049754439550016 0.46544803312529131611 1.1029405402382237433 0.42621985235495524202 -0.17676417999414018034 -0.80974045813496853263;-0.34464768633873837933 0.3907550806228475837 -1.818408035583541027 1.0305392824168702059 0.89931863064957751064 -0.023803196731943230408 0.40152818095730224979 -0.10539582353417545346;-0.037855392974921409088 -0.45079470081005529325 -1.2696391092289061042 -0.78450130006555529594 1.9323829638532474817 -1.8678938127371860212 0.81830581920419809272 -0.96596438545329088576;0.1289802673106713593 -0.18896993303729522906 -0.53779680613875457951 -0.71333972819902646822 -3.139254425515348057 -1.435641732130589876 -0.10773600528972956947 0.43740174430666478589;0.084556691069634426494 -0.036661978772897951817 -0.079619746930639279103 0.3444307683705103984 1.0383343888004432642 0.028686705080620114594 1.6344034115118224459 -0.61953447791653215582;-0.16575279595051883863 -0.017472113150363070139 0.13312478872425012599 -0.93090554186198415376 0.39201625179851229408 0.65675867613323690719 -0.87649004624866466706 -0.66915761447583854693];

% Layer 2
b2 = [1.5261281093603393799;-0.64555152405335003429;-2.7979529212199913957;5.1156189904496374155;-2.2718826247734389412;2.5813739428241726515;1.0756074452453914247;0.88171075475094107254;1.8923694501698165027;-1.9575870386168918724;3.2699945603399354077;3.6783953185091720073;4.5024735397415573956;-2.9209031659853619267;-3.8690789738346489379];
LW2_1 = [0.42494275831039823865 -0.01952525888681702626 -0.87665346060034354814 -0.32962667818555252408 -0.0052788259109663717913 -0.30369161467231842177 -0.34241021137877181157 -0.32051164678781962625 -0.0052278464813274990938 0.89364806967840859642 0.56892145902297652782 0.19664698685066051009 -0.075292202714581146772 0.045942416739487526578 0.026949210771743458553 0.070623680353276208077 0.16050058085103099947 -0.87581063874342635067 0.023088244421979767923 1.6715033779979904249 0.29906413044952284697 -0.064400194373949073978 0.2347506665428248096 1.2266821546731496717 0.98081137331838341709;-0.58437726258517008748 -0.506007382916544457 -0.57563264086489407401 0.29361470553763385771 -0.41204825981280912561 1.5541797291280730064 -0.31335178541627095594 -0.5189429975974721998 -0.071106543816875297792 -0.71777658252232123814 0.44232831469296174376 -0.15635219594846760938 0.22378529273922428966 0.53973836613809234031 0.83104095018104873205 -0.10822183711023183883 -0.1556424823591101847 -0.62903376212387218658 -0.018056034191717348913 0.80529775003011483392 -0.19938374685465995451 0.81307913536586395864 -0.53599337920201117758 0.46373263678096432461 -0.01484669259260973026;-1.8527541499075133036 0.089860075388173032973 0.99787340583124661819 -0.69088734724137002008 0.94011845929305026814 -0.6718103380073615627 -0.6485682622863542246 1.5691411021877068599 -0.40613191749630006555 -0.91607455262188375311 -0.56652137921977308022 0.53004484618381098038 0.06192918941541995792 -0.44654011701738738349 -2.0323395586235935006 -0.11472060703866084064 0.0058867721803312128212 0.72553464048587057444 -0.30117144687221197863 -2.1822377646918722505 -0.43090147218441815502 -0.12158240647782711519 0.01140162672956422131 0.025054214465371118375 0.74383922794127399492;2.6747236552106934404 0.021461084702065154339 0.35221643635265353911 0.3152319759651556752 -0.53650088004812113152 -1.9222098029813934339 0.15375641044259924506 1.2422299636122575706 -0.69071580288616130083 -0.45290868943421624948 0.58686286446993773946 0.79846290760783733553 0.60618379019375412398 0.33138329263738819375 -0.060283820058482727455 0.52917989039041135602 -0.25802107013584579809 0.81503558657295094303 -0.36371802599369368325 -0.76087068911036925822 -0.57232730982423329902 0.78078245828928249939 -1.1018418224853290255 2.0620121317812931672 -0.5396995351257215523;-6.9885257855714035813 -0.16234971677608869833 -2.4408490500807142887 0.8029508832943790253 2.3922777042728564645 0.98387740172565674168 -0.29238861271537142095 -1.1710583396354090624 -0.54932777235639673741 -0.90662242056721831229 0.86608997235803231174 0.85151531332812879782 -0.32794253498305264705 -0.28489017038399749859 -0.67027314284388561649 -0.29081160404628741967 0.38616808812626129965 -1.353609285756150582 -0.10343111925550402264 3.1847688739864845608 -0.60923281626181347015 -1.3433651097799343965 -0.64077228996905699621 4.1935067176762252927 0.96731199570653225717;-2.1508658019434538922 -1.6544690062563556232 -0.59100652235274642177 2.6193188123995443384 0.85401608163929521034 -2.0725495750592508593 2.2216830583874673799 2.3350568070640473906 0.3231875219597489246 -2.0280555355248846006 -0.18620039228008133936 0.25711214659946568517 0.85754599558511868107 -0.82380344154196361028 2.3313557750712061711 -2.6445047104510774894 0.73520834356531916676 -1.4058738011327105433 0.97336496571592889815 5.0020513933911336935 -2.2631276034146727838 -1.5890079346427288254 -0.047122650106257220337 -2.6091069682533940544 -0.36353830891248806179;1.573332524843223057 -0.14011857356667772878 -0.69650742214518102102 0.40164742507440071018 0.12675137349852924618 0.42732928647495999064 0.50223643520265504048 -0.26457433386461243208 -0.14470916012110690563 -1.1219127692837518318 0.71745913853360043788 0.93853894345253852105 -0.79964100874414378861 -0.82277879620648108894 -0.55993360532046776434 -0.0023461613795860616211 0.26636270509001969131 -0.25216317152309158933 -0.32165138974440637964 0.27473063511665218872 -0.37593805951844444513 -1.6586962885249152411 1.217478000086482659 -0.59319460583341165894 1.3280940427452807562;2.7234551158042763497 0.24231311933630270139 1.6981854524504893522 -0.73898720028340625365 -0.878041114677266199 -1.809132697562518155 0.13401081365701547621 0.82441570146187959089 -0.040530324315065438534 0.72159578043051653129 -0.51002779843907486423 -0.10034424500186510387 -0.20683294507562846598 -0.1679053049170158296 -0.5431916514545533925 0.5024820350830278004 -0.12732109290605081453 0.17011710903160207065 0.24149579080479602333 -2.2287223543049119456 0.31016069682615854397 0.52242440897580766368 0.06249675096894363735 -1.6059489084215450294 0.64289100132136423049;0.44304978232748654099 0.093913331793642446099 -0.81028720056657477233 -0.31384010356638386563 0.093981785397825909345 -0.36541543060090619921 -0.28790650490116836968 -0.29737276639102200315 0.14997486320770991863 0.96377892106901441593 0.53249055310149140396 0.12044941211918760993 -0.0123203056570791291 0.11557129667527336891 0.16353204733764903356 0.11678863212271624328 0.11012243247659676826 -0.88463030235190864836 0.0018213194033408022267 1.5811018798727360934 0.46880423569660828464 -0.11807941735488983181 0.17967010277840381183 1.2023398836524004096 0.95560381478773315234;-3.7946876202228594366 0.19302716511317891124 -1.5691164709970841784 1.8198896168075795465 1.2107858460706519832 0.74745356384363981217 0.96287342992628965899 -0.086827305648286529083 -0.02000914115605365029 -0.78031768624524822098 0.93493150952203152304 0.34691063462896243319 0.22603978794898740845 0.70291260338837002752 1.6287900213023436269 -0.47034075047769835498 0.4384211583822315017 0.68123774609022003101 -0.27504072938520413016 2.2318963572635888148 -0.34128038101946434058 -1.3022967808795999911 -0.5162600706156694752 -1.1196310196238563606 -0.89141879609513174909;-0.37581918296995098228 -1.1449560837491528087 4.1642023256744726822 -1.2520748347952912471 -3.0508793878757520801 2.0389286485368987023 -4.6721236309057347924 2.6269219792511888656 1.6416456884658436888 1.4950758178097367512 -0.27132331819389793326 -3.7707070133592224259 9.8409274506745472166 0.41024698179020813171 -2.4497579987988182815 -3.2438146311253928111 -0.11422527730085781095 0.063517988608688388807 -0.49379354193745528256 -1.7736607845597283095 0.68625847722545385565 4.3321811996008570134 1.2786948479241926346 -9.0216728223990791946 0.91797788957464010462;2.9460875870927925035 0.24578439563799600598 0.71538614621060236942 -0.13764340121017673568 -1.1453690371687861305 -0.35914186265644992835 -0.97044032758937626149 -0.59733367763227784852 -0.88239995099255719158 0.34464724621553283201 0.27105704855746315385 -0.37799807449211719756 -0.16788542848903609261 -1.8233890201987623669 -2.4447492220302797783 0.20928435677086162925 -0.20075829024957003721 0.651797943888492326 -0.13212660344799082024 -0.35625348926744304645 -1.2576949855846959458 0.84816356830033279302 -0.23198097506234177656 0.67141715876808705232 -0.79989251307652531331;2.198903311274563066 1.9795945831128598513 0.78482603908426240569 -1.1739407180654597429 0.04049041349259379502 -0.62124479763438378299 -1.2515125492978804722 0.15664663513385884697 -1.1718438556211836499 0.36393154431133595272 -0.50353717996839697513 0.02854020127041131949 -0.10572515102441898405 1.2321990482763280905 -0.04228496568466331984 0.59393430563536475653 -0.51463416025582775504 5.6581562064657351385 -1.0213802735949055744 0.85879802128158100327 -1.2744597557455239478 -0.81026358917437613094 -0.8874794627825032256 3.1808350446971069836 -5.7883008421648929342;-0.83052266292548093674 0.1087301566797147373 -0.076409245622039462242 0.47945324873865236004 1.646186193532495734 2.7714037684994501376 -0.97991239994355894272 -1.7152265923402956194 0.61801401378657006003 1.5480284236009707044 -0.58664814736301118536 -0.03250735062568742284 0.19745656712810263911 -0.38017023683393785793 0.46822293631379247048 -0.72665328648246485521 0.44401320718821951639 -2.9040084821864136089 -0.64029124687693139695 1.7340700935458364107 -0.86627083622551293018 -1.8468993108868330122 0.79129002137198500844 1.783211888006749346 -0.2844775953729000717;-2.5238700070074684412 1.1138784307115410233 -2.3678329257441594891 5.9167608849080775002 2.2848638483395116161 -1.1733189567389941121 6.4950359289831594012 7.9422561106643367879 -2.5171916411963000115 0.3835467840154255792 3.1629017850284131264 1.2807703513601043088 -1.3196171571873245298 1.0193205955929733886 2.0471284981430004635 -7.708229922549984181 0.16436523522535398389 -5.7771464670860890323 1.2367622383478680259 2.5594447748633744055 -2.8393817756249544182 -1.4406112866627822644 3.2926801802310561307 -1.6862874040573998524 5.7637178727568967318];

% Layer 3
b3 = -1.7790493074054722911;
LW3_2 = [4.0852633821516510793 -1.4354495767337713819 -0.79138144805101473533 0.39273908643601623814 -0.955647988198920606 -0.15085897807406414839 -0.31901490956427397405 -1.7414431811680364248 -3.9674549228314575089 -0.88916516363123110356 0.6101391083664955417 -0.41227551937113371405 0.21144370051758162865 -0.15161557567402428326 -1.4804391365638243094];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.0415843642790311;
y1_step1.xoffset = -0.9;

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