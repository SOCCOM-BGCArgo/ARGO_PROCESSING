function [Y,Xf,Af] = ESPER_TA_5_Other_3(X,~,~)
%ESPER_TA_5_OTHER_3 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:01.
% 
% [Y] = ESPER_TA_5_Other_3(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999540494077341;-75.751;0;29.715;-2.05132667284751;-0.9;-133.803853032615];
x1_step1.gain = [1.00000001035699;1.00022983623563;0.0140762793578401;0.000302617642608564;0.251952632905014;0.0584431850247228;0.0422119037568594;0.00432440760873659];
x1_step1.ymin = -1;

% Layer 1
b1 = [-0.60331502356586796942;-2.7886483711786840267;-0.80992289160045582985;2.1761513107314800486;-2.0611254535426897228;-1.0478286785638288858;0.93494239581291305186;-2.4864800923435796953;1.0871246335802142724;4.9950129630244139989;-0.008739324530173456923;0.32136692397680172162;-0.92590775386667545366;-2.2572530720812578053;0.41577843292330368685;-1.1107779217282369721;-4.1055399526424221079;-2.8572515558936948388;0.28957193891380295092;1.2895035188620818101;-0.051134315345334808989;1.5918972745810251457;-0.26536904785558118602;3.0199725077865506861;-0.55024121589926766607];
IW1_1 = [-0.0079692248562014782187 0.13447057862672004425 -0.60332325829569410924 -0.47723960079488536579 -0.29030176985423356006 1.1455518034852987252 0.52234229991791314429 -0.36400318711777901459;1.0944832705542724138 0.057265675097556145712 -2.4772654930598680956 0.18989552428952741203 1.3065039231603741143 -0.62768117368190734506 -1.3975854017562840959 -0.697668686092335677;-0.60009879960284417955 -1.0406451226240385211 3.1994638416399121539 -0.57685971251323164832 1.7534583014940061663 -2.8663090528629813214 0.77980252673680483522 -0.94869036356221370099;-0.14317433246935906688 -0.34672492166277563141 1.8847031024053082504 0.17614648488793033243 -1.4777392569713023729 0.64958992931480474819 -0.22887991637257942834 1.1893885063197471297;0.32556124626369348229 -0.0099338145459921811664 0.10028080177542418594 -0.22302388057664515597 2.611664721032604497 -1.212235970727357115 1.0160364307482905311 -0.46834106958059740222;-0.23177339104438385298 -0.022044867000346385355 -0.086264003447523565637 0.23691211558750754906 -0.53737604766582025562 -2.9362512439632109107 -1.0773211718062136732 0.44199665816426303433;0.73400076365962529934 0.31273063579390841804 0.97916774145085405667 0.087998067774474650338 -2.7433426364251882745 1.4401687071242703642 -0.45002775069536865571 0.2116886012998139277;1.4649926962793013985 -0.15737566263272165235 -2.0859014976013776632 0.24399044157706967217 0.42759849945486977374 -1.5092561919477578414 -0.0540090812072210949 -0.30418791656536653401;-1.5233378696785164763 -0.3670956986983093806 2.6934333735510964658 -1.3895726773464309733 0.080986175407011445904 0.78133075959814213274 -0.56285099888154399572 0.23987164517596670787;0.034022315575544935518 -0.12543653008875776944 0.13705161051188372023 -0.85449711040677323037 -5.3949452775318471609 3.9072641244412675832 1.3946488233112384059 -1.9723762126869537337;-0.05047076945660003422 0.17029633743387209699 -0.0011368295149648667639 -0.46935427086060077517 -0.14934156780288004862 -0.01233027138518407069 0.68813813212797192431 0.25126336412765504225;-0.46061663566118604685 0.54807036572949696662 0.42267753404851754473 -0.54190865955114075891 -2.2899998354546733381 -0.38889077195170429491 1.358208133803657347 -2.6982880778982654135;0.50178743448918261549 0.86785867737563981983 -0.53376647799054544041 0.49352640935002017875 1.7876105174702239253 0.048579175552722023146 0.3695074382903979493 -0.17786304259046764598;-0.078669801796543620265 1.0980604419412114403 0.5269746619735955484 1.6052808458943224057 0.44631764751229058286 -3.7615887921541664163 -2.4114316656253209636 2.9362885648222785129;0.15511159628509524344 0.43137432937012831768 0.082181127464284159378 -0.17485022745423181578 -0.18211594157859509191 0.14076639993623227265 0.38977226121800451208 0.60095405105830335124;-0.15031324839706941554 0.050332269731436680982 -0.42445115148166429053 -0.48679966862261220273 0.050276932403599206012 0.079319816713066385039 -0.34753332197723668662 -0.19915678452316018521;-0.35669404075329586012 0.8863181120398785362 0.99081595047755521843 0.92778063264208687766 1.8859583825639145704 -2.7917922784038360362 -1.9231209765011190171 3.4941253224493871166;-0.42819564178739322191 -0.025822521585563470298 0.12014227651323723667 0.31190849344669041576 -4.6595524558356702371 -2.9625355245614248823 0.011142579732478079604 0.94028847051896358167;-0.2771029832249297109 -0.070037991466883511849 -2.5417727469431428133 -0.56928680390295194336 -1.4539962521319589328 1.3928808205104388218 0.0007314378306018535017 0.02836717045237400342;0.1028764676844218795 0.64443821941665713204 0.25927448827132132392 0.5223331991986907008 -0.93243999500466878683 1.2418554186906620718 0.55701010173521570135 -0.45564012079329690996;-0.015258848123692275375 -0.46109481588881034986 0.1517808925719325297 0.61103271239689727068 0.38639919049927889949 0.023797547962340709349 -0.41139918592505764217 -0.066098998115325130476;0.66986895942449153907 0.62725964590342919713 -0.30012639146225017495 0.72513694626763358642 2.7601284435674360296 1.0503194779967159889 0.9124468047849704222 -0.8347468772763675382;0.0073130097342290045076 0.21138657235075120067 -0.24203249887788619299 -0.38975025669626478519 1.2478530828832021804 -0.68795471627466109332 -0.033393746309130133143 -1.0264761515294378302;-0.52302201409967463697 -0.60248152244963593116 -1.2637504629048721405 1.0996773722062485046 1.0151553933292820631 1.1747205095766775074 1.3893753901465493161 -0.2166350041418154404;-0.11780663959414605646 -0.07293667255374146341 -0.17879526252611799286 -0.23527666436723559951 -0.36239659763709580753 -0.28875422223969882118 -0.37624028621965655006 -0.15400520366869732691];

% Layer 2
b2 = [-0.15560908568259057549;4.6337449687508120988;-0.38675973704754235127;1.5232785264175208173;-2.042577780875593163;-3.1436931260129790111;6.0482488208300209109;-3.3670861749564249266;4.5684510729428309261;-3.393145103866284984;0.58572868344689954512;2.8476095315833118704;-0.5643735822947970604;1.195462727834528005;-6.0663557600869637554];
LW2_1 = [5.3805628336925934718 0.74641352880237343825 -1.5020508736735596056 7.7301382472247288646 -5.5003612429729917821 1.6598498172053572919 -2.0431429200327553275 3.1840415066020635315 1.9927510011747782137 -10.499731132936950928 4.9865732630579131524 -0.35687107884274765768 1.8817782274400005615 -2.9954158879937100579 -4.5014635799422260121 4.7016377449652555498 -2.9596307765853584471 0.42514021816984565172 -6.0791445554790604078 -2.8303040814579234485 -4.923168547195698963 4.485500736045667125 -0.52144922133960136534 0.93834196444353024091 6.3470079370021332821;4.7281953344924163574 0.6924926920953736742 2.1611884292724132806 5.1246632064004353779 -0.64419996172121418709 7.1370446983017465925 5.1693185998820130322 3.0269418090701005752 2.4780984970450057325 1.7424387302462258642 -3.2592792490369655312 -3.8307261905203624508 5.3976328589516837653 -2.0632560008262910145 -3.2057038668329829179 4.8065652827906202305 2.5601172162097571849 -0.44664167572164570341 -0.0062585605178196177314 1.2038344123596806945 -7.3665668020982577602 -5.0591014658647726066 -2.7760816288734120327 12.083101655295980237 9.2919001425528140459;-3.5971658810123288141 0.32522686047016358124 1.5371778155285089618 -1.2745766901034947249 -0.92357089417511595286 1.0240203369133582534 1.9067319135688860499 1.0623946478658625736 -0.091755019877334226264 2.7832471404145389116 0.73852289612835819721 0.28081310827561611898 1.4829801312888941656 0.071529132277268522855 0.77449063903539350218 2.9947015542542523292 -0.20322516922504693659 1.2334245974423094783 2.540895647542366298 1.481667871207431153 2.1521444927501409694 -0.26271769005313666545 -2.7414110296138662903 1.8181628710379169611 -0.32388836150438216066;-1.3291122909955337406 -2.8089809716605236289 -0.60844945509950343343 -1.9066542026814454758 -2.7361378622444325615 -1.2743112376440428601 -0.15969581291478829321 1.787371013702143685 0.21839085104522762415 1.9794761864671417406 1.3258472944461803245 1.5588452008653832515 0.90891671343316360332 0.7056400632236455861 4.4309901913323468392 3.3529746629465502217 1.1214578603347917785 -1.0357507211600762709 1.787901022179142263 -1.0097818283660005001 6.6221667581830443439 0.60713054964357215937 -2.0907235624528910378 -2.8465389007882762051 -3.6440324683021070307;-6.1331581015371670063 1.1923184215145274933 0.14370335139414652792 -0.091827647321328406682 1.7829867876373823865 0.67017882286645191758 0.94014176675882432654 0.38636246807470930875 -0.30451322761106447423 -2.1353589751215902304 5.8526515447779354773 -0.24514886107543726323 -1.8195437589932845768 -0.032131701754099406543 -4.1228654809919609292 12.116663028799978363 -1.4473967697798819643 -2.5021600757097384182 1.0627869559378460007 2.77908039186801048 -0.33741342148775504972 -1.4989888492207705806 0.63883169031041275865 0.80820835917530953196 -11.712975981019546623;1.6670152891967895403 0.79028849549361790405 -0.11247537643187451029 -0.50579627326577603075 0.46514467884849169943 -0.7079391949652036331 0.72917161431606747346 -1.9266819948079265679 -2.2445175064825182432 2.9742457290072938925 -4.989617188212258192 0.53960034148150337607 -1.0015244555845919106 1.8236437154467128874 5.2817570570402470054 -8.1559338572015143143 -0.382287894936598871 -0.12739783287250194821 1.1011845872439991467 -1.2102626788411985359 0.078373660978330011995 0.46446204231333626522 1.8200945408045154927 -1.3230805690503755212 5.7170503219945469908;-12.338851430277021137 -1.998964434242025856 -1.4401018277387724886 0.12195298534988355255 4.9414223228495277596 -1.2810928352962958066 3.1300495230868885521 -0.29784425804373254598 -0.91597423542970091148 5.4734741801552715401 10.451422166910782963 3.9468135432009106012 0.97870849536447546146 -0.51071583122954278267 -12.435485400018539082 -6.0787881959267329179 9.485027651015013106 6.6766582176628181955 2.2886600797492695492 11.603128663178903679 -12.694433636180717073 -10.031979319491172831 -2.3963411317769458542 -2.788358971707113998 7.1297988156068781151;-2.3201646872415655487 -0.28913620971217918187 -0.77275296344568822349 1.0705170365436047497 -4.184662623036022211 0.81240597935449665457 0.95406012395046257968 0.27957272803566912689 0.27922449683439382051 3.6852261572914670751 2.9581789262033830745 0.81961939373103531725 3.4252356767630458556 -1.0012441645926069 -0.89211937860758794283 6.2055606988246179156 2.0949000701010391801 -2.3380915434284585963 -0.34812926934273125124 -3.7392277479734405965 2.138820285986731573 0.098077618115155043577 -0.20020069580048813207 -1.1885647337872369267 -9.9411681787336316773;1.0429372921370356764 -0.52173271598388870629 0.87332606231954357234 -2.0247090611844296149 0.33244944982628926722 -1.9305495662086744613 -0.05719654046836517014 -0.27580207863565431303 -0.079284774262522714205 1.0425207216530569809 -5.0450657784579924225 -0.31702147806806102448 -0.12159556515705641222 0.25662998959298832791 2.6272237753034715091 -1.501996337271692461 -0.36094090605737455535 0.69212252896217141274 0.14502862037089525593 0.37335000873687290701 -0.096922441411087301155 -0.88603389042238334117 1.2588851411004413627 -1.509745689025749682 1.6297949078531692724;0.61783039614731016798 0.53720394729580533966 -0.070780648383001895652 -2.6054045744606120216 4.1352551263953305494 -0.036616016192396995399 -3.7472501243315305963 -3.9148340554715055362 -0.56223069903718170259 -0.53279225001751073787 -12.764940283613292138 -0.10633022762916005677 -4.0949421358083997191 -0.54439861385877619249 5.4260705933771848919 -7.4151736001326291614 -2.2085291339396473376 -0.093033312617898669683 0.75858521805529255388 4.3217958598372971935 -8.2138188980598343392 -0.23740753573894241013 -2.7191477269202692923 0.61876531361653730468 5.9551469273966155527;2.6394124902635405228 1.0699740278902989399 -1.4283537326834561121 1.0283690294408585064 -1.2945114685789722575 0.58449187829344939615 -0.5403218533369752441 0.15222151849106399513 0.83343707507717923466 -2.0856677810046932997 2.9044400192257682214 -0.80127215577085819653 1.9776485558132026465 -0.082591784478653104684 -12.157991909251864016 -8.1057559533825340736 -1.1668833350193681397 1.3218140865470124901 -2.3951604558159256264 -0.1243056598794682055 -9.4629237731315640758 5.2073194512685105906 -6.0609452284424136792 2.1130086629050990332 10.827272969132458513;5.5241049011824125969 -2.6644896480612754175 0.5150377126278141704 -3.0382720404347640653 -1.6950288505377537618 0.82563668395750278162 -0.77430480224692810332 -4.1296834102888979601 -1.5952968478930791285 -1.8006597835633402216 -8.139443618777454148 -1.2743582085456097719 1.863184735826441063 -1.5327830467743970289 3.0888915604241269364 -3.336514080621055367 3.9806591491206431144 0.57796978504995255577 -1.6928506226241266663 -0.74027890522838102338 -4.6268043744447879817 -0.39983681225981387852 3.6201207993699866527 1.803660065852974892 -0.84590610659942888638;-0.75491847377892828952 0.93916147129414884187 0.43548263651506441407 2.135119484042152127 -0.30530417572500245793 -0.30506201724650444662 0.16215022461530886888 1.5016481403069665834 2.3388054645592730729 1.0107095165544026649 -0.1510987005332613542 0.80228505180190934354 -0.0027784127903626616718 -0.22507120805054309232 -1.8354773609831216419 2.517200123790007904 0.63178748929123018652 -0.11064622363661130722 0.42784198488003083449 -1.4272960290406659922 -1.5131868523937348492 0.5953155885284499016 -1.2052281888935809651 -3.6484517827781743904 -6.5402331577029650944;-0.84358089498331334521 0.069898045211241516261 0.010856376547025110152 -0.41452056155951066962 -0.48475324768701810729 -0.0022666774743772847775 0.12477655303838444678 -0.34333771655194023076 -0.29369481065873664916 0.43428555487082687314 0.29492427753188810158 -0.031946617780192344827 0.10903950431537255683 0.14001158183952672198 -1.3133751555624884055 1.4620475716551479817 0.08987199038159002229 -0.040275762504796246144 0.12782226371211569127 -0.33582554791200036615 -0.55992488117294392236 0.0030389751033902154091 -0.26449446572932899802 -0.46377399733088531253 -3.3234879299756214799;-1.5418209280595691624 -0.10825476427093476406 0.25538043264620535755 -0.34577271689470961968 0.29226081052992775033 -0.28372404170827569825 -0.09438223993009917101 1.111020152059330357 1.4692301703681338587 1.4763475312055327482 1.2870529203819203179 0.64436544515984861548 0.11613575585510549359 0.35985882285313941509 1.2958994897238813593 4.6835012409831966451 -0.13640622014134634599 -5.7580847069240688185 0.66274421243878456345 0.023115073539509544548 3.8838431922491385428 -0.15679796814561985929 0.054578608417066847014 -1.0985284362490252974 -6.3251169668514641842];

% Layer 3
b3 = -0.29547854412101431931;
LW3_2 = [0.039371874878479422299 0.01714422740471350029 -0.058553299724902883872 -0.042285370495346689956 0.12249120699641090348 0.070883193559476365131 -0.0080729677408183358328 -0.052299577204008268483 -0.41080265456647990785 -0.098818265896243748236 -0.043141020003954803907 0.047942638288724324924 -0.14105714117689166498 1.2892073090626783305 0.26367532537229232403];

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