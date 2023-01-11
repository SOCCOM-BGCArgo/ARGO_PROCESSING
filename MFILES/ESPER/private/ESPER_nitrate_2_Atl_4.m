function [Y,Xf,Af] = ESPER_nitrate_2_Atl_4(X,~,~)
%ESPER_NITRATE_2_ATL_4 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:32.
% 
% [Y] = ESPER_nitrate_2_Atl_4(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999994516886;-40;0;9.1619;-1.99149636107278;-0.03;-0.52];
x1_step1.gain = [1.0000000017546;1.00003626070655;0.0153846153846154;0.000309885342423303;0.0656325902039552;0.0604107502911417;0.604229607250755;0.0145232735458572];
x1_step1.ymin = -1;

% Layer 1
b1 = [5.7997349789571677192;4.0654349136737657489;-0.99180510299562818499;-3.261285042592863892;-2.2563216024963117512;-1.8654584894198484868;2.884109598316376033;1.7280178888121573433;-1.4111843821460592352;-2.6540187494494453802;-0.64752839816865093159;1.0814780185361496656;1.4183870459170642153;-1.0504621004518570615;-1.144546979342887072;1.0182563653170861873;-0.81875475178871637283;-0.024220818442758796502;-0.5477099813582487764;3.4865720235334936028;-1.2207192216044973065;-1.4273161047314064653;0.4307278048353196942;0.87799232691792494787;-2.4338999420546065799;2.3209478129867169294;3.152572007618659633;2.4676889472973853046;2.6267098973158260478;0.11410248935836497797];
IW1_1 = [-0.019633454280961550736 0.1283776703935667618 -0.85179901465383534021 1.1046343390243860139 1.2219551880900116902 0.34703537682788310192 0.8310990003977013485 4.3455563107888073304;0.053697234095742195359 -0.59925000479928447472 1.0121755567433010636 1.8407001600589991686 -0.27712858784496646747 -0.22570587263104105458 0.5736075363616650602 1.378822731659881784;0.13164889809402019472 -1.990258241032088371 0.68083106542171378806 -1.4419246831714951362 -2.4567093907803476682 1.0217115262981020241 -1.4514655228268924514 0.97889220312595692786;0.35003782696479057712 -0.4235308447961146916 1.7171963522251372236 -1.2987796781382319189 -1.781586562612020419 0.27066532784231978059 -2.1432660789756452502 -1.7032801475621670573;0.61613023618171880447 -0.76627312085871390757 0.23759553162951169147 -0.13837591553938125011 2.2216327363250214511 1.664663322920477917 -0.34936055073493871292 -0.81302471365132422676;-0.34135977904216913137 0.14994521557950044044 -0.37165265968891614268 -1.8159967248180699517 1.5044656418036106427 0.56237904443941688637 1.0668160878761874244 0.82707872827167716601;-0.7514186911279674419 -0.20392067654962978129 -0.78100305887211163292 -0.98088989802491921566 -2.1322835464332579036 -1.5256567536420062581 0.73653546689900228905 2.0716476514849122559;-0.47765555098969575409 -0.64819122446871713095 0.79953162412124012537 -0.53339468900500164317 -2.0047237423893591313 -1.2219983295268017631 -0.53971595492884494316 1.3692324048720450058;-0.12941903439748497773 0.24305174652321986928 -0.60425306565601788478 0.20868219291779854796 1.8329337282725619573 0.39641375891049601687 1.0989888191142069829 -0.87834374767199829837;0.053133304330292792372 0.23164319072245537212 0.049669228229204193115 0.98670265609555707353 4.6954279209079095025 0.12982497341184012773 0.31287627784162053146 -1.0084522614634874316;-0.5233483398412279719 -0.44593123022104302899 -0.10543056204408902954 0.24588667429283855848 2.6207004385825491966 -2.4692844250484253799 -0.93514830710963559035 -0.099076606230988517598;-1.2141522526429149309 0.20225688815985001789 2.6718002892095076284 -0.54786896440502230732 -2.6028213620683575691 -1.3763712200726243573 1.8981748905693081042 1.2835871615506180277;-0.067593503704259941256 -0.16754760773514773553 0.021910532069776104958 -0.42462075044942981794 0.41926298896957547857 -1.2719351774276910216 2.4768001987725796198 -0.70993749275319006919;0.062997754066672559392 -0.03321602914018786934 0.65100319301347098211 -0.60565956918326713243 -0.0026455070120165909142 -2.1102513633556925221 0.17031220335266208932 0.098573723529220141915;1.1640889230783879604 -0.23815600942149148245 0.027716663157952538959 0.27096026861177680001 1.8146869394743776471 -1.4095619121367852777 1.1671349744473262788 0.20817470308177751792;-0.52025911714613615189 -0.051102431712097197425 2.1518908852647609109 0.90767281264473331781 -1.0981437595755136893 1.3036184299372761686 1.1101746822838916007 -2.2621075005121533685;0.50257824809334816774 -2.065898906795390122 -2.3618253429128683862 -0.66751611353288864414 -0.19712284150014097728 -0.69683649327097341164 -1.564212945014845868 1.0091123139628237482;-0.39369384497360626618 0.76760281345657055407 -0.13718247488816900925 -0.18044255872543457642 -0.77144780484452779579 -0.60564021887185826021 0.90448881697301397597 0.79737821177278367379;-0.74458018928985691964 0.42148381404820778329 -1.1439406548155548471 0.33721386988629986625 -1.7033140894518878916 -1.3893822347403448347 0.0026771479706548806377 0.81635511266137883446;-0.18629219857078058165 -0.22332201513407839988 -0.33418349258692775505 2.0840829260957156599 0.23252770238943620829 0.44022363148205284622 2.9766292320280252071 -0.42594575520330602014;-0.10810763340093899887 1.7618638096611307642 -1.6115793349191032124 -0.0078163894671814705684 0.95731387132245981952 -0.98991310973428636721 -0.89743827822032851138 0.78565703491250549195;-0.95101514518932106501 0.55427996294665315968 0.56055125540861372446 -0.73863779886383118534 0.20592419464332997747 -2.1759961442797113307 -0.37103610968603534248 0.32444525050880329564;-0.37209212255922874402 0.42928816129264601997 0.088821114031377143561 -0.13059283996927048555 -3.8262299244054966962 -0.74405058339543506118 1.2571722019262625558 -2.1700429430050540702;1.1136078236643311978 1.1961650863035322434 0.38207825204254652895 2.2788972510863101384 0.71188784089373424724 -1.4210574757622054243 -0.71173092364893220108 -0.68467380517494735859;0.042865578937702306095 0.22852601221396851017 -0.26493852431282732818 0.70373354945473831634 3.9019972053064773121 -0.071407184153266120563 -0.59646599401550537412 0.10524656411384254728;0.27687378847623961287 0.62218112094899913345 1.3688123728868175633 -0.62082148492908673099 0.99290812465118438723 1.951217732381157699 1.0715765357701780847 -0.47276025831450507741;1.8520839797040249408 -1.0882995370455967876 1.052354693803448038 2.6629894237115405708 -0.015173701527928139948 1.292420984736890599 0.064092516499043766798 -2.4653507878206939452;0.019324015046207475288 -0.75679486292158104632 0.59838547681059806038 4.9662310095111470432 3.2937275450179015479 -1.2239013608358468854 0.31612809171974148859 -0.4260310717827476279;0.98363773352497885227 -0.0586775852168974002 -1.716865188177623569 0.12962102442427189675 3.3213058148802296898 1.6091322688861724455 -0.18345963926439848302 1.8212447072627813149;-0.58586116521502162868 -1.0315795760037269613 -0.36827364513005805824 -1.2827809509921888065 -0.87366512756122260974 1.1210526750078897429 0.72308616288363825042 -0.33373739658941348196];

% Layer 2
b2 = [-1.9338675486044383245;-2.1159981731509156511;-4.2964835314468130534;-1.1648298350628119024;-1.6173919622200207957;0.63104101022144687416;0.58760188095806920039;-1.9515788474471884584;1.0420173321107055742;-1.8388120559427316714];
LW2_1 = [-1.2473914793461367001 1.2031646549907666355 -0.86999852398366317452 0.33957401351427968228 0.3580916887671327653 0.039103830771079331474 0.25351488438398184222 0.16845396293461514792 2.1277721998583989382 -0.010744696671625909162 0.4774551025739244281 0.31104383612254632441 -0.74370361783597072058 2.116869153841275164 -0.40238531001306138135 -1.1799904629025959046 -0.061343220212717139428 0.57513849493845281557 -0.63591786118863213861 0.61664460718045643883 0.35422714531315280251 0.2322203222041537285 -1.0544462111173187324 -0.60229186015881708283 0.98144201404982844839 -0.012890533892312892872 1.1696947610198265011 -1.0248020890638258606 1.5039327013220957774 -0.29105997169444536299;-1.0900191042324487434 0.75163323285482808167 -0.4871688183279357931 1.3001616904383335172 -0.7985158522605004805 -0.63267790139380564174 -0.24144384297290158448 0.038068324641786416285 -1.6782142911835573607 -0.16464407502835295971 0.81524147221213738579 0.042766114311584342567 -0.39058717788397701032 0.7512830257051248406 0.37349506648048730328 -1.3292931727774228712 0.22146263586558481706 -0.5703204670055170844 -0.76391432210383447554 0.22242175365444910584 0.84351215447741023645 -1.8332321609424837572 0.33929679661299705362 -0.46017597539619886637 -0.77509086325449050303 0.07052639307080031672 -1.0954653243447709521 -0.49648524173808838977 -0.036754081301174870711 0.86079101256106860252;1.5201400984923489901 -0.63766500717414542621 -1.0812864306681206905 0.48616941976766070432 -0.21611221615061199097 1.0493838872718015409 0.60240065946410104036 -2.6344077798153420211 -2.2705047854598765511 -1.0504239902300060994 -1.0984624280494263981 0.10190814489207747351 2.0534709301230731349 1.408059497563377338 -0.62488694596291960437 1.0640439481164882807 1.0542212619461723833 1.4107616792484001245 -1.0521954889190496818 0.5099626332318919486 -1.5588830532620010505 -1.322803387433266753 1.162122619307037974 -0.78542437280899546348 1.2586460648804658202 -0.53718986155279391248 1.2401471903707579347 0.30236091789757918447 0.53876792713635768273 -0.79016941763252468434;0.35677583282328834136 0.083625518477293131481 0.76241425393628903251 -1.0615368981084727373 0.59028827813877526864 2.2742846686614202056 -0.52849880917327407559 1.0927390970681267746 -0.26269884329531761891 0.53411849569942859706 -0.39539472602113290289 0.7999081999915236274 0.46599190878000429361 -1.3247422112136171268 0.4053594671729454646 -0.37865272060592486403 -0.42795931271187709344 -0.55389288702827965682 -1.1263064103450777242 0.10480310990273560001 0.034837971720304089507 0.53780838324870716693 0.13167770903764669521 1.7504149964314168475 0.35852810288383657467 -0.44358119455481820337 1.4950888093208858365 0.42575042788315109465 -0.61710248829030844853 -0.24507329401823210735;0.7624403079980938136 0.36361580232642498611 -0.52518966238546238579 0.3577871165821887578 -0.77323427589824056394 -0.76596329100793070221 -0.043183995480906157971 -0.15154400056111924644 -0.78728620930435588843 0.1186169221647740829 -0.97514235002856308743 -1.776431332522037998 0.98698999756862348942 -1.6376564360450507341 1.1521047379667674893 -0.2750174663017754817 -0.71595249064637145331 0.373599304909625618 -0.24301044091602766462 -1.1523513485953329472 1.0212080353422838819 1.1731165691770668591 -0.55029264143365741457 -0.59726347872582075116 0.55204961207515945709 0.9992643809353357609 -1.2665201266435490712 1.6016531703737439951 -1.1538750612086798863 1.8115506858973939863;-1.3275352278717793553 0.023159537068601611376 -0.11073546070136039421 0.32789407319314012312 0.25636249631569091312 -0.75614353688521263042 0.083034329948895022055 -0.46275825652960517198 0.92967275182435271663 -1.2424068775478003257 0.77194339526523758632 0.11775095904014935833 -0.39561220453373707739 1.3113382165714408956 -0.65215701756606592188 0.35139682500689262135 0.48566945974298159072 0.063243341315044554474 -0.91518565787559635893 0.57243191028159257083 -0.1362491942015427171 -0.50184889173363411707 0.37600739869681848448 -0.46655631515334511361 0.86026070627154149761 -0.60259409420927667611 -0.37934359902164310752 -0.74664839692679441097 0.96802568818630541081 -0.22733898186399201236;-0.98760418355387158407 1.3299604203131534241 -0.25789252497710840606 0.40760429101471612467 1.3741204554229884405 -0.014231552639688121661 0.044236475556292544076 0.21369408241820034289 1.7596818203406825454 -3.3278328464497892725 0.52116185580010609968 -0.048885854245968723975 -1.5291613769536258616 -0.44079674938771806758 -0.08090485589614374029 -0.70212786026462348232 -0.93816274238886854242 -1.2627640650170528502 -1.3283507017128479699 0.75086246229252373574 -0.2297440630946568918 0.026684023734299497366 0.45778821457457719024 0.61117371700958100789 1.3857805812054893035 1.4184861801175561169 0.042046331008524832906 -1.927297906738885791 -0.5058485412238664658 -0.37638425003517717027;1.8678119601264158156 -0.28059071579658356743 -0.34418875970560092181 -0.67406535832931102625 0.059068987008376719339 -1.3034115497190785771 -0.20663516158704051118 0.46816209687709120502 0.43362073575242432444 0.38473428276948545834 -0.78292520058299619468 0.16936603013384668448 1.0765987539539716433 1.5739111166655268548 1.1268123940014254103 -0.67136787389071683751 -0.11155231100477043527 0.74527202037265938284 0.099973064468775865721 -1.1729396946838763771 -0.19226365148788487902 0.89740952228415360459 -0.64804100247478324537 -0.93057872284393183815 1.0789137349259421139 -0.35886432941307094824 -0.073319723316603496732 0.61312458337643027662 -1.1824028581105310298 -0.082760119213535515703;0.16368533604724194719 1.3549359380806897946 1.0239491267212337267 -1.4027071150624506046 0.51191129745947516749 -0.038980193927437928325 -1.0695235628591217569 1.4477445806690021346 0.40250391464203899572 -1.3034899213769635118 -1.5420194087156053886 -0.26427870378428758302 1.2461565560780738515 0.30798270234246349908 0.25466713216722436375 -0.27642365114325412323 -1.39797775508941724 -0.87595735906124716497 -1.8583746496180102437 -1.7124288079078779923 -0.13647132650421603395 -0.4303732794708987619 2.7959658366609962954 0.9307867720763087771 2.1256234645223042179 0.64410589989365085728 -0.62650168461429500688 -0.47674994026602712793 -0.23157606730464594391 -0.22879433605964477616;-1.357561932304286012 -2.2548120971782705801 0.40378979182547514171 0.39482055567439233723 -0.94157633087532088823 0.36645834797497939261 0.60321306016633247093 -1.0814524782700509764 -0.29787037224120022083 0.63701398615508897372 0.88228074902695619031 0.40821542397950949521 -1.410504371286687153 -0.11977562279328007877 1.762122970566160518 0.75344551614367982761 -0.4607276598291673575 -0.39219427498462938741 -1.1043797246705457926 0.17496243236419659062 -0.61339082900120345254 -1.0450450541223095247 -1.1148989022698525542 0.46942054395378163267 -1.4651245977681790489 -1.1868008545836941092 -0.96927226050893655618 1.0698651609953762076 1.0630823620512404037 0.53405123185996694612];

% Layer 3
b3 = 0.087122072652120324809;
LW3_2 = [0.19046039866802594309 0.20791349835853947803 0.22040348088712766184 -0.27995461101619584232 -0.2442284990078418816 -0.92080446942782612219 -0.26047699060477236266 -0.3212494774290896582 -0.15510356631996885524 0.38471229598347628231];

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
