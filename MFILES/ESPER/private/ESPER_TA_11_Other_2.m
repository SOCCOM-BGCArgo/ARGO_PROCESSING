function [Y,Xf,Af] = ESPER_TA_11_Other_2(X,~,~)
%ESPER_TA_11_OTHER_2 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:03.
% 
% [Y] = ESPER_TA_11_Other_2(X,~,~) takes these arguments:
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
b1 = [1.7197131155151237358;-1.4478824360884423594;0.43278260124505496531;-1.7679901724287365727;-0.89962350824913572822;-0.68174436851128084047;0.88471478643628709282;0.44925115868372161865;0.24284998335099045041;-1.5802771478421511198;-0.42498703616698685348;0.64937058343985099018;1.1647884495014804163;3.2547962329367980949;-0.43151606191285146741;-0.49615312866565319805;-0.90342186499553334134;-0.47639733289221308787;-0.6536460358463758169;1.5314597856897651518];
IW1_1 = [0.17249154273534103576 0.6827876673140388375 0.63376831205864248009 0.0085500054806486539688 -1.4128575555956381926 -0.023186526335547387456 0.6345189831182725948;0.3579929223556910789 0.40939595954735069627 0.65746618374759879799 0.55235201216822882309 2.7548460624883523984 -0.13958551031085758565 1.773929845848745579;0.094439889263914833517 0.010985178404866392105 1.2813476206269054369 -0.54700111413440288821 0.23361360277391918228 0.052921163217383861144 -0.30758765621999828443;1.2260049701320261395 -1.2592280271795412983 -1.6022927059697285035 0.68593474413357369279 -1.8166380847441194746 0.64352489652179045709 0.49207447877220444887;0.36150917110865948834 0.63749639162456328556 0.95457514301435153481 0.43417371356298567475 0.47016521071298755796 1.4936736860611030675 -0.020221591071278317148;0.14811767830879085661 0.061930357538570617038 -1.2784254663740630686 0.084758564733615121556 -0.23915562354851932492 -0.16286284076013082811 -0.10098630364166741313;-0.29872604550214282293 -0.86073092730696121322 -2.5382789078150809203 -0.49160371648021911062 0.19875885145947363997 1.4033763678343755821 -1.2030928766000870134;0.96308074975625779235 -0.017088247872616646461 0.25420800508501922854 0.61530166614054182528 -2.0289338358151036701 1.1601128337151505043 0.75774176679251881161;-0.028597821043888114328 0.0401991733205676674 -0.081829949112543529099 0.041197169587828139092 0.10326325321920347411 0.1467343148573339584 0.11394562640673359943;-0.56287445962799331767 0.035148052218562764182 -0.87578026907546224766 -0.62067837864571917539 1.4387185411978684613 0.546728358841461759 -1.5216187804394891714;0.057646845329264485847 -0.077225960873733956458 0.27660607571303197361 -0.35477070097020385608 -0.85745423243656648804 0.55568310039988533244 -0.12105954670724740352;0.69956134815281756367 -0.1109352634378486524 0.10848369533944564325 -0.55042627695126866705 -0.32013902652316161834 -0.090172580154017148435 0.51712921560518987363;0.75974541032902342952 -0.74261430814618334395 1.706539129929066867 -0.60039278708416066177 1.9112881649698600928 -0.51261692107945078245 1.0970208359219852934;0.19296878644927609181 0.16282453380220668349 0.031976801067523465127 0.88953573866554880301 -0.77575681378089655915 0.58939356962901490178 1.7801879819683688044;-0.45718398766826523483 0.42878797309560301443 -0.5151288979050434591 -0.27266502341257448094 0.037204364024986426307 -0.30831718282209913395 -0.34706499890819980658;-0.87999846356342081855 0.01440949312984029021 0.23408942919538577865 0.73958969526198725664 -1.3010355106244493317 0.88298800146813261058 -0.26388718523167925545;-0.090151409399690754509 0.27423091636814106442 0.43893814447809786916 0.021437023256191400922 4.6176458827533162221 0.63127557585128069029 -2.8477762454090442823;-0.090679459892808622623 0.064730145724437845867 -1.4753625129857295128 0.58245496106570926464 -0.51004663183218867939 0.24747271501703888008 0.48187649029325291261;-0.14685625750624625718 0.15116152391126053134 -1.2021875490073783421 -0.13284471847746304185 -0.30828114394696370937 -0.098413629237455266252 -0.54547059695512789546;0.21009722321839993664 0.38987863751848750304 0.318044857246503887 0.094491706596729896361 -1.0805740331087250183 0.27887452094829812133 0.6668196979307076111];

% Layer 2
b2 = [7.8996207335025587781;-15.507258057721307054;0.71379717873683745655;3.2961982508602329922;4.5207995140532375444;4.8044189722405983289;-0.62661281918831024385;14.165418875290930245;6.1881326766695465835;1.1657366220034079962;-0.42930153357977918205;3.8852256027507197089;1.2728743197347414107;0.032159692168113387933;-3.7248454743600367678;0.25230437453540111026;4.9652601545315997456;0.25187112635030989294;-0.32586938084263483706;-4.0739023906590174917];
LW2_1 = [-4.3857527400826770148 9.2797763984432695139 -3.0763369588993629833 2.1685721912143649526 0.91171797317898617408 1.3771172760169361204 3.9823963041736640811 -1.473376961652134165 -3.8245610321348104677 -2.7413847095758896444 5.2783381935403008711 3.1708967561205945174 0.034746365260213366022 -4.1203628132242089421 11.935068219183641247 2.8364370389231803316 1.7710977657322859447 0.59449046702486707527 -3.8074101536059039574 4.5359819414531301618;4.5230226806579398868 2.6083942237803716679 -8.2259630714604963231 3.0798268155762218434 -2.3544389036678512461 -3.6114143816258401642 13.807249492950143477 0.4698303131580667813 6.7144893102577238864 -2.8626990382093775445 4.8296340833643958135 4.5338430559988314172 1.5100278021244835092 0.87045020207033896309 8.0308582404544033295 -0.88158152412148127475 1.6572289634312458162 -4.9971247619587746414 -1.1257805074710436344 -6.5124372002784278735;-0.9344587818949062763 0.56268528290887087628 -0.35145226807250484757 -0.51549399417829477876 -0.19641198836728157584 0.19646205760987608957 -0.12122807988415790215 -0.010202413149659564709 -1.3888201023312658489 -0.81407764267662141577 1.5468235256290541191 -0.63234422977439375124 0.12518259158555941624 -1.1305380507603215445 0.25588240880998180948 -0.61214550554046598396 0.20921516683950092119 0.39600301690559375212 -0.78568586211051549828 1.6964484958441001794;-4.9128001326116370606 -0.59307804933235064304 9.3194625219358364632 0.6772600805784143585 0.19763380167831781642 -2.0181295607361589006 -0.2708572897934436452 2.3542727990528207194 -11.882421312857067264 0.97840818937682039191 -2.1280380249068819865 1.0178138250442423374 -1.5405429363152289834 -1.1858068892183535326 0.51067923908127854649 -1.2916640370146519334 -0.00064843732047012263287 9.4576182649147089165 -0.53621977154993372139 7.8127775349268793192;-5.0926928633162953375 -4.0931023392821357376 8.8690180467239603246 0.20248661544760085396 -1.8837938015281134874 21.207541121961369868 -2.1076601671175567709 0.90634022012222936748 -8.4916510781080472015 0.78271904877379205878 0.39193870945013342677 -6.1783757781545762811 -0.60889428423191971351 8.0148366723437316494 6.8077125445949100424 -0.4873546547012198582 -5.2906392669375330584 -1.1829772989630047419 -16.061986635250018907 0.95173042620968295591;1.1872681934741398102 0.041426963534164404357 -0.10810904206528679228 -0.018280310789819710637 -0.66722662644469732385 2.4166295436304250366 0.20460041365678927061 0.30934654021479601749 -3.4552348323071422342 0.11171424428266271689 2.9226302216463806438 -1.5837378838014641858 1.4376779697528068169 -0.62364995079080953744 3.6163526785226407512 -0.49379202472829719639 0.29708601588241728697 -0.089093678952582702757 -4.2052850199674391973 -3.4869232663761606794;2.9536463955070422216 -0.92824319753497186447 2.4968457696787025313 0.059037810782074359728 0.96245804285146607793 0.10260194376762053758 -0.67761911053931311866 -1.7123517919459618764 2.4290806120640198174 -0.74521692308052056575 -0.23964767337919251622 -0.055105954274389711001 0.3519628480665462078 -0.058156924902607008698 -2.2355734578079338171 0.25182917575484842221 -0.33761047733126559045 2.4614312209592847935 1.5077649760693878633 -3.7733376980393140698;-2.2218756245259716486 -2.7797529217282184177 10.306480838072770823 -2.9297633078031184439 2.0226749832926973482 -0.86369286213862150881 -14.941826409279473253 -0.77672504383840113018 -6.6012963407407632133 2.9102695132555447266 -6.7912031128060110774 -3.3808674280421002756 -1.9084680412108807612 0.65872646876971863161 -10.558139409976266876 1.3368739834927052534 -2.0568970170632945305 7.2479740038656430556 5.179402670105368145 3.6319477211395265748;-4.2817269270027500383 3.6948016967751593498 -5.2615328553257887378 4.4182122767708547073 0.46714319452031710656 -32.638927706620570746 -2.5569745860642596469 0.36054219515953345176 -12.743685982384898026 -0.51275629565757741002 6.5177632526358415532 8.3181547884776581014 -4.8467264439411650301 0.28566288160245212335 -4.3553564541883034167 0.84879916781085984478 5.4512399772511157536 12.696416221272830427 14.373108822156369868 1.938010054990953801;1.5214784384930370997 0.77879413674018915792 -5.8419466730164826274 0.0026564453930772802606 -1.374383473789840604 -6.964678891922208237 0.97283508153931630424 2.558919446470824699 -2.5514283979216880738 1.9968897242823808735 0.68267804985317259714 0.81030671927845188129 -0.34149977349282228944 0.3918329508775454384 4.3739244149836293829 -0.41092542517695873627 -0.083731967527808023877 -1.003583684438654533 -1.2114139781728641143 -3.0070129140947572388;0.96144650081934346364 -0.57559514406753942151 0.32608904677570599295 0.54362118237791934305 0.23753200962777612504 -0.40615262770148202742 0.00074672972206447463077 0.041520098715895667474 1.3473132970472063707 0.85029499529129359825 -1.6980959431367532719 0.72652184961870980295 -0.093681328923202342174 1.3442772950148367261 -0.17254330603572975988 0.79172405113073052174 -0.22363207803902371285 -0.35242787296871919134 1.0303228525992338405 -1.9664366205575858348;3.1028816375426053753 0.11315491183863604641 -6.4904375291486440247 0.40392946816604013982 -0.021680192533583823866 3.7723107990321795846 0.27447971179619451432 0.96447308989878222096 -6.452385288912400263 1.0938031377157138024 1.6829722935931401562 -0.2237559113532729127 -0.36206203687911076017 3.6961071659469961759 1.9030200695012495782 -0.17346066616676325545 -0.9036008110071576116 -6.9515785723418908049 -3.2850434106209234209 -5.0908069610303572361;-4.0638818635294660098 1.6216175243650978732 -3.9553964393418858947 -0.25933784801042958357 -1.1662122009632684971 0.43610031576366314887 1.1062258088062839612 2.4151097525397955401 -2.2170155599365801713 0.67216533184473736817 0.76582194438412032333 -0.67782395703115394525 -0.26785564277428453028 0.16852115981452986393 3.2120709558102800152 -0.55306252015954060486 0.45491181490489340788 -4.478664126391803002 -1.9981498337738765603 4.2705823642815055052;-1.9668789755835913713 1.2568783479877905229 7.9688216673054688144 -1.3479939692173481536 1.5822713428442474193 6.7335874312841594858 0.67537835361659870337 -3.042376497164225313 -5.0842124839673488168 -0.98044675352002441659 -2.5116657287226384909 -3.7005555890391184093 -1.7887075796230758051 -4.2551500094309231415 -8.9927383315490860838 -1.8423680088352571982 1.9430004721673006518 4.5563297887067273351 -2.344941366437070851 7.7636040871191349666;10.019658188367152718 -4.2820097595400659074 5.3781500813382114856 0.79149148990120110625 -0.43733723218363296237 -3.8951371623964696767 -3.2545543920104167412 -0.93650433301402824515 8.3694847441571251778 -0.77508809986287074967 -0.36644262028491902949 4.9960397766585984414 -0.49410953982512101357 4.486329828969463307 -1.6748983783261570668 1.4096973955974052561 -0.21192516640652594995 7.1910627495794336994 0.49614686495051707471 -16.701279762125039952;-12.066172979793005737 2.2513382667366919065 8.5635125701619223548 0.36282355998677284781 -1.1679001290402786228 -0.72180803959489803212 1.3204302954270152881 2.798884303572847454 -6.2346853266948274452 0.58005550655670745641 -2.121334824134380348 1.1767430377682910336 1.564868477119950585 -2.643658940922654299 5.6007417145235818268 3.3875492630559969065 -0.74906113236667903887 4.4561610964797431222 1.7406374083262190489 16.808441604159213512;-3.3754264643707636573 -0.60254953455970816645 7.2115777645812908503 0.6614953626502079187 2.2199230529181743776 8.0107849929995502691 -0.16690084377638053637 -0.72388870831648433057 -7.7857285626712169346 -0.80887440084731032641 -4.4815383540707260934 -5.2277033935369754403 -1.6925331804951548875 1.5921731622681802865 -2.3326361478378809799 -1.5595302310787495514 0.60969374102235573964 3.6555229910664395199 -5.4457665000473989281 3.0409167840490853862;-0.018054235645727730042 -0.82037027857660005381 4.8603400870616360052 0.07926789377416418314 1.000872940960290336 0.4544425396833195463 0.45377477035484675705 -0.67437236753766260922 -4.8354552459398458808 -0.067920383660082572774 -0.86769572254024096569 -0.58637602795558785296 -0.33760620837328209065 0.063730527976886511277 -1.4477149897816004032 0.071517551427808828679 -0.40586959189884092014 3.7671127571764246866 0.99258854954414144078 0.32618440149180022436;5.9549658772165248166 -2.7772685206405856029 1.3382500903624789945 1.6515653393661071657 0.23552567455358816573 2.2354758719876652862 -1.4105510406066277262 -0.40345087388956502039 -5.5433156548799455976 2.2871435837575919692 5.0367863643942696328 -3.9229532217002280881 3.2993193173061707846 2.1600125219741261873 1.9351901680085672464 -3.7925872474903270515 -0.50636661825130957482 3.5687134648343845633 -6.5617077469854416805 -3.4887341485019724985;-0.61272591216618965682 -2.814365415980408347 -1.2044963240876298549 0.67702500683865263209 2.4781229567777218747 -7.629638110218592395 0.89524496174951140315 -1.8991723670092277576 7.6186822960189051201 -0.42210746601048532156 -2.8458781222934521615 -0.4910941103208508185 -1.1602440279215402441 -0.19230830635086221259 -1.9583607049145115298 -0.40645784249787469955 -1.4755078980030342795 -1.9624895916876101687 7.8334221592015067515 3.1414271396897048838];

% Layer 3
b3 = 0.82405223518005610295;
LW3_2 = [0.067635471206596645133 -0.28072135537952352946 -3.511533630111157489 0.084519340157331285246 -0.030162874510414456469 -0.2078012869232216675 1.3948355496703450651 -0.24681061553326721913 -0.013439469271266260217 0.081132781445467649917 -3.2185896620259333289 -0.29107166163138037396 0.73233376820872586599 0.034838206636440376129 0.078655464036015629303 0.071017182406856907417 0.10624143052070936233 -0.56568210193222001614 0.03247510441274553683 0.10455391272388410429];

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
