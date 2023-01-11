function [Y,Xf,Af] = ESPER_TA_1_Other_2(X,~,~)
%ESPER_TA_1_OTHER_2 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:29:59.
% 
% [Y] = ESPER_TA_1_Other_2(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timesteps
%   Each X{1,ts} = 9xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 1xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-0.999999984769129;-0.999540494077341;-75.751;0;29.715;-2.05132667284751;-0.9;-133.803853032615;-0.1386];
x1_step1.gain = [1.00000001035699;1.00022983623563;0.0140762793578401;0.000302617642608564;0.251952632905014;0.0584431850247228;0.0422119037568594;0.00432440760873659;0.00826143816767911];
x1_step1.ymin = -1;

% Layer 1
b1 = [-3.9211180601327266793;-0.069923837183868173262;2.9415251109907289973;1.4212654319596325081;0.070136214531825821772;-0.24992684940528317461;2.5933092043394831983;0.18615327271862339154;0.33856106563427162959;0.77258271141977674645;0.53330689020244437248;-1.1126598194981169865;0.39713814863604190997;-0.70250920628908708476;0.10070658109308242012;0.14738382332725746826;-0.35768768315144455761;-0.83077791680767654192;-2.0424772411137572803;0.43420791693142579692];
IW1_1 = [0.57936717689201000514 0.25350736049186362031 0.59031608108868993678 0.69139281582890410593 0.64242405510954192494 0.11256573398035961409 -0.52381166053531669213 0.30868434335170680249 -3.2038864511344815611;-0.77889961398325058273 0.14723114993453612076 0.1310816391259732594 -0.82124896533724500713 2.0393667928970660164 -2.1288242581279832599 0.0054547539406692921449 -0.8552698129459870291 -1.1886269085078196017;0.12737116940176954061 -1.3625022500255965507 -1.8514181629465220613 0.18964202513278921591 0.79821796203626438349 0.082740386591090633872 1.0592196039614218961 -1.1826392109154448828 -0.38463326342898046528;-1.379420118056993827 -0.81007682878753783928 -1.2497880121827802391 -0.28428532070253537123 -1.3895198402714072916 1.3070464953512037631 0.49113041184028272124 0.10504496338341717809 -0.27775965994562074046;-0.49726531852647309906 -0.23751352381224508092 -0.36452453076813662314 -0.64606352658397547817 2.537228203249684455 -0.81097699370736164859 -0.13842327374927299033 -1.8631974792862504575 -0.89628629766763212583;1.142876823567457123 -0.33203430957876672869 0.4722724684501376502 -0.48187748170143290816 0.58886454480732619565 -0.7940505235849707466 0.55618091750364495596 -1.5156465487721286767 0.074630339541248577606;-1.0911049317514340817 0.3162085257416324402 -1.1512209035864828799 -1.4913155414571839241 -3.596751380413408139 0.70619918093241496049 -1.3022104529001363726 1.1211112923368729 0.70094620361868187342;0.94162800626533049719 -0.9908209012401872906 1.7547737704942272252 1.461517796043880324 -1.4256568306135257718 1.0460025435655440074 -2.0082333116967641118 1.8471345071312581521 1.26987724907345223;-0.32584170218397706664 -0.0030355820853320220412 -0.099294678950503617587 0.057747196351239191148 -1.2716734519665244818 0.10604503403550426255 0.12876698545920223227 0.14898265451377612489 0.74369738214415603839;-0.20954681798113353186 -0.036509187634205075845 -0.0096940067441765778028 0.20904946781469607719 -1.2959747144922644235 -0.25896956699191542306 -0.033437322835052539494 0.014397585244109598057 0.52905994618395901785;-0.32728125381748846001 0.17508799581330267148 0.057965851917532562609 0.0748122082940294586 -0.49899433263274017847 -0.60735474125098154463 -0.14564335206893308516 0.045562515247111620709 0.55788986989266997618;0.72432027347002991835 0.35189513505319769404 0.15063118211914913736 -0.0296055956908027329 1.4278602602829413293 0.35255150918325500742 -0.17262194587701037984 0.13688239922467509979 -1.0275889812274745516;0.60340262474545180993 -0.15874095909111168856 0.56673823827619307369 0.36481227901296364724 -0.85018057672396574187 -0.40470332681833520727 -0.13604078522977763877 -0.31548875452561564448 -0.37183729599944842503;0.0056709974208898790426 0.02662368513987945981 -0.2647784326084808848 0.016248882844869850378 0.10273793447264255307 -1.2266286566361923605 -0.49592858865945127489 -1.0630299709731765301 -0.37550994804919712866;0.17833772147068835401 0.0060078474642343178413 -0.12563391992028427602 -0.0089215620119109523278 -1.0595673005473178474 0.98872487316238377719 0.2010708961012162288 0.85064195499170269787 -0.21850971724843087918;-0.18253741400028361541 -0.48267284593657899805 0.31851356624412296981 0.18836896616767367929 0.012875850881982831703 -1.0757622863297879867 -0.41635307528966675195 -0.91491822719543114406 -0.21952239001488380032;0.42965263701157691001 -0.10169136761417466031 0.15010164471611800452 0.1942958031406842867 2.9236120226074633344 0.33118267341074397736 -0.18410762446551237614 0.74270895155624805373 -0.8127239173787219606;0.1854987774997760297 0.12893102930613833945 0.24389677183216545986 0.26730708854329077173 -0.092593571641302857556 -0.38842441466145088969 -0.15682529478430129455 -0.33634549076902880982 -0.97191331767989186385;-0.72192442507416709763 -0.54694804812328101651 -0.32634465605254731058 -0.44381366933553845211 1.3096591184379839934 -0.13923955301746435143 -0.30101973546388527403 0.59222238974284202584 0.22493967533095321487;-0.5863391850071431044 0.093997261804957182862 -0.27400166714113477484 0.22808830485274741995 -1.552199120364858631 0.1852403526037407766 0.20406692234473119973 0.083855909734607866901 1.0138569771996834845];

% Layer 2
b2 = [-23.437499332163291399;-12.968871359433926216;11.417417990440577924;9.291454025176710374;-3.7656438317041560637;1.1265169403831030159;1.4285933908903827305;-3.060369675490830943;6.0816851529791531306;0.51247259866377015136;38.025893513063429907;18.287556045474090638;-5.2608195581600352142;-0.12896052553207670854;-5.9422180630220644559;-6.4380745249783766226;-17.8260122046037921;-6.1624984440402048591;6.0638224730361400816;16.986835462310477141];
LW2_1 = [-10.219673844022391762 5.1838344358864292261 7.3850516532286984983 1.5110727834374362288 6.698627127052652952 30.077768022256677227 -34.619949110278433579 -17.502065046610926657 -17.645695622112359757 14.177631137755831503 -20.044270452568579088 -24.948675999903432654 -9.8113356737445620581 -16.221573152571931331 50.565842717962027564 36.169140860173058627 36.824306637463614322 5.8694398934602993734 -7.0320342041309356063 21.851997229130446954;3.1066806998986193911 16.300364469707101733 1.1133228994347506013 -0.21684604877065810635 1.9780654134620929696 -3.7007530279690392661 -1.956952167419286237 -1.2915850202168617233 9.8711137926736860493 1.6555060477816381237 -1.3056875151927438417 -1.9068882469971633054 2.0114899136184245521 8.6305086103096080308 2.54846553607602333 -4.8565836353559665994 0.6792912035056896114 -4.0586325505611666742 2.4690863504389670702 -6.447045408822818402;11.222165055827799662 -9.3835595939482292493 5.0152032852180674638 3.9903799708948115565 1.5410872550398744263 5.7053500548326576691 -0.56191281553944694149 0.86010174532123584434 0.54900616936613977348 3.0531105459325185159 8.0989046596995599714 -2.1867859108330200968 -0.31633262083530000464 4.8596994155514030567 5.7908383258752538225 -0.15210759951772717669 3.5885174463694506919 -0.67231375586328279148 2.7904983420648914461 -1.53401072158117735;-0.24180317569429926605 2.7920367280219329231 3.5942626635077190578 -6.2911908849576443359 8.3153205132593637217 -0.87874908796549633383 -3.8426934200774880424 2.7207307414956951597 0.32786791452089952825 3.8553227314769435985 -0.82625525651358933121 0.5817332939788822932 -5.9753032942751902468 -6.6604063643115765103 -0.026207550166044417395 -3.3283976414414944145 -3.1849463973695950614 12.915934277369437666 4.7261756806251868923 6.0442443128276437392;2.1813596692299466184 1.4650905160204350519 3.2252384295654068502 -0.80008408242894624163 -3.4219739083832458881 -0.14886776622384895186 -2.6941864901295691226 -2.5742083924765788439 -2.6209213301145593 6.2239345242383778967 -1.2348909486665777813 3.8057511638816361454 -3.0754112714927033245 -2.2587761213065200572 -1.4440540387267806199 2.0157987615413013316 -2.6383617592418024955 -7.7735612587251692318 -0.30796330364819762826 -6.4024171042771191864;-1.595534978313035479 7.7003521038517597219 0.68109015792849414428 4.2564173315731892444 -3.3035557848262486758 6.2379383137922470581 -2.235252513878647207 -3.4833847649065776686 6.9648353318613702712 -6.5204085987616986486 1.3074430699668226907 -2.581018216281339317 0.18942537887809093866 -1.9788340060157010303 1.4723580824992525962 3.2155647019348703175 2.6078084489176562855 8.009001469587655464 6.0421928942679272723 3.1403452294544864642;3.9394949045806360033 4.3408505841958371363 0.17966046849831984744 0.048908129561752257397 -0.62342800911650897433 -2.3450237518232643907 -0.017322138174165474522 0.33865094998884182065 4.4907226166863569716 2.4903557716414947976 -3.2939794060573053258 1.0120647777759894659 2.0937786214057800827 -0.15928074976497927362 -3.4757408439219239504 -2.1875875738967489603 -1.1562285399191361002 -0.48569458853363894901 0.1454929139122268511 -4.3700022975912720113;-1.5510843812425947341 -0.5152596783292824334 0.19335275744249083574 -0.65412161517680555978 0.17835280278451418057 -1.1414171787021607507 0.96298998170455851753 1.1565143649108369228 -11.83709358657996269 -0.95816226936550308313 1.2510294249605604655 -0.46926471025170890528 2.3336538893994882571 -1.5492414683739674786 -1.5646514805879672139 -0.43958913675006061617 -2.9495255956936499153 -0.94949332316898316098 0.33120558507103903256 4.85906418271928775;0.47628988457399740186 1.0799657434216889396 -1.7490546186931297257 0.95319984090801823662 0.31178732488103039211 1.6749153940276897412 -1.1665161228789300285 0.71488633970412474028 -5.7983620238810171799 8.4937401895825939135 -1.562626175592359612 0.33336581339849047367 -3.7167693238487982121 -0.66107841721645721478 2.0236698902813992085 -0.44376580526183012632 -5.6283278330076500495 5.5190064740936541909 0.39831309755041416354 -2.4159474187798686273;2.8350300842493751929 -4.7696466585672974148 1.3459994883058898907 -1.8016205633987729673 2.1857971383047645197 -0.87174400045890043653 -1.4439687723739191849 -2.045143453959380242 5.5328897851928111606 3.3865662390469157828 -1.5541455500878644802 1.4430549898632787365 -2.3893972276938924004 1.4851919812860616688 -0.1096938501970749491 1.6375126117847289553 0.52006009277258913315 -3.2023864651590367991 -3.2359389034885985836 -6.285470228635571921;25.51848743160918076 23.752758540731605308 -21.796127301374884411 -10.995176340877822696 51.834845917181482378 7.401093001995777243 -8.4353585389173684206 -3.8280338015318613465 33.552970106059383681 2.0833899807193096976 -69.806429839961538164 -39.187327995212285714 -32.288177459739692665 19.342815249458535476 28.890230913310134042 -28.779758804046494447 -13.021056821295895034 -27.374575849776043412 41.085174636921792057 -6.3833656938733565056;4.0009899002762630005 -0.56833771602798988098 -1.3128480121862848851 -6.7367896285643062981 -2.977107711142354507 1.8865953952005669603 2.0607912905376473134 -2.4231841629261947624 14.217477107142972415 -6.6120871256196744881 -0.052554033092710478181 -4.8140617173240176641 -4.129867386835601728 -9.851913269221761027 -5.0760403219761531801 -0.33883467461387600794 -1.833514727125195698 17.331342230440657914 2.0770578373085957224 -6.7780438011546415566;5.9866425165429699717 -1.3198694144613281054 0.2226853501203751029 13.721624748704444485 0.50549072241929660088 -4.8547506199740313448 -3.5936499515185960085 -2.026127144601348018 9.9813216947253540212 -7.9009887163996177506 -1.1268290099112308056 0.25258768063273595583 4.7220656920975256909 -11.898446518424018237 4.4454693775769857922 9.4873079425178552526 9.1461798659280937329 13.710541542984683971 4.8465243712998935521 14.213652956644374115;0.4162833095959653229 -0.13267374183787231101 -0.10338690973877606716 0.069543753759833071282 -0.048997414000488244412 0.27212676814486835841 -0.071093776154342300244 0.017983723024558070924 -0.80372970405812360628 1.3696289027702619467 -0.29085264410920669276 0.27965774364226925197 -0.6083311993243581961 -0.052107581505744612893 -0.28865093302562294664 0.069734714530277422395 0.34344224773628906355 -0.46749858620293449629 -0.15383545726876463045 0.17154691137122948796;0.22824272351011665183 3.2704111553697785553 -1.5761493431553215139 0.76434554877062432787 -6.874608778949220067 0.49302324477351205934 2.7498518382695928963 -3.4229923306175700226 -3.4216979295713150222 0.33295625932624484111 3.0833200766728205977 3.4679845806262226127 -5.1660270432904162874 -6.0813947442104749896 1.3719583889485860162 7.9660144845060374408 -0.69523785887118572102 2.5642777101207583357 -5.4514562817734280031 -0.45677642127752970946;9.4508309484481305418 27.348269459630774492 -28.864227493193123308 -1.2559813985168721828 43.685391082415492292 1.738017152198791937 -8.0854976234727313766 -10.521076191078897466 32.623160860441544173 1.3154564346444932443 -51.125553830268358979 -27.546533749839209548 -5.8822689637343428259 28.99129027434333139 2.7546826410099876092 0.38400940833356889126 38.452955436962383828 14.529623343373557987 -11.678330978856003242 56.652135272852930825;3.0287549843869552291 14.164105364484338878 0.75303768549195693183 -1.0612942861492071778 5.5886105950645381668 -2.3510326893738611176 -0.26783387614551701361 -0.08501283402245132792 -1.9334389856848441802 -6.7482402562431449411 3.7663673697404185248 -2.5468260669039173294 4.9592552448237645635 6.2055704743400550427 5.1379760931888425191 -0.53251157667441539889 -1.8213757395056646704 -5.8754264682119741181 0.67259454531872608918 0.53884054388761792875;-0.30331415565554642466 23.75721102538928875 0.91541033520629355724 0.22613188091616248654 12.354302630232959714 -0.73552965415946514049 4.1865888727131430613 2.8262598299895493881 2.9856209031259588471 33.393082275721063468 -1.4286341213929310001 10.701218315815930637 -0.58600963701707309284 9.8175733282072492614 12.120158116748196875 1.8526216546962899212 1.6873219691282899113 21.070524570277676446 22.888015758345790829 2.0121895326710377283;-4.2331120562362256621 -4.1401892002171791773 -11.23576505059164532 0.83879598513304920449 -4.2373365170051480177 -2.642170815768301928 7.6592202925921535339 4.3024755542224344396 2.0105141569212481301 -0.4718769264143719 10.106444391831356455 -4.208780882379664412 -2.1605471025969720245 -1.212867554350914423 7.4524361256837652334 7.4789807688646785522 0.76493068016544296661 10.225905460831407723 -4.1786717965765332039 -9.3566411210450848301;6.437547107788466505 6.2427725007284085024 5.6428336157073299972 2.0066852506384083021 -4.7534691327309550601 -1.3058822977670663246 -3.162607771887214092 -0.28564052345310070313 10.473488630552342471 -1.774810704448275489 2.945747390192467563 5.2003230323537534829 6.942895317485263007 12.534962335524120647 8.317636020637058536 15.655495451498852688 -0.73324774556891514354 -7.8651621667031124119 1.6596838491206700184 -3.6964821498744013795];

% Layer 3
b3 = 0.4186328628177523048;
LW3_2 = [0.0025428022780115081686 -0.019492361663563273488 -0.020272913475136735556 0.0092259662564575933841 0.046507848089880461651 0.018965995763976100513 0.19994454957680263263 0.12860619547007495767 -0.05193807148281434638 -0.11445509281264239221 -0.0027846282925847610362 0.021439275541070640374 0.0074325493074973523339 1.5251698333063061774 0.051023985528987070293 0.0045841031125121124365 0.035939900122293193252 0.0052133022368010314279 -0.011805868766707385692 -0.033026470481652310529];

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