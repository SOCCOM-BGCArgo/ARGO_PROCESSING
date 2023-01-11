function [Y,Xf,Af] = ESPER_pH_1_Other_3(X,~,~)
%ESPER_PH_1_OTHER_3 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:20.
% 
% [Y] = ESPER_pH_1_Other_3(X,~,~) takes these arguments:
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
x1_step1.gain = [1.0000770825255;1.00022983623563;0.0140762793578401;0.000302617642608564;0.251952632905014;0.0588888685539285;0.0422119037568594;0.00432440760873659;0.00826143816767911];
x1_step1.ymin = -1;

% Layer 1
b1 = [-2.3988326439376699284;-2.1144534882876260262;-1.216338603029536225;0.72631272696474757922;0.16352976003220898171;0.7375964595815637237;1.9006698501641119403;0.022561025098709535514;1.2952298504648656063;0.27062562371487286494;0.73743269622462248591;-0.37563464965922688243;-0.99986989122448499678;0.16159174855093771939;1.6762576900857304452;0.63180078707873554844;-1.1296254986122862274;-4.9918544432650922005;-1.0970390294527452291;-0.0054803351857748029002;-0.0066158310607637551073;-0.75549912912600425408;0.27764866201710247662;-1.3651055016392996677;-3.5196604792576224874];
IW1_1 = [-1.098329124611966634 0.70680703717584225654 -2.6980036539877887414 0.56936200918099821866 0.60918820067864132284 1.9879908929214504809 0.9724266597394247702 -3.661859069096827568 0.09534633173266060524;0.16566045571805884662 -0.56172409794918842785 -0.49407846526730003767 0.23263537657252064683 2.4359030886253427717 -1.1517939182933230047 0.044129526603666605067 -0.24605798973489037551 -0.72275785032175265954;0.98464678696031415761 0.24004461096560023026 0.13737760973576768064 0.26624438569807412636 -0.24395416447506512725 -0.27592854382363846133 -0.38186493016297734515 -0.96813562574948019979 0.88454136656595605537;0.074689692802209703415 -0.0031768439162245401888 0.70501777339351212603 0.24610985803834825592 -0.55181014343523115961 -1.0048480313939949227 -0.098119300888000346439 -1.7648027914821315942 0.22213464240194283339;0.040825498935618931284 0.12456419122294570823 1.0497421654117251588 0.11100123092047126105 -0.75992445532079322401 0.68267745061566076359 0.150401549098747922 -0.1910693027379369302 0.727713702839921317;-0.59047703360893377678 0.47655708698452292627 -0.22636616925797040101 -1.2181350837492879169 -1.4991656729144930971 -0.11256078832841270865 -0.65274280051891686583 1.5068128448123494501 1.1847675744121830554;-0.61865078327863987084 -0.053857709454495049373 0.05813979268129305239 0.41544457600700923638 -0.75562262059337848008 1.664751032325125335 0.94202663492594518324 -0.24990914984111572972 0.62984486412470486183;0.30360729294124594313 -0.51156268711279284656 0.048732798404053132013 -0.26516005585219448726 -0.54034000990169206435 -1.2673900332660790191 -0.40407793472001707347 0.19864110703876000641 0.13585321468169497816;-0.87638716755586265617 -0.97081517678657913706 0.35557766566839782429 0.39955355617854187988 -0.17444040196762242156 -0.36835934344153276054 0.60221741330641120893 -0.76387439558271941209 0.35453626751271233308;0.59981346965884152489 -0.34720875425597919062 0.15561791540529393285 0.2442077729432588773 -1.0869101643551237757 -0.16886326319904121362 -0.7443232478398683627 1.2832346374828687008 -0.24506807373109784631;-0.063603732534680809674 0.31833108412876887083 0.28324850306554882895 0.222016680843684866 0.74742422831849097964 1.8707967766526778508 0.32935166064718446322 0.26839050444002254814 0.3716613466694784762;-0.8470030245599732277 -0.30896900381210284925 -0.035476303519588560376 0.22525847705172163948 1.5766747414847035458 -0.034567209026870333866 0.77327146839987537241 0.65738300603876087536 -0.29180855458305954286;0.22024009196577964964 -0.24059024029129696953 0.70164804173235351126 -1.2889071299572647522 -0.15792919896499968413 2.2232194617792009339 1.7709002890303309208 -1.0270368129093909726 0.18538551834581815103;-0.36645956518191430407 -0.72793305630962912289 0.87733506698072072361 -1.3851582172895442469 0.092859689531268677087 0.043377596541293347854 -0.0019207582791977463599 -0.68685047347439853738 1.6422116622624944871;-3.2971514097283578515 -0.61320189025644500358 0.25932040230588732088 0.042982840364261136468 1.7099625017535649008 2.3800106688865212412 1.529518209853839128 0.20838849745201734609 -0.052852432419222043769;0.16823904093916094982 0.24816542230822657977 0.92097547792779999032 0.1576883972910535503 -0.11433727766787776803 -1.3016557488372295648 0.033480857753502679675 -1.1227865262325347206 0.090249163028842963041;-0.5464967805492457753 1.2594395908535838124 2.2594919181479777848 0.17313482589392165112 0.75117457949853294608 0.10549111332573497812 -0.9656221102457699379 0.28456649206110745665 0.48824853780691085392;0.096395989616109956089 -0.73239745158063473962 -0.075832203231363096152 -3.3485447182175533953 -1.1286842073906442341 0.79446808570421523221 -2.735241452692222186 -2.3521894034617187863 2.8533799934029961953;-0.71190009235783879848 1.0324408278024355123 -0.36380004378310520918 0.18753049801161250643 -0.14256762080739421306 0.33191511975951681901 -0.037598443756518398762 0.068424607127823791619 -0.21590315296712389581;0.48448589557250865134 -0.32759416741836583364 1.0099552837510739067 -0.43445548906419928503 3.2894037560114171015 0.098847281468151679262 -0.4396625444381642156 -0.27529182912479399636 -0.58189530508763076533;0.95423068468893590399 -0.69555667602793658233 -0.083660583270370153519 0.18974477816410717512 1.1886515643590409574 -1.3923563180740701206 -0.66613520282917670912 -0.74126202145700437196 0.57856921294036156578;0.50247510955012930634 0.082563882406724298235 -0.5258597657913133272 -0.40026853938551204548 0.011907132358558525712 -1.790378481618494666 -1.1604890829440990263 0.63892679873467850271 0.43391738410623470479;0.092027019206315646693 0.14444094538699805974 -0.84998215208058647274 -0.46900242568183969638 0.043432908165991865324 -0.33969595080915077068 -0.17771426811906898546 -0.6381879338183837902 -0.16235108613295168301;-1.1104306897563043233 -0.28859165076554665896 0.25805335388384181838 -0.098917652188016444437 -0.068817318830542886787 -0.82437586462099177176 -0.42882711645019683244 -1.5164460644694699454 0.92630000097279308058;0.043176931550997985076 -0.62069669660389914512 -2.1797628803225213012 -0.98217315698798535184 -1.8438266760702448632 -0.71535669663233203419 -0.98832003581671179493 -0.93756713771518085387 0.0276249642480242176];

% Layer 2
b2 = [3.158076912193076069;3.1292136245400898886;1.1062436100795707272;-1.1741407461145421109;0.80484424593144054949;1.7357979211385106133;-0.17987004818265067696;-0.67951786452903684133;0.93909632611274340697;0.76821095228690083889;-0.78017901218233431937;0.19427567012109692168;-0.70450216071298055187;2.5987746383393344196;0.80694380008565991247];
LW2_1 = [-0.31696745779559648559 0.30154717886421056328 0.60861207360054503379 -0.38057484818188741515 -0.45178981699619641388 -0.21819179608389843716 -0.22825770230650327397 -0.13355806637007897053 0.19712579948957995035 -0.26140986940795468696 0.56390433279380736131 -0.37289767826806047291 0.036887603865171099404 0.45732611654039601046 -0.22440981262362563742 0.68714102218605110917 -0.61551012203165167413 0.14460895852572247522 1.0554929525677634317 0.23147449252441504308 0.49615237341323548126 -0.22112880095489612087 -0.67721021565645178608 0.30496498814170741598 -0.017620555653859035911;-1.1740889444841209333 4.7784482779971586552 -0.89736240191266369859 1.2668479584668841387 2.9725991201755248561 1.0459655570823340742 3.9393898659318637989 -1.8741926541769691195 -2.3370357073587966212 0.66427183221313657047 -2.6405078470769485222 -1.5569635403225567938 0.93783597215402680902 3.0457696602185300172 1.4955179025386526881 0.4367197350652087251 -0.3210050387689924456 -2.9526005799038412825 -0.20135369711254241798 1.0984942225107985347 0.47138672640717393936 0.56399671204440970929 -0.86907852463754942107 -2.5086477196070111617 1.1385630169510736476;-0.35563504556263036971 -1.509465219954780224 1.3134370854703070197 -0.39150430632377958284 0.47227570866791845905 -1.1288139643513424026 -1.5726720864547638623 2.1935417305910247521 -0.22881165040090484242 -0.77534099023332492262 -0.71212736366236384367 -0.06689534346345452076 -0.6997669659668405151 0.3617399350611986697 0.53505192596664441496 0.34049546519335838202 0.13455024650229871486 -0.021749178745501078491 -0.46600717500473298749 0.50107690511932478916 -0.073480831357356626854 -1.3035383705515326547 -0.21763701022935183493 1.5484793580952356251 -1.5018610974970154626;-0.1670063787267550004 0.51639086067302086835 0.89198643758172446727 -0.093680298852538201881 -2.0353763796911161776 0.39765684956089975 1.3810185541146056121 -1.3309480677637108847 -0.19691842636140535761 0.49501489177967777922 -0.91277369094781335424 -0.19436518794415208466 0.43939768463177292235 0.25383764100965000576 -0.012254172178201083038 1.4771423301560311359 0.23875254376187746153 -0.06667305044786633883 -0.38186207280308964718 -0.059936825676887052483 -0.31502730205407880604 0.25978901256532738184 0.46822505861773833225 0.5422040601003167426 0.5907936759223844092;-0.22814991778318152726 0.35978683185323073745 0.62134498729157172647 0.75321198876886386042 1.2129619601640304616 1.2665916103949603766 -0.069057400933631302165 0.33855076524197025289 -1.4546994432441062717 0.85091312167265176214 -1.4590755864986368895 -0.51604560028201074129 -0.52931895732491440754 0.59444808677027205501 0.89761357545512132639 0.7390830724105919014 -0.39179535196089820825 -0.7679529420272050988 -0.077671360550483223295 0.1952579428637479797 -1.7829822428994668293 -0.90351375580936044152 0.46726033346481415931 3.4379472469223153475 -1.2154659958150442201;0.14155580570003625795 0.044312565789360602864 -0.77390616896357000698 0.63807686941945351844 0.5947207090548681796 -0.19736279934394568292 1.0769528395041756408 1.1331559884503559221 -0.027636299547018069317 -0.22723321136694896172 -0.29666223281924702926 -0.0094771351759814700944 -0.15537710339175586638 -0.77973137121188484944 -0.002326144920821626022 -0.73476080183750946961 0.15670748063685813189 0.315526776042899193 0.50202960113662298269 -0.73974464700260655903 0.21643705897202433763 0.40914361939789356537 -0.14050616899983653374 0.090491255601308343004 -0.11431389528388626042;-0.06862334349785313703 -0.1260875177479164333 -1.2252438663845972577 -1.0990829491889222425 -1.4425510773101457573 0.52150216373386237834 2.0868155132367598448 -0.072380023417521366369 0.089698193636682341245 -0.61225860269973231276 0.50247401713304506998 -0.096736455873854765297 -0.24212943265429259787 -0.19138993454943944994 -0.85119939726459736828 1.3961835983928887472 -0.047342947843250736406 0.80970958511009616387 0.18978635531453599827 0.8704161151575747768 1.0661965490595186612 0.91834909380493823239 -0.83816086145657309192 0.30328402739674420463 0.17196070088867543291;-0.32093223671613291259 0.52654499501559304253 -0.71405631872460673115 -0.48375208708454126016 0.6851902293575723446 0.51877441353710196381 -0.08932851864085372684 -1.4172083358678524423 0.8790077966628125905 -0.47569021113450749016 -0.2922910279387330168 -0.8150780451142199956 0.92059956205991788636 -0.44773685314995698148 0.45868025794592548472 1.0012039951743485489 0.081792577486836534617 -0.042075314011120225344 -1.1066569390635139225 0.15724078100998134522 0.37843793119688340365 0.66377687723004130671 -0.91109484261820139483 -0.63077803987588987322 -0.25429729072925622013;0.040154708785213648337 -0.71037669804989700406 -0.17012902605513996468 -0.10128455813290941634 -0.9739096578427933748 -0.34496691861341988172 -0.13656358016222377993 -0.71570188141219626132 -0.39934609062235537635 -0.070424060371416144499 -0.49713765141270138681 0.40740911684136821291 -0.62898596276630946988 0.030738408373691707837 -0.83922720243624626679 -0.67917743910962180554 -0.17307789965531875098 0.86052717104526066372 0.12972420724667385605 -0.32154399199747257798 -0.090007310789469371914 0.30838384854038347616 -1.8502566527047765721 -0.3268415241680289185 -0.29536789361734194781;0.22463200903958735477 -0.058309742105335612972 -1.1458150465694814102 -0.42404396770335500699 -0.77211418918858254479 -0.31914162649813931916 -0.033584448304265344942 -1.7722410793614720159 0.069745835366523623033 0.063700216018519995043 -0.80483868584473217123 0.32049628627478854437 -0.15867253255312696592 -0.071872057661211186108 -0.78252857017833210573 -0.81578906583357602145 0.76363983765260923242 2.7416998458139532069 0.32743146192307359676 0.46255643747495739371 1.6384921310830793395 -0.62958249929695742075 -1.1137868135673216763 -0.191878609889868601 -1.1914147489896094179;0.6399493488935420471 -0.56313660002754195588 -0.98518902254077811431 0.64416767293462451338 0.015942887949338336445 0.18421235173202601954 -1.3129520456578376297 2.6997976254260369622 -0.037669322205922756153 -0.88506955927093811098 1.4209522528460094559 -0.42808448149511119096 0.030703599190708211197 -0.44664273541476562368 -0.26354250061857031895 -0.88972792181873372463 -0.24362002129525345273 0.7377260074080488339 0.66519060358267623201 0.21983085129308421735 0.56309932249596761356 -2.0888857224341750296 0.88375301334639699835 0.39456136696319538126 0.22453722979199058907;-0.065064278408231229767 0.1152925078037846246 0.42221082376065055541 -0.18968483771902019952 0.20030976547685028222 0.20542886594142054091 0.55884685133833722492 -0.18804206135622192342 0.41819024060122844677 -0.041354087704204359199 -0.25921838655418738551 -0.18204129348261902099 -0.28966886449649331681 -0.48295170733474018787 0.076640148778020417142 0.61338603001236668977 0.12838741202937015151 -0.099598844199511010755 -0.18671434334362976837 -0.14480507187909649747 0.13147311756294280394 0.54778538123757947176 -0.072559906768325688198 0.15330690630686655274 -0.013884643849116683548;-0.76043413938913306804 0.67361844798765846409 -0.080445515455990268538 0.78576576123996022538 -0.45958538464941950386 0.68212768168779835953 -0.022393050053909076269 -0.43677378375283454348 0.034274987107762978911 -1.190789862965278445 -0.1063769497933467062 -0.59521582745098133849 1.2471725466561143225 -0.17543542674378220791 -0.020576020783247921803 -0.33053983779801809906 -0.080033719724674740048 -0.023200028467618045991 -0.77849206550141747396 0.11744052913344912858 -0.15773039924545553814 0.60709565430637191241 -1.13272562948964306 -1.1862382429917661675 0.091170370438998171725;0.10186262731999337627 1.657868044159793941 -1.1476777849369486084 -0.2413070936992929294 -1.9495692049703727822 -0.76419113791641857247 3.245712427316893578 -0.002605942298634020271 -0.83704016951648729794 -0.06622852067405810006 -0.85595433846625079699 -0.68846092854644946879 -1.9285422038786821375 0.032358849049767442196 -1.1960509102418943606 -1.209907742008832443 0.41059334606626479713 0.25190486281620710907 0.28801222555488459331 0.72019484075772810527 -0.57934752597273631025 0.50440438689648581239 -1.2440946937184209631 0.88245199130104878815 0.90137701710603046479;1.0516875738199165102 0.040410102295392602567 -0.29304569495361665332 -1.2308145458169548192 1.8705109312821024403 1.0535670930418801206 -0.92824784914169888861 -1.114702044448244278 0.38021777262184985879 -1.8656498661097689062 -1.1171852371603507859 -1.08654320060773113 -2.3138427757965871479 -0.78775358942913031868 -0.4386116635935190522 -1.6903804119649101789 -0.72660209405196207744 -0.52118210180790591757 -1.8778668436100398598 -1.2899989855979538955 0.88400519140931621553 -1.1563681773387008977 -0.91985487853040570627 -1.5409693458036384239 -1.1983607185986340227];

% Layer 3
b3 = -3.1286609900518786986;
LW3_2 = [1.8614250524452635727 2.1566965349076427927 -0.11368519951531658174 -0.33797841698305336466 0.068008270933458556118 -0.56178334925572759317 -0.56089958665612871247 -0.59575121317023116418 -0.38817081127634672777 0.18131769155485200073 0.25489922873348053178 0.87545837346634924891 0.4220406286256633277 0.39991948019649764534 -0.25121494904417007721];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 1.90239167587468;
y1_step1.xoffset = 7.36300175998845;

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
