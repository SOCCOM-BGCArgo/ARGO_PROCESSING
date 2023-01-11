function [Y,Xf,Af] = ESPER_phosphate_12_Atl_3(X,~,~)
%ESPER_PHOSPHATE_12_ATL_3 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:30.
% 
% [Y] = ESPER_phosphate_12_Atl_3(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timesteps
%   Each X{1,ts} = 6xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 1xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-0.999999996490807;-0.999999994516886;-40;0;9.1619;-0.52];
x1_step1.gain = [1.0000000017546;1.00003626070655;0.0153846153846154;0.000309885342423303;0.0656325902039552;0.0145232735458572];
x1_step1.ymin = -1;

% Layer 1
b1 = [-2.3047411961521517831;-2.2717949502678855822;-1.9175692918075555315;-3.6514982750885112139;-4.2445063553285748981;8.6556256446810611038;1.4267058935958878418;-1.8113434506951344716;-6.2860732141493702585;-0.87971530510193518548;-0.13604385693458054263;-1.0082746416435890424;0.16535046864784308518;-4.913129588860481789;2.3299578949338428657;2.4038880940627880278;-0.5944740098075824708;1.7070416264465573253;-10.295374754202924805;4.8956762030350295944;6.1896420325468870516;-3.3286450093544290851;4.8081723427372624613;3.3966644634749143528;3.5456747792338427772];
IW1_1 = [0.78789228556420598526 0.11027278647940184919 -1.4181120603730625884 -3.8157772992043130778 -2.1555864413320189321 1.7363988258190814484;1.710041989736449386 3.0532671378202684842 3.6171420382803693983 -1.1318560588416122226 -1.9220014409112009535 -0.4586628489016515986;-0.27324283608389909883 -0.36367168091498180971 3.0174569748350466014 -1.4753552703394756662 -2.2278902802652358339 1.4937121019799548538;5.5425203037052153121 1.4896404601475310425 -1.8765846252114046777 -0.73070235977764608215 -3.7626695780064167707 -0.00011140683346440999335;0.95736928819548572367 -0.41124849495922777276 -0.92594960705064599527 -2.8493035362710816827 8.6121042853424611963 4.7151410337953754137;-3.2182344053421130248 1.3755296665762946073 6.2349679192822575757 1.1049482591560559896 2.7090464755131882235 1.7676407093881734145;0.51601899527492467445 -1.6121727368058771379 3.3041397140877122318 2.9486412106425854418 -0.0017114099070988670476 -0.088108903282133194113;-0.35531713427250616322 0.49045833328980348309 -1.0743050316499813679 -4.8724096453489726954 -5.2237668075797332179 -2.3282382742159262179;2.2974378830240778093 -2.9686337399773270462 -2.2673160449032510044 -1.9828006167577727492 8.7590542460274178893 7.5556347751319945516;1.1373429236434577927 0.29482139134361923238 1.7326386043361798883 0.75650552913586366355 2.8957197714600675376 -0.1716063687944510141;0.67816670245664134598 1.1263152407587417869 0.75327394105977685257 0.74943459546953083983 0.8315086180835266072 0.59654053147080277064;-0.1885001155920605076 -0.585659570070164337 -0.22580547416029039809 -0.69326287541968911743 -0.52202211478407578582 -1.0085789423210400617;1.7051577497130552263 0.21517979354730712438 -0.26153071316892234632 -0.36164294163514260072 0.19164820954574823819 -2.2812885914604192905;0.058200195915580857364 -2.1939718064015809595 -0.84929594868817803732 1.369469986746689294 8.968199423074167953 1.6522202863601573863;0.0049853235263881152256 -0.7584091380374766711 -0.49087750057568896311 -1.0831811062933431877 -5.2910659141147817763 -0.2223014610640978439;0.18735448055371420328 1.1700705390386996196 -1.6136605375307853105 2.2753982641810326371 -2.0550296530886336122 -2.1042950563251920926;-0.11507915806675098058 -2.6346538630486375787 -3.158096080307243092 -0.83177444257560773089 -0.3985668904407482227 2.0832861387583312407;0.14958759009835309595 -0.36492325546772602651 -1.1691938850508942771 -0.71510679278846611684 -1.1863285340536426915 1.7932857191419653287;-1.24397000997203655 1.0791964359038725974 1.3833004526562240333 -1.0205388018015846274 7.6856029627412807415 -5.7532277464590801941;0.87157187978040839837 -0.54316960029280447753 -0.44996277309083243523 5.5084133133553097039 4.4482418305186968155 3.2296488369329581225;1.2902772977799152887 1.8697139572721934719 -0.44420268562970255388 -0.41094618713569958102 -3.7768616787364202025 2.7100342300259119277;0.3899351677336178601 -1.813023584967490498 -0.69639556123477752703 1.2164890655418760801 5.8410012748732569321 1.1677529873996963516;0.79901992639451435707 -0.46647421928252463807 1.4352616909213726792 1.1419520321006648711 -5.6027951336946379968 -0.33266430986646361445;1.472508015537358661 1.6048882562288884213 -1.0227725684057784594 0.64028323005962917147 -2.5521345342910266574 -1.1495626597600148155;0.5050288552716224455 -0.84844246316470672831 0.36242692629997386167 2.1146853008175541255 0.1564067941480777324 1.8790663406538572477];

% Layer 2
b2 = [-1.5922769443616391349;1.0864972985983110121;-5.1235858811939474222;-2.343686834774737715;-0.50522602764639112927;25.215204000365396553;0.76470209499190411329;-2.0593232547677224886;1.5097219082080344243;1.6757002137351697524;0.13810721796446731591;-2.3886368079285120736;-0.8571825491937467012;-11.654230947755626957;2.3847191250006085639];
LW2_1 = [1.2182313117196814645 0.45250001175996334979 0.40006283050675350843 -0.18334317697681526416 0.21886292451376149204 0.92647688858743282303 0.43263008409050374148 0.43847610730248220978 -0.66015857618807960439 1.1649448015851202598 -0.83303626210458947554 -0.97213421769780672399 0.2949469204072420192 -0.38191511026576896448 0.19047322793456189505 0.072216498370212101054 0.25254995459425916282 -1.283100323328335568 -0.077334972893406919714 -0.048367110318607853259 -0.037568047987225176854 0.14102350857538117901 -1.2219328889298186613 0.055853211847075805163 0.62084396506645767033;1.3423069781965280001 -0.068063613108441028965 -0.69661543444309859119 -0.69768592504256310427 0.21177127138765294845 2.3290155266592327266 -0.68716276541564436098 0.44889462404175445309 0.8724741421583777079 -1.5784746825714508489 0.53905947480397653493 -0.046298499544963608865 0.73488268630853625929 -2.1782798057260555069 1.4060910422872168191 -0.053505801103981599776 0.40460060584525753269 -0.80372476345248822227 0.78582774886943529413 1.1484609569341546198 0.40295935576408875578 3.0259408591149443701 -1.297573886812779298 0.29129275019004879921 0.27354983235393387497;0.38725498193067453556 0.38236599661054365695 0.49096580114873555889 -0.10537941855854428275 0.062385867860992372247 -0.37757521117632208973 -0.81834278185063935585 -0.11782546165493740231 -0.074137194327399680294 -0.55830679707498098718 0.41060325147407988888 -0.077022709954469345539 0.31826843113533248575 -0.17571903629394261537 -0.43826199091315565237 0.024974634156930294177 0.55278598273680379371 0.098753600031217983468 -5.2705147858968324215 -0.59279767404643901596 -0.044226453657586831114 0.23614467779917555634 0.15896490664609438936 -1.104886979450788731 0.89684183876068379604;-0.079381841204175379589 -2.2602473261171698304 1.0971532927992668771 0.86759699876023899101 -0.054344045598080675807 -2.1234453088442926116 -0.39812714653107023066 0.015068047067869273781 -0.93409119320666234021 -1.082115260584803007 -0.7398761136554798945 -1.2416722760717955332 1.0335431434588588928 -2.4082482859897469396 1.1836316083855722425 -1.0761301104470730738 2.8804897088007583861 -1.6973821299481703928 0.32112353436668095163 -2.1722173389440495583 -0.33827641055665497172 4.8399869438882658912 0.76196925437554474669 1.1997794593444646161 2.0556252861336563598;-1.0102080531387676032 -0.2344050615712654495 1.4852844446098452114 -0.075956414067350680464 0.097186642744939433225 0.95031568174912861569 -0.37193193874179542036 -0.0024645815587142193975 -0.51757216937650729616 0.10146098933138689158 0.19594092426803116913 -0.51381967651736004665 0.25354973957731274936 -0.010738247937961529002 2.0333638357778713868 0.53651570777731216744 0.48981830565719286508 -1.2748075438977140017 -0.81650669280603238587 0.71991340255195579445 -0.09226317297193459277 -0.89041634974323935481 0.87927389410716338958 -1.2078946553031371103 -1.0919813933170254838;2.0805366147252901676 0.89965939513161496954 12.194295611319262207 -0.21791228712004145907 -0.11521397865340644862 0.29390840467655898749 -1.0759828070288004476 0.13045749027708508638 -0.47642237497641976018 -3.0022713566385812456 -0.77013322697932495853 -0.33954695116325139814 2.2153666064072270458 3.9942298321081053381 -6.9614296604221754805 -0.28542549415689760783 3.3980100566371840287 -4.9863071421896965063 5.8701296393880344482 -5.5677809917362486303 1.4129968109404187349 -4.5071698432564444303 2.7238907607298887825 -3.3304915937077237409 2.7838528278731788923;-0.93600244488028994105 -2.0087736532332352013 0.86092939424437986418 -0.29126416586654901852 0.76276828869052115678 -1.8236542630621563887 -0.19072666223863377066 -1.2922121480698680607 -0.93974003174388320847 -2.4381528080728425145 1.7730985070992342223 -0.75541579563587168966 0.7721664933494216676 1.2489745663142204357 0.034927247804321777391 -1.2711085894441460642 4.2012530299161063851 -4.703593899986806548 0.19223251945909738958 -1.1179861465274256727 -0.025795841349238613632 -1.3008020838899059246 0.11270928689505108067 2.0443734320846109753 -1.500885044128968504;-1.2316105546385764935 0.71560047938169601967 -0.042131393874360045793 0.89935223109346762449 0.32175012973876526701 0.034911488590675242294 0.30571500892946146255 0.83812567104927848671 1.0082367787765749156 1.5430098974458257288 -1.353243235198749117 -1.7559507814030752559 0.30162259140799463353 0.94435802912263933084 0.8954368669478307563 0.091552631807148995846 -0.18837208665061880297 -2.8824490410961192133 -0.41621479058776233995 0.32090002180902232887 0.64461738328225726136 -1.3599251632910089871 -1.1017103374055317033 -0.96339271249179891932 0.70277354138744341228;0.19824988201672075205 -0.039273530866308345444 0.94643408705127363145 0.30931677274516022891 -0.00096109458236877777851 0.93247969766225247135 -0.17006681789976249575 0.403875947619491249 0.39825877444124085924 -0.034894378032025843983 0.15195650776491606559 -0.9303033783476812868 -0.12015919700740278209 0.18266609383606038919 1.5138150499064262355 0.31638487253941544042 0.18178418324581657739 -0.92403860725115172237 -0.83808435657238955496 1.4789022140381011816 -0.3133607090923845595 -0.68891663848049067287 1.8583231557582540194 -1.3233178199663058194 -1.8442992434943543678;0.61168029656363565039 -0.33176442580152776252 0.23684079168301483409 0.94298200963552558651 -0.11192914867816750046 0.34127362688705398241 -0.028401039904054462687 0.6325926603891911526 0.59703837665166981097 -1.0805766928727813525 0.25325305624074001543 -0.74886960384329193108 0.030292526200283477478 1.0459634317552908289 0.76113970985483847542 -0.082981595161903079494 -0.096156353896825891292 0.43086264179648048334 0.60306000393259362458 1.8840270378996650802 -0.2653960887312556971 -0.4696374814816368648 0.87497075747620745023 0.70067638086591377267 -0.88693951496196721251;-0.13858362488305067672 0.48579935953646552482 -0.45127518393741405944 -0.068329616193640002608 0.0035841129817557335416 0.18668747268048702104 0.70154441480688478627 -0.0045281476513937234107 -0.23540715371749143525 0.47925428080798809782 -0.59584022115105983453 0.0069659398580865583908 -0.30046333605034036829 1.0185358917984999039 0.072450080373606201567 0.50865005738978141814 -0.93807441656287338105 1.5738231362778429823 0.037370896082639569469 0.080542609672194967474 0.10300280523918059483 -1.0128000945803892208 -0.68881099059050665545 -0.40252268077653202694 0.27949711193725901293;-0.58517301787180431294 -0.54275082157640908154 0.36492627896074192329 0.34648735942416974964 0.62635496792325762172 2.0309492454010813489 -0.70649969995770756093 -0.12693332797122000577 -0.79827253000724240906 -0.02430021665921979282 0.92593927800652509053 -0.21479605795318071348 0.018949171135454436626 -0.73514350481277057625 1.2669663068015062368 -0.59826746595465774003 0.84131328382247161368 -0.021524630142303921576 0.53164595561316185801 0.0030911226199914708899 0.42771579849349644631 1.1926307212381570544 0.99293818403100653214 -0.11894640846566474635 -1.672815413889046976;2.6237291871712260516 -0.098745266328340061679 -1.0466760822815381893 -0.27341804395771684977 2.4224636435910227306 3.8471902915686682078 -0.83398349062018262146 -1.9354726490623326285 -1.3391859730772064374 -0.10537373056178192532 4.0115741902051267331 8.8931318898318423294 -5.3227198300917493157 -0.46697181219238292904 1.3114220911832343841 1.2671700340879219482 0.36956246711897855484 0.53078434871954971452 1.7985296877770304835 -4.1696697322505524852 1.2534103273347476648 0.98295735175264675743 0.94938194714987433898 0.44190192820896601056 0.20124610344325230549;0.065878617456972057842 0.050611044208813328604 1.0268327178098624053 0.69017134308740757387 -0.083498812245581463531 10.037507392335777823 0.51413974945786700133 0.078858062826881597518 -0.9888000804872422167 1.1899305080608346685 -0.87326847617890746012 -1.7204168277206677296 0.63050992393039972939 -0.83089367388974633144 -0.67241794243426200595 -0.37433524234302134603 -0.054797035597144751196 0.44639782986419795918 0.056716650400332566107 -0.77312499406717272166 -0.88233209320770578099 0.17981621540617181987 -1.7198365182012063102 1.1906504618697379971 0.79814267111368975005;0.42578191091985229111 -0.43682106268701959007 0.7936350239924053751 -0.21509598030436885563 -0.10857693580685820467 -2.4300196957760880068 -0.12060561507836867201 0.72670444425672953148 -0.2464983297357966896 -0.98550942514577699605 0.38238499159777483705 -0.3164565872659437229 -0.086588442870166343335 0.60559634487976554684 0.52484592156391951523 -0.053968133068278342779 -0.005255533203383629115 0.62597914461123249641 0.78378615574111787367 2.2698122797381947713 -0.56441138438685645173 0.049658174478985857292 1.3132258345950467149 1.1648592907198669355 -0.97863814787289471475];

% Layer 3
b3 = -0.045798730871322121394;
LW3_2 = [-1.7656839333055072139 0.45490684379215745592 1.1675577598650586086 -0.18773033589453907988 -0.65822871564708573811 -0.12700752219200392434 0.20449764568462694414 0.84069363469354851937 0.74065307480131437057 -1.5854596095420927515 1.3258600129984985383 0.56703700275774748807 -0.33537368952232077257 0.93151126581789545078 1.0920862726018265576];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.604229607250755;
y1_step1.xoffset = -0.03;

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
