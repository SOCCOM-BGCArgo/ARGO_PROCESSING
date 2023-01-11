function [Y,Xf,Af] = ESPER_pH_2_Other_2(X,~,~)
%ESPER_PH_2_OTHER_2 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:20.
% 
% [Y] = ESPER_pH_2_Other_2(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999540494077341;-75.751;0;29.715;-2.05132667284751;-0.9;-0.1386];
x1_step1.gain = [1.0000770825255;1.00022983623563;0.0140762793578401;0.000302617642608564;0.251952632905014;0.0588888685539285;0.0422119037568594;0.00826143816767911];
x1_step1.ymin = -1;

% Layer 1
b1 = [2.0511053094745683723;-2.7643812674086873926;1.7177484167424008632;1.3738002927948254062;-0.70771660941636338471;1.4895870162691737804;0.62954872316151166967;0.19836503866403878971;1.1155064653511159722;-0.36589865433705881514;0.70694668005640015629;0.24559638788295626299;0.95454034579949276207;1.1010094452595078351;0.46667423066348284744;-0.54210870966154922623;2.6295065821535890116;0.75947027632810570896;1.5014020017185669698;-1.0682903761473732107];
IW1_1 = [-0.25136426694977448415 0.0068570470850916401262 -1.4977099376995506308 0.95756277070267581397 2.0436149460841428471 -0.31416599490847618714 0.63640830564290129523 1.0416638354435168168;0.53783011055021823221 -1.0664885553320775369 2.0499645533724111246 -0.77933934386326952914 0.19919510776972940147 -1.1381441927766209332 -1.0599664865849389717 0.68321925013984019515;-0.091935915647500282555 0.53603181422934287337 -1.4337572146087813607 0.27418339125774138232 -2.7493194781782377412 0.41710204499789627075 0.47781518579202875152 -0.32688600158160080467;-0.054510045853684969963 -0.1351930079411350627 -0.73402804446979252884 0.59251035853983291535 0.14601920995559761196 0.47111239113332092909 0.48938167106434715681 -0.49372770926564973903;0.13885665165978047786 0.59675275272625138978 -1.0214313522109865762 0.046743960307589052516 -0.94748230921626153123 0.35922909944981123997 0.58927491711668267271 0.84988796251846321361;-0.82387281186385608045 -0.093960784224898549755 0.77148809101054249826 0.010533616185552331299 0.93160165768242286788 0.98162254799991555121 0.65643359168824066963 0.98267694709725350588;-1.4325212022267870271 -0.13780623455028709334 -0.39730290000499379754 0.026611312157689328423 0.70271427739375846855 1.1591252880077298482 0.76062065982692772526 0.196197632682440809;-0.25565610365288804484 0.84407062904765040035 0.41473209697776225457 -0.22289578508435889392 0.71548284579557264173 -0.64803165018509356621 -0.36218627132116554268 -0.0080033984499265926199;-0.70405508521238557851 0.70311673379050976251 0.11029523392059999154 0.023006548936780564507 -1.7984164968269187312 1.5714308337914273395 -0.61363640975642808062 -0.078005142866759732967;0.055486697839047777192 -0.68775082702896039866 0.30482495548396204565 0.16648917554153028209 2.0794748512002021457 0.53467380894328153662 -0.049692083258944634572 0.75060614383308532638;0.27828389057619673519 -0.49968049240409451173 2.278060203477637824 -0.61761447836991290039 -3.7120118806716990711 1.1422771150455499622 -0.049851098183165784561 0.83092464177536062842;0.81475582612327746013 0.043919039832328379824 0.12287994353326081587 0.21107849374867890258 0.92310732298744735402 0.32350848384636110566 0.14636722061437326681 -1.2916948817548767892;0.64526920352556749183 0.15468842246701616783 -1.3194248727793402853 -0.25368113772425665564 0.65278508840521198842 1.3137061530039795265 -0.12250150343520908869 1.0830390170014265738;0.86653735063452130838 0.94896427395836480123 0.3967352980586307809 -0.15420886794557858712 -2.1747330382525484893 0.14684688294546341392 0.82087292048155946489 0.21567185310494696449;0.19819328291168794576 -0.36873759984164816705 0.58925178716711279492 0.26160585945879499814 -1.2591829213649055053 -1.0583640557270426186 -1.0537982408532593492 0.36977682324451965901;0.12954849855095340594 -0.35998595629621449321 1.2440896567538601847 -0.8355999385017288672 -3.0622253795221388017 1.0999629290858243813 0.49116829963248859192 -0.68080306719918337599;0.85387485900078841095 0.032515050036180305482 0.99774638723429864839 -0.31073328301064589985 0.22508173306343584597 1.1087548458810227991 1.1684049601913029104 1.1202613241609937766;-0.088972791005864443337 0.55757598791499973778 -1.2389813940323040065 -0.059575916886963890817 0.060555251643034012365 -0.34544521282961704323 0.03243589857189868586 0.055864192502392877315;0.69973542658272425143 -0.95096047542720496537 -1.5568158271783560931 0.48498675607400798526 -0.43980151048392596147 -0.42126580989141887423 -0.84179982155940746846 0.079649929285644655019;-0.52434830438091994953 0.084035877789886692657 -1.2896473605106379967 -0.086597990500198424035 -0.25843246785613788186 -0.64034048291977696543 -0.62803857867098467072 0.27348313619736253077];

% Layer 2
b2 = [-0.40266710259428373231;0.67038114204065979429;-0.70733026580118785187;1.3451256755476797;-1.7071750464846471296;0.84556725895933115833;0.53314142209775228221;0.17107113408347052186;-0.86571204322463357617;-0.77172596409856974997;-2.4719097638346196888;-0.63687361175106937949;-0.082204134176392676392;0.73542841317343199403;-1.0671513552714593231;1.4518677959703729652;1.3844675014508605049;2.0292387803196310081;0.35833830078493850158;3.074514469548750295];
LW2_1 = [-0.3663592553641041305 0.82053566302008429378 0.12381493639087846892 -0.20319488467266225151 1.376329509123679129 -0.41354125485979659205 0.87459092894356194048 0.060748415020221914473 0.11794797922954552982 -0.51626142781508210788 0.0041264473159044658818 0.27847920541593634836 -0.14746994605452770633 -0.30370495296729904666 0.038376268675404126562 0.4455626886690259858 -0.36439182586396184904 -0.66777932941210693674 0.42069374400674852676 -2.3628664258122440422;-0.46153012512716040661 0.070876155532520579428 0.024823669258942755722 -0.78441538252557330146 -0.092993365907803049453 -0.19545330652845849251 0.12651827147547209385 -0.10071363982262510062 0.45477041099480153274 -0.60664471325315261208 -0.14972903180412097668 0.17287523262937845581 0.32684431779734424151 -0.29140394079631126711 0.053505081857093059194 -0.14984125077543564153 -0.43966174505740718281 0.35297133249038348257 -0.10012192599458685782 -0.14110935744924196777;0.33440783559373049583 -0.66982280383991654737 -0.28030760245190089464 -0.78518220217256584625 -0.46810490626151013416 0.77603449478398200956 -0.70970504980402238093 0.067985123092213814444 0.41584146253474280064 0.13620905599303409494 0.14333298066164340545 0.081286504124866557075 0.31017739780583492015 0.33618898532411184421 0.20573011327966453377 -0.13678335865210958233 -0.26798160898514439276 -0.28883736236608770209 -0.057364459675778911585 0.26460461947319424691;-0.99759332627750907374 0.38117490260876690789 -0.49956898633978713553 0.19564003403737348363 -0.46092571163679407764 -0.51855400201574997077 0.48509687599431727811 0.42143105959528959215 -0.78160031950278907598 0.55932585025023040526 -0.20519366755510298761 0.42266653042493806947 0.1033185126556544764 0.43970750029337068998 -0.23068317733660745339 0.60315244656207589458 -0.35839351987057926907 0.8629627180646197715 0.10496716966923325121 0.2845225610072059097;1.2053072825236406107 1.3060973150341623761 -0.77647458820490289355 0.8633915914547509729 0.1663281189203562882 -0.089301097630946874029 1.3163578941010349865 -0.10416764888328725047 -0.61256731352168825033 -0.78588542824621421268 0.39940974088921016305 0.31471800326101956324 0.82750896907264481772 -0.44380221470117248161 0.83242791837177410041 0.70980625332213642675 0.35717902809495699623 0.61449236831732878539 0.09929424535047898237 -0.3671218011003134496;-0.10341129252087283286 -1.9434795641970687985 -1.1131675044393749552 0.15140094424875918944 -1.0947057613440533963 0.057412099301987082334 -0.39722433335363233065 -0.26409551032652617275 -0.27450619313283414202 -0.42792920734256678639 0.89913634107848028343 -0.58080372293291659958 -0.39396973286460101882 0.57537277150259658054 -0.20675805553443116214 -0.56613628021125061984 -0.38884237266779347886 0.1661423074540373801 -0.98849055256089757293 -0.3842006398412662338;-0.94812827242966601915 -0.044741457067759470012 0.21596467002878039221 1.4113937379466603428 0.20838725095919385932 0.58423537120295554459 -0.042635403879954306139 -0.26828440370496720035 0.22394106319205611677 0.33953558383693699385 -0.33060139990737613669 -0.69494919477724304002 -0.079505248282329854526 0.10727948436097262119 -0.81142987188330173431 0.40506303589932773912 -1.3250644638707471668 0.41157427540727259396 -0.14106840461101549655 -0.73863224497991253337;0.98391909374064756566 0.37962195275791582594 1.0650862760613299951 -1.1901353235220630822 0.11207962708290777898 1.1474297145284380051 -0.60727456366462051474 -0.72691760589472176335 0.031142878886989470294 -0.57983522421497202259 -0.12347193805679217149 -0.10005889917443190618 -0.47908077394891729339 0.56804432006090088692 -0.021324424451626231014 -1.0945235461774722729 -0.49775655577059729762 -1.7471112857969142151 -0.35002770631023599757 0.72475576968746191486;0.62703583344484092876 1.1945749659482129967 0.46505895781733164185 0.8138572514511535072 0.48093783454277982958 0.13835870915883988208 -0.87426371306547534523 -0.35620173486862383161 0.067130381723656470938 -0.42704914919974662491 -1.4615335826110855688 0.47208447163223477006 0.31577735979485138662 -0.20366036320930175352 -0.023692848261683351063 1.0220528175023895212 -0.49539811877552614172 -2.6182393656364721934 1.0459004732272523253 1.0861919137007545189;1.2592824110088740586 -0.14863208922472614937 0.83863999612498163483 -0.4307343714119429845 0.90002331341342045601 0.73038151175524002046 -0.41150310228531167267 -0.47183464713963463311 1.0648987368571658507 -0.13562201216251904001 0.44257911271533023623 -0.26340227710954000617 0.062575159771550217802 -0.61513834147331192259 0.53137894020482301372 -0.24831209976727230959 0.23954072695188086128 -0.79646589086967189974 -0.11855385866089129809 -0.099384034281797969945;0.17060258087448576525 -1.746443489584527109 -0.25796893086673994278 -1.3324022286049508335 0.12900044829860060625 0.55181615727919874903 0.048319219631238932711 -0.20035849006411501905 0.80737970287288662874 0.35361911083101199571 -0.42667762553096955092 0.0097528428127002329134 -0.18960822643880359539 0.20504551647725000052 -0.72520526665168028035 0.55788690800007434412 0.12116389416737304019 -0.39894403668918154704 0.60030721301563316761 -1.4638586319321891249;0.35438948145085230612 -0.40733739565085569367 -0.59312198027085039076 1.4377255288850676607 0.15144305209251576749 -0.040328864481461910996 -0.79278121062063000579 -0.74507940997618993251 -0.28290304662767207233 -0.61205215551191816115 0.096677906453547401977 0.18683504223585670201 -0.1995743424171104452 -0.026749576643458724406 0.32043761830854711192 0.088974199514831153746 -0.418860546528926736 0.32873480067589627707 -0.030978050583194263745 0.79335368525034388743;0.76903865479138189265 -1.0144814455118471574 -0.024739367793798372475 0.38190476170011289359 0.1777528776490447382 0.53616658453784082905 -0.26647609040120751045 -0.5706919893196804594 0.63875858113309458286 -0.62098805203959306365 0.84788664958421644524 -0.087900758855239041756 0.22558397992499118478 -0.12004720160417757524 0.04376372733374594981 -0.02067290935831656018 -0.78238934694769934897 -0.74230437281303374153 0.044378790474528945209 -0.28285235553519899598;-0.38710286303014518339 -0.53523819592634291897 -0.27273376455765518944 0.93899403559505956185 -0.3102040878590274775 -0.14786874956855972085 0.93323931158978223888 0.59751865387466029933 -0.78603268900689826637 1.1538012094536647112 0.89842056654966018225 -0.20888686709853310219 0.77742431195106165287 -0.23779866481912512177 0.71484230061325537697 0.69426146265918786682 -0.39673131642400155661 1.9802873026318592586 0.53103020530142219791 -0.45686591061021991633;-0.088888221251998369987 -0.48830296197085509835 0.9489494812872263374 0.28691276369254298251 0.20751431860990060541 -0.82551775594224241495 0.85905683233665863785 -0.16306939563209235655 -0.19202778469809467232 0.040918072841108585336 0.76047641551870726051 -1.0368920163782846533 -0.43626884129125304002 -0.47691347918505511139 -0.82931792921686331432 -0.94697125893728939516 0.88770030810687505785 0.46605565688155770054 0.33892735630784692757 -0.84773718620174787208;-0.77250159807840401704 -2.6593923297181016885 -0.71806554656803101011 -0.83957885544507004916 -0.63819411025921646274 0.072165891721326549946 -0.390491852128372674 -0.088822859372739773609 0.097279624488640983881 0.13722130162376544593 0.72145750232004723923 -0.50198441918266833994 -0.3663930137897343231 0.47182347986832517961 -0.22640978071536063676 -0.75065096339625025923 -0.23907002872083985001 -0.045108157232934167902 -0.92303409111196466075 -0.28277320773908931617;-0.13243786075774224842 0.62677106769818746734 0.023381520909931133301 2.760965908365323962 0.73885271434358457121 0.28169893779464744155 -0.28281858065920673084 -0.62474704564938088325 -0.62470226958879926382 -0.42820945812049709955 1.9793531449524164323 -0.38292027783989340817 -0.21665041294616466105 0.24332637549373536956 0.6668986690507087145 -2.0771643481241355289 -0.052096305813492359282 1.4765036476289967737 -0.88520993097833089447 1.0457207700989565868;0.56628837863703029853 -0.88371180301893026421 -0.11687903209986437447 -0.2380224868505426139 0.66679664699200513756 -0.1018316929801508508 0.27636754987049488852 -0.20490206810345368282 -0.00680572968539491438 -0.027876648801187377258 0.78668195166230125981 0.57155816396493863696 -0.085355243556031940022 0.25703910235696764142 0.1101153389882304906 0.52116075815370865865 -0.45929910279104779747 -1.0477957575974770066 -0.12301733664288419523 0.73942140075029161306;0.33354597952299530617 0.47097841289682762511 0.44905477600945481464 -0.61773609090069758754 0.80040725334708295868 0.0099195847888628105316 -1.1620191195433062781 -0.61858566046068452771 0.90391128698560696542 -0.58938594103973729332 -0.79104925715872287384 0.27345470621863665617 -0.65279418301425862214 0.0023508043596093663313 -0.49563088572290620037 -0.63319243393491719818 0.29852189696490361026 -1.4266826352121646515 -0.83909923171905353989 0.58081688060328795498;1.4054139783382180173 0.68498994059108930799 2.5409159008425254989 -2.3121486104450044152 0.44706428690623678524 0.3320587314170240778 1.2913553085139259924 -0.8185719337341899049 -0.87187845538398800471 1.1843351890334947729 0.98755996772992460375 3.0896082675114175942 -1.1916670987623307898 0.19984523234323098295 0.60847279437466850194 -2.3464205189775295679 0.58884981970645311478 0.58985217328338312992 1.2569844560130556399 -0.26769000875523291105];

% Layer 3
b3 = 1.1334391988928791406;
LW3_2 = [0.90772416628854768472 1.1230303866769277832 2.1262822446239288965 -3.252652808569977072 0.94780763290717684466 -2.0704942929971279497 1.0196253398371093457 -1.0419549586607308456 1.9793366283110385151 -2.1242143549407210834 -0.84792232382231103038 1.1658638732105770508 -2.0852413383061088048 -2.7930607978487858034 0.94242731253404066205 3.5894862565147320765 -0.76007906415010295653 1.3189295139188499384 -1.9048990177822777703 2.1607685673121439507];

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
