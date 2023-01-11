function [Y,Xf,Af] = ESPER_TA_2_Other_4(X,~,~)
%ESPER_TA_2_OTHER_4 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:00.
% 
% [Y] = ESPER_TA_2_Other_4(X,~,~) takes these arguments:
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
x1_step1.gain = [1.00000001035699;1.00022983623563;0.0140762793578401;0.000302617642608564;0.251952632905014;0.0584431850247228;0.0422119037568594;0.00826143816767911];
x1_step1.ymin = -1;

% Layer 1
b1 = [2.1476927923382604924;1.319125844752564447;-2.4967799123871015077;-5.6629209070918768987;1.7080246827353118277;0.18034423941328339702;0.9102581383714289931;-2.4480214183213764301;0.50525977813979472408;-0.75819854089708993428;0.31502304111199042058;1.5370892048153006471;-0.11438844935853942353;-0.72423480399829343046;-1.2108924494571762231;-0.55092855571471888165;-2.9103575440355733406;-0.58393382404659655016;0.57438910665639553166;-1.0787172930014365946;-0.24988173510583860937;1.0073449116105099854;-0.37302123788974977936;-0.77403350971367956834;-2.5035414583919122222;-1.2549791884445453682;-0.26753701436829280258;-0.81122872096676146292;-1.7724210117779624074;2.2311906374462004976];
IW1_1 = [-4.2800255883538245172 -1.5000412462372099132 -3.3563230257549596658 1.1933200819736908826 -4.0717846192820266182 4.8991363863687746871 1.6441888216553395274 0.82892573530752866606;-1.0093889941801907906 -0.016318261530049830987 2.3542226085865349816 -1.083201812286310739 -2.1681317276252887183 -2.0869264728763563532 -1.6551820792556566353 -1.0570558651321873;0.55965717238235090925 0.29363133375498828848 0.052342979965341485116 -0.58446571810797365298 1.110852409536200458 0.90603974998962455434 0.17281269085528286333 -1.5892403732201765276;2.382875060702535297 4.7484623899458169305 2.2626056864795600632 -0.25839592966355323655 -2.1637047113648999463 0.050470496086672443636 0.065403145759592037911 1.4455656140487971317;-0.55091067059483622703 -0.14134010238896838052 -0.11823650223729582909 0.50709931303539990122 0.49768025905964752287 0.7862646690353961576 0.45445404983203463445 1.0504650383778475486;2.2555075500823500434 -1.1094076142475237656 0.514155431015170028 -0.14064300135037219319 0.56068549725076022305 0.056109228023860369117 -0.89082459467310026735 1.7212339494995536615;2.6260321234867047835 -0.95846029970406687593 1.4702179800288455169 -0.30414962080295071534 -0.34533987296387969046 1.1202317896360414196 -0.59638963023615620962 1.0488716723186606306;1.0886258759841764832 0.24760516142260879935 0.78348104039257049358 -1.1144108009238333779 -0.62269521982144315775 0.8429287664080575837 -0.078897071871582347136 -1.7638055976123812041;2.7345886142199717916 -1.0112160288203424496 -1.6598949633083255328 0.19859614111387979829 -2.2916256693431367175 0.93680516942595637442 -1.4080418172970374791 1.7138720482534048895;0.46526578051283945525 0.61688140175528460851 -0.010871487819010674375 0.85340451711053810779 0.18674050957771085035 -1.7073421321826989239 -0.23765796342646203221 -0.82764330935068741724;0.10967999208134789069 2.6783089095532748303 1.0550516057528207536 0.6347754543272783323 -1.2260052527133418199 0.5617452303163791294 0.045533898026221437949 -0.48753256461035232583;0.52157492535721239335 1.4649822009376662901 0.21744116013735428905 0.97147031650856929819 -3.484307287863689151 0.48730346448847167151 0.18287453167965700418 -0.056858076711564897732;2.077767831304973889 -0.27980006864709000558 0.40393802660530364612 -0.52435558266794801696 -2.2395029447856500404 -1.6370764939346702516 -0.47589098464663942556 -1.4483585987660041461;-0.017793946118065812056 -0.0089126592993088761541 0.58215107131526222517 -0.36195715401308842241 0.58906461395741205322 -1.3250899267474665066 -1.1187897298339966046 -0.13654977203506560302;0.61968690456355590701 -0.39155752941340776774 -1.3302632749438458859 0.81088168300761309659 -0.65863362592766694359 -0.80722095217866751682 -0.3217491940629824132 -0.7207653109030294214;-0.23712763900612904289 0.41434545999137606387 -0.34440547735866272605 0.18992601347099069553 0.93771135852478060269 -0.48480508051710341499 0.10059255594853938609 0.090417527355619101148;0.24601873482107283753 -0.32402828206908151909 -1.6204358401911587162 -0.61759730150928737746 0.11336251973934877157 -4.8196019075469997972 -4.1437737884062659788 -0.9444919080937993483;0.40900156706278367746 0.11528885216773238076 1.3838975889014384357 0.29848615662414834659 -1.1801268793028694137 -1.1818382593291900662 -0.20572122301198220162 -0.12640709965790286606;0.32856741470494227331 -0.29335949107468650698 0.9647103438188588953 -0.24267439979547744677 -0.4264295071006386495 1.247130590378624948 0.2813475498857572199 0.28759443204321860632;1.3666578983193029551 0.95017895437778410184 -1.8056518381473321 -0.58896714380278403222 0.010290844933638729208 -2.2369055625066893711 0.21857557222070722602 0.079376253330649140794;-0.088782169357478540883 -0.070369728136675108177 0.37982476310532370389 -0.053205205939271735249 1.7897528352287750408 0.78050942173553616854 -0.11101414169481384209 -0.9571507751584307222;0.070225389018844439071 -0.054886929270326942354 -0.10455411706708321595 -0.038100534275124643691 1.8158747038570532961 0.12987424535847974516 -0.10337023510102495705 -0.24239878193926675531;-0.15158133940268181394 0.087861299735198558625 -0.81644421413745682514 -0.064390280348374615005 -0.43881194374097370892 -0.58341332232196141394 -0.24087775807233249847 -0.12811152727154276332;0.057333805959679955133 0.4716562039698833364 -0.65999668665300670334 -0.86848057981477788658 -1.0310465170510276778 0.19953886058351394395 0.67056190300837026363 1.0402234587287964018;0.43256023392107112446 0.060908993959588861311 3.0207548460349515906 0.54755518932969093893 1.7084977885287297994 0.44682798485537866817 1.8632760224820921202 -0.77831213542271193706;-0.34602599543157502859 -1.2798340424940386484 -0.68979055066529715301 0.3991854071812863114 0.049619185344636440593 -0.023554752985954569705 -0.85031400144800439733 -0.40686061876275292359;-0.19604008654469673445 1.9714146503006892974 -1.0233409719066253185 -0.42385734279427866866 -3.1536208843064996721 -0.22365359256408476329 -2.2874842214486714731 1.8999887963618700137;-0.35888785932921390964 0.13367469288110792047 -1.7693733221265595112 0.31232594006196678915 0.41532067649004505983 1.0548059008654073754 0.50003433988352741224 -0.36894557427213842127;0.16095156569944193969 0.11504064525159507504 1.1702942571125485305 -0.28121601938731072279 -1.0669755030719663136 0.60103725712943811832 0.084607404632218788376 -1.5687330285994458556;1.9898475600094109872 2.7769042498129028473 1.1672050963845277405 0.13989914339381764874 1.5054653517428788767 0.88561040541529312264 0.017963204552627554805 -0.46726366713088152149];

% Layer 2
b2 = [0.74456032013819128945;-3.405737168491147937;-3.5784341542006550263;-1.9029637755446426883;-1.9716596961474701644;4.6667415999178594177;-2.9279233889016169634;3.4693648445695779436;-1.1711361925171948073;-2.6024568325506036715];
LW2_1 = [-0.092622961499836148347 0.033511854091698808911 -1.0252325897573228985 -2.0427690672458815868 -1.9218790900703615065 -1.486654730455726936 1.2258401050026030976 -2.920979135869438359 -0.28116972903948395102 -0.92040029504325482357 -0.75559565822281338221 0.33093935145578445844 0.15173437185003335603 0.53176894663587037648 -0.35113506155425355804 -2.7559683715660465708 -0.30179731156058875463 -2.8820039416768894469 -3.5407497481461942535 0.14897676488642300208 -0.12199252946533697961 -6.4928792407045028412 -4.9449573786737293091 0.8939626017006264691 0.30001497685585987174 0.30005636356904197282 0.53970625757779677745 1.7662782632486841994 2.072427620551815064 0.57966979937570928261;0.65771614115014975788 0.86831444404663893355 2.4863015974447777445 -2.4821233306709724609 0.69252244297898235548 -0.26534224572780928941 0.38896874644104195706 -1.082009489805118907 0.3618350642752056201 0.56685648843291769339 -0.97027163930714999118 -0.15837682210927139792 0.32126739519267727418 0.40434070876259375904 1.6866209277785559895 0.18374391795822669904 -0.22314487840681337949 -2.2014273066605762885 -1.6867392044199103207 -1.1601792681345841629 -0.18245233650083450549 -0.76494746842833183376 -2.946926435149095802 0.15615382755463638742 0.86297884226405752184 -0.9839340663161701972 -0.23554297609335775321 0.99080618718204105377 0.028776390062217500898 -0.97481656589627707632;-1.9874936143673527233 -0.64824463991662317763 0.69005425324277180898 0.72804536814324782856 -0.63846224818924679489 0.78472186968519375139 0.013089763995812554731 0.24066939439310328086 -0.87180721838764252407 2.0265728517596932612 -0.34685622727267212406 -0.22678428460856497884 -0.67575795747559119775 0.43413131480158761999 0.65660213848420567739 -3.1278835037200565239 0.43487497275927716744 0.071543894446726305492 0.77964905352066549149 -1.2115365067326950843 0.082320381673500364617 3.7615561015674954248 -1.4611432486437792022 1.4730285530578837161 -2.0884457316211020306 -0.19117318070580102685 -0.52262982434645555152 0.62735517479789260076 -1.7417057129371462665 -2.5199437922340504059;-0.033766142117540628997 0.10190493887936505346 -1.7500228709641501013 0.41480201213049927578 -1.4000759650097611697 -0.16244886313486683882 0.10280882346314688947 0.76244588413214731126 0.20988575230221478973 -0.81563928518080308638 0.16538550150471598155 -0.3118258502862968351 0.22454598830009228627 0.24649906647295416473 -0.77107641657265879598 -0.39133093725334527901 -0.22094577970895579178 -0.23550968907354599691 -2.9691009823459735273 0.26346518889048475831 -1.4549584185676005532 1.9355862135206070018 -4.4557863163756037395 -0.50114733476733164252 0.1264196556891583334 0.21155552965026691581 0.21840650686127618951 1.5898702416839727292 0.3728678780493888123 0.1505431379898940214;0.60582795993592275519 -0.58489581457172146184 0.29158698362101964641 1.317742848290763602 0.38447945046940151803 -0.55399288535534818578 -0.10424708546372753182 0.74021540951134423558 -1.5787954909558634675 -1.7967765107430522242 0.4752864537024893421 0.5117462753321434521 1.0655483219551915575 1.5609002787035464177 -0.044824917806936724618 1.383328007221282574 0.13746078277934631329 1.2610805447613802066 -0.51955845327169902781 0.68568868412249184274 -0.43125720639662629141 3.3800289264732970018 -0.87236740619321895274 -0.10171584288563463583 -0.14310916337128773423 -0.26647417532306483245 0.55484242887803836286 -0.066391276637637053337 -1.8222437802635176318 1.9740088569412017616;-0.018479791313761279731 0.72829114589602328422 -0.28983891784532722635 0.0063510455480102409598 -0.14449645046494002942 -0.1470302477822213183 0.35317384732225276522 -0.5613935235733843454 -0.27152713057228866633 0.17152819625755943989 0.0035109495328092604413 -0.12553435273951485396 0.18861768329663194943 0.090965367593161233772 0.41829027999814233363 -0.78783680216435825194 0.016161108612947752 0.5500149861675996954 -0.42353890158672236055 -0.12533716092275162812 1.1249112968292647174 -5.2063162058641774266 -1.0297522125197431375 0.92394444049196278179 -0.54831675238963550889 -0.31338627894878462454 0.096584104884353469411 1.1321306466051941353 0.45296591006060898943 -0.24593436484172656598;-0.11985328598674485634 0.48125507973990810928 -0.20363425442226140705 0.6100272337297697467 2.6303401517331717052 0.14494163728089920484 0.060568315779992117076 -0.36947916744255976385 -0.40250179417806658533 0.70413788313724667844 -0.44938828213482473473 0.73827024037688726565 -0.097086709377932248199 0.94577936129357809136 1.3889461738161483773 -1.9498320148008820851 0.29216193717687272136 -0.42724864028178460851 1.1421597156797438366 -0.077434954769137115393 0.58609573883439047837 0.85154867049768889764 -0.33005249223774918432 0.9450245568877903013 -0.26859234252696079226 -0.8607018213502588555 0.12008329696337897985 0.23477003968547072166 -0.99259527405327785399 0.098986304900988211775;0.059642372854670548665 0.77281225737373493434 2.0729752970873782658 -1.0517503963180983018 0.52514859288654269776 -0.28837858628123225202 0.61283116328509013027 -1.9729811503269680806 -0.6556696519797194922 0.82672085596575539679 -0.10915028221110133211 -0.032541307936852385818 0.050808410658661283532 0.29621057663426342277 1.1900347009588017055 0.59128196558126222548 0.14999127049115659749 0.55521792866786034981 2.1179568818699641142 -0.51886555909686071786 2.0935808644791746858 -4.512806604524359777 1.7671342351752534761 1.7196463809634427999 -0.62113505915740252483 -0.048990252539078085747 -0.19284358814750796052 -0.11486681309527957406 -0.33232130919606761577 -0.52016307479903134414;0.021027420430010507896 0.49309202233355109435 -1.699718562877791328 -0.46516886366429616162 -0.17123254596251499815 -2.0186883065348166788 1.3300535059092695178 0.52458123617600793676 1.1245734821823178962 -0.66649492450424618895 1.4858796161789882273 -2.1942896754965675044 1.0036780975453176712 -0.029376237076671724013 -0.55782754481220542253 2.5616438384558817276 -0.1534517166494328011 1.2680065212704454414 -2.0121860476468373946 0.34366475800003271068 -0.87792221058804453815 -0.60236817564216360577 -2.8608050985825359547 -0.78701927003326588306 0.06186160757950866923 0.89391506810314735887 0.37535046867205168475 0.60977973128187512586 1.138354050059839917 -0.58585944958133961968;-0.56228736181619010104 -1.5395766678912510272 1.0216649473181016461 -1.1369261788392932555 -2.685141922118341018 0.044431545369470880347 0.76671078817456983501 -2.7233805625789906379 -1.0826319524026983832 -5.8512666952652558905 -0.50292076008999997416 -1.0118816800977323833 0.34924717914618058634 0.80112996342935693939 1.4992668321508617613 2.2303994243252431851 -0.055401084944049105996 -0.50851961532360046014 -0.10998161761720454122 -0.67527614395547475024 -0.24829739791532437487 -0.71380296131588460984 -2.6345462467596556522 -0.25660775950048303473 0.12495691757944621392 -0.080561085080358418864 0.64961637879512490468 1.8725033816064828862 2.6459064158096463615 -0.92418329966319956448];

% Layer 3
b3 = 0.66691042400884525954;
LW3_2 = [-0.078060258362952300581 0.098541772582423858839 0.054681641551844580029 0.64603304185548227068 0.070512550170984386599 -0.45897680961372078201 0.15446067167256993802 0.2614101453931220731 -0.2545526537037812731 0.26789820937094427356];

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
