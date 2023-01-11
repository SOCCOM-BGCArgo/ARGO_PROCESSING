function [Y,Xf,Af] = ESPER_silicate_8_Atl_2(X,~,~)
%ESPER_SILICATE_8_ATL_2 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:40.
% 
% [Y] = ESPER_silicate_8_Atl_2(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999994516886;-40;0;9.1619;-1.99149636107278];
x1_step1.gain = [1.0000000017546;1.00003626070655;0.0153846153846154;0.000309885342423303;0.0656325902039552;0.0604107502911417];
x1_step1.ymin = -1;

% Layer 1
b1 = [-4.0515730314671181134;-0.22685689014983859146;-1.1042844528573061424;0.33936468426089017036;0.95295088946005335195;0.77661595966143548431;1.1338108625828042797;16.043652319762525593;-1.0497315550675865303;4.5060314746993492818;-3.1474432919244845763;-2.2003289231103284607;0.19425376891171805549;5.3324660294615080502;-1.6404522433688542016;-0.41344469209103401619;0.096475920626408223635;-0.4178571463168718747;-0.4385281602290406755;-2.0187376941278278863];
IW1_1 = [-1.2896790092246743775 -1.1502329002519431356 1.5154805541329126051 0.54378570936730197261 3.4315092681886563142 -0.89510602543064321512;-0.0061845681819300574872 -0.067385409056122308469 -0.90021686258424982174 -0.16214862842035576707 0.25113483620190857337 -0.17460123595237253546;-0.86429405360918265533 -0.23414438841326434826 -1.2373154033692477238 -1.2402755003276715318 0.33756062050630114557 -1.3817749232532468007;-0.56078381290871026632 0.28020538240191295021 -1.1416341263069009759 0.57491806231135755834 1.3980171316305307627 0.29164606972552237218;0.054171796558712202729 0.1310787163894450158 0.11392509867282030189 -0.16146801047324066403 -0.24275975453114922997 0.53854564362803336586;4.8505232602664358765 -2.4671994514796824483 0.24559158010240392245 -0.80312788606770202371 -0.6420637506390273197 1.7869951458866610849;-0.04574900806990396529 -0.042807284744612109084 -0.19770809090512400652 -1.1782309546607596662 -3.0089961691194648807 -1.3833644543935146309;-0.37867363663574721278 0.71143624661652449248 1.0988630682340505551 3.2113184801667600254 -13.52471013584495374 7.797423243711505414;0.62140254938793282324 -0.89630645644113682291 1.4325101254390923966 -0.74262709609730248594 -1.3312250377416607883 -1.062977871421658671;-0.01118583606340214226 0.28719415780016388506 1.179703192727548311 -0.69035805844976427803 1.1310758520177943698 4.57792866748985805;0.32043386666225870396 -1.0904628274335146365 0.75727681002530655707 -1.4682647252903935353 0.89756949432405075395 -1.7674684048626807531;-1.8843812848194900766 -1.5116557430274866203 2.6232551150558771624 1.6884758803554844597 -1.0323547551032363767 -1.9940842043197590794;-0.30337161430293835274 0.11224753196104456099 -0.13494124056235926035 -0.3054806583511839313 0.29081268256190873078 0.046487903883151625306;-0.67463282610961516905 0.20706224780915211126 -0.49012145089296571543 0.61391349414887641522 -6.6038623206684334832 0.99218080703084199357;0.55697091289322864238 -0.59071002500212232977 -0.51108831509212260258 0.67677487895132515394 -0.42017233156428374086 -2.6267724588829883814;0.014652295554442645517 -0.063687075042779206502 -0.69895543651481784053 -0.35969942051141745587 0.18929936241534264796 -0.11431733874182148369;-0.059317285043684442325 0.14803923159657894626 0.73950356138701145436 0.27584909671300700262 0.25227584854602302933 1.6513208633620048982;0.11396322017026611073 -0.4427956573554339692 -2.1552296440475497974 -6.6221106982250885764 -10.97113015630130306 -0.94481654067031928257;0.22225135103181328788 -0.27594181191975997169 -0.25810823265438542773 0.17884499198638253992 0.45986251008392470263 -0.4608843552532362331;0.39636812149774380831 0.34360997226872597166 1.4655483865332419313 1.6778643076068882056 7.4257907915957419576 1.6021591852833052538];

% Layer 2
b2 = [16.742329004852894769;-30.232828799882636162;-0.26836035691887211163;73.243946330020733626;23.555861958734521266;-19.84292477291456791;-3.7645626708320656384;4.6689688804389701815;-15.020448433127068455;3.7767488296809603732;-29.591387355264792802;-69.486556788524737271;-26.906315836395286567;-29.103601484499428409;-5.5241700404536580038;-2.5668540168259297296;19.769191116338742376;-14.304295466184839825;0.35671741383005672876;-11.762633184351988902];
LW2_1 = [-8.5514433754896117534 42.657687477856200076 -40.131177568899843777 -32.807524087229843701 -42.042146692669248864 30.077022158181542011 3.9134074592037060292 -1.093426552624292114 -7.6745642496328505544 -66.216734678042769247 -21.868217567705695359 13.214739343543689287 29.22593758142909337 -39.256776519402841075 -5.1195708690827848031 -29.775917234946792433 -12.541567000620378991 -7.0323134632879247263 27.641105925743278959 -34.058291485418081379;16.467913509712598596 2.6329739051522396132 -4.7351703584990945117 16.37873503728530622 16.96833709315326999 -4.5206481116721954905 15.051798128874612104 -7.8022606136897643125 -13.944266618332751051 47.499757734976029155 5.9016837449514367719 8.6183568610324154946 8.4856321408279846708 -8.2555345796945704961 -0.96378058813695000051 1.3049549856044029283 -6.0436882658096724441 12.098703925956852956 7.3380938996016356057 4.2993442661848364494;0.11353816219283847933 -0.17678103832540967288 1.8500258126257129732 -4.6038318822996506796 0.64853281088037162405 -1.0804678507362990292 -3.0870202289020420849 0.22965574907299654561 -2.7642253513485850469 2.8168207825131656996 0.67602901894942835082 -0.76289262039917260516 -1.1955860609462720667 1.5463855846441665332 0.52990708882988679651 1.6254416488456966405 -0.49929028005155484671 0.29508735995854956657 0.38704152752280290839 0.61601273327558014881;5.6818821533090098086 -14.354631956747509136 7.4680948812060057307 -8.7780176144121533355 35.126271490111840023 2.6175861659362085199 -8.0845254615236594731 -1.9599648704238943342 -9.9465130167270583428 -37.365243151931295529 2.7528552299167849959 -0.36904303536619281756 -47.225481003777112221 36.402483863168249911 10.745952009003101324 0.40724954527953849404 9.5487976217575223359 6.1174660798898692349 -9.029934513688125719 11.64606727547807985;29.826205365870993802 -24.336923150009994288 -251.47196806422914506 62.007979824529016355 -129.78469603653965692 53.527486178856378274 103.27913707505241803 -43.540066990750382558 -99.086839417473640879 -13.202912831387388692 3.3614266847381273351 28.59131805994797304 66.667509434992084039 -68.827393259743161025 -76.606349060207264756 -26.294528142446722541 -93.320783981268050411 -4.1660676246557875047 -78.12247244428819215 -42.677659661238024569;-0.12800232093452429405 -3.9601686800087341211 -4.2519067535881163167 6.2433685476239348944 19.033244861611969156 -0.80253397515373348092 5.9357282815088208139 3.2142517084716795139 -14.857747440490959079 8.5227828089825550251 12.783431087527246106 -2.6217254304563550882 0.071977227797552500999 2.1900530656370418292 4.0602275224897841355 -5.5357095617077529681 -5.6697762520287957599 6.5856528911516738489 -3.0773532393436813237 4.628309732269928034;-32.256533701494817024 9.8451204308489224104 8.2227570014151893929 -5.720257313305438629 45.856722097563924478 -13.648867883853116112 -11.88426300584144002 3.135387585145694711 23.552795155057172138 -46.258725394878261739 -26.170135626494950287 27.575368602564886089 75.066455089766975561 17.882971996895413724 -28.072604499145718648 -50.064732247996708736 -22.684878957797682375 -10.27775955923347162 51.301960346070224261 -24.328077440225243322;-31.632225743179141375 -7.6168840554189918279 49.735746797987545165 -61.992222959008870475 -0.77433609028116512363 26.473445784492298571 -16.189537120238153278 3.2920618236114047583 -26.856912443565690296 10.548847725573509138 -4.2813875462475436962 55.136358163053991177 -22.444186103484842221 -5.5141557477531177867 2.6476604318079082212 -34.608606526799107428 33.436722903262754869 -0.75847579121090713539 37.122995657249056478 20.868141169653927847;-6.6940300559672056835 -54.888866170843897407 22.604852972695969271 26.166082961009390573 59.707153732055950002 8.0713622300526548514 -1.6531474147546265652 1.8884363259816119296 6.3473146934239936812 9.6970413881521917432 2.8299406001255613674 13.833396545025838265 -55.929642073518984091 1.0273267953545932585 12.312477126746587786 33.115768042145752759 10.461594645441094187 5.4041957381099559043 -0.78491196762938086895 8.7770303482263170025;-37.480342089239130132 17.727035823437894635 73.829081026728360371 -57.968142557453511188 65.459305709085725766 -26.923139511196900742 -46.740552041564122021 6.1209548472941310848 -31.898235932067937171 10.967813815753018503 -19.440174069609959417 68.490876388915609141 4.7303640033386278319 89.821118766991475013 -64.466582197366079754 -69.57912890707815734 -64.937518228368176665 11.619477916309332244 51.208852232819836559 34.551566294249127509;-8.8749868073185051998 21.677138470859890873 -2.2139775614463128584 2.8804846313569467675 28.181935130786193611 1.4667017108845818996 -5.1528742970449332716 -2.5858693752380488284 -4.9470382053757608887 14.662143995829721277 1.2111551917304181902 11.238173778733200336 -4.4906714873563071677 8.610804602616806136 -9.1808483560099247711 -45.89313485440757745 -4.9506137810893040552 1.3979530323434641303 11.394857652487214494 -1.4129174923114875551;-21.204169458913021629 9.6517414997022576983 19.772015682177535467 -39.470657411610211796 50.425530503261335014 1.7316982446005213525 -4.438271818923920975 11.633489051414382942 -33.613206813733484069 48.585736586556031114 1.1214771839216561844 15.315290540137210584 -79.221707239971109971 52.031099468304574884 9.4713294467571564894 27.94554330976819756 -7.3397717653328022891 -40.610587877298094384 47.81830568895705369 0.85863084731204819544;1.1601558549530011089 -80.736290049280086123 -6.2398470696581771477 21.302951372354819171 20.060840083397817324 0.6837258490559409152 7.0615792691851781981 -1.8094820331180372364 3.1723088452965551021 11.695088717616389573 0.632317221984450617 0.42796257995561781673 8.0993449248767550586 10.816020204807022864 3.2052607174852427718 74.425900717898628045 8.1285670451881433252 0.60010368942686143079 37.861958165276682564 8.2910843781121936757;20.155357842245265942 -52.131642514264811439 8.3333196551375294803 25.742619076358803198 42.538308426283300889 -4.2834944163743253753 1.1368488461638086573 1.7469172967264792184 30.302965058691349043 -6.3663249135007902169 -6.9563940817206697531 3.7620214968245178078 67.939003785929926948 -0.92522367005545314012 -0.8364703231090940827 36.40255575050729675 -4.4033216683530591595 -16.7290459790105146 24.127331265196087173 -15.073930398639934225;-4.3875802952048932326 13.5178282123466591 -9.2947034465096862732 -17.530765860569470505 10.303470306700045711 -11.633236272319255278 9.8843066052096411056 8.1056371883241986609 -7.190177031454738632 -9.6775018392084124486 -12.854538862399479626 19.355411784202260606 -4.2861127456138863323 -46.242217226722992507 24.020695269627243107 42.846396316162731921 -20.886291829579995749 -0.68496745256667834756 -36.61306108566955686 0.68799488532340691282;-0.55920409162404038916 7.9988452455388188156 -1.7913166383247371982 1.405880837708357145 -6.4358937927179757921 1.6654182259744858552 3.6268916447008590076 -0.16898017603463943992 0.1955736458385551324 2.2408536446912568074 0.12056182111706680027 1.5099390720207710448 1.908944255143508073 -2.4908548937570156312 -1.3456051027075610627 -9.3165958387031668764 1.3241006544361260122 0.21627591232567344726 -1.9292814390466490604 -0.14249882134598992889;18.963191214618639435 27.71092997641064315 23.432111422756381813 18.197521892206854233 34.650530943460026378 25.705208123723611635 -26.503076520038636232 -8.3997531986978941632 56.262624196802391907 -26.144093979753243673 -19.879004110446373943 24.655697774651976317 -86.774780448706437141 -22.083963350183775987 0.87737667518578399406 11.990691852387831418 7.37307743151150774 -34.995008341218898806 10.517108246722264298 -19.414913031504909213;15.343281278755972252 -53.843994211940675143 17.750629124127147662 2.0632623322191463444 -12.090874507903633273 1.2622766432647469781 -3.8132506293046737333 23.307611476466838951 -20.023373040240461052 -21.878305279712876086 -18.069001088382957931 -6.2108215486551241469 22.040729279380201433 -42.91059294238432642 20.338421367005963702 9.0286469227102852386 17.306560712043427941 -15.654196603922652997 4.5899555151415540877 -17.853537981529953527;-0.061190648374762880712 -3.2434361227776795644 -0.060071937846482209578 1.0539221647988707975 -2.7225627430145409313 -0.033552815749223412201 0.053850556324345551229 -0.029340096292079863305 0.9681375258448586818 0.54898940376477478598 -0.28089112756777090407 -0.14672000039587276832 -1.5562409308100777494 -0.13530191070942679255 -0.15470716199581782679 4.216814626187998627 0.11129645942089233523 -0.049670987999047823414 -1.2372395336086858819 -0.058644462716831331772;66.15283370911474492 66.036194505878995642 -62.933636745677141278 -14.997614500993828202 10.577362163907523041 12.88386296305460732 -15.747059481127186231 9.9279071038581339081 -0.030365427152917106451 13.566018110508890615 77.00175871114687709 -2.4847729037261028218 -26.628400803739051383 -63.761760814588576807 -13.456661148839012654 -34.950764416012880531 -79.263291949279931714 48.294531746221586843 -40.90619404203346221 -9.9381676343499840698];

% Layer 3
b3 = -3.7673711106600777931;
LW3_2 = [0.014593846721862965832 -0.01645732393416158082 -0.68225637689633256144 -0.027308939378094667511 0.0069169509625968232658 -0.035570774334558234731 -0.012195110581262037092 -0.013915624718463488058 -0.020874453710065445089 0.0041130816363646842404 -0.01838687449339057936 0.014913395048318041705 -0.061380322722419319859 0.0075128717802649888635 -0.0097903740841899858371 -0.56837186781863280327 0.0061728894024765364537 -0.0091334681130172418878 -3.4986198240155732542 -0.010647265563123214233];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.0145232735458572;
y1_step1.xoffset = -0.52;

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
