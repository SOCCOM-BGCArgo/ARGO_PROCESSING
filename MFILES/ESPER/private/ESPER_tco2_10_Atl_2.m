function [Y,Xf,Af] = ESPER_tco2_10_Atl_2(X,~,~)
%ESPER_TCO2_10_ATL_2 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 12-Jan-2021 17:03:57.
% 
% [Y] = ESPER_tco2_10_Atl_2(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999999390765;-40;0;9.1619;-0.2178;-0.12];
x1_step1.gain = [1.0000000017546;1.00004398834998;0.0153846153846154;0.000309885342423303;0.0656325902039552;0.0461352500991908;0.0145655815308426];
x1_step1.ymin = -1;

% Layer 1
b1 = [-7.106892583887789;-2.2839026927111852;2.1616552661292503;-3.7407515795726169;3.0464052954340537;-0.4265849049245532;-0.085675533816980717;1.21306783343751;-3.3156472326361413;1.405291452683022;0.34841727650623971;0.72296714062049905;0.2530498594114563;4.1141790294999669;1.964904465486019;-2.2215552357792325;2.1632281880531736;0.39325687134132015;5.686029543899088;0.060645328719802993];
IW1_1 = [2.2547899833490144 1.6838871695528952 3.3656751061706975 -0.15055001837957996 -0.15539407121572427 -0.1026555313036199 -0.29932979322987041;0.20070025537580785 -0.35096806761860228 0.36393247405298351 -0.69912414232261122 1.1090826584126181 -0.2656201607166343 -0.52543667977423769;0.64786443455972165 -1.324805887811759 -0.57325750326486347 0.53821356021237488 -2.1547448977296169 0.38320151795950613 -0.34897202005800243;-0.16983257490986581 -0.4022555877231479 0.75336661685860562 -0.10071596679912782 2.7061066173568529 1.2060744079265671 -1.3103786167780751;0.8627988818365443 -1.5213476083307269 3.2028057686761824 3.5822273817320736 2.1192325102555438 5.0374096824563539 -4.8205727948365853;-0.13065873808097478 -0.38816995096508594 -1.5782407184109963 0.29433766571666536 1.2712509767928097 1.134722988996812 -0.80825867831236742;-2.0527711269563 1.2528283958280249 -0.63979854697892102 0.096938222864280069 0.48466742676366464 -1.384880152026231 -0.79835258513831742;-0.32001543709722963 -0.26685620337161631 1.0516482774146318 -0.1969601809678829 -4.1562251203302392 0.79192756126312269 -1.2919118176079323;0.107974625195184 2.1386612447057103 1.1290393230750928 -0.069628687412637044 2.3469621666867666 0.8495056511501814 -0.17187462840374063;-1.7151422247975896 -1.4953966205982916 -0.36333720481137693 -0.22706911794615836 2.9138477389538497 -0.25709474169211455 0.76902290234631954;0.41376027228045631 -0.063534352925932416 -0.069897649690902272 1.0484795531530016 -0.11503318587971428 0.90539814635814053 -0.206995291545506;0.06552391437473909 -0.32415213115910596 -0.98756567961917197 0.94000820497919046 1.1809473006066336 1.3902678875951053 -0.26965469367266637;-1.6599589958982186 -1.5091938943792857 1.1816437184391411 0.63223453675410468 2.1243221953494045 0.96016621797359147 -0.42356349976789076;1.2857239232062698 -3.1055740176680331 1.9111677369714022 -0.63531022831813577 -2.8322778455710798 0.12591764555442161 2.1196440279742514;-1.6461612303925703 -0.28298616484981431 1.1765110315138345 1.1394518515507228 -1.0384656329641724 3.0129312503750114 -2.7077724086389052;-0.0064372635358857918 0.014970553195755271 0.41700487910397649 -1.2077364919714446 -0.017625682309118037 -1.2446374964904876 -0.23173005316882925;1.461583313500288 -0.76355036770847629 0.52269317557888451 0.18215514546148828 1.028080678495648 1.0642739401326564 3.6952243214403131;0.81346988511818974 -1.2421545766985589 -0.24917757138023808 0.12243148078378549 0.87473253546755836 -0.29147275769608871 0.10249375402006843;-0.7473239969364347 0.22467240848614231 0.30321620405311495 0.067702775184982938 -5.5394948296830018 -0.16306898799284536 -1.3848366897070679;0.34880056580983565 0.1147510039243365 0.094952249874811337 0.27798524382819756 1.6562640695767339 0.32211795954530958 0.079665323155660731];

% Layer 2
b2 = [1.4651422468844439;-6.4146733307502775;-8.7133686244282451;-0.35812139446898966;0.74003281710007174;1.4270139927238361;-2.8464152876804203;10.260203319111474;18.209807819230313;1.1051522793314921;-1.1597252281654542;-3.2958338806700991;19.082658948928295;10.65051977643745;1.7545881158476122;0.92342700832370506;0.0024829241540882535;2.3438868058830002;12.758803663381086;1.7597680257091612];
LW2_1 = [-1.8212502697976489 1.7214550783949816 -9.320834430160577 6.4523766243281129 0.092363056589830253 6.7445530664312843 0.23139514466321814 2.8587849598213828 0.77230575189899209 -2.1389530958139877 7.6858580222396462 1.5981195742263246 2.702578429788256 0.91806028716790411 1.541669395226563 -1.4591747214089288 2.7605537827140023 -0.047891440477007642 0.22868953404632866 -1.4580181651749264;-8.2257114570481367 -1.46224620663625 2.018280987240872 -3.7605450409407255 4.6318131304776919 6.6175581813775013 2.2125525384676004 11.039103268552404 4.2928576916525145 6.1059929382026414 3.8213820416143278 0.10431248918792677 3.6206906278755064 -0.66038440067020987 -7.0949350065058949 -2.0307801959169316 3.2257303774731838 6.1853405906357262 -6.4348031170615796 -0.9461289572676338;-3.1009208454610153 -2.0136847361380892 2.6155673837787305 0.74685104081810028 2.9012228716544639 -1.8243094305268535 -1.7433825641933034 2.6270807250598938 2.8843299862599641 0.81318552626024343 -4.3632245923292192 2.0472903102233504 -0.53248489480443628 4.3197644155995665 -1.3693970856510496 0.88603661571944214 -2.5315963296324679 -1.258921729468687 -1.410147948771691 1.7008780251584106;0.84395157449364078 -0.12314745472312727 -14.355832831342699 5.6274683299511175 4.0476167375898955 1.1168539362720726 -1.6220602633484231 -4.0663502497993917 0.64184169209822239 -6.2728565664315603 1.0580202827321703 1.7862386185647741 6.0977845137686488 6.0179069562482548 -2.9979230639929999 5.0639417224848637 4.3263207094314238 -3.3500951256143043 4.7912147019832503 0.78579307439656754;-0.61816795962545323 -0.5349672678242311 5.5823747046101042 -1.0092856768545864 0.56622226115737051 0.61492996046095116 -0.83049380389289451 3.9519505780648836 1.394290014283716 -5.2030027608479674 -0.20199766419062065 -2.6101849743680949 1.2976202125312715 -5.8217034144180655 -1.2994577719693343 -2.202690102333869 -1.5789150676806454 4.6898116770837612 4.769302316704267 -1.556407153608693;0.3909602208081564 5.6072739237085623 0.24668524985300602 -1.6174253078374483 -0.012346242265205871 1.0694579890202647 0.49704094353120076 -0.53433569033928141 0.71109215508730261 -0.069529367231461969 2.5794712340610992 -2.1451702235646182 -0.64474057676472019 -0.834363468214676 -0.14438504877376768 -2.51704338696906 1.0607202029647289 -1.0889613467346386 2.3523135724031241 -4.5746975376685839;5.0142841043979951 1.748680119835567 4.1679356877553868 -1.0468938547981583 6.4356006929654876 -2.1330117757833578 3.8490994886367362 -2.6941381933231243 -0.20647166741402362 2.7811929205979449 0.099524610425424154 -0.16642949868255741 -1.2847256813188503 1.290377864244252 2.9085356282820687 0.59047836100420859 4.2589894927403433 -3.065995588980774 -2.9341860238797053 1.3567359726856558;8.6264015469507456 20.684775046541205 -3.5483276287842136 -9.4209973928667967 -0.4570182505352528 6.5623037985198831 0.40432393408516193 1.4458053856304001 -4.4620330198464462 -0.81014712930178201 7.1614768850387343 -10.653336149925121 -0.95075482383134746 -6.0798273195448465 1.1629544987318992 -12.441564376033025 3.9378891822938864 -4.9734453703144608 4.9098762321774068 -5.2483731612438813;22.014045925883234 0.85705197786517762 0.012038289824544354 2.4924608565418449 5.6549794519570549 3.0007670705133473 1.4858309337347579 -2.3062976597062455 -0.30083525620290064 5.4665290062359322 4.392484303911278 4.0407961655819786 -1.2474075706110901 0.51121877364328194 2.499107526624829 8.5653302445878321 -2.9630061012709952 -2.0899761431534505 -1.811908321358997 1.8675999116752047;0.061241001122773732 -1.5210322472006963 1.8650425647186721 0.014338811000345092 0.41358363723438774 -5.3658782463350052 -1.5665937012488731 -0.48776734800354543 -2.2322929479125482 -1.4086851056579794 -2.2383926116411459 -1.641423878037001 0.44231384586429368 0.60901837488831623 0.7582491932861094 -2.8740620958884513 -1.1497078897176165 -0.62795189270768492 -1.032935484481242 -0.92351028158239634;0.16844394334427443 -1.1002388041971789 -1.5800561329667615 -1.0880273901543827 3.8435135111216483 0.74081694303414269 4.4794973470394934 -0.16892685517548212 -0.16730582506017672 6.7630729129132607 0.61685858184408959 -3.0272790976964261 1.0042930417960587 -2.5949725228932361 -1.0468422733201448 0.84271170928731942 4.2808299352917079 0.18495185316172513 -6.1600196789658126 0.43911208706736965;0.42189008787983928 4.6407139379094504 1.2849885351541872 -0.59365135299356153 2.1209898576867854 -1.6011634679839561 0.57436509389684887 0.74068387087767895 -1.3142107783621508 -4.7772100017037573 -0.55116171794711333 2.6771357380521108 0.64273965161787472 3.1380804867167775 0.54169351490705364 0.20448446666844769 -1.4409982648590989 -5.6438485486892418 2.0584286549311246 2.786893630930082;27.961109277267894 -1.7470619471332061 -4.9874600288171251 1.943572401015238 -2.4280491856089594 -1.3647316912099177 -2.6273617407969376 -1.1290010230488889 -5.4442905810538242 -3.7651358168643103 -4.2299780228285 5.1187745231777955 1.7381952033961485 6.607630826583768 2.0658923981431867 -0.37649728716640057 0.9841214174578724 -1.9525420272564726 -3.4435334997633098 3.5184773168900807;7.1841477078714009 -2.018191041706586 -1.567280122013095 8.020719090727896 0.60184969495758855 -2.7622848951134724 -1.5220768859294247 -4.853133403213886 -6.5784237826758796 -4.0933186345169448 -1.7972214129823605 -0.46724555123969352 -4.7437826560349103 -2.5239249879633858 -0.60231684559138565 -0.90860365921910224 -0.135056747906387 -6.7148210924029224 0.30300135334451683 1.0721811532073626;4.7826907208574321 -3.1610190792333781 3.9244610922474865 7.4529908244592766 8.0464305293313672 -6.352581800928319 -3.2207359176353458 -4.7213357559563747 1.4179399705510938 2.5256975930726444 -3.0630485744409097 1.047076677251755 0.22912857375827866 -1.0759366939677919 0.56877109791215219 -0.9809881359821967 -2.0502597824230131 -2.0451854053622394 -1.52974505894794 -1.5462341059391034;1.8993171793234185 -0.6296998609030805 -7.6103437235335534 2.19945588064502 -2.0092127617857689 3.5658064562643879 0.030022470893662386 -4.7450847926320074 -1.8262084880095601 2.3735331724966695 -3.3574595645507381 1.8345421510457176 -2.0458466952499741 -6.4961820693633703 1.7690673880898256 -3.484985454196655 0.74216164628346848 -2.5148771572619126 -0.74608852344590126 3.9581133525501961;-2.2910561590763203 -0.61888715125803473 1.473178294759411 -5.3424947672456158 -0.71712887049287288 1.014412999405109 2.3381377542468433 4.6547239248001402 0.51027148959404345 2.0134937855945343 0.22013951726894454 -1.164642338902014 1.7680891531371077 -1.4709962435414952 -3.1601573056938985 -0.44218506108405359 0.34316547492375588 -1.0390518772061943 -1.5648509822111727 -0.31633076607792499;5.3877322639944145 5.3503765509101893 2.5094769392183012 -2.1400530662675057 -0.82916340180892023 -3.2302470824318625 -1.5191076378708797 -1.6999142360263451 2.0032780894436781 3.8534116530397537 0.61018681741731562 1.5824145215549255 2.9943407637093449 -0.49080627228279405 2.6567345118436396 -0.95029961912455263 -5.8154093272048453 -1.0712727529326169 -1.5580292734736421 -0.87343450949902035;0.064965153413632529 11.066146048474756 3.4650677343358227 -3.6010748432681834 -0.4238461578741316 -6.4233598906445257 -1.0086237649197412 -2.5635935474089488 2.2083523262333866 2.8727182778745957 1.3745547978844013 -2.5697946654130597 0.94925469655031358 1.420546533383038 1.1425429128631535 0.41199511729683819 -1.496068174439698 0.88831203507886336 1.9203115009938052 2.2345940474571959;-1.3648969375548525 -3.6457718647740065 3.1404202830325993 1.4042142844920038 0.19244288043344482 0.71185965049208377 -2.3947491394270934 -0.27246604432064253 0.52871833505085175 -3.2337736141445577 -1.0554511931803241 -1.2331961718004376 1.14640457813838 2.5829339199746277 -0.080591934296616208 -1.1813586037192703 -1.1178730416338887 -2.9214923115489198 -2.8255583764883214 -1.824803626837632];

% Layer 3
b3 = -0.39067319592212901;
LW3_2 = [-0.0024166202330624846 0.011729547988607371 0.028787792040284036 0.014466320312981163 -0.034876273267306947 -1.7249542646333957 0.070801712724510121 0.01911392787634339 0.012159008700887987 -0.05699988368444877 0.027577328023615108 -0.073868579059833669 0.023452588057469116 0.087787672888078885 0.061220146088110446 0.010244308202019733 -0.014624070373197023 -0.0062632201793324785 -0.51550941837671183 -0.060358251555823146];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.00115315130997576;
y1_step1.xoffset = 678.286531932822;

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
