function [Y,Xf,Af] = ESPER_phosphate_16_Other_4(X,~,~)
%ESPER_PHOSPHATE_16_OTHER_4 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:32.
% 
% [Y] = ESPER_phosphate_16_Other_4(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timesteps
%   Each X{1,ts} = 5xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 1xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-0.999999984769129;-0.999998793562952;-77.961;0;29.715];
x1_step1.gain = [1.00003298043486;1.00000060321889;0.0138606862425759;0.000266418009857466;0.193068829037552];
x1_step1.ymin = -1;

% Layer 1
b1 = [5.5262943563131745606;-4.5548200564743970276;2.8209113485689796264;9.279545395111650663;-3.5118457917538643542;2.7843852579804067027;-7.5727611175915710007;-0.33467177325255215159;-2.837171128614571991;-1.4406196227773622187;-6.0583247591940931898;1.5136574021977997617;0.68417533320435730548;-1.5199500879431377953;0.24961532597002758505;-0.86620000158586085703;-0.82480491015291246271;-1.0940903490815330201;4.0475979978648171098;0.60278994734685609203;3.7652295754849016696;1.421287907149870966;1.0108026146803532619;-3.7514139413118265942;-1.772303749939170503;-0.71046733876198342461;-6.5770053355864019906;-0.16380199188536495591;2.7628213295832084739;11.4501577259315237];
IW1_1 = [-0.62935312791250519471 0.30300647371153005549 7.5624976751328167879 0.38254991817905253937 -3.1256031551144562997;1.1694330496118414686 0.024551693875409136975 1.7690350088286019847 -3.2809301974116342748 0.99350936507529163944;-1.1728681535638381117 1.297468706633769564 -0.76871974034827017519 1.2925461374066944398 -4.6448661720024251309;-4.9142152217124674962 -5.0997910831689861055 -5.4355635263594246354 3.9892489918668649906 0.044340703587343606451;-0.43799746207780054785 -0.81759944734221878804 -0.18747500685068849147 -3.409873202316980656 0.7692791870944469812;-0.19606850850962548316 0.64407021279264642377 -2.1781029622339467267 -1.0257159996085873921 3.0419595054319947813;2.7832839145985994733 0.35764636532671167934 7.9774849834001244631 -3.8225893170217406514 -2.2266238556024497974;-0.61353719991954491775 -0.38636546250316561313 -0.099633121936780144612 -1.1999707949810689556 2.0509632345170794032;2.7514009886208823374 0.26026951325156905126 -1.5888603359792012881 -1.6300410839981802624 -0.28711370236599020878;-0.057203224682170233928 -0.28628139600414120869 -1.11489936664028777 -2.1424529588014911141 -1.1525095429330383112;1.3212678550409626688 1.0846847138987663062 -2.7162619928248452617 -3.9116617900657391083 -0.86122527346382504732;-1.3865530022208576444 -1.1393497958223590683 -0.66895244782569251996 1.3488269139315085532 -1.0522460463756657667;1.0743872087703796314 -0.11037734392819796569 -3.2974489427657571738 0.11478854239278102334 -2.495420761358252193;0.31323164150464438693 -0.16941255818673209288 2.0470689267148940438 -2.7775297044623394882 4.2238727266120434933;-1.2197292681602089903 1.7805059630519706193 -1.8264624409563103846 -1.1152284700862249434 -2.9518137864633366974;1.2923965664716270396 1.2912117140859102893 3.6546349759311076433 -3.3954240057150295584 -3.3802828572691634434;-0.24478104157428509646 -0.3012525193725243966 -0.57380305163763933862 -1.0039482063343483365 -7.3196136721058175922;-1.0582279981416704562 -1.3598771487921472723 5.2772473270074087282 -3.1349613393246902326 0.80784003229066647656;-3.6986044726006306327 -3.0265132339831715136 -5.7562774110860015497 -0.53696990487510132262 3.4108376043698847013;0.85058367868970230496 0.81402937340978620551 -1.4637789795468618692 -0.6887016370887943939 1.1051770979364612391;0.77842377294690678902 -1.5765674040579298332 4.9221830631596050765 -1.6617321801996323138 2.1466561787403861672;1.1780010815963968263 1.4226605763606983146 -6.1856451740641107406 3.846575651569473564 -0.56242791977726847019;0.17067296470800341623 -0.32535489322764216125 2.8871898522044152635 2.0048292662274369569 -3.3641759012615981561;-1.9416321302910146684 1.428565828738077359 1.0038592859130042179 -1.100847426456160516 -2.354549794438874688;-0.99908383077911677717 -0.99143416095438319235 -0.11538326655881379901 -0.81635372546845430541 -2.9184373696441165968;1.0521359651935902058 0.84956434213941189881 1.5305328303468452855 -2.4450821629304315685 0.37948834427861938279;0.2758359935566863097 -0.032359024091884791563 11.518316186997834549 -3.5425271574811696418 -0.43261915963578662092;-0.30865395724979882708 -0.61725638008730210338 2.0756979455527861944 0.4890681213224224444 -1.0254597263191504997;0.70365298135552434289 0.50891232346917547602 -0.94229384048998698908 2.0161682168433010531 -1.7760739887124841196;0.03306220691894643654 -0.059901194157208584312 -0.6495134259297579149 10.453756176039302517 -0.63046917266160595883];

% Layer 2
b2 = [3.8444358577022450518;-5.5473669077851006648;-0.52562752576058591192;15.581669801283506871;-1.8291127466312502303;2.5523611451558227969;21.875317832659060713;0.43105359449467112176;47.791758786507948287;-10.19589836440511732];
LW2_1 = [-1.9295779150575347316 -0.14398320003397666045 -1.1070019942643951616 -0.13457877492211417492 -3.8485975476750353685 -7.4419387025864889296 1.2362630752059844408 -1.4112216291697827764 -0.10675489500780169227 -1.2828270742044531882 -1.7748017877668860187 1.4950686981326397618 1.2833576881196582242 1.3229714725941315034 1.0282886293239965614 0.48316823921853341295 0.79221199079909632168 2.5574160513275052686 0.042217962324355602644 2.5320248122394239942 -0.16637692263773887213 1.9914297861350991337 0.12992873484193234734 -1.6139514725194290534 1.5847642638470729182 -0.069294535750635544158 -0.59227658984797393327 2.9469369061437657997 -1.5316506063616543987 2.7556624841101635681;-1.4659660766174409385 3.4252055282534046299 0.60858467325363918565 1.9100637644525255698 -0.79728028771848846734 2.64917386813211575 0.39624945082040435862 -0.968327079908355981 -0.23641860507351164511 -0.14936850983185914754 -0.21448096786260262281 0.30635466596199179001 0.043236279703191729529 0.27692497661833898048 -0.030437490936757426174 0.097411271877174465672 -1.2193655162974252892 1.631772618421611476 -0.033064737959910321663 -0.47392061921947092884 -0.48000096230978389622 1.2779189278514981876 -0.13364820128806612076 0.4550608479664626782 1.6664363988516333048 -0.13469920046630462496 -0.12305858254777954364 -1.1324581910712379074 0.96481955893285131776 5.4565055365028873879;-0.25492300089647063066 6.8917154052101716744 0.40983709663828227221 -0.65128389706845468687 2.6107875585552351438 -1.2732920775081921771 -0.029030661712362026278 -0.38839906874514501522 -0.84359174707243167468 -2.4079706951177644036 -0.45796094958139454967 -0.60979926307058907042 0.62377462586813359202 0.22695183867388504795 2.3783396111516657001 -0.45625055695988764315 3.013988004501116702 -2.3392526319419761904 0.54067072887193390329 2.7263235166992241254 -0.815689394703785986 -1.1285292349266622391 -1.3294430298398007828 -0.76188213041954067073 1.0392075000114233685 3.5348827913501312814 -1.4187267238452523355 7.1439199968248621175 3.7840861224924253747 0.46921154690322192105;4.8552326800451393041 -1.3280261410334790462 0.99603764805320882392 0.49349924891554253215 0.47680195789497858971 -21.750819851679210615 -1.7489133823957478953 4.3476463283904802637 0.85880469779001922248 1.2518400924574648592 -1.8132752735615198514 -2.1343957046750712259 1.6074846308796362226 -1.5815014697921869757 1.3615976893737777864 0.42387881804372479211 0.13348712328421588613 1.3592607095916471316 -0.75869812814781711907 4.2899095190552687384 -0.18525679622602611585 1.4757805022947774187 0.43214503355959432707 0.22999990017594523772 -1.154709015312589182 -0.51800493534141645036 0.34557142354266734419 5.847219189596227551 -0.71732265536482964219 -6.4763078439078674009;0.84726085951509488847 -2.8452000254536291557 -4.7643627568336794909 0.41558387641414867852 -2.5532328196894957806 -3.7252763739960448319 -2.5998619674481484587 2.6370418582605870128 0.98356394011667869215 4.2873912545594290435 -2.5648860645620388254 -2.3441912401925457665 2.6131756955165501566 0.38663086973217963838 6.1121098085060499017 -1.3530205533789496641 0.23241669093884173569 11.914551582218710024 -0.44055661097970710838 4.9520328392545156504 -0.45882444521422599548 10.994863037181689336 -5.3778473916541393862 0.58374842293653772529 1.2330224516350603015 3.0416835395871428815 3.8839507356389582071 6.6064017576882028138 -3.4794342916384590758 -6.0076114774315314548;1.8550262586027237433 -2.6974653809710602559 -1.0706274035163678704 -4.3776997243961437434 3.9800477081608334196 6.7199829268124720016 -4.4509948586082845878 0.78136459230680865851 -0.52882863297253612789 3.5562877151623886185 -1.5438778223797520184 -3.8898057644891106399 0.59807989849189413523 -4.4397763530135367915 0.487639251370475868 -2.2957374733406221701 -1.1012360955059352641 1.900136501538578937 -0.98624585824940558254 -1.6523669128060269884 0.55023820443594806395 1.0767691953262035209 -1.9429964602269622631 2.5093272652652678367 1.0426235715908767077 6.4835702537960777647 1.2976498445942590276 3.1297149706583367035 2.984081048532926328 -12.655980458453399251;0.10220286643498746604 -11.833698690483418048 -2.8610609902170645036 -1.2514877371962704888 2.9808008594909098932 7.0206864253722063296 1.5666031908919229032 -7.5190274270148576008 1.066127879259504363 -5.8414347496211789235 0.59848199997108864068 1.6927947496198894051 -2.0310634037672676833 -2.0191998093553751303 2.3720012307962021758 3.1886635021742013585 -2.470463916848517627 5.2600741589985178237 0.95345184253993453538 -11.028996058639989997 -3.3985685440914576994 3.245741705414871614 -0.54428672356382357123 -2.3036528406466141305 -8.3157561601818237307 -10.334843908005245439 -3.0518041270714508784 -2.8174642641721119318 0.68949575434993370937 -19.620498992009576256;-1.7364897031372545921 1.333969629025947734 -0.36845731679887799848 -0.20044922538393250555 -3.4568252710584848053 -3.6635109784159154422 1.1337646208220191202 -3.4066416633585694385 -0.18973709573382929738 -0.22664120760208869565 -1.3939023776829615464 3.0745383238596430964 -0.24322979608618788827 1.9489124178993368908 -0.078726986040323271299 1.2455867671093541205 -0.16135051467573896522 2.1330862051812897739 0.36387155393805875603 -1.9531706991142581131 -0.23541995527994763471 1.2508785207007426798 1.3220673432504130851 -0.8565992916656489653 1.1268939359470464545 -0.93135978447655576318 -0.97143265397475753176 -4.7270157543146291701 -1.734672084706872397 5.1061657408862002328;-0.36906295543490402755 -5.8699012351283998967 -1.2436991880794390219 1.5684248065325159249 -1.3343022963257340496 -37.430437497406444436 -11.524396946056041813 4.3782655730350841239 0.78946832244014253366 -6.8012386547995440367 -0.43456078493994304557 -2.1601647284815115313 8.9836831306182975965 -1.3630116761389081681 1.0314630689174832945 0.64711672352960913912 0.14229750252736456106 -4.3232582185196637425 -5.614821918883698082 -9.813496827222346397 0.73049073053692314961 -6.9133550405213490819 13.200582109719448454 4.6294911095547357505 -2.58820560823141399 -3.2931684359965980313 5.1747638378513993729 -21.918910229435606851 -18.298213418136484165 1.4366326781144183045;-7.2075274873751871141 -2.2689780532494250309 -4.5473967726250767996 -1.1191993799870985082 -3.3566896058340303099 19.68152772000722095 2.4102197341921582918 -3.3277581900017891847 -0.42149232664533792247 -1.8367885166686728216 2.4963554864713497494 3.5950334006785653607 0.62922993821963779659 0.60997727075842633759 0.58039284267572932574 -0.28295078045342286766 -1.2935291604517389352 -9.3934927062123243502 2.1110627607815679596 5.854981022759214504 -1.3316477317299548755 -9.9681867661126251079 4.7387036377167115475 -2.8549180269932081799 3.2697878702193707845 0.19535911183271281732 -1.1704970706294677285 0.18207608170555339422 -2.5184615759040003269 4.186696591289946312];

% Layer 3
b3 = -0.80131737569441485736;
LW3_2 = [0.6311905757580703602 -0.33449420451910799601 -0.5957302858045343541 -0.55081785831904472861 -0.6267214285176432309 -0.1573909019692305844 -0.34393299594300913435 -0.59293051033724342158 0.21874983911960729777 -0.40673446875849794679];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.568038194888224;
y1_step1.xoffset = -0.04979;

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
