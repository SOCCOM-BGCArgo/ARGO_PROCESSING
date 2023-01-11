function [Y,Xf,Af] = ESPER_tco2_7_Other_3(X,~,~)
%ESPER_TCO2_7_OTHER_3 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 12-Jan-2021 17:03:55.
% 
% [Y] = ESPER_tco2_7_Other_3(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-1;-0.999998793562952;-77.961;0;29.715;-2.08184650291172;-133.803853032615];
x1_step1.gain = [1.00000000274156;1.00000063368067;0.0138606862425759;0.0002963841138115;0.248478071810163;0.0588359962254402;0.00432440760873659];
x1_step1.ymin = -1;

% Layer 1
b1 = [-6.6138443349704845;0.68702227850740005;2.0020637532913739;-0.61284400467912359;4.6908228392180531;-3.0339242227738272;-4.3950428568513757;-0.58679276799825864;2.9297571044023818;-0.11500156475831327;-0.17347068579696537;-1.6148786838128415;-0.0011135059296847336;-0.42094737634625717;-0.40096158473477483;0.099967300747467638;2.2959097940374797;-1.02239313927263;-0.69650943078023475;0.89329290060856248;-0.14661952586241789;-0.045994123510268707;-5.7948921763405687;0.47271436789349947;-1.8744514277966482];
IW1_1 = [0.18902068726081492 -0.0046128293885170438 -0.27608098955877913 -1.4296946082439164 -3.7706674574715224 -5.5432469361163577 0.92548003949184576;-0.2104141001957231 0.42882921157376325 -0.6959718102885154 0.40239931038812132 0.81879985670215494 0.23651462434050791 -0.97256292537289835;0.2589285030614677 -0.13290859488426982 1.0178641697032338 0.5172824407365143 -1.3587794073861132 -0.70025654935691073 0.17639539469075763;0.44200088672894283 -0.054999862918255663 0.68809840731209992 0.65554241985340755 0.73994499938623248 0.49136543700072854 -0.46665732599888549;-0.60443760272610525 1.1627956604988718 -2.002371218857232 3.0604721111642408 -0.3498780246434498 0.093162851803257313 1.5745962812729204;0.33418576261552496 -1.1561112070772352 -2.3522545622820163 -0.0088402947338884711 6.6703900262466771 -1.4167694223437659 -0.71876131126292209;-0.1701097161494817 0.44102184653880278 -4.1330937299854575 -0.53090814112695728 0.097947677090675284 2.7219274560622311 -2.2970868320347657;0.42874093564845256 -0.59240707653652025 -1.1100756766317055 -0.17406404378887669 -0.51444030558156217 0.33086890091431653 0.58301347216780441;-1.6053151500891838 -0.37482860691534209 -3.5178199922681803 -0.14289389262323154 1.4116108660680016 0.069121677731060052 1.1068531668987718;0.46766451471097314 -0.66385408406746182 -1.2041432729332366 0.44370599683546913 -2.0321553491304716 1.8417083989305174 1.2036677130811948;0.15134017334534333 0.42283891466670737 0.22717758594604867 -0.61723397540742242 -1.7074578159203238 1.3678224316493488 -1.2412950670824145;0.40982501264282167 0.43821816481218545 -2.0357901935877623 -1.8586865102564623 0.23001326828959606 -0.43638567527802108 -0.98320123778124591;0.087111512189850276 -0.41230963339478799 0.13908525532725668 0.5605390193686095 0.45330783228043159 1.0111325181106667 0.17703963891947952;0.59098953024808443 0.41703185126271741 0.61874142090308681 -0.8223577377752489 -1.0697623665813309 0.47110445280480773 -1.0879407365777392;-0.4354077242242686 1.4910382630688352 0.0271620830963584 -0.56591600293227762 1.7629086712462025 -2.2884289025707942 -1.9456902465053147;-0.42511439889852348 1.4367497655778811 -0.010289118515614062 -0.23811581107109972 1.3126667384890951 -1.2653122863326689 -1.4758623850774597;0.1236622982439684 0.095172434759084198 -0.37275772077846492 -0.28868561182209845 -0.67764233928107731 1.8735135178627953 -0.98138796774360737;-0.50555907135670786 0.97311199143241822 1.9893756020911009 0.28355553025765273 1.4111767455662476 -0.51081083809969585 0.917176067958418;-0.1506742672425975 0.69299107526812176 -0.92736852819005 0.059247436519226683 2.0308092329810381 -0.25963283661737674 -0.98725589279246884;0.23078312364262399 0.10244321911550863 -1.2316440959312638 0.67496179240966514 0.5837706419470593 0.56808438362665303 -0.041398353021032765;-0.14839257522329105 0.37418850783038099 -0.58865060231818434 0.077260708654104046 0.93382318630726147 -0.1493061577295951 -0.64152280086904234;0.29838560546588488 -0.0028072277914052979 0.3638913222206564 -0.74786446801489337 -0.077312518825059989 -0.73387987403136845 1.1543147404530905;-1.1044729414137664 0.80032986541312823 0.13024286446564701 -1.8937125669516988 0.60257381535428944 2.4589615583807514 -2.5423695351265931;0.61554308936404623 -0.026807581040028269 -1.0610372097484373 -0.91431635056460736 2.3188498037057563 0.28971795512379878 -1.3366918152822127;-0.56330606111332704 0.31393261078892803 1.0948466518525088 0.55080890572527597 0.47569675760259938 -1.7510872539682951 0.54611118721326279];

% Layer 2
b2 = [7.2529001386943754;-0.30905830869316508;-1.3900834344043327;-3.317802134186262;0.47913556213453706;0.74427033141626131;-1.4001315976059314;2.6996078614633419;2.5864157943022854;3.3939255373732178;3.6666298976039067;4.0144855084257456;0.067829005844834908;0.61567141098422817;4.7273450147783631];
LW2_1 = [3.4693184386220635 -0.75226228498904368 -2.8224546976018816 4.6309864410103101 -3.1670759928783698 -0.56185555079638794 0.69652507707424316 1.3603939616570018 0.89700078778400205 -0.57782489611144183 2.0565360628887315 -0.80259266380756689 -6.4815547154844273 -2.885698989611087 -0.27430264924114223 -0.057843579537357422 -1.158693808402498 0.81714927579445673 -2.0096037385882695 -0.30948557347018196 3.4513333182027064 0.34475235560779666 -1.2310782877942794 -3.7682812562097849 -6.2278572667128715;0.12884937770031152 -1.1159028179377974 0.87368021184328493 -0.39353846702308337 0.0313195203452779 0.055911343897418621 0.089253936355750602 -0.36355037925747569 -0.010871224731642329 -0.068261242860378768 0.24276316516292826 0.07738554342920638 1.1444565193034262 -0.23377014834727908 -0.031687153128836562 0.12712493897897792 -0.20116256750055445 -0.092334815424701816 -0.37486778746114879 0.00031439703399797686 1.9977360183355417 0.44350423559022289 0.015889466198398771 -0.0056138654354668865 -0.14715797227226288;-3.0454849035569476 0.39985779523892384 7.6783014358466781 2.0568330353867696 0.47503426647538438 -2.5603213743514313 2.4001100576862346 3.963210632224567 -6.9444987357715071 -1.1460423414921503 1.6450889209955131 -0.67132687090473586 0.61424195643912682 -2.039871704006381 -1.5464383811832285 2.5575271715977821 1.9955066408926303 0.76151591263387208 0.23107507915737863 -1.3585595450519572 1.6283779078242737 1.2315973689734454 2.2249400925412388 -0.80520348027662836 0.043459720632341971;-0.25445221572037663 -0.33832468963527779 0.19777649139671627 0.60481852351259424 -0.41725223761985181 -0.39046416769244596 0.50274307561949505 1.2478380129625657 0.95269491297892106 0.0051496979100395402 0.39951384616040581 0.28000456198501855 1.7481611452437777 -1.0468668864570918 0.25140837584625431 0.2952409643753951 1.3628699544520406 0.067837906600790587 -2.1222043201782439 0.92865955955865864 2.9331786266283535 0.86318633027243163 -1.6641314640307752 0.88175678886069542 2.3115525052548898;-4.1285564738108382 -6.1792306387790168 3.7826050476621904 2.3782927002857868 4.1340549258563604 -3.6177607979996695 6.6031245315028659 7.6283075732930641 -6.5390819146373236 -1.305219364982634 -0.51503479640894312 -1.247071877382157 5.1007728524680882 1.88505993497789 6.67853711892006 -2.5933827013717763 3.6099369019005838 -5.3390963520012429 7.3988014322651079 -3.3006672286017662 0.33168598348511757 0.53406353420555397 -5.5154426609150748 1.4909029932506805 0.91268147909082531;1.3955017491589916 2.4247669911596148 -2.5769632502237161 -5.5810049495532787 -1.4186621393665251 -0.11180696819102762 0.61611333209797214 3.0553098438745159 0.98117474103707814 0.010375222320159156 -0.85303600080766839 -1.661402343697604 -1.0917578950429914 2.2564289294723552 1.7895914373524084 1.2398529571936197 -3.7654243140249184 -1.5387793193972721 1.4308817639297875 -1.7875297354901083 -4.1810065619322634 0.36228872381667199 -3.1123781337188876 -3.9486057624710109 -2.3839284210522118;0.30683528161079537 -0.70736537373534103 2.0609377255370913 0.76075396922503502 0.41810775260953892 0.35745667877953369 1.2220178392116183 -0.66572490074758106 1.2214300068505404 -1.0592563163062878 1.2978986550538179 0.65037835027680269 0.81442926559605333 -1.3105169944870454 -0.84519138081809797 0.33083547150591014 0.77226255854229431 0.041293699771908159 -0.77550157343437798 -0.009996261014217649 1.8059010584835755 0.164106994271281 0.083623249489588089 0.1937344359964813 0.62240880430039991;0.27433220558546006 -1.8295989360863687 -0.74599765826163844 -0.2242675186130208 -1.3111399916897495 -0.47789320024023796 0.2525911006863123 0.43639418750467712 0.58603150229765122 1.3752538900011742 -1.3233755447889208 -0.36971254168079237 0.69886323617753221 0.94155073238519971 1.2714371863417155 -0.95615085487277118 -0.84118666168704548 -0.63112771813449997 -1.3756625783982197 -1.4483055815441339 5.6882918544703847 0.35011155881757711 -0.19594095701957975 -0.25509538878898408 -0.64206481563525897;-0.30698828356890562 0.36717506923522125 -2.541704595400784 1.2938231232147472 -1.1568172882054373 -0.53519911663127839 0.10514289775087146 1.2844019394204214 1.4250198835889187 0.15356914268086763 -1.7534354181406722 -0.20855730173434078 -2.2012439353267679 1.1014071846244848 0.023509669300589835 -0.49453274624693416 -0.43294867667192832 -0.25916888075223732 -1.3575914197997276 -1.0918317881426289 2.2147039268538946 -1.0353514341854309 0.077820649671683145 0.070124838735303521 0.41213949246347631;-0.47302931357804373 1.2235717129675201 1.5783372938545783 0.37533208434131843 1.1396856895322642 1.4149809926420442 1.2932993592670639 -1.2325169987405862 -3.6327287036881613 0.032111648897784083 -0.0070129684671145619 1.3583214059989044 1.6174800038121513 -0.54925627220529594 2.5408775233106415 -3.66572970695686 -0.40405627455377552 1.1281072216112433 0.78030180846182262 0.97059910677150807 -1.1905449231766099 2.2701168161124716 1.5085807044952133 -0.63668118207473778 -0.8996739349128392;0.99427412259242143 1.6466059081019513 -3.8082625011947129 -2.4924211170846706 -0.4490475002116035 -0.72321753884036244 -3.2651560714775778 0.14739342809484618 -5.2738775306982273 -1.2469208655986768 0.3424365805123718 -0.71740519478895803 0.72607794431243 -0.82167240217737403 -0.56476286948157273 -1.0482063710943628 2.5800016234305088 -1.6634075687948053 0.82400482487244475 1.784678918541923 -0.64041998957494339 1.4697685610121443 1.43121623963819 -2.0245124843536595 1.531489334082109;2.6637793926940141 5.2313966355709915 0.12022715840881537 3.3540143200617387 1.0976939411031128 2.0893847913273444 -0.91715122922276904 -0.82058607083953894 -1.9655033429094413 1.754345850074476 -0.23504038440720432 3.0640201298856207 -0.063715560813307051 -0.49371606676760427 -1.4152116371196255 0.29203408819195997 2.7914899376240503 2.6707768344902409 0.48910133722993832 -3.4910516191388798 -5.0275290169141229 -2.6799400288016546 1.8588290947408781 -3.5711802710546916 -3.7414048716188919;0.35754865510177009 -0.63520220339808298 1.7327479184476449 0.91260774419395874 -1.4804541289171349 0.55957484249176037 3.941481421477937 0.12405165194256183 3.2196860465059052 0.18049503099510933 1.4963788490779018 0.12239529750562005 -2.551128471372873 -2.6131411116987184 -2.1130985723128264 1.7013313272562938 2.0738085071224068 -0.016129452564957805 1.0368010557365641 -1.0269827133240319 -0.56177196960473408 0.094111061947385263 -2.2073712311695948 -1.9608775681902448 1.1042375587373763;-2.9932697561231367 1.2704828057327826 3.3720727977923644 -1.202217219114867 -1.0489545216355984 0.48165883481860589 2.2658082759331153 1.7153556667190741 -6.2689317783886729 -2.5632585738997036 0.38989663575926209 -0.028540989312733297 0.97294904013099259 -1.4817116271235196 -1.9296610193825057 2.9557554791148171 1.4367994582532009 1.3937246561196328 -1.947446508641532 -0.66154314752188592 0.79667169697975249 -0.68944253117295229 1.6580351092043828 0.006186341136376834 -0.12697039460896986;-2.5754233037345462 0.08129823141661556 -3.2595536509094636 -3.0347053128893715 1.4002273163186816 5.4925340602343251 1.8202532219168324 1.4907313786633176 -0.45770002689194017 -4.4512337471159187 1.0366811599693855 -0.35519669388587166 3.5678957204590329 1.6116735108217386 1.7726573456142221 -0.99640700687374517 -2.0461533954783651 -1.431227144078421 1.1659648930623208 2.7036723572631396 -1.7193502604881568 -0.21907950791387307 -3.9749457897901515 -0.25523365925122488 -0.89554352081523692];

% Layer 3
b3 = 0.35441427513958312;
LW3_2 = [0.049460743438616443 1.9443392414290772 0.12439193899174085 0.28594705759766664 -0.014818307167287244 0.037014506911186511 -0.5507671153234468 -0.49822561478168426 0.43990209156519877 0.04855676282001406 -0.058079386845430452 -0.045710962420128003 0.078658449566114924 -0.11422445970831407 -0.04390011120708754];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.00314743155299998;
y1_step1.xoffset = 1774.26765577007;

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
