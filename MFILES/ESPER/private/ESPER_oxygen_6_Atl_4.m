function [Y,Xf,Af] = ESPER_oxygen_6_Atl_4(X,~,~)
%ESPER_OXYGEN_6_ATL_4 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:46.
% 
% [Y] = ESPER_oxygen_6_Atl_4(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999994516886;-40;0;9.1619;-1.99149636107278;-0.03];
x1_step1.gain = [1.0000000017546;1.00003626070655;0.0153846153846154;0.000309885342423303;0.0656325902039552;0.0604107502911417;0.604229607250755];
x1_step1.ymin = -1;

% Layer 1
b1 = [-2.2929624680051921004;17.133549254377896887;-4.4224280872764651917;-2.348351529214806277;-3.7245669052646062092;-2.213876911321633667;5.0128330422893068175;-1.2734288563464453148;-1.220921824487038565;-2.4897208411682640872;-2.8069970388228098201;-4.3295794915349556931;-1.1246101061214792693;-2.2579938030441915053;-3.6265705137206416531;-9.0106333442546180379;2.4727818516541439919;0.89628793426215025431;2.356298529296100952;0.65367833351626636329;-0.3585770338060012552;-9.4508038759143833119;-1.5655708608535612303;-0.59671477614662837841;5.7722476112153513483;-4.8530410939444212914;-5.0764763832453940751;2.9158935883066416928;-1.0522333565619275131;0.92512997397602259042];
IW1_1 = [-0.72626275905250481379 -0.71532756852449874696 2.5624044838183461792 -1.5676351767251881508 -0.97424166668227951682 0.06067898561776607913 -1.9071182596215832117;0.25090676527326383205 -0.62361764061643754076 -0.3800101414871081551 19.611952002657776717 3.7740734859976461379 -2.2736758657223501423 1.1593755852958917973;-0.0010146455422012372501 0.19118594326423907215 0.073659814284389488193 -0.50654152185244583073 5.7038356642076388781 -0.051366390919311534502 -0.39428758932709867535;0.94700623374996895087 -0.34524128421631383734 -0.30495870054127927062 -0.16514704836525337384 2.4604754894557268052 -1.7319025038159223673 -1.1793631423064001584;0.55428546000757139556 0.2171504675514435434 0.16613921004535067705 0.30858541374727566087 5.7122852463556608171 -0.96205004022868911839 0.57389776887621424262;-0.16125942758463193694 -0.036675519949659812768 3.385466094798250225 0.072029109399395177538 2.4895817824966490761 -0.41189987009646428229 -1.2797110246440079884;0.040613554703793096345 -0.17346885042650594944 0.36052925720260403653 -0.48956562958277310971 -9.5382763690311040961 0.53630628808755687054 -0.54282325815056509111;-1.0848585731014572797 0.18899444137002721633 1.3922613451341467883 -0.19162652194500864877 2.5623722539244013774 1.0008149630026741406 -0.72431381398905470004;-0.83170262447331477773 0.68145786914245853882 1.0907535418967704288 -4.1254898206042049935 -4.6413503130938327246 -0.41651928150884542923 1.7030347601109558919;-0.55337723118331927363 0.16405879944307041884 1.0989565449859024504 -0.052725437999039781456 1.8930215082785899483 -0.24784896987877172436 -0.65813682198628109532;0.31518237416318023447 -0.0014686092231913658399 0.64964585472323077564 -0.10269599577613881813 1.5783320545562395942 -0.25811958071796614611 -1.8105248368106447998;0.89418048909997205342 -1.156830962509047378 0.11526493536689877417 -0.14070102636839337817 4.7885346326044899712 -0.3007883670724421088 1.6342252819845906053;-0.03012581717566993264 0.2586690036121693792 -0.42491719686585377458 0.21300601841506106027 1.9887338707023258788 1.9580616608587064764 0.545522288062848415;0.23517523678295670919 -0.29796625983490698797 0.7501280336955881145 -0.57851796795714027777 0.83854968365877224912 -0.94317984931124010384 0.2699682389247175518;1.5067170279266313138 1.1694006245951777956 -2.6381327017618465547 0.90413950364323603104 6.1149475873004544013 -1.8637071586658198186 1.8048960530797046875;0.46592719136233434751 -1.6120567691011464806 -0.94558137024885735578 1.042454970144533366 9.1496917970824878097 -3.1985116381236520233 -0.46281528743131411519;0.08239814420701296882 0.40715046811863858656 -0.101993478262261697 -0.80127161647780331677 -4.3140851308835923561 -0.4922697529076949241 -1.4866812769258366256;-0.77450645932405459693 0.52712393968866788629 1.7234281204788659547 -0.41032929389129729758 1.4837284864000255258 1.3748837834692759774 0.30471152965982750693;1.4519678199277663566 -0.94215587557028368515 0.48241134731188439755 -1.2265855307821240672 1.9009766100377323284 -2.1734740250007700268 -1.5985040796909635308;0.47088301022219769543 0.20372950120224264658 -0.61341611981466737102 0.70981756791631622772 -0.5017665717621822008 -0.63474588596734871082 -1.156427539973021501;-2.1590097262162442071 0.3523247439722195784 -0.12476051711651843512 -0.27767268340076994848 -0.4913702738232094136 0.25133391390456299996 0.085412189272992572464;0.63492593828877708084 1.3036365516688288579 -2.6421372093940456161 -3.4808212456911458155 0.27240651322512882127 4.1127373210995230579 -2.142879139251301801;-1.000585159882536157 -0.73601562726223590527 1.2756130833680967562 -0.91834951751445326185 -0.45853999653864968566 -0.80253017096422407661 0.53300632499277666998;0.27405419962213017904 -0.22076534027485403278 -0.067231585905607210707 1.9754240199448349635 4.7292397991958736725 1.4966361944860824273 -2.4690693387558106942;-0.0067056732337382992173 -0.024841091116378600373 -0.17534919448021302202 0.70812960826437720829 -7.7546081128906898527 0.23673554256793932771 0.097066004459563620665;-0.69940419009004695106 -0.025079767426673633263 2.58333400999049001 -0.25113159536147416917 2.8809034277032838922 -0.77416446269135263947 -0.7198073491495547982;3.5580326170122478224 1.1585108083656991962 -0.57713720892842923504 -0.093084278786573146358 -2.3581472519403261501 -3.1747074092116824851 -1.990730148234690855;0.039196683933562480318 -1.4053576646683914042 -0.23199142867346783659 1.1438902016431535813 -1.6107684328099671767 -0.70014801312516339316 -0.87207158383174221417;-1.5335502556412492492 -0.66893215472275047162 1.400472283261710249 -0.95500205816982763896 -2.717712530946402083 -1.2115948428940019888 -0.26956206874124477579;1.1895906919551992864 1.1103015506793734168 -1.4148316162314342748 0.95534127234306231991 0.87941163253145626655 0.75269461278909000068 -1.2534734964942038093];

% Layer 2
b2 = [-0.83040701179846654334;-2.4249973290399013948;-0.19476453999351103774;-1.1670611623635551712;-0.17117062486116285669;-0.98182255674429519399;-1.8702120539227913287;-4.6100298352043012073;-2.6249499854148723443;-0.7105773286946697187];
LW2_1 = [1.7483561551213493424 -0.0718530476546730551 -9.5398831006734514659 1.4932685762125283979 3.9815716325659589181 0.030964808913781366462 0.24542178657337920811 0.9894274755163743329 1.5164050988623283267 -0.10603534631277565514 1.4863282159364317181 -1.5696959028303127948 2.6191988851138670036 -3.8392715094542539944 0.96629450183669141072 0.85315411829113019238 1.4621850830215168404 1.1389196670095425024 -2.1769270951765782307 -2.3823515142642306408 0.84431302746483161137 -0.25686909458433254683 -1.0980297775883993694 -0.10655374855515772903 -5.0576898929980105279 -1.0700862827739694438 0.0028659293263216499559 -0.28041698257915692594 -0.91413329228499407009 -0.64689126043845768255;1.2634283127712342409 0.29490295676113609247 2.0897639011751079074 0.45965435764818479702 -0.17849405184565469984 -0.84838208532066905487 -1.1310982546089309864 0.20766212775530390466 -0.46846760149643623672 -0.10552085175637443415 -0.29771739868361657422 -0.50239693801604545786 -0.34873872167883812612 2.0007852325733361809 0.87191309562314733839 0.39931243212957934219 -0.023807113575203462524 1.3650118223454881239 -0.40208900527787250212 -0.75592722922382438178 -0.18466840898495456846 -0.12673173082693531266 1.4726445212206720203 -0.050549929910230786967 3.5474946410623555693 -1.7474261043471832622 -0.2491635215260197822 1.6066832013480374375 0.36153165603836523445 1.3437650994915975033;-0.42633557655878162596 -0.33937116726504884534 -1.7648296090638528888 -0.21073489836724249025 0.68229288903058138427 0.32399737756827651225 0.70923099975046199894 -0.065442315474670290087 0.47612043917301760665 1.7448588926899535334 -0.25837380084661537571 -0.65901664250004043932 0.43751516118382893783 -0.87348859796434030756 -0.27187375270410935046 0.010615102262922091608 -0.092217271789488713951 -0.91729031278076189615 0.64055209073644536755 0.64858603378189705335 -0.37490256074964362076 -0.12749575139768981602 -1.0468492091366434771 -0.026877709905955764896 -2.4270810068834984818 -0.52319401566457335306 -0.0010394801873063409716 0.27135966258307936405 -0.2981033626078001264 -1.3134403009903665716;0.62300573645067292006 -0.83125865236874696418 4.8636191211484174346 -2.3778640493864555694 0.7039512073922307156 3.2416058520789983177 -2.9620236571878839626 5.5902041753213058684 2.3155573532060174102 -4.203479492140258067 -2.4523582704801940579 1.8586747347393055563 2.4217994986383382106 4.6485693007975106283 -0.40378733632792140895 -1.4563659236567416677 2.3828678874908648666 -5.0041506182635373534 -0.95444973983995262667 5.7551646544815167061 0.40960001879536866998 -1.2195105555230776329 8.1261165846595240225 2.79887484961582933 7.9254834226947572517 -0.91255519506142612762 0.71813011151461536574 -0.80511353276205410712 -2.9127294254410394103 3.628428231843388474;-7.2970355790130350115 -2.6123476186257077458 -4.8065509566003461472 0.6787477369856182996 -1.5613905474477109969 2.3661191156361653931 -5.0663927023356682966 -2.6893132936866890326 -0.73268320430451316927 -1.9870359100360714955 7.8405290128845850361 -2.632512422189500878 0.89768216346513063364 8.1437526275974612133 -1.8851649834072496947 -1.642486088976578662 -0.80203160292830810718 -6.4028226404400427896 -0.10630430714065114417 3.3888325435548227738 4.1956169590957506088 -9.4415674592507450313 4.916920080732993803 -1.5833603852568196579 -3.277318614523146767 1.0609768264099770363 -2.4369831072346141987 -2.6692523686215494649 -0.20797864270596300162 -2.2810313733663480917;-1.9847258655437536312 -3.4541902999830358389 3.288201585633524715 -7.8921043970981310522 1.5524008021223059117 -0.36389919026394573764 1.388753537514420433 -2.2880273032238980235 2.2190162746261554716 4.240087807543469367 -7.8568065151306827332 2.0716079877957604971 1.4587701145080376186 5.9078955716222170125 -4.3079662529141629079 3.7438178641901411758 4.11555262694112578 -1.2309343453695840154 0.014849055594402624236 -5.7634629644706283358 -4.4074858178863944858 -0.11816651546082240776 -0.41112346310162017904 5.358971207544308335 6.0212748972021348237 4.0209210174646816682 -5.1856714835611885661 5.03430948961936231 0.16834466983047946531 1.3761043091548086359;0.44007520595587196199 -1.8737222781895814805 -7.4144744879086950107 2.4818374553334137111 3.2739485922393281037 0.19699863490325569826 4.3220747867666373665 2.2823315530178183153 0.40605073040466449497 -1.3856454181817488003 -1.792649168197410825 -3.7879222075715066609 1.8625085606714775199 -0.25603985151364688777 1.4480037322093923802 0.24748656175174907346 -2.7478616965545543671 -2.1916525338668204625 -1.2379198626670715466 0.73621471755667533543 -0.27682918323205296662 0.15269427496160867253 4.8948717930908429707 0.33782462520388562011 -10.296367895471806975 -1.29074116157851293 -0.12070295671040001673 3.4119419358789544461 -3.8426515862850383343 -0.17672724301077707976;-0.70072548029913761525 -1.3776573392252839945 3.2667446495137077989 1.5584526419511823825 2.381084326149573549 4.5021026984679455296 0.93696984549798945086 0.81290325554539466513 0.72888578186726216135 -4.5145314735669401784 -4.150982822391773297 0.11535627553723443173 -3.2193760868680936227 -2.4372250613022794496 -1.6603023659742330853 -2.5722642059378308943 -0.098192690563499676615 -4.2318625794104036686 -4.5070069010549707045 -0.03219350855186836613 -0.38734493913381562447 1.6977830344789210937 -7.9440128775499818659 0.25101071670000063563 1.3294137512697710157 0.04155518112149150789 -0.26364527273179505507 4.7865450790610628573 0.61461228690819058595 -5.0702605150790756028;-1.5364313841157752893 0.269326747297228275 7.2627355268169937474 1.2959551440575627357 -0.77362443542767500571 0.57389485247250615618 1.8164843717587599858 -0.38680358231283668013 1.1436489231441571146 2.0923947854308084793 -1.0496813775232913812 -1.1684612136132697113 -1.6889583122386979674 3.3839302932401875701 -0.55437300541418499389 0.60003616793528202766 -1.8210159402273602414 0.45875716795344240184 -1.3896891418545063157 -0.90120258107672079895 -0.51889995089499119096 0.66897685695912667558 -0.8538448754809231378 -0.51943300453305341069 3.691219175551443854 -2.2002698301183967189 -1.3209409865165970821 5.3426182463242035681 0.18947839431051952119 1.8802130069141500535;-0.80090110868519115428 1.2860678568358412388 0.1776789213158921632 -4.654847062310216721 1.3866650962141087167 -3.04613092312346323 0.17972326091208079557 -3.1141591932392245212 -1.9657615240122543465 3.1041421586096742402 2.0749816132706264682 -3.0994468860005857103 -0.6294881770694957801 1.9980372225225182525 -1.9188508466091058757 0.2475572227858759411 -2.9427012821676341581 2.0110569536158453552 -2.2168874936671323361 0.24047464068899276901 -0.10430951117834978692 -2.0941726476964550585 3.6209347874120947353 -0.38806561344070400477 -0.0037413531633077312666 -3.1158705210252293227 -2.1994906496236787952 0.46818184934408996201 -2.832675918830387829 0.61538382466976526697];

% Layer 3
b3 = -0.2716431993590759375;
LW3_2 = [-0.19281680598536116156 0.36556279849095485446 0.69070810208835020649 0.17642682695859024378 0.032981701009492582488 0.02063808408082170795 -0.14485181158988694183 -0.074389407872620147422 0.10352603908462526539 -0.068816526642872910546];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.00368391969055075;
y1_step1.xoffset = 13.5;

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
