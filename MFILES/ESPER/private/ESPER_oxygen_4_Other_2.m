function [Y,Xf,Af] = ESPER_oxygen_4_Other_2(X,~,~)
%ESPER_OXYGEN_4_OTHER_2 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:45.
% 
% [Y] = ESPER_oxygen_4_Other_2(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999998793562952;-77.961;0;29.715;-2.18999082656054;-2.02];
x1_step1.gain = [1.00003298043486;1.00000060321889;0.0138606862425759;0.000266418009857466;0.193068829037552;0.0582073296585987;0.00819772922900357];
x1_step1.ymin = -1;

% Layer 1
b1 = [6.3333038005115795954;-1.9626920039988755562;2.5581586600165278078;-2.3233269215855547607;1.9208856562932710155;-0.3497736536475370972;-0.54863870880071996528;-1.662291242181861417;0.76995242744232439414;1.2154004216046916387;-1.8202513628095746689;0.76986861303027498415;2.1864894979443660894;-0.056081111871410363467;-1.3004120828940011201;1.1576074454275830927;0.0247437017304816681;-29.397213654718822085;2.6639320819011800268;-0.43980975147038192663];
IW1_1 = [-0.014101556846252479016 -0.052324409607501599917 1.2812921081465886708 1.8698766685954033573 0.60380852432219789261 -1.4295970685518148979 3.9079127551877239632;0.15764527564669983395 -0.12236422829170662263 -0.36430443346646340652 0.17071178855859506229 -2.1178580714237260629 -1.5581516359444342879 -0.98273106080766625681;0.19039147788392629668 0.13543535572327081162 -2.687448459212884444 0.5666459667541567935 1.791835415091292516 0.36870490017010598693 1.2946413718103608748;-0.048462882069191609569 0.70722257898097107276 0.36282789281627864098 -0.10223322909594305397 -0.71305336025095356067 -1.0810746386264920016 0.33446922102095749674;-0.019738291648437544801 -0.5479121213037762006 -1.7966787660258847215 0.77798825477056809952 0.0077195400568231282296 0.45808468236168653487 -0.56608315729074965006;-0.033207485934555837359 0.078379211098184756179 -0.25029404358530382702 -0.12700107127465040668 -2.1362269979540493559 -0.54874566160683035942 0.84661996795057725684;0.092951479306904505329 0.21602473662224805362 -0.59041028098310954775 -0.12808282242485558022 1.9355805280319637873 -0.29898002179880506235 0.65437869827296046843;0.17688479013770400106 0.72249553340508310129 0.7169072608146154213 -0.34932332832032819692 0.089338363825532188245 -1.0747157212290383566 0.52000607980484236936;0.017638394011263124062 0.23064145492310930696 -0.40438987772873563609 1.2480621199725021597 0.51048689659884671777 -1.4050213677295375625 0.056339229958385791952;0.52238893101080408066 0.18942164435428479607 -0.67714740851707999614 0.12026929897256367907 -2.1043958750044997608 -0.17759093520526009202 1.4411576653141224824;-0.036456470736934712407 -0.30370580002262281827 0.65788627160437840313 0.15438820839531855467 -2.1739957596443999499 -1.0488613345637991525 -1.3407252428004579059;-0.029752572338284637848 0.12804659924796257808 0.068877135709147371823 0.14771966834488364784 1.9713267793081084012 1.2228717397700514802 -0.35814879927720599806;0.2390930134242722771 -0.37565278530520140876 -1.873998517327470692 0.59951214344225600694 -1.3797838483755793249 -0.28949984339239354103 0.81242407945221328802;-0.39652385986079735636 -0.46591982310577834125 -1.011206905212485907 0.61429077650542396327 -1.1930104794493101572 -0.9612905432478210388 -1.0310645062274059747;-1.235863531711359764 0.34283409019478322532 1.011515420053026082 0.0023800353512290358358 -1.2893766415487115662 -1.4067479763773371726 0.34428445796257461931;0.002242396881071070068 0.14423012031194573646 -0.16605125152963762236 0.094652954989601356672 -0.24184503313059677043 0.068647950214551750614 0.82382008725929789517;-0.17145491272560683949 -0.15731242111358870184 0.55516361480143394669 0.0041075888975430309077 -0.84961656654518891152 -0.14059761124519484787 -0.14265116101814651062;-0.01088547707817202892 0.18989242796386113232 -0.27708108397906117126 -29.957850860358426814 -1.787030120080862261 0.97916063236425388006 1.2395906087098906845;0.11525165988383158366 0.16747772182699732779 -0.34000644993958872853 0.0040390415055396771196 -5.4513458658127253997 3.3765991772768688683 0.16511497364593299975;-0.11093446573925207221 -0.023101673032358334414 0.78977596404651695394 -0.24689056266797154704 -1.4718620174577368154 -1.1372491682097405352 0.79674448299248290528];

% Layer 2
b2 = [-6.3475469540018103487;-5.2702882986284658173;2.9690258793137354054;9.4970962385134694017;2.2460308734698561928;5.6957749157947690222;14.135671794457566719;1.4144376075728180009;-9.5429156570438209428;-10.21075227388268658;-15.613759823765025914;-10.487966762188143832;-4.7943297232733295132;13.82536471265319733;12.41964425261328131;-2.712899637060476099;-7.2431232274668024829;-0.41627751131537993334;4.7370834743520351395;-7.8480593573548667052];
LW2_1 = [1.6429229666739488813 0.90899955239527030493 5.9611840588483531178 -1.2258881943907262713 -1.1040801749652082719 -5.9574652237043377312 -8.8549540783669051791 9.0055852528516666666 -5.1388335541938907625 -4.458482131790535874 -3.0143351361551675716 -10.408664976707624916 -7.2869650317529695016 9.5095027221127512718 0.81219085146986402801 9.1562545116066811346 -1.4503269418467852336 -0.29274987927563989043 -0.26380818885079015956 1.0442853726464136344;0.32758490226656855615 -2.4953677169889862952 -0.13356979472577659718 -5.1380477650311835447 -2.3140646875988362119 -0.92583816044187394745 -6.0972141778056894168 3.9825596880297315927 -0.42334869640419081627 -0.12943551203793385307 0.013124063272054998369 -0.97523578111141107083 -1.8009404034745533085 3.9257939915645465589 0.55958277239466913677 1.7013153501596209871 -9.0164628312606627247 -0.23326334487243186033 1.9484582609135070452 1.8783603101891499509;-1.5629631784814610462 -4.2981367260406146968 -1.9318920924921900717 -7.7575572385058890745 2.7449189351680551852 4.4873136615998623355 -7.0824201882187525214 5.7137940423152508984 3.3943724001004822455 -0.029587254071357073359 0.28345687512228384675 -2.2806118701322151132 -2.2371013571661211294 1.7142994505860815746 1.1992160603704136168 -10.174825399737507325 -8.9291270959342643465 2.1738804900616504057 -1.7102743253364969256 -2.7152251231018986743;-0.5769752004612129781 -9.3748423659139170638 -7.4294006170646369469 30.7946986472903248 -0.16374661728778217018 2.6903510033622271003 -21.407480302748730594 -11.732763656695199828 5.1297064755068504383 -4.0944200555233702943 -3.5868020567616349936 -7.7017292053629011761 11.78279942122929036 -2.2743882516029083618 -1.9646913999340946955 -9.0646258275125273229 -19.043641240039114848 -1.0265187165001108394 -5.279458838655428643 -1.8901763701295022368;0.82970598780042248421 0.021443779008745333936 -4.6096041480977021365 1.6533277119984579606 6.9271375856162409335 -10.894123595465703858 6.6392050863171956365 4.2410375797991317626 -0.76851860180585351845 1.4452996833351445805 -3.7038302766536950195 -11.374090642311182719 -1.8809871386328882448 -4.341571981452746698 1.2344647602454883462 0.55905059757340147542 10.514050855064718704 0.32646675817284609433 2.8485715276703169607 -4.8790833560558537485;-6.2387901140210377804 7.3842112726232516096 -1.6582592737775423419 1.5094256476629088759 -1.0282861918965766623 -10.848553766103300688 15.468719530915294058 -1.2058042055444666563 1.7530579421011389663 2.1169014831195669402 -5.8726450367356584081 0.15508991740453428876 3.513742544749916874 -3.3036727739144646776 0.021876212843581772927 5.1559042051830390463 21.217545523656635709 -0.079087157824563672492 1.6717579227750987148 0.58185233329066909924;-2.7664547672732284767 0.64581223545716681667 -2.9380987484422917966 13.480236739799078549 -0.65105648686169093065 -9.6015121520501622854 18.121308069269836949 -12.259697298818711175 2.1547960740540310631 1.529851956812524616 -3.3852636578806110101 -4.4069319789264342901 -0.4475232145488355906 2.6150397639407461092 -2.0998013511732240843 5.2340233922680967638 2.565094616361668578 2.3787000917757232799 4.1256955884945938351 5.2228662329308912504;-0.99895751911739205653 3.1317848063909305978 1.6234408553132173569 0.74408219069856518235 -1.401733327929191919 3.0613871402579562186 5.2958743812678381602 -1.8603260880430301238 0.14921117882140130795 -1.5430252548134388846 1.1728444085345548498 8.0967629143508155209 2.0166363050140598823 -0.048669233831244697319 -0.62918340811122819289 5.1876907596439529868 6.5109855833937437808 0.12282476332852455825 -1.285636623669829115 -1.2381321677864127739;-0.68975182710954063836 2.9993942301039493437 2.0266867746587378107 -9.4132983811342203495 -7.3046300358417051157 5.6112786468403186291 -0.32823690587281101338 2.481458904692588785 -2.0222427438826500179 -0.23777631103216689823 -0.79602532816800486426 11.672078594646924543 -1.0005255076165004979 9.2267409242562692384 -2.4915582171465224803 2.7565634268982441846 -1.3642352119292902035 -2.1042379026024682354 -0.084016000557877085209 9.1706326147336234555;3.1390902889519276542 3.526236622660961828 5.1642347400015990999 -3.5245681289510435441 9.5770770862208411955 1.5841780770123787825 9.1248700492992256272 -0.96074301333191614027 -2.909137709327568988 3.9801538469409916132 3.1655769289711668257 7.7298888102617846485 -3.673582201959188609 0.063354437595989993026 -0.32541278244113858165 -2.9106939544161956057 9.5874779097035300879 0.15084997182448842135 1.2798486661800432262 1.6462710518325041065;4.2379975737263970714 -11.173855655121524677 1.0547351902179482597 -9.8321486617489046722 4.3137953942306301869 22.407578238119405967 -19.442832947518137843 7.043192614125199924 2.9843714944943551259 -1.5836177886914928337 8.9811509324417020395 3.8584288928192425061 -6.7908077053276612034 0.34256651783363745167 0.88632173942197778338 9.2973436310523318582 -12.101586878684617687 1.0036298883351422173 -3.2215104799865863505 -11.145373589334408138;2.859887594322338078 -4.3762745093700097243 -0.61747128978861254289 -2.5431983756513676376 0.64199189083496244201 -3.3465362848558211084 -1.381004590052416825 2.8498994430602513361 -2.1569406923519314212 0.1692634558318386051 1.4794991563678276947 -7.0542459740259513268 1.478837958542876807 -0.06521572196478167549 1.1826706318455859623 5.5657857969333912607 -2.5836287962768125226 -3.5035011217653182491 -0.71165490585413071489 -4.5422646688848233865;-0.55432000277235615737 -4.3101254486663247434 0.22999742835333716884 -0.38899977819560560555 -0.065399034048286286414 -5.7747560496692509346 3.1421935852762659103 1.1974328536325973271 0.18141219226345553506 0.52512415962821612503 -1.0124121175350335733 -9.007596162357311087 0.05789708603258031383 -0.40458300168851357714 -0.23343232432576091484 0.55220370234564253309 5.5612613304153395433 0.26906389024237115093 2.7490509793037500863 -3.5723742164468412774;-1.9266394436789757716 -2.9594036480868086691 -1.5027128973072583218 1.9509004737256594453 -0.62607201118261468498 3.3164448859318365059 0.24439517044185679606 0.59381812435688730023 4.2838419433928605073 1.6955984135886463804 -1.1459012096518479407 0.6638147598772746738 -4.9493256962980503033 1.4466840422517646125 0.30007594835927220434 -13.297931012505079806 -0.37115531309387084224 1.3397201747001643568 -1.7080365158260872871 -1.2460019587451325318;-1.7313079542377682873 4.7818467181273032196 -6.2612517273726639999 17.249198880312171411 1.2514334708813503738 1.7118903646665413198 4.547522322863462918 -13.703993547861324842 3.4264085999320728071 0.84591984501049632961 -7.9556845686769026571 7.4603737695448693046 -0.16290926411980888311 -4.0881824454103936617 -2.1143690860105732732 -6.7296115509754708128 5.6981352639816638117 0.17655724320115478987 1.2954943131806235801 8.533795311203817846;0.29704102703387536 1.0433080412042259333 3.0279234399216914397 -4.5121709592469443351 1.3593261397712408378 -8.0283144576534617443 10.028310275224297854 4.9607475350231844402 -2.9268321380851896052 3.2629250933466535933 0.59144903149357419547 -1.9422699602553801235 -2.2993635820546347581 -0.096912535373391645033 -0.50979980308711492221 3.9902018499222693393 12.894727151106897267 -0.72007307380868001978 -0.48149095112051915057 -2.6281919871585834869;1.6618203990217887966 -3.954097972686229312 1.8307231141856081802 -10.418543616051389122 3.1473763934104264095 -3.8815085641981794673 5.5093502329290808817 8.3873217168227007789 -1.5689986642532440797 5.1990965213250133203 4.3139409364553742421 -7.2852415006092181571 -5.2801021059135546665 -0.27168787668770238986 1.178438639206673022 2.6631585830691224537 8.6970940945998709282 -1.4458055596814121113 -0.89705809766500399505 -8.1060043097761020903;-1.6213311299463537551 -5.5972222676844172184 0.011697727590338168752 -6.0805384069465819863 -5.5729315116712667688 -4.2233037810798377265 1.7964079223236280036 -4.6796945638000302381 3.671825016720235535 -1.8776748164380239192 -1.2346469277917146989 -15.440210202001253492 2.6000244788974877785 -2.0335973550922554764 5.7309633457662050304 0.088071532530425877816 3.7999613176142021942 0.5082983062302839361 2.3681205648377101625 -6.0108857784185030226;-1.4105783751643243829 0.41970373890515833004 0.049870459936994210315 -0.39311432046513089533 -0.32468627949834116819 -0.46664639505733551683 -3.4942134584363935268 0.10669432403622550187 1.7796942517917009319 -0.56916848238383044301 -3.4929235459457776969 -0.51086785029488945842 -1.2929548080791894993 2.0169642207973734749 -0.44787606338273600048 -6.3305902984318613846 -1.966718325548014068 4.1928371042298540061 1.0568907431503253846 3.0488962593646826704;5.5953336140119409592 4.1146755999238635582 -0.54597142108860696741 -1.8465170044539775951 -2.4329121036068070971 -1.3226015714670669166 2.6631366019236235587 1.9155175466722618172 -3.5282799660610959513 -0.0047997229878981723106 2.3693338155364962461 7.068193459603649309 -0.54477834752353548886 4.3577120386777377092 -1.2854391366431294763 16.611825237791819632 0.85964410536642721361 -2.4693116489252289192 -5.4100599938397166966 0.074682561312849199409];

% Layer 3
b3 = -0.43154986017767466011;
LW3_2 = [0.04037179887717737975 -0.20515813024680285004 0.26127994977536767029 0.041881776411460069909 -0.084247642277779963282 -0.090617905302059292838 -0.15597387650791896818 -0.18565886111353371857 0.095578858995248469621 -0.12275938773453039243 -0.033000870422246829372 -0.29028618731655725371 -0.34776359247265703489 0.088954406796297388893 -0.12908112095210758086 0.16561941857779183263 -0.17740201795902926785 0.39247484060644893811 -0.31237110507929155556 -0.076137066686940688132];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.00454442172233583;
y1_step1.xoffset = 0.3;

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
