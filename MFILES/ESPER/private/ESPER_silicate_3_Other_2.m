function [Y,Xf,Af] = ESPER_silicate_3_Other_2(X,~,~)
%ESPER_SILICATE_3_OTHER_2 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:39.
% 
% [Y] = ESPER_silicate_3_Other_2(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999998793562952;-77.961;0;29.715;-2.18999082656054;-138.904784718693;-0.9];
x1_step1.gain = [1.00003298043486;1.00000060321889;0.0138606862425759;0.000266418009857466;0.193068829037552;0.0582073296585987;0.0042772329208737;0.0415843642790311];
x1_step1.ymin = -1;

% Layer 1
b1 = [-3.2452434717086959814;3.901488021786686744;-0.48628768203414474147;3.1183810558867164353;-0.40292663425890379303;-1.2235085779438608089;-0.93258238029989792839;-3.3903125164514027468;-0.57512368283635406652;-1.2299430244592579609;0.39935548294103989653;0.46682407745618953454;-2.1685991575785452667;-1.5294322237096051875;-1.0637292525582340286;2.129609379978921968;-2.5954598066852740956;1.6578641784506318313;6.0334458413579428182;2.5030040217239779565];
IW1_1 = [0.13962214462746430632 0.064380097621532045293 -0.24308547266493985006 -1.0495165168461133565 -2.2587312551887803735 -1.2876517547194619162 1.2828537308887781876 -0.3130999655580181118;-0.31288873201739814123 0.06089131643291238738 1.1877009179738990952 0.30432328216035714474 4.393775654818012022 2.7561054449789748055 -0.49239172725098112604 -0.33931962452061159574;0.034001263135755195599 0.33591323285841306889 1.6204090762214700128 -1.9444871275614588058 -0.36829170213858181127 -0.30676076047906941646 -0.61208181815715856366 -0.26085379231849986637;-0.011883255454123824169 0.043054607773231549916 -0.072125805457418701305 -0.17046472874394721919 -3.3953849246170628717 3.7032449447146102983 -1.1590579816345185638 0.086765155106147476283;0.86192346285175203136 0.12379823676172513336 -4.0623087071428773953 -1.0564150215582839287 3.8401024744588738713 3.3259742458903276408 2.1691985534471145414 1.1843182689732780766;0.12598074010250095989 1.1074302741245540815 -0.63388826865990643977 0.52006693346063970207 -1.6320725823848822866 -0.45202819809752203817 -0.18068615160080181425 -0.61377438031698772214;-0.0098550027620704847009 0.086061118332758645977 0.012688598895918360918 -0.79908150664830757126 1.1818184499154262834 -0.1775837950675688004 0.82951576702940488239 0.21594800171801337152;0.44333671008852898465 0.33470178430602254283 -2.3853074847125217417 0.42866191441240453219 -0.4507754928141653461 -1.4698880201071966578 1.8862773281818541182 -0.54600318628901156615;1.8102693538790706196 1.5515671837145206258 -1.985848434153065778 1.05992072171376428 4.730877889816420101 -2.0812263356330782393 0.25886237048701338237 -1.5947274400050801546;-0.10325631646084273974 -0.026800655711492222394 0.27815396248554946634 -0.11081821157626214891 2.8332318624766812398 0.27329861757998402672 0.80047569355841474703 0.48688087316190387099;0.46435774399615503683 1.0112440362683463579 -1.8025399691611907471 -0.28391079928108559072 -0.80756778215577762481 -1.1166469147585618149 -0.27979586545837398237 -0.59249805079466855151;0.0030982230343697811023 0.10068296209588760426 0.68672618724175904514 1.1415839556610838379 -0.075106121344454157551 -0.25734476548438345311 -0.70746553874271789208 1.1541916654918944474;-0.07855062286245398373 0.024707876075409539868 -0.11414700787381469593 0.055699724837422831536 3.2060032701833689472 -3.4791135728383517645 0.9597574325459005129 -0.76667182383602261453;-0.77927889576412867623 -0.39713834836116068683 2.1541506528184828362 -0.8579021937844806267 -3.2264062715554193694 -0.28493793529731181735 -0.52367824132216156396 -0.14619970861895234537;-0.0084235571841447046126 -0.74661863472345046944 1.1609407557492483765 -0.38614437943541568554 1.3997660531966888975 0.54548179002130592341 1.3446583336574380318 0.29492660732630149312;1.2236667991300733682 0.44862033099525255864 0.8649577785071772551 -0.33347001902080419811 0.87671867293517680153 1.9671800934285068596 1.6150181493829449675 0.4874244829419114855;0.41925834156730190649 0.51402130720075811521 -2.2529802920829880009 -0.0036797205469089495911 3.4057624240588078024 -2.5455411807287005033 0.59073983949587560716 -0.7788115713742905788;0.65954242171927679461 0.29256971005621829063 -0.24894355985961882927 0.17774482767366078795 2.7535740459088762222 0.4619140104529737556 0.22455656773528026582 -0.92386328313548626845;-0.18238884533910504349 -0.20181083651164197224 3.0507427942472458149 0.067627928249940988392 0.49436847492102325985 4.744412192288907093 -0.53411289432591957382 1.1363684420625288496;1.7693274475143156987 -2.0247974802822104046 -0.42162387885954760725 -0.46163512780674315072 -0.2898119397310375267 0.82974589726763048336 -0.28833003419845426585 0.29605672833118451548];

% Layer 2
b2 = [3.7010348036647759962;-8.1615727616230184793;-4.1426592161236106548;-1.7997433036788514915;2.62837659038016902;4.449750182504359941;-5.4704601615634214085;-6.4529490326667735545;-1.4252976310962928164;-1.8652493224338593247;1.2182211634000210676;5.2699456027223945398;-5.1100558087976901689;2.0318719512061518273;-4.0157890479122064775;-3.6594454694127653838;3.8480036683465757008;-2.726461948214381259;3.1310814326753537706;-1.4281620366984368875];
LW2_1 = [1.980194679800023394 0.19901466129766867685 -0.44517000803461403979 4.2222500084988423552 0.61307764781096774076 1.9298694535286535778 0.006707313677714148159 0.062804723332986756024 0.74880889335785971728 1.0242011940689925709 0.2788422419594145385 0.21023774109516951691 4.8725818852427140726 1.9421253173461376917 -0.57399554470021130204 -0.52948029259511841804 -2.7115608423305310737 1.6832491178087192196 0.78361850175550173958 1.1416150199294363698;7.9951718935869395821 5.3777790495912318747 1.5651495587632771045 -7.1645479121166486536 7.9582798914900472198 1.9881269520611453139 -4.5358991331681037806 -4.3307526442509880482 -1.4189768762302321647 -3.8371386224877954696 -0.64121544362732840749 6.0338620363375223832 4.9055010566848054054 -3.1462427926201224082 4.299798014563154247 6.1611347591021283776 0.58151587827452710933 -4.5568286022607793484 -4.9146158996306281708 -2.9690257876004575444;0.1190754954364850321 1.9151114087661746943 -0.96121681192969898877 -1.2173098418669030529 -0.44810992691297796275 -0.3450991092624736134 1.5652453552393903458 -1.4033458528832289947 -0.38452416054097904308 -3.4925071298697094591 0.10460722207437676834 0.31805670058719520688 -0.21453870797410146132 0.84490840548565349089 -0.27846533048969734914 0.41132917064579754829 -0.7736252429978371925 0.43613160987026389215 -4.3443012420793021278 0.00041900248010238788296;-0.35433075268145025616 1.8803203644495223301 -0.67125072175044253608 -1.3544671015416400728 -0.5459861625275939101 -0.37805988336260476501 0.65549116429498954783 -0.85195207637647996091 -0.043845471240930307399 -0.85534652285934953397 0.34749124066745051831 -0.39882899686447625953 -0.035303217818967898289 0.63817664844925392487 -0.049244391541421919656 0.53763722985497419682 -0.56652114804566577888 -0.71869863185096183322 -3.8650604209651286958 -0.071533856313661214887;1.2089228546454244828 -1.5855715127584601998 -0.25041725243543355717 -2.0869573939097918469 0.11719870843180645037 0.24327849591963454645 0.045895521535501360155 -2.1662498493767365915 0.96731698783452824131 0.49314315528015018053 0.14479836566884748961 -0.32815369429608898244 -1.906778222420043134 0.22156381361790519757 -0.47851734293776476692 -0.89740128619596326587 -0.72418197860701927215 0.26295151116062254193 -0.99165196822967827117 -0.32622124023823306782;-1.0064231568901904623 -1.2760732249075268374 0.45552217050671206522 -3.1701781888882325511 -0.26183057555515942827 0.069296495881029263053 -1.6202360693829747085 -1.5617676704430321077 0.0041395744747870898833 4.7606228366987899392 1.4123341150713675241 -0.47957971012171790592 -3.0861625680833069119 1.0709716823254553653 1.6416757211747070588 0.31587908912991913279 0.14737250023595904747 0.92604486141901043439 -0.32431155796290545013 -0.048074219028473801074;-2.0076637387885734043 -0.029939067841280104371 0.042925833521604589427 -3.4190156791560069038 -1.0048255248007780249 -2.0104965783732597906 0.83738462336225572002 0.26466897025528557386 -0.3797666111538809397 -2.8835193226360438246 0.29540209327252597848 0.14556332182142900367 -2.5486800635210591359 -0.99430665695934961068 0.52851729624495513704 0.84016259445315699672 2.2602105762847779502 -2.0607803491693452713 -0.93438027517144717216 -1.1734342339290884105;-1.9978090730739015424 -2.1008568103132065907 0.24783282994317915038 8.3282682033816080036 4.7991297174702323147 4.5990618003734500263 4.8473327864574082824 6.4462135506186815093 6.2022820213256988353 -0.58742039805053081203 -3.8706724154709779562 3.8113962780877446157 0.14011903043453546869 -0.36214199466782293069 0.19955387861064569077 -8.5794869526108978874 -1.4465609413685862616 -1.1490019739565173751 0.98437046126345362218 3.550056161007049127;1.3493001714467618424 -1.0738074542274658185 0.0947601784919273743 -0.21655734116749272844 0.19044276085534825316 0.074385358345390450996 0.22875168997812458938 -1.8514026357941777423 -0.69388235163029188257 -0.69895527297483439622 -0.77619149733914172273 1.0245534066839894205 1.8377345944712657477 0.20398795334637165322 -0.78869813640132069299 0.26200975135844639663 0.087613941475641290979 0.36281861349755256674 -0.32889684838350974339 0.21388360326298033742;0.71050531058508359372 0.84067385442010778007 0.072621904489141284045 -0.60357399639992959095 -0.24214950279481467499 0.0074721461954712232317 1.4126964388631997327 2.5312331557082554667 -0.87955257550278609369 -2.1314620313596375212 -0.98221034059730483179 0.47538832459863267221 -0.85501297379581708835 -1.426141202899235072 -0.45129550474587520892 0.41175506165044362117 -0.63133781057616145116 0.15825041043124449258 1.0770091804138923752 -0.66302983078608979106;1.7610321511963038521 -1.2129439947366269514 0.048593389348548281237 1.4390604434898746078 -0.086956181474596147551 0.17910222968979097602 -1.3069764975168689514 -1.5260448769474139752 -0.066958413775464076245 -0.075143043846063004021 -0.29220930399410061096 -0.1863868688511140792 -0.11360156896433269702 -0.74760769090940248915 0.072986771129443683837 -0.11802643753601479992 0.42482930369774096757 -0.20628450951903709587 -0.20713686186319682503 -0.30613502018217603196;5.6753129052639224028 0.29008663136635859381 -1.6316972964763041265 8.5554582759437298733 0.25879934856654895903 -0.79743745520131259497 2.1698861674192189142 2.0533365671656813589 2.227773040743183941 1.641530426158112288 -1.5692318969415635088 1.0919412082081345616 2.3845439299409383338 1.2202462326515748359 -0.31707553545524308491 -1.0106758357166956142 -3.1135235262774569875 2.7726705800081861142 -3.8750033669410144377 -0.30659123498702001154;-6.6266041459518731926 -1.7115693186066671672 0.44667582539756944815 -5.1673255076689050824 -0.0095457309728619321731 1.9000131532941537671 0.82707834645229172388 -1.8765667147942057813 0.37249200695603162936 -2.6972096486430743312 0.0016206285772366185727 -1.4523541516275231267 -2.2998005099799976314 -0.91291885873908862337 2.2948207439909085181 -0.28267294598009912177 0.26305648635301026594 -1.1632275594179131151 1.7976932650847534667 -0.040829621926117852515;0.10287246222207394042 0.17211498123661533866 0.033071866227957762152 1.358293358911244253 -0.035644948523021198483 -0.019361724135260518359 -0.51034135643973599361 1.2013178803663606686 0.49490712071599257671 2.0904383604507938976 0.075577704557030556121 -0.34694721685011908896 0.64129599659840053061 -0.078146090596456302402 0.13467321869309384352 -0.15767479971170025865 -0.4613129192229247777 0.01047248603463788591 0.65031928107173542219 0.015436744269476745245;-3.9092254869930411765 3.5391179160510750634 -2.2373628864846346609 3.0715507126823227146 -1.4659342521170524165 3.5011133443556983202 10.118996407760047163 -0.19289777502672769693 3.151619668447292355 -4.0183732002499983338 -2.6736493715459737075 -0.89513287448820244574 4.4357329007969434898 1.7293925240669987264 -7.5735544224490629617 1.1081744864919058635 -7.7918079896218346292 -4.4679370365930202169 -1.6513851484322759067 6.6292950466189939362;-3.5830027647563951199 -1.9333687174953335663 0.093513307296623757181 0.16262174276666477302 -0.41807438117414952305 0.65826331993118636365 1.8090594314741870186 -0.33395271412671190259 0.37656289783102547819 -0.52243097713242392377 -0.50721162388386986652 -1.2317357958016819097 -1.7469348907584152375 -0.75867371384051285332 -0.84905690748544226931 0.010283396880168986509 -0.60850088968025917158 -0.96888936521477098118 0.3310820456967255021 -0.40839458544669476892;3.9357247769645273827 -3.4498669827500112461 2.1995264626985333578 -2.8938604583534988102 1.4375986882216302831 -3.4923514682549683386 -9.8099321924791471616 0.23376172600304490068 -3.108043439004311459 3.4629803001419525721 3.1260362532896994736 1.0311038839686488267 -4.3957978305085090298 -1.5825906849585593683 8.1088407768261490105 -1.1378343325539359565 7.6070635265904895661 4.5285893194735713507 1.5771737207697198446 -6.5506186007720028641;-6.0071645642030357948 -0.73948444778491140283 -0.20637572069627577176 -0.57328953381111247278 0.26933578138295705129 -0.32192884450183062439 2.185329744731668189 -0.12268934372314739545 -0.21814095824807222535 0.67511319721464446708 0.95658045680789505205 0.37402316362889653778 -0.76449550798477128311 1.3973702222313917343 -0.28368479514688477172 -0.13402581026273693854 -0.21691295815487160459 -0.56381507783086692864 -1.5169454473567671737 0.056956393706439290003;1.2191672905149075756 0.17655357557484530062 -0.30186843080769570902 2.3505052975789380021 0.44736747902296797319 0.34569911614233195252 -1.4050800942377810188 0.81159540368329963211 0.95554186835167653769 3.67597762897513336 0.68748407248006515591 -0.15338259803201453879 2.186139759055884646 2.0144510465189875426 -0.36774192728414190068 -0.5400884420842152478 0.082095298821411047396 0.8969827254838148578 1.7240273797399872091 1.0906984372990859988;-6.0209563350855264119 -1.2161759130379601324 0.69742168502402002161 0.15261720307618276138 -0.53607100490081682764 -0.21034204839933617892 0.70002249778712233308 3.2332431529003926585 0.98367471697334707414 2.1894845880475508615 0.38479800168657696258 -0.94362849546417060242 1.2922709019068767677 -0.51659480661102141763 -0.21862019759571588162 0.12851676353469171366 -1.0738528095994035461 -1.4343502869471225214 2.4804301213825228345 -0.013491792878401982075];

% Layer 3
b3 = 2.6609226826784531106;
LW3_2 = [-4.9277860277664480293 0.25810844833061546977 -1.1185081669381680136 1.0326898094606278278 0.1756114301840969727 -0.10422976625008285867 -2.7261210397547772644 0.16751661734500553069 -0.32004531135748753856 0.44265812875877663668 -0.40485698110908868719 0.030000406154552208132 0.061358807666467805875 -1.1889798621639098286 6.145240580880426684 -0.18587890964429704121 6.2101692289300780558 -0.12580619936961348615 0.3270895034122228151 0.1915851596849624805];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.00819772922900357;
y1_step1.xoffset = -2.02;

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
