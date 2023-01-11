function [Y,Xf,Af] = ESPER_silicate_4_Other_3(X,~,~)
%ESPER_SILICATE_4_OTHER_3 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:39.
% 
% [Y] = ESPER_silicate_4_Other_3(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999998793562952;-77.961;0;29.715;-2.18999082656054;-0.9];
x1_step1.gain = [1.00003298043486;1.00000060321889;0.0138606862425759;0.000266418009857466;0.193068829037552;0.0582073296585987;0.0415843642790311];
x1_step1.ymin = -1;

% Layer 1
b1 = [1.7188902581100344857;-1.857862002201525975;1.6314051953863302202;-1.8717329485144602597;21.594034872866650687;2.1558870189538574813;1.2840430619155973613;0.12299013275542113044;2.3181895007487152149;1.0106914149247832047;-1.4566159530112723353;-1.7545787601664937227;-0.6020794602800014772;-1.4829675957169852207;0.9623605731089072135;0.78689710509263555149;2.2020158000424521205;-1.9483766219754812266;-0.92749647238175580988;3.3007191558400079678;-1.5872341680231318861;3.0205132082080115019;-4.6169926762660820074;-1.1606705047287746435;-0.20584291613848015867];
IW1_1 = [0.90298367518707423063 1.494555764935715958 1.3863691814250509449 0.28918929187725556096 -1.4275251623760354391 -0.37039157366063563126 0.37821977351258956324;-0.06980965015137920171 0.068696698104971171128 0.26435968123895026149 -0.035881512016477261351 -2.5771290530229729043 -1.5795908941126277636 -0.7291649293145801014;0.26363878154822295574 -0.26263839805123839755 -0.046007331576324222988 -0.25898946636112746278 -0.51923797275099081716 2.2817547259646615565 0.23829485300802344394;-0.11083990045406258429 0.14491552798376325373 0.17892719847345261996 0.099336801050016326764 -0.88119029624346389973 -2.3946230505415413425 -0.39840505378273949333;1.4579766609326982874 18.563107691429951984 -2.5240951004887861586 2.1727998011287601443 0.090306192351799124962 2.1908567465292088094 1.1254980204688800338;-0.42063849574819323074 0.095804062764579478451 0.70724932905835158348 0.59108517190229625804 4.2829021376172971003 1.6677982026813709115 0.30188505459027775712;-0.7117325885711638378 -2.2025997658035691984 2.8989070754436063737 1.1757637747476330503 -0.4475908056810606328 4.3406397429939032051 1.2029551810383816157;-0.86902320956556500242 -0.39016008926542872226 1.1214234563723242122 -0.064082967584188119115 -0.58983901957520790127 -0.1433313532693580028 0.20203108485357212931;0.7792558960084837727 -0.10932773312026118628 -1.0770326544179207318 -0.35391612791408133543 -1.1881898446741001596 1.1028165619963441468 0.13973458113665071512;-0.10597328141454510908 -0.43028043917919323036 1.3763560529209077909 0.52729309563500781355 4.3117645172324756686 0.14501322733139554444 -0.9196073148603058911;0.052031920179794771686 0.081149017569503073544 -0.33874305477634247818 -0.33898913792658796895 -3.9467143924934697807 -2.3102629584033977928 0.38639544005875037236;-0.19187709325514654135 0.62198243145988496572 0.23718881752596995827 0.21993269165686601951 -0.63095053879493578908 -1.9091803485580354227 -0.45605792647642090509;0.39668236787264776755 0.10819297360541385578 -0.53743071352883131642 0.66212716779449987126 -0.10181951251790519186 -1.5531625747279431682 -0.12102639296561532034;-0.22807520901279731973 0.5513407029735251097 -0.20913275474071613136 0.16615057954005990615 -0.19931057955800510206 -1.015329997017838215 -0.37325302541139898649;0.035381416723173818306 0.18632866051077601188 2.1027472700241003345 0.039158355289459254034 1.7644255095178749126 0.57739882788730467755 0.69377383363910694047;0.32570807133776014108 0.45565095300642916021 0.29670769940999852654 0.054159326907775820481 -0.31227161971684219699 -0.25962358306394783813 0.18654088591804529607;1.2182504321445519757 -0.15007254010894374718 -1.6913100780595591299 -0.40951386836319630946 -2.888151263622924958 0.90420509306276719563 0.24126887727296300601;0.057810371415629599468 -0.025850302612776415029 0.31186440431457512767 0.17921227619543109166 -0.93956351674618998615 -3.3194124821874924081 -0.7399768301049506336;-1.4758046477711583844 -0.75880285404714697961 3.2882896423638690031 -1.6758063037838817344 -2.4419413564670819916 -0.58896083322993841858 1.1752475114308202908;0.32794853938775903046 -0.33601333939691396813 0.29877852784317660007 1.4284411033641570832 3.8470578797677825733 0.97194541796956512858 -0.74117992065674209723;0.17527904215185483028 -0.026093867794946663807 0.05607321074114606646 0.25487614240196648785 -1.8508543129795445026 -2.757861107355225716 -0.90488346985497269959;0.36524788034741134668 0.28090620024997819959 1.7880320315939421949 0.068920035595845238752 1.3655996787267519199 1.0784335391766037837 0.53564773395389686517;0.094508616459352629002 -0.43215404134009322012 0.63692597245788618032 -3.3538354513030133042 -2.9176543081194341589 -0.53148186193176194525 -0.85243627871019311826;-0.3940142878398483206 -0.18579878603785349167 4.0712667707472816048 0.37780138308482108522 0.65366660640499718582 -1.614072296113521654 -0.76548142787226158834;0.24165446692521863614 0.0082899925999373319518 0.19887873290420932104 0.53025040816834256852 0.86960835975868422398 0.016564029806701848629 -0.88944352368493451166];

% Layer 2
b2 = [-8.2159524598013948804;-17.623133238519599786;-31.393884820702457006;1.9176517682250016072;-4.7450678477110024289;-16.917520827529916261;-0.44645324790119056413;-12.519845859341254268;1.5228360178205317954;-5.859576492916600543;-3.2257299727857544624;1.2024649076295714778;1.2056250030785482519;6.6085200106530406217;8.4146470176144205766];
LW2_1 = [-0.090121965679320611975 -7.1866261619435158892 -4.7786557591563871839 -6.7086532282115403092 0.24167267471064543538 1.5147775537369279064 -0.26765448607170089224 -2.9370768916103373591 3.8830305731346737019 0.97180412038753549719 0.0060859902195469443159 -0.50555495359617819151 -1.0178534638445251304 1.6878644362472536855 1.6971134576980475561 -1.9940060218720139229 -2.3748344273270469706 -0.7626773702673982358 0.48428732076012037666 -0.98746366375828675821 6.825244726149345631 0.1258762759632118533 -0.16683051703659104725 0.2363180319925574091 -4.1859269399889607044;1.7994909170614019622 1.4760874351182828867 -1.6037381617733734718 -6.0153104122908587215 8.8423729995862796471 -3.6115327126294891613 0.90871293263514918426 4.3361932902501543197 18.757099673286216301 1.7178553676629320091 0.21471557757905207486 9.6323387294501774392 0.036483064405200868463 -6.6177946948102217561 -2.8662283596469175784 -2.1163670860061434276 -7.7558456626151448532 -1.5660690598127200435 0.60817899474139147653 1.7188675022821151561 -1.6142356550478882404 1.9906979242356337423 -0.1501135542738521178 -1.5067879079642072515 3.7984950855172430195;-3.2569380454154326365 -8.1162706482798210317 3.9015301261846730974 3.590629668498364957 4.8196220030104406362 2.8429524846492855161 3.2525547078327718964 -1.8852590939617326349 -35.47928581918764479 11.051869713626276237 16.275718813332392898 17.69272270065624042 -9.4043949901140244663 -21.721349069253896857 -9.1908933270265986692 31.473633454521756647 28.596350976422296952 2.900604315967847846 -0.5756773655564288994 3.4466127201364158417 -2.7272506024856495443 -13.157320823704656831 1.6009013632593263043 -4.1778151820702831287 -7.023694916659056986;-0.052052352768183401954 -8.5933470627996957347 5.9037078251423533715 12.805406525108478988 -0.27919164578197053217 -1.697288300814883133 0.3549889701893145233 2.2587692741441203204 -0.40937498388257731463 -1.2558490627572955756 -1.0241232110086484486 3.9227791185908387206 -2.8032606775180468617 -3.7907842617845082245 -1.6591227744381535913 0.040170179506741829356 -1.2993614374801392142 -6.70298839008223446 -0.22567275798569630396 0.92720448350358652956 4.0155496931471548194 1.7266454571457272582 3.4283497317316262709 -0.062660626119850570959 2.8074843018337793943;-0.5392576202585609435 3.5495147290940649576 4.0121480453130731192 12.751887822470509448 -0.47182619815319815082 -0.66337085288812791539 -0.41889730751664788766 -0.19647352306221960161 -0.217548865761893917 -0.94640566775469014971 -0.82120758418741579554 -9.5971700335658631786 2.0014097590635828716 6.4081720443128844522 1.2635391942422458111 1.5620918003584456368 -0.3288386435111004813 5.2583099796191925535 -0.40921587270454612417 2.4691373454302545376 -9.7403042101202768066 5.6878607387989212896 -0.54461761289822940846 -0.30881830350500288906 -1.8178554771303898629;-0.52601585367953496153 -0.50047860215218398228 3.4859766668140750312 15.694807471892119111 0.035609519439640385796 -1.0346555458089858881 0.13750296168597470059 2.4054995329183284092 23.074559508232322713 1.919342803650102347 0.87201890019814920496 0.21257459822780941372 -2.2029596711420804667 -3.3635145729951694271 -0.9032632672149905062 5.6574834952833779056 -7.7985525637098955798 -6.9701358714710694997 -1.223479147054660654 2.571754268044780023 -0.30563801761778786048 0.13008900845075521691 0.35651181912150664566 -0.99280739641791770911 0.68406530192901759957;-6.2754311057340395053 6.1815134087495273718 -6.0406618220599570535 5.1919788107092559315 7.0379054768147515375 -2.3164989624613965269 -0.19775739760638366693 3.2775070713340732631 -0.49154052483843380106 6.8749910254115773611 -4.1990438579790048124 -6.4703381818974357031 3.5414717192382112465 1.642452208594194607 -0.38891596310397968939 14.808185205766907444 4.8702059784225983918 -10.307860774197830978 -3.2508473541323326117 -2.7909710160482643637 -0.20137372216022925153 6.3684800413370012961 0.817571878346798675 -0.42243228494493206826 6.6114320954246066364;1.0991288664599996316 -7.005898268585407429 8.0156577149230532342 33.268093840894216839 10.018541760665415907 -7.0359443658730373983 1.8198992586054283382 -1.1047477680808799327 -4.5563697672740257616 13.624015756878169014 -5.4695034584736719196 7.7448102621193291029 -2.531523540699480268 -11.147172978995836701 -5.6454292322733641285 11.701912643093162814 9.4586217760162014656 -21.707773581006534869 -1.6829817551254333541 -4.1888380772358431159 4.1871002331040800826 -3.0257292070586863098 -1.5570199006986884527 -0.75659575898812370109 -9.848772416900226645;-2.0894901068312781156 -1.2363061194182742852 -3.0104949096414514997 -7.2106115113478637113 -3.3310971543387037919 1.4795835648748063562 -0.10792403147555321152 -4.0306172762564482781 8.8050585031634156508 1.3781149234389098446 0.43714324616029731629 -3.5591189329751968806 2.5146130312402128482 7.6401482177895809755 2.3317168298210457777 1.3280933301686470571 -3.1926088017887521708 0.18537413014474582451 1.9442075823509064048 -3.5174247907191493745 3.1370870295506483316 -0.62476904686648204734 0.77314279482248671194 0.71674022878985055307 -3.8405732727447130515;-0.38911327276514157125 0.2646727861030908624 1.126505128211530371 10.731646751621561009 8.0088131612502646561 -0.73651680081652803889 -0.20359920792660071709 -1.1909954151788832988 2.2398592923644256025 1.4331105377613078122 -5.8776767898257924472 -1.662638956459550954 -1.632522037439194218 -0.46071110709530410654 0.29202912014079945857 1.6083185244094035937 0.33867596146679557201 -11.747164443896368269 0.54955303105870489944 0.57246950133137253225 5.7631755275149751938 0.95348791433129664874 0.62015093946058319485 -0.43809584119329170182 -3.4998343986419908624;-0.0025215490990342197264 3.4411296374380722085 -8.8596260908045554316 -17.881282211399465609 -2.0740507621428534435 0.81360925941013106932 0.048445546857278420916 -1.6034332007788087981 0.97828542506805293844 0.71735022074506138434 0.68303835013130709264 1.398670732980807685 0.6974409091752387102 -1.6336090756266776758 0.36237347612403741914 -0.68790673227321263994 0.44327660700175125497 5.9198440278833714956 0.23635701523610647135 -0.89164995261096491319 -1.2653801484986351156 0.84493859551034189437 -1.4067778618439541471 0.51029714268162063018 -1.9792567719707057083;7.8625714015126533241 -6.3510732621640126894 1.5991860117424774312 -4.585214039262271335 4.4313346049168824692 -2.2649451421632704751 -0.68953641514138341861 1.8122127980886293397 12.586581897344295555 0.1260686078991355985 -4.3275114789563859929 15.864147116993255082 9.3033089772817323393 -22.735022044321489432 -2.443279160084713908 -16.73662228383900441 -7.5457868605522042671 -12.49674455492968228 3.4539811422012700604 -2.8949010684771825375 11.745413222067140779 -2.8153262317356468891 4.1705165018236929342 0.21794669696850074714 2.8107102521789015803;-0.30604586394194654986 -0.54060765597172610875 -2.0109052124082791302 1.9293116216991348377 -1.0084461047546016133 -1.7020274128640440203 0.31958199290297445438 2.4204741337347175367 -2.3899270046414273772 -0.1894072342737415926 -0.049521169295859319981 0.58950679877579736932 -0.53621836126608379836 -1.3611764295887049947 -1.6066095147312515756 2.0526849856721929299 1.7876870474988908466 -2.1910698049754468997 -0.85568251022041108556 1.8299771893835450332 -0.85125528373938241788 1.7102390846674850877 0.44317252920367744551 -0.25993089311094563776 2.7474122279605035502;-1.2026744581945569834 7.1630983303627022707 -5.7609842175468628511 -10.769320528349522803 1.6499678683471370455 -0.64419640929995958256 -0.81742897667334790679 -6.8051742399967389829 -13.488618980867364883 -0.85784761439136270056 -0.077316417716485966172 -9.5471610730194367278 0.66460386638512713997 11.658090696169912448 5.3246867046617545327 -4.1785824838485510213 4.1535413906124425409 5.6566661361953887166 1.0440628401090588095 -0.92931772648038046825 -1.8491412584409439024 1.3519918825193091294 -6.4341274599936557621 1.3250507364987809122 -2.7269488715462149031;0.23362871462963774594 -1.1702865356855729484 4.2683697387141670276 -2.6900002698386655808 -2.5238063190621287113 2.2542440476089180557 -0.42632069093106134128 -3.8767064781173914412 -3.0690333624609551499 -0.24507318571350933389 0.42049674790411872838 -2.423275056498881419 0.39393051159670156025 5.5583501736143299254 1.9532038746530342443 -3.1392384395735710712 -1.1528374584628582333 2.626085308462789758 1.729686409386196333 -2.6052281890543595999 2.8123046103536877283 -2.0294391579345534282 -0.64697482571489950409 0.64803399954767115787 -2.5109467273456536063];

% Layer 3
b3 = 0.044078786041187541234;
LW3_2 = [-0.35662578791693944069 0.16623926410038747381 0.25594679910191509009 0.23118087470540768513 0.14757431150923308594 -0.26935659330808003897 0.3301015299910391998 -0.032802081550327814863 0.084761319279585162212 -0.29577866631880861403 0.53661917401406789008 -0.099234460417189265802 -0.46367843464127861397 0.21063133528720720533 -0.28356126600147013583];

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
