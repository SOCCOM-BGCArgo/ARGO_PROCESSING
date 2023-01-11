function [Y,Xf,Af] = ESPER_talk_11_Atl_3(X,~,~)
%ESPER_TALK_11_ATL_3 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 12-Jan-2021 17:03:49.
% 
% [Y] = ESPER_talk_11_Atl_3(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999999390765;-40;0;11.29;-143.241523403244;-0.12];
x1_step1.gain = [1.0000000017546;1.00031699124092;0.0153846153846154;0.000310077519379845;0.0705602618858232;0.00458552991639557;0.0145655815308426];
x1_step1.ymin = -1;

% Layer 1
b1 = [-3.1778808064454687;-1.9894688830310374;-2.0590235503452523;0.2480764091507747;1.4615910883288823;1.2235098210410871;-2.1564409282616057;1.5281208061034335;-0.49900554076087522;0.6485087037104208;1.6140787927902316;-3.244441144632662;-4.1768608840990247;-0.10291950426183429;1.9675946306141652;0.22188927997323732;-0.6977745579361665;-1.5778481597918936;-2.0908606480193992;-2.913547446348562;-0.83130240172847225;0.68418727549578884;-3.5514336594440064;3.2893803999086813;-0.56878508419757146];
IW1_1 = [1.6292904013262994 0.70698200753031415 -0.13599214756784248 0.20662224667751275 1.4066834483443775 -1.125150794152816 0.70796321108239224;4.005534969236809 2.2477470407982145 0.51397266097499905 -0.28205066987886379 -1.3138881025879348 -2.1632725125232368 2.2067051518464225;1.0296531900958861 -1.9427970109547748 0.31748189884056816 -0.21639747995028721 -0.90298367790498335 -1.2984974194623951 0.64564827553048287;-1.565514004210546 1.3271580612051264 0.85351027232502397 -0.38376068281939901 2.812599740716089 -0.50109300640982257 -0.20945693998694495;0.066612946614965265 0.39951781026299144 -0.050994530581660619 -0.17796354328685679 -1.7115315141125498 0.44419683201269444 0.49340510611866062;0.17901482385938622 0.79665681170576752 -0.79402152088188294 0.05349057866524666 -0.46862949712832741 0.23122875497171472 0.01624646231385498;1.5814497316004481 -2.9043570285015061 2.052624882210778 0.12877564893493817 2.3610452993150828 -1.3135609102398831 -0.82481512786245492;-1.3898646221681392 1.5511251793355547 -0.6518659903092624 0.3266791500306287 -0.92056869808231734 0.458007459033555 0.54755816985783734;-0.71359845896197049 0.67328031403315536 -0.33574611237012297 0.18211913569371263 2.0553450102740944 0.1438170872647086 0.98364022944666851;1.0884610145764346 -1.5865504583063539 0.83045526629217836 0.60098236394100457 -3.7570556947688147 1.3383516001126412 -1.7192643819302094;-0.48803079874552219 -0.2460106320230053 -2.2093873574707947 0.60366425213487662 1.4575525678001646 -0.031324299960061142 0.94760478517969537;0.2246680495030263 -0.015552354347103317 -0.18656795455058003 0.33892809993344603 2.8841895613296504 0.17968060938208427 -0.49309932740641876;0.28096413689646449 0.20555405137920024 -0.86221619825185258 0.54331132291605788 4.1713460974568708 0.59273694519813624 -0.38419980610508531;-0.32163244432393273 -0.079996752675074084 0.19234419916747494 -0.038245766502682042 0.53707061103555298 -0.33203204870736425 0.85268992074070848;-0.027578541896256754 -0.49750627029310079 -3.0213439276518788 0.43763788393961067 0.91677236019342367 -0.22268740007664808 0.43437040745214728;-0.34854779316596696 0.92978765353294424 -0.63400561427943547 0.19209893787725868 -0.19712273892140453 0.6993191212250105 -0.20360372026979512;-0.90811109944045909 0.95041209058610843 -1.1206074938744213 0.23629726207672813 2.1207401254555651 1.1660411726106619 -0.14967406566219435;-0.20866669615503186 -0.56311276587027881 0.72615771913679206 -0.12363100622623772 -0.85556295421646811 -0.30013659875287191 -1.1830344873018108;-0.45574012168447386 0.045838586091579321 0.67853408258428127 -1.174350351795034 6.4347936020402567 0.067055421994326569 0.3723844664998493;-0.50037707391628805 -1.2642588730878168 1.3920791789718945 -0.4229192366208791 0.011770730386882138 -0.91208847374857516 -1.1775246416261365;-0.19950690374787944 -1.7251034061652517 3.6643367868411096 0.38192043335426773 -1.4300214734662686 -0.41535068442913037 -1.341139630593903;0.79110145000101706 -0.50985436339876167 0.43623598732392921 -0.014459287622457552 -1.5914017078319453 0.80422621050799559 -0.90662339997775521;-0.68579233060736366 0.058855530074537893 0.88295035942903866 -0.13310826199112194 1.4966291324658136 -0.82849677112419495 -1.0479852205605622;1.2850355896918531 0.27735917384380093 0.32863740074942149 0.52474635596654484 0.66207612169766739 0.44405121222427912 2.1295787615559512;-0.26226303013195967 -2.9285035029617434 -1.508021323257914 -1.725624553220886 1.4326128606981496 0.19895537414976489 -4.1202154024638267];

% Layer 2
b2 = [-4.0904541066325653;4.7983906000695606;-3.6749874323039564;-7.4779414870006198;6.8589565899923501;-1.259529027648661;5.6566519674558009;1.43613357081481;-10.672354269105943;-2.0039769326197496;4.8897650277507259;-1.6666920136045;4.3047621379392584;-0.45812582604414392;-1.8496292443865361];
LW2_1 = [0.90281401383934023 2.4062127461170144 -1.2338924236807871 -1.0917495100959931 5.3657316628092024 -8.9570894500593177 -0.39220219859288541 2.0241901434268406 -0.87249353718269174 -4.3206613620324283 -2.5132557147043397 4.3787841710814419 -1.4921064037225997 -3.1187403972891579 3.0134924381467538 -6.8808931499034616 -1.2886533513655374 1.6133765270186669 4.5793220905189624 -3.1828881865131882 -1.3388980635156651 -4.1753730773664861 -1.1820521807652151 4.85382372081689 0.91006674266032817;5.2176264743280454 0.24771385956923236 0.71937994222532786 1.7702330829473163 1.5691551262309438 -1.0845597534811198 -0.27854630391273655 0.17822828866397256 3.4468274896474793 0.087826894297189484 0.19106632004779084 -1.3353243232209926 -0.97858957969842186 -2.2867210358723078 0.38259150947828263 -0.92705723094818082 -1.6134259321348214 1.0126261006982824 -2.7787471789840157 -0.81284712951772908 -0.76788306112754545 1.2765283094408175 0.46653610080411678 -1.261262346516983 1.1715171741238921;-2.2837361985218836 1.4974094656193229 -1.0196463590829348 -0.24781987801493258 0.94412640199910591 -0.61703809624821038 0.60650795281768954 1.4599362165058893 -0.77095019426772293 -0.23991709734496869 -1.4182507201142853 0.84143062595186613 1.7254928116272223 0.7744983402732013 2.6700454598649945 -1.0857213349513986 -0.55521822371711094 -2.2761222659137426 -0.13314379222676212 0.72723955854144273 0.079314485262509046 0.36942323840945279 -1.5116486660011883 -1.6863907003163789 0.29224962054448883;-1.7515366573651767 1.4751153008321398 0.25777217700393734 5.6026306542597339 -3.7773873835756722 0.23131651501180461 1.1915486793473427 1.7615901433180234 1.4474414747760458 -4.1199798033478388 -4.6686419021541639 3.1570812682206406 -2.021828704336341 2.6889436317684576 5.6389672770783195 -1.8611478152920549 0.5890275619131532 -3.4969818776061423 2.48453676052804 -2.4924399630211327 1.0712856735286351 5.8631513044660215 -0.5109039075729116 -4.0867013330473609 -1.007942221142432;-3.0107643967022324 1.9203761244542583 1.2909033212857468 -1.7770770581154711 -3.3938352991325167 3.0070440732888502 0.91834150159961658 0.77996838579457228 -0.70018082725859621 0.20739199317313645 -1.5532093350090153 -3.8877118159501558 1.6349263082338839 3.6419769655007217 -2.5777253721916176 1.3279417899768109 1.745087175873127 5.000877529530019 -1.9557509356962115 -3.0162336263623244 1.0346245803216065 -0.45818273821075528 1.0103773423814975 -0.24484430278549696 1.0404948175724384;-0.35184830424324925 -1.2666662039403715 4.9995154653767262 4.2919865508335056 -5.3080713843178868 -2.0044131177274949 1.888974189809091 -5.2868790708288005 0.7549224479784501 -2.3997702294154659 -2.6501123912339408 -1.3222718309245896 4.7319989002696543 8.1081102189084593 -2.5457557340640165 3.4526286729356239 0.15468499086413132 -0.19602818057012542 -0.11158442332023977 -1.392110767571503 1.0313492547004284 -1.5636210612891099 0.24410765512965707 -3.1820640216264642 3.8071005482487275;5.5121718501527077 -0.16114050905159991 3.7724788887131764 -3.2302873080004719 -0.29014141023813866 -1.151465746695669 1.264427555816181 1.7891654058491937 0.73194134853172665 1.2670951324471436 -1.9885647488071094 -1.1335427948122367 -3.7906920843395522 -3.1587130021542862 1.5528337470074529 -1.8716638413236866 1.3383119151488037 -3.0651063594338668 -1.2765789204503066 1.3743707925497166 -0.99125868015220575 -1.378187831073922 1.3474852369950898 1.4607043212244515 4.5885335172964616;5.8379197634271547 -3.3548430959019813 -0.87202269221370798 4.8721368241854055 1.8558702593779339 4.0204143811897231 -1.5594758590845699 0.28343829110911434 -0.2290639533369451 0.93815738137460203 2.8559332391998016 0.52988886799668222 1.8678274663578103 -1.7839386828287938 -1.4064978211128383 -7.9762219620795056 2.1771141791321966 4.669687171608496 1.0377458997637978 -2.1401264332720604 0.15803813384914983 -1.0865116324718642 4.1674968841830173 6.5863886381975139 -0.19382939726820134;13.426465016121133 -2.8089356150714111 5.4450826275893185 -3.3035164309379992 -3.6228608879882116 1.068952666335329 -2.30636142760385 -3.5397762034910505 14.671838750070753 0.83398888293441853 6.0526695635467336 11.421966827239528 -1.9047022566468561 -24.778005754681384 -10.959429734261334 -10.584814101517251 4.9773168254014939 -10.978069535577903 14.049556207166907 10.11545660967554 5.8719844665800505 -7.7414155294593376 -12.132366539270643 11.260782002305808 -2.0568324516340906;-0.98108860436359002 5.6197760362484388 0.16666237582321133 0.96087561338968819 2.5908724703587884 -0.47081446714823605 -3.8550230247835424 -2.1558066275943397 1.2877010730603737 -4.985973762792967 3.0576663891302513 2.6539164145710235 1.9968554727492966 0.64020846255679464 0.66012640733537276 0.80353139125209805 1.8063956558453098 -3.8518721445696351 -1.1437703006961004 -1.2748662913651996 0.10277188733793606 0.15737205714606481 2.6192139072100824 1.1742358995482702 1.6635394753451107;2.8463946481416347 5.9555937716472531 1.1605776905246601 -0.020475207627937728 7.362050006166732 0.63459637319484641 4.2654896160381668 -7.2469429447715878 6.670524720224047 1.5080799636454019 3.7434866337819459 9.3346748384157383 -1.9457061272844316 -4.1261496031244782 -7.4888473487967886 -10.642822291907795 1.9712371864465159 -0.46971091476429189 -1.7074918648961839 -2.9378745212197868 4.4421998478983911 0.9600786239669491 3.6438268377282976 7.3720325529547779 -0.34919806690469646;-0.42253378056534924 0.23807179473736736 0.40217834472647856 0.68297686183769579 0.55756928939946049 0.82450022821981139 -0.32424307555165194 -0.98139144613059126 1.6874770327871564 -0.24438758443741945 0.18423459351017601 3.8541585658062081 -1.0450329528737103 0.45858608896455888 -0.31906530451145376 -0.54673654659368121 -0.33690934614175388 -0.059797303959504867 3.5049269691207323 0.54557510948332288 -0.13296877022830195 2.2311499372868338 0.53979162195794272 -1.5640066517505895 0.58034825816794344;-0.67162253064129507 -2.0760832883355809 -3.7681562984271806 4.5580675488597375 3.1719118737193845 2.5955029898250785 1.3616088499682053 6.3764499202670883 -0.91849172548442515 -0.16229462009726117 4.6203836108855825 -0.81570680672051898 1.7055952598535813 -0.56005587375327914 0.029009772227819211 5.615502319221056 -5.3123690012120273 2.3286936019098734 -4.4498249856314356 1.4808436212987039 -3.15184768866358 -1.6123493369004596 3.5462719456011365 -1.34923600794405 3.2145966499054848;-3.4023393575618295 0.12567172046799949 -1.46045740280428 3.0837541398040984 -2.0576458066896595 2.4365725496885351 -2.7620848983788489 0.24463680230271015 2.7088565931764892 -1.135637192612899 2.226293353231573 1.7715289821065971 5.6259350485513115 6.0660527084128573 1.6589704016658802 -0.95830172665412872 -0.61462064434392238 3.1659643896524678 3.986928008618948 0.19324429589442871 -0.61302668476814337 6.6020931661572533 1.8140253593548181 -0.41666019119551378 -0.40399794755740642;-2.2083739238385798 -1.9669547058931383 1.3435041162010108 -0.11410469224281435 0.28513294678381679 4.7240598721331306 1.3154204052796139 -2.1047113046948969 -0.32489724326706099 0.85645364556579129 -0.95004073312738435 4.4307080229817162 -1.8514996722869983 4.23740850448453 -1.5408131387277899 0.39304142705862888 0.81123079492708794 4.6175453390759236 0.41481820647767959 -0.73496772253702991 1.7921859152031339 -1.1913044067350356 -4.0284732838028487 -2.4808052257895841 0.85835004373022317];

% Layer 3
b3 = 0.8355011868875114;
LW3_2 = [0.0077930937945823991 -0.13332902842848782 0.78553299674239019 0.0085882062457330304 0.019098182032069696 -0.0057548647870085947 0.011130421675925088 0.2666831420516651 0.0060950032298995142 0.0021776334851771571 0.0038831667496206593 0.097719896743025739 -0.010359814180409618 0.011583921015348657 -0.039818468252587252];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.00111018595614765;
y1_step1.xoffset = 1025.5;

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
