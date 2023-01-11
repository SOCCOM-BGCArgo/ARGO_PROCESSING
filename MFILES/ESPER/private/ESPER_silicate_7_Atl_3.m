function [Y,Xf,Af] = ESPER_silicate_7_Atl_3(X,~,~)
%ESPER_SILICATE_7_ATL_3 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:40.
% 
% [Y] = ESPER_silicate_7_Atl_3(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999994516886;-40;0;9.1619;-1.99149636107278;-177.232103709389];
x1_step1.gain = [1.0000000017546;1.00003626070655;0.0153846153846154;0.000309885342423303;0.0656325902039552;0.0604107502911417;0.00425400504498441];
x1_step1.ymin = -1;

% Layer 1
b1 = [1.0348901154469924979;-7.3528399064010185882;-7.6856373851851778056;-8.3350287028800540412;2.4411236288744615486;5.4321086864627377011;-6.3469999099297540468;2.5913656077019062707;6.2905596249463044956;2.0794035381576390087;-0.035054946601465526712;3.4353068626834128985;-0.35897914889866855415;9.4505545115766231845;-0.55480768490887300004;1.0892741704325545893;1.2684370357882825964;1.9367562649788381268;-0.67928679004412162445;-0.92556902840187138537;3.2569157362740317474;1.5780985869050636605;6.6452150194030812358;2.1550527724525427864;-2.8712054171004348646];
IW1_1 = [1.0734273273117511671 0.033805940056888195144 -1.475375173559514419 0.31755597825362935627 -3.5755990464484521674 0.7827586588546934232 -1.2421937023406948164;-0.97570345421578419387 0.58910446775252234275 -1.7696579480558796149 -1.7855117664931019483 0.80359221114669732167 -5.7102931702764330169 -0.41891822618990481875;0.55759510659286559697 0.14458628796973599373 1.4418248094150596383 2.3167915058343915291 2.6042074612674350753 -0.3219422946705214339 -0.30751367943169077845;0.99565568547065563809 -0.59489966403059635347 0.07979183737480161076 -0.1675646887145314512 2.6171377547096517091 -6.2852032121484926463 1.0406457860158562223;-0.091048083289004808072 -0.016087210282664701766 -1.7023975820032182238 0.21414243910156249839 -0.83085411973178391332 1.6455230621551943937 1.6054413167560308562;-1.9530353088713694465 1.677166668714152431 -4.925114340561367321 -1.206709816270933322 2.4771394178393784102 2.5228462784083243164 -0.71738729511645760173;0.17538811886321326039 0.17175460188757418756 0.84247984365147576824 0.12112972255150775136 5.9720986334660333128 -3.2782260234304780866 -1.920793581464708355;-0.22373812382190916281 0.25813174453440790312 0.44179917627237214051 -0.19405118566014481485 3.3246928807426536778 4.5493025592975495641 -0.88140736918580508963;0.010469000660203002517 -0.049530570117294975274 -0.10848797946388366586 -0.4245329101331091759 -8.1312632091049756866 0.58001640522376929354 0.50933309556831063425;-0.42607893754796205554 0.16024535131900344287 1.5081347309315806804 1.1528312665305027185 -1.0842336785846753422 1.9790439469740250811 1.4532604322893329929;0.90353387148854535926 -0.4186809797393400423 -2.0106557905950022125 -0.9554228380142033572 -1.2345398961275024785 -2.9150456698080371432 2.0304427356055496645;0.025896224872917775195 -0.77272920051219562421 -0.25673749691700509246 0.78822683989095188029 -3.8646582611762059223 1.8463934875063092722 0.85258330500861378987;-1.8030771777259786415 -0.10954896414437516328 1.1250245295531438039 0.86161811625479522991 0.78570980500617848552 -0.028850895166465030661 0.58371932980190976981;-0.44136923073634021231 0.27412889610395657636 0.68984040752690289899 -0.78877083952090276675 -11.779748066568735965 2.0318459480402553119 0.040771324963104739469;0.86968451394836843171 -1.2327108267103863426 -3.7991709240658044777 -0.92833374793454226381 -3.7085723375218786479 0.17719099121457496104 2.137972459840710826;0.013421026363145025298 0.39720168802805583264 0.85318497432891926557 -0.57050310317405095173 1.1947712937052330506 1.4204802614911438141 0.015011877270779377946;-0.16126966408275172915 0.28850007110665037979 -0.68245127813329631294 0.38777939673453887259 -5.7638876235161422557 -2.3034658836284287453 0.8978120542671637061;-0.47916231646969975611 1.1561986210474313186 -0.39511719211526780127 -1.0160465761541872531 -1.9141517176912268283 0.15991852113563548299 -0.30530482651902357549;2.2440205356236391232 1.0769862840644575108 2.4025673856354536895 -0.8487237144116897003 4.9976009413351061283 1.8169102558759602228 -2.5813733425559250989;-0.35578672800797761511 -0.098463742222592992581 0.7165059938796939365 0.16081050110965394717 2.1932581508748048016 -0.21751782727704549458 1.4384620978786444301;0.4862638110214082543 0.17913908922610471963 3.3620227846297496832 0.17243890796796415565 1.6701273523261235265 4.4148909099163331504 2.8403015028579816637;1.3477846895478295952 -0.60143675932493456671 -2.4649523974926270675 -0.77318312948024847664 -3.3115387972974756892 -3.6710820900770957387 1.8852050395440480557;-0.67882872146345862863 0.016209379561207371695 -4.4813994045170248981 -0.91845828739981405509 -2.8813560696864755073 2.3795321411392991173 1.4315157875000044108;-0.46205647424138163171 -0.55304002294466136913 1.8554470476113271715 4.6362793232493837436 2.3626407714731123555 -1.2086339902881262365 -0.43367506426030738576;-0.74674499516401049437 1.6871641592951580613 -0.58566850749400023446 -1.0047011030194132886 3.0284095770373755663 -0.70587264720748765789 -0.59348857660574760864];

% Layer 2
b2 = [-0.23109478300011679353;-2.0151139153181465424;2.3689377227257715219;-0.36451063258514965204;0.011321083605662328464;1.7435513559947490858;0.84198126180543808594;4.674207746728969326;3.0801313098311626781;2.7814523129732120488;-3.4754727547744774796;1.0278423650059254069;2.1228165904843123712;5.8669563615793229872;1.1186160352466669732];
LW2_1 = [4.3223774767389882001 -0.5168371359100794038 5.0246051481289457641 -5.1218722636474725007 2.643895457227063428 -1.4955713627761340057 1.9962516512319674611 -10.903370962937675159 4.5049053444570201776 -0.62488136974203078999 -7.4100559636365765925 -1.2620187500610398068 -0.563063994391767908 -2.3163769654786463548 -0.67193738773896427574 0.77809906559429842243 -5.7249715896873203391 -0.95880939679842736023 0.3593637203160739424 6.8057122553098574613 -0.080072054790680108338 6.7658835848959357762 2.1065671006872475068 -0.3019871882892319892 1.894192283371317842;-2.8072583947695406437 -0.1627749436763729074 3.4247510456071386109 0.59574384293177451433 -0.48715443874640085298 1.6460148710784783255 -0.68045743755462839708 0.73374273164471037756 2.2220146236346951518 0.8159997734138687564 -0.62004008988911629707 0.36055557726626807913 0.11349810531950903225 -1.8916821018084166717 -2.3926456116355012682 -3.0557206661512950063 0.77694466960464692029 1.5944053372183741857 1.5888934323495580436 -2.8965384424467108815 0.93427813340324650238 1.3314179855932826158 -0.93179796956291516263 0.5461030997719795721 -1.1351527522220259403;1.2135669227055276842 1.739413144561983593 -5.5364848552725920428 -1.5229113557676137169 -0.81492300374880677349 -0.21768228298541961174 -2.03927775148796675 0.16240175085964150825 0.48049740355474912512 -2.1132070031799496945 -3.7796979981647873537 3.0361895436039612051 -0.54591566043359252003 0.55030991850385824193 2.1985474560486739648 -1.179529583654184588 -3.4966039114092426487 2.9997789444916405088 -0.050854407245206609267 1.5099925261626339079 2.9500281945733570765 1.6885508201590613275 -0.8240283761899429571 -6.5130321131274913071 -2.4077094872645106172;-4.1128355998151420181 3.6677190444654286239 -0.12187555987759779441 -0.095461615790904577605 -1.6782949333470347231 -0.38623508035168230279 -2.3010853300344162875 4.0750619416142281892 1.495429112262335547 1.3137999351147802152 4.4258026230991953653 -3.6300192976107648057 -0.18584518579354911139 -1.6418583464599356692 2.3821696544084551839 -0.59066196040676877121 3.2237923932864616106 -1.4053295345674827743 -0.23933135624645851536 -2.5460666127257494118 2.624166431251074183 1.1308807993448040108 -1.8472702490076859938 1.7025781256408321074 3.8612306079990301733;1.2391735932310419255 -1.0163024448495832086 -0.76667070657862157468 0.85946747461070560625 -1.565298783990274778 -1.0806886107005921893 -0.97435622767573937253 -1.1020964913478028535 3.574195826855634639 -1.3656606794997054877 -0.28646640803876610715 -2.1219183949799269584 1.0730827782023115713 -1.0369635443905784733 -0.40574563761007664819 -2.8735945832541727185 -4.3188198457553532705 -1.6441650644103105172 -2.0874849519669980147 2.1214609948075753287 2.3449857707987362687 0.065745535390470843939 -1.7355505042916701264 -2.8379821868910188876 -0.021923073142188144458;0.28096588358794516438 1.2966833719393295876 -1.4389428194721372023 1.5830284044622160433 -1.6719299245037071433 -0.042707124396136469935 -0.22548100670150664215 0.62989600532258027688 -1.6548300810930227644 0.1057646591348933296 0.039193818560430561393 0.4161915082514701103 -1.5194056823817925217 -0.027766352058978462092 -0.85748598052484203969 0.10975673218425233135 -0.18452227197873902198 0.79082589157966409754 0.10833223066377170585 0.75391527750848996714 0.33206385646264202638 -0.61206693970651515624 2.1581459641069389299 0.18950957431971965361 0.070700085454047142952;-0.79106502905241649071 2.1213660928538122796 -1.8443862514493916116 0.19277337291230636773 -0.20798609766089590223 -1.6001784405601533479 -0.37216676910886220542 -2.1450868266836602416 -3.2884838585815692191 0.044756760216625397775 0.23049014374064188959 0.95136402508370343512 0.39382388597224077254 1.8450239467434605967 0.38547810996025905217 -2.8647253107170116415 -5.0635199363398344374 0.19865666178041394008 0.59430275888470740497 0.15139424846879009912 1.6420871567157488968 -0.036475590518217235192 0.81494900811982484701 -1.3244650836525486692 1.6066277516335178976;1.5628444396792819937 4.1118667102595960472 -6.232503112504716114 0.63860365626151083251 -2.8435213652895487968 0.30873883346594072075 -3.8138869398867867311 -1.9909695143072707602 -2.0821222986978584757 -0.85363432431123642541 -3.7496225076197648995 6.8722271643798018559 -2.3239880559568879548 4.279324264133048672 0.49518516382593813496 1.3547150905857141545 -0.19815445026696615338 0.85868121712297984516 0.50612515697006543736 2.4035281855106820359 2.3630567391104260189 -0.22550833526239957849 -0.073928005194610596496 2.3292559831517198887 2.7751538939832984809;-1.0726809014924629793 4.0971709863112817374 -3.9406339580612192997 -3.3307575204950228276 1.5366186306871774114 3.6272783929358225485 4.2356998441611226625 -8.0312697638832535318 1.796413856073908466 1.936312759410411477 -2.0045943402321189097 4.274365644479980908 -1.6464664097197405646 0.25220188924081371473 1.6604815790854026147 4.3192124549744974615 -3.2222160082880986209 -5.0483564573600032688 -3.6959270655155793683 0.38099816699985361579 -0.70420226973575927865 0.21162487094034238 2.7465143064770018988 -0.18452015392075038491 2.3495667192495179876;1.2298682736658179682 1.5112173175636367528 1.0750804381172995861 -0.076354206836569127326 4.2535452062376908344 -0.41420372203028166247 -2.6846864897893469859 0.41332634745898277906 0.96248292298167781311 -0.42944708038313550436 -0.51166705110528598599 1.7066639059505279707 -1.1603910391868297314 1.2358528143215163908 -7.4234031514582277111 -1.0033414437006717979 -2.761887924574137454 -2.0494652930134602542 1.0957599458804996306 1.0805574779888618142 -1.6720409951854888764 -0.85615592195810852338 0.081827073782515857836 0.48817838260080892443 0.6006470896164711748;-0.81209907599128516864 13.514815985683595656 3.3078624612502189883 -1.63363563914803267 9.5356995274315234923 2.4807890858866832318 4.3713434623067577434 -5.4593062188498544529 1.832036737838010021 0.016668372904145061647 1.8352300981970797178 1.4749572577231151449 0.71455895932004398485 -0.44799746876122076289 1.0550932607425398402 0.56516174331824420296 -1.5338523846780220961 1.9211550714913301352 -0.029285969414590027476 -1.100077025715832546 -0.7420921082260495405 -1.873711590589754783 11.084210294219978721 0.37670769225779066058 0.40053112471133706096;-0.16133453307108261421 1.4802559097654932607 0.13929995181081281097 0.2793459576032096825 -1.1080747009812186299 -0.23555215622124792141 0.40891802954842904061 0.93363728696882153368 -1.9920402066117162132 -0.043467865797023622076 0.92319492210271192345 0.81826398922617660237 -0.67658826696821083857 -0.32062324814143500218 -0.29213452006451601228 1.8036841684762088267 -0.082078662492284432117 0.94349669840184657854 0.058320775202050111063 -0.41249613914814781124 -0.19845584123546472499 -0.95152494674547694054 1.4557586448836079551 0.45615205112418516498 -0.058412784771496702174;-2.0035420516612560959 0.99708488315042320682 0.28406047080107804659 1.0273899765626679503 0.91388016577024389075 0.51945952772273595865 -0.83419369981432200944 0.11862590874877525227 -3.9260652115297265397 0.3606729617726387449 0.90410270737706832644 1.0565848076996799954 0.046294013390188140655 -0.75717630537598468177 -0.11866226004989563481 1.4538010811236219411 1.9474622041384457294 2.8198854228424705681 0.99993417870290890459 -2.3206974639158675089 -0.44801756911033879804 -0.54362635041300511052 0.65680683851996313205 2.0475641940048889822 -0.34705979065525294214;1.6663096796687422163 0.48025372406140043324 -2.6202830520551314564 2.8387879610836259658 -3.5366672270904024877 3.5940385752870249725 -4.408143339723004317 -0.090644680476414135617 -13.063554520395960878 -4.8602555657276074186 -1.8503557199437399294 8.0894901913458081566 1.1124890045665571936 2.7443658492323965703 -1.1083485362732405921 8.857407340555571551 -1.2213397999803512484 0.99708559567362331499 -0.77245091863075776573 -2.7161058938944582231 -0.72479053621283640041 4.4962756582661977944 -4.0875563951095372772 -0.1275554156410505513 3.1292358593661551147;-6.2756923699834228358 -3.7104054121894614049 0.03459934686229152595 1.9844640770831758481 0.14649829015551010802 -2.5094563416748729701 -0.8612168743581248842 -8.9095264818817678787 1.5598304824293904858 2.6735838655316652712 4.7341105294674648007 -0.37660014874310460131 0.4442918868783551356 -2.0266908624659474292 -3.0317939757591010341 4.9688873054739444157 2.1591199604588213745 0.96485192321428825402 4.0205264815668622091 -7.6224180238580077074 -2.4550066782937509657 -1.9885463609548685859 -1.1739320861985893529 0.89373788061233894187 -0.13524357125409697256];

% Layer 3
b3 = -0.56437915012268247139;
LW3_2 = [0.44666563962378080799 -0.072202643429184337753 -0.081574745401214815432 0.019208717339667231994 0.038242262236783473861 0.58424611202334064686 -0.054412934845333250922 0.053638932141695482014 0.023506058929161657151 0.057857132620929202504 0.12153727193362609538 -0.31270064014431731003 0.068351910501597196168 -0.17863055748620243479 0.028856519838192408456];

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