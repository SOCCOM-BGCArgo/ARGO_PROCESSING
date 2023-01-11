function [Y,Xf,Af] = ESPER_tco2_8_Other_3(X,~,~)
%ESPER_TCO2_8_OTHER_3 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 12-Jan-2021 17:03:56.
% 
% [Y] = ESPER_tco2_8_Other_3(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-1;-0.999998793562952;-77.961;0;29.715;-2.08184650291172];
x1_step1.gain = [1.00000000274156;1.00000063368067;0.0138606862425759;0.0002963841138115;0.248478071810163;0.0588359962254402];
x1_step1.ymin = -1;

% Layer 1
b1 = [1.674592665992318;-0.51078342816916278;-18.213668633195862;-0.52704869755734196;0.67770735278235084;1.7480461645945069;-0.75666652671786905;-0.35912055081210442;-0.41233702275456624;-1.6767701241696522;-0.44621044670568727;-0.011242375541938927;-3.2379161935255576;-0.31724927210409182;0.19689864555173145;0.14739413768586881;1.2549675532169309;3.1065133695758629;-12.203803551054866;-4.6468874385174237;-4.5340996899119732;0.82718234787128231;-0.49740302006141457;0.28576282869622294;-0.65096042515819663];
IW1_1 = [-0.18252270396693282 0.4961765988486741 1.247228564386647 0.21048317809488187 -1.0713818131144197 2.1631920375745288;0.013137696649406696 -0.083196785253345576 0.059202160380384151 0.13783208139613634 1.2912739210739328 0.44445420493386972;4.6227685557269993 -10.323717178490813 -22.088116730361232 30.850638359431045 -10.18752501736309 22.720330369895539;1.1846000966771963 0.64092995302530409 -0.46982131762282364 -0.37823704542600223 -1.9379162903666949 -0.12948184275350133;-0.32354942504937667 0.66309411531118356 0.083399493746508305 -0.43058218264638337 -2.3219086471076742 -0.94728335081187909;0.11164239404652024 -0.29205040924570791 2.1084398101562387 -0.17679913240180875 -2.4098936337647139 2.9651258902441384;0.37690303708260936 0.200039060647658 -0.3534835397223961 -0.61169706359891551 0.36654498434989308 -0.3253976150843198;0.55389869759803279 0.15413231919262571 -0.64082838001869269 -0.18224356787353083 -0.17003187800163541 -0.011236122637741692;-0.087390483405705335 0.15342961678717065 -0.11250268640584857 0.24219623473125329 1.7512504009688483 0.098388955053219554;0.34447708281084005 0.2997471554885538 -0.065300253806837572 -0.95813532159562742 0.55765576461145994 -1.6290196744638545;-0.30394912328849688 0.33059749021542195 0.27059017673163654 0.056144617848519936 1.4765794972451294 -0.32027009663205674;0.44490172471792278 0.044607034455946555 -1.8377717419638357 0.35302243302573733 0.76503918392317327 0.026556891582285937;-0.79220035140340117 -0.060806958760962565 0.46231990243634691 0.32336956806536715 5.1516694127106852 -3.0691368857643622;0.69915010962948954 0.34991937391823591 -0.34492397301706562 -0.54288815194436213 -1.1987425312243971 0.3328891729625777;0.57696278413587521 0.11984393904183309 -1.5101075098586449 0.40289712957914536 1.128166403227195 0.48274559863917865;0.35485057278824639 0.28057279546798009 -0.41952320607218724 -0.21424979493179325 0.34383297192322637 -1.2325237750650375;0.38340532585293191 0.11182282959315937 -0.4743104928689425 -0.62127592851347335 -0.81570613201306708 1.5977951651297189;2.3305917158996396 -0.52613903876016743 -4.2725934410151849 -0.68256353303805106 4.5066708592812841 -0.9618628952137277;0.34005012820302533 -0.5569082430873743 4.3139855468543145 -9.3253541591004403 2.4491866656808501 -2.2906084533354374;6.4736309802974207e-06 0.16016745212468533 0.2366250469563157 0.098398096910814858 2.6196426682269331 -4.1984012395599928;-0.19321101658869957 -0.09204028970086163 0.94117369848913301 0.018271711322831116 -1.8973523794531717 -3.6227272534610857;0.29606724306342252 -0.35873296497766743 0.60213662582106164 0.13455061068665411 -0.70230891049046429 -0.089546018886596149;-0.44627722645602286 0.23730069865537248 0.36407311296246625 0.025287192882264387 0.81841923740528744 -0.64016832505382681;0.2685672711747939 0.049655147045805592 -1.3852454249245287 0.20015364561782581 1.0862407272456172 0.36819843116046286;-0.30044550809154458 -0.15872966864996099 0.23880918169401191 0.0068686820144256717 -0.26137901092411542 1.1319281691896861];

% Layer 2
b2 = [12.162369412794831;10.791952887355636;-7.1969591201598906;47.0093715021995;-6.5402762612686178;15.400234910869791;-6.7085177905047741;14.643128007303437;-13.659790131557459;-22.585664056014473;5.4296740917531867;7.7836196983805959;-6.9014978079674627;4.0897265566339724;-3.2021591185378724];
LW2_1 = [-3.2076829968337592 7.0236483233658618 -23.47120470993088 0.46892850131812897 -27.662812133408025 2.0498242138543028 20.399790511796237 16.581414705215536 -8.9538169505517757 -6.3898879870836414 4.9796531257546714 -1.6157717518194601 -7.2885080811861442 -1.3987940340601057 -13.870454168354684 -15.497235780642857 -7.627125701760419 -7.7631505044227183 11.336480223659319 2.8657897818537057 -14.088700788310524 16.491457417148595 31.043601935112843 5.2654297303809683 -4.9284146417007131;-1.1604401536042308 -17.526429222589787 4.5709864307906862 2.043930360257844 -2.0259126088344792 1.4380756207569225 11.211622308293315 -6.8228607918990614 7.2639855360614511 -4.453640606685827 2.2037103172984116 -0.22961136413412248 0.6082837320670369 -2.5501632693579754 -0.20198302729929254 2.9630109581435096 -0.5253351710728511 0.90944828021351343 1.0705913617080272 -2.511670592371364 12.105715585805937 3.4175371802517232 -1.9485914612196455 2.5161990655306101 8.3895331663206036;-4.1335595657698336 -54.671881643851137 8.0901114062001387 7.2835315636239786 -5.6300669222413822 3.8717395223709277 56.563001575750633 -33.808905762265077 13.689707548547908 -15.740125055425302 19.914445193132675 -3.0567221984554243 -3.7170465028442417 -15.594317394194251 6.5225848036105791 4.8911212272582807 -16.461070635792627 0.90415761518972193 -0.65755525674777959 -5.8834838005790804 -0.45179971878489922 8.1656093580962708 -23.248361280675596 0.1057020857325176 22.502202480956004;-1.1554023532825384 11.818228813551929 -0.92361256514120393 -3.4348437576249329 0.11912883765395481 -0.095338881615600959 -22.885780424521229 13.634574985226182 16.366706803249045 5.5941081573637765 -27.735683436151547 -3.3760360483648149 4.073824273860903 13.189838789745545 1.872656207523828 0.52445541087050529 11.406624978859201 -5.5958004111908428 1.7260155688921925 7.0141996020182829 25.073720872786158 -4.1859090597067636 22.358482278069967 -5.3100456904269375 -4.4710354693329801;0.81457289161535262 2.7929294281341761 13.991088135734833 -1.7639395367174802 3.2776031846654652 -5.5853984915877515 11.658148825284256 -14.303948472546038 -26.112534050903975 -3.9550247804078063 4.1706660171279752 -26.104003736139514 -2.8923879951283733 -5.0613865775450257 13.792208324744548 -1.2999528250544188 2.2378605948840349 -5.6132905393252441 -2.6918142007523729 8.1977362750566272 -3.0531712555137926 7.6633370217078349 6.5092962966236012 28.063339009227807 -2.5270478461835557;1.6839131870271895 19.953776038465467 -1.8892259392852202 9.6189809469103924 -4.3014095598007511 -1.5601793877920562 -31.543132829683071 -17.478869934392996 -10.058145458494492 12.146395965716868 6.0962924315250708 6.0328279911622618 6.5896272866188959 -2.7753221470862157 3.7114938296844353 -23.763502371799479 18.881993636725674 2.6021021444452472 9.959862253278331 9.9183210611742147 14.29335695228524 -14.118700343447671 -25.17475721660076 -11.937888190980784 -53.325162137903732;2.9255815066847677 6.787380758374451 2.6114538443006463 2.4772217741182971 3.4549306544078018 0.10248231235730534 -2.4898866595461517 0.3118573040640204 -6.9468281843656587 4.7260263625864845 7.3338404452115693 -3.5137631008785681 -1.4709899260858805 -7.2465402955484768 -4.5115223136015663 -8.9333238712463547 -6.5078298337111162 -1.0931854512436552 0.50677249400282776 -2.3188620814263459 -3.6355852220560352 -7.0672475809605668 -17.880197581127945 9.5715688012271656 -17.44716808710718;-5.0037906024238863 -37.29639712066578 -2.0232660780016647 -2.6853752457318039 -5.7900538845391045 4.9043188965549627 21.201095080514758 -0.28518778676857343 28.179381308526228 -7.2799511780774635 -12.314331652325507 -7.0107462345038174 -1.2462889148853202 0.3744350563500784 7.1414148064532155 10.763060806619725 -1.614842807074089 0.86763416758738876 31.585963578590334 -7.3847404173287847 -10.404905591792911 15.628283407972601 21.608341501384302 -2.5558608233651747 26.730934763145463;-6.2095818135492191 -0.90891954916997419 -9.0558617313140974 -1.1150834076739271 2.5143610842605653 1.2481329268416952 -9.4112744701528754 -3.2565842667499409 7.2707720360515387 -0.85150677017539433 -2.9277682403232417 -0.74558192530831513 1.3202890353589372 5.4061141799476342 -1.1865627099696714 12.306871311925663 11.076315699702395 2.6853404403918231 2.3602631531062097 3.1364537210293304 2.7005104728619793 3.5234729698068756 1.2481207992075771 -4.5127454028813219 9.9767735045851182;1.5802850710716754 -17.394746575560422 2.0339578926841759 5.735155774647823 1.4372752961745143 2.696431594191421 8.1210348314390526 -9.086431777123968 5.9865881961876912 2.4150792517232427 -15.34305801841847 -12.765632925203363 -2.229239045546997 -9.1146731375914793 7.8096261774571785 -0.41526204375138576 -6.6684197740117845 20.287937665587883 8.6374961504243419 -12.076435368065178 11.751779614511781 12.019542111669857 21.127614394861112 25.226877268606312 7.9853632386167597;-1.7310077960934909 -15.268310702671757 1.5511372361623796 0.5124763119676361 -2.8767352208173378 0.25889035651304687 7.6804591460259637 -5.2458843119906122 5.1487137705533259 -3.1520578464991966 7.410495108069048 0.95960557451112427 0.82323913899132606 -0.3842659936172424 1.9599409447286842 -9.2950980421879983 1.5853809196289752 1.6553716799075162 0.69401400292036719 -0.70856948277716869 14.841115119331919 3.6679349973073547 -4.1887078196049741 -2.4749662657879155 -11.343818560558169;2.5081579150905773 5.4876832066031405 -29.561447452735941 -4.2820299211485171 17.615880185262665 -2.8989434827945884 -6.8723848843575874 28.980654573476965 7.9166145919368933 -4.0268899018582287 2.1631385766192328 4.4764302154412583 0.87886437909101756 11.945099258582786 -8.7517690915357811 4.422520028863004 -6.2168907391578907 -0.094581314063278812 2.4595953699282482 -18.470939901765906 36.221231275979541 7.2833491739633427 0.2006519406584428 3.1686160256631455 18.784264373685449;0.50666493613204111 -0.47008390886445689 3.9057019161541042 0.53333095880630488 -2.4139454852691498 -1.8432577715780285 -0.033308403020864796 2.7392427137897561 -2.3210749434190219 5.4801405275464568 8.1891244618776984 -5.9861984374392581 -5.5387273605260363 -6.3000673847885444 2.0616283772915613 -9.5322115967578007 -11.880551261969009 -3.6689865833338828 -4.5447371004588142 -5.3353746488941551 -10.154419539402875 -3.2555496820935779 -13.6101969654234 1.9789294914940194 -9.3732396030579821;-0.77448788942050339 -9.1095759320528895 -2.7457018004317044 0.51743719587438597 4.3708892956458092 5.2005981970610327 -16.469528017489981 11.251809227980733 21.915800859123348 2.9194758979036974 -19.68194033004087 8.4443590195323743 5.062218349536324 5.2601639103691546 -12.930618082232732 7.8115037620238637 -1.9270161814579365 2.7468027308739931 -18.852984630678588 -5.3892871966245162 32.688844646594774 0.97037594369416158 9.5661159523030186 9.4285032115213898 7.42175794066178;1.0541475685356414 -9.8211569880952982 -8.3806450756807127 0.99899997243360084 -2.8903864021340762 -0.079411973189704332 13.349570281898696 4.8501023448083984 11.392993984980309 -8.2761084761867778 -8.9911073936385204 4.3644946274224088 -0.68546592124143979 -2.9403579047641855 -6.1038686103853292 13.561159004253161 -4.1415398681851352 1.32467564684642 1.7751282178813996 -4.871920801723876 -3.5042930572622049 0.65981292430041627 2.3569243103253958 -6.13072472771099 22.044109187499249];

% Layer 3
b3 = 0.13936090325989331;
LW3_2 = [-0.01673814464056796 -0.35368730616540911 0.14260432194658271 -0.082457731801840503 -0.040788352685022837 0.036752099848464621 0.22532353051216844 0.14975706239981945 -0.12017102415503059 0.051427418250180351 0.32179276601174317 -0.024567133148654191 -0.08142037656520483 0.020826844879420797 -0.1255940793866494];

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
