function [Y,Xf,Af] = ESPER_TA_12_Other_3(X,~,~)
%ESPER_TA_12_OTHER_3 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:04.
% 
% [Y] = ESPER_TA_12_Other_3(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999540494077341;-75.751;0;29.715;-0.1386];
x1_step1.gain = [1.00000001035699;1.00022983623563;0.0140762793578401;0.000302617642608564;0.251952632905014;0.00826143816767911];
x1_step1.ymin = -1;

% Layer 1
b1 = [-0.20473474424329354826;10.011752607505613355;0.22102124605081549524;-4.486914259581242348;-0.43458549713545691073;-1.5299331164033567187;0.50329016197581288861;-0.43931178691817140969;0.59411970752172515997;-0.078885445245652038859;-0.30101202433724805196;-0.76362711687957562923;0.45690427177022752137;0.13735841294477676988;-0.63354759116456227375;-0.82219474171269957274;-1.2531454485491178019;-1.0586488169770031931;-0.25491914487603706263;-0.18788673035843170278;0.71120612829295026813;0.68035342533755105432;-0.30061732273505425983;0.89178785008766470188;0.44262129959654644695];
IW1_1 = [0.018273394377742124134 0.15239115782217244055 0.12230140365156258464 0.09829683189341054117 0.22495363166315640591 -0.55256272201865352933;1.2195548779859346489 0.98075366063261415484 -0.87376301309580128862 -1.4463234618755629324 2.2877639713289981138 3.0723466259056362304;0.12517286776600150544 0.17211982704236306052 -0.66392237523523522924 0.72236810661984640625 0.28114033813945771278 0.21866072069943515155;0.97465540432911568214 -0.79816310556014291588 0.85492483418517539739 -1.6966707268823089105 1.0766713085172654463 -3.0411100686660019754;-0.12185595068870147539 -0.53529753124048706781 0.19382454160022927203 0.040548405304537467331 -0.48519113744238373531 -0.13875198374290470538;1.0165258824636600465 -0.54634380628993906015 0.10963647173120513734 0.22447301910593836505 1.2378334400593626885 0.43950195761116789273;0.096403321960290708215 0.20699659935634781571 -0.50330704099252732142 0.48300633048153601701 0.071670137428345506714 0.55113744834040667708;-0.2811463772567410202 -0.22351763765884868906 0.57128648550525362726 -0.3745115288202285897 0.26240158930748147537 -0.38274117527190326671;0.013757919672120477728 0.03921244173183574544 -0.011172263645274512628 -0.82525402074329845536 -0.95723905893703153946 0.82384292687253024567;-0.055565639847794806583 0.43856166507324034765 -0.039601764152937833641 -0.18704497936503391986 1.940583077664205458 -0.96090286362388188657;0.10930221008902502033 -0.20256025509176120436 0.028233747347393027549 -0.080861031806020999824 0.51328998716541651692 -0.29881294358689769464;-0.068996957050001000988 -0.0042227410708452362692 -0.013253463160574845137 -0.037268268718088568003 0.069013114306633246975 -0.21878928296798808018;0.50738507470203886296 0.0011537864241384094676 -0.14295536334825736935 -0.63936756285375817566 -0.79454316546109815533 0.56344386211388564067;0.11048103386145960436 0.0014597500074368181984 -0.10152584504876337179 0.4443641678128130712 0.52480241301711405555 0.16471742732798952402;-0.38815492247207966647 0.09232846145703695051 0.034153344872525755183 -0.72926450041795809476 -0.55450296255602715689 -0.01553128750254706425;-0.66349790809111741652 -0.14493446091762143824 -0.18934041215954813819 0.2885992558772532135 0.60456428790221261593 -0.14795733793593235372;-0.34936781670696898106 -0.15057656612241818173 0.167369299371606578 -0.076672647104822688746 -0.43143928479278137456 -0.30162748215144424835;-0.74861001494096579556 0.058122819505663073092 0.22389137797930894203 -0.18016281271929587837 -0.15690550148546103526 0.056463789686524670786;0.048031885114489400457 0.0038345701446195973663 0.71790063579337015653 -0.20939062541849484678 -0.56711913373194089427 -0.25636876923582829768;0.26921440023437614997 -0.01311130537840135557 1.1538696182715384442 0.046309995687898314332 -0.32758510807467489112 -0.23771671098041405323;-0.20145074488973141458 0.038619594123832659549 0.040293052725456642016 0.11381885237397501576 -0.147549063315824458 -0.096585905564618254604;-0.18287652656560723918 0.04022091733721745721 -0.84108568802102046114 -0.2492014171818854873 -0.042863590142786199555 0.16033756878591515216;0.083638884483195757857 0.41833090150820845254 -0.029381043894374750391 0.4054482575999938665 0.31557587747082405105 -0.26465479718897716443;0.032350012282307105682 0.39383130399539278965 0.057715516476293013404 -0.055965732493180118945 -0.34232770864734746175 0.36937271993108367063;0.42817510262170815905 0.24037455382055769593 -0.49132110667482997712 0.25232790986465153393 -0.44865285423011180832 0.52926462284573982942];

% Layer 2
b2 = [7.3419519294469521853;-10.051764379430160545;0.8377516759533694124;-5.7441453634665267103;31.490617622026110922;2.5162342533233976738;9.4445266848027902284;5.320056769357772275;3.9153739539057101915;0.56505383123919694643;-3.6832886077304438643;4.0157850377262231945;4.9724504589653726327;17.559726800120273538;-0.87597345064715415042];
LW2_1 = [10.409233872552155376 18.025760799911843435 -2.7846995010510409152 -6.0300505929303547248 8.0337497359460776636 -1.5258748986766550448 -9.6386654561890399151 -1.5482120152428502369 6.7425150258137094283 -4.7536661902094268228 -12.802287978154813786 -10.132494632961563497 -10.027903231586359212 6.0926268733000243927 9.9258283990520066453 3.8699172339693497236 -2.8624202077714975267 11.288088754877190212 -1.5648799662368917396 -0.74040366643078359132 7.6292960213469109831 4.8948483606319594941 -9.760731656694112246 -3.1639937897356573693 -10.813256322860191005;2.9191517389531727034 22.543598780945007576 5.4552208590934858634 -8.2833931744095963268 -5.6159608191914403008 14.552451024648769717 -4.8753645569708208285 -4.6775679002791807903 -3.4359615782043322341 -12.134672741565388776 9.7823072697241233442 -19.598528313615297236 -4.6758473333228556257 4.2270729915468230331 2.124871930497393091 -1.7206850422799850886 -13.265301245166183364 6.8763094351922156022 8.0030859295225322825 2.2427202509808386033 -11.004186263144893232 -7.880856510170246132 19.069507728280235881 -1.250988546502317611 -4.2147750637640877258;-0.55952692069900700211 4.8835614690157109763 -1.273717023720794872 -1.5336118808663847801 0.55850402668398990702 -2.3336834435813242905 -2.1571397631281978846 -12.462256503936849583 -0.75217401343715906137 0.28032357182414241459 6.2990483059943542798 -3.695027921886447686 -1.5454572600108951175 -3.1706708837103514398 7.4704717264386610864 -3.484725723964056332 -2.2145246662966102669 5.7973506348310470671 5.4148114002813425572 7.3730396355008087639 4.052793648057263276 7.0222009801192069034 -7.5629158093535577834 7.5136273802455182391 0.78927011037363659085;13.107783107016341972 21.522828122911437987 -3.5475564490646194216 -0.25580803493055537245 -8.6478830760291671709 -2.9276753787005951324 -3.5950626222752632977 -4.3964047163837438603 7.6207291529779341843 -6.3459509831977420902 -11.823310550572498911 -17.568227848040251615 3.0545396661958563733 -2.6456919362580029187 -4.5380338042072025928 2.3660127535612982541 -10.46954553053483572 9.668846788668330916 1.9555778090712656603 -13.104838490539696849 -1.0550045456139045541 -29.520106536881048243 0.55183722496189258067 -30.536451771653108977 1.6183313839878714635;7.9346117876399908653 -15.13671443397563543 -18.760732450491214962 0.8535122063546860538 -6.05288267201575092 4.6703768640068430429 -8.1109552123696548875 -4.4936765961363667188 -4.7794274842101822642 1.6170728488016532509 1.3649382677161776556 -1.432113287334632723 5.8528849342198911998 13.563901691801458327 -12.145427107123868637 7.0657483958007807345 3.3603632559120542034 0.267534369542861028 3.8602669260945203789 -0.26954253407367856354 -6.6903398381916581172 1.0633091980661009135 6.1228009689770033219 11.34507682935215378 -2.4790589091242685171;0.27581327860505211458 1.3185473507038212304 5.9208762244232309158 0.27591517564719264533 -5.4737405981436006996 -3.7664528148090568038 -6.0351530762418867937 2.7032788762040937947 -0.056558466622171753302 -4.0382106551888332646 -6.1613579774369782172 0.043706939778597617141 2.6796814724145181508 0.47425547397265144234 -4.4927407239645171799 -6.9197624342394803421 1.9555433422347809902 7.734620493819445386 2.1173186290386141195 -1.6162690139042903414 -1.8941085519648501734 0.59670760553423229933 -0.056441401883878285217 -8.4746785119554939314 -1.3969939506748285662;-13.565799728555491654 -4.3138537908396559217 -0.35343857296187281003 0.13261080844727707317 -1.7318885561358896563 7.3972276464194886003 6.2890998447392938431 -10.956563083907989054 -7.4204370723677008925 2.7848820939675373864 7.813509734182495059 -2.2711020002601740586 9.5519103128791780222 -9.0783524710311578332 -2.9819635467060372491 7.6298978025946464498 16.790677551753184105 -4.1790771922911158498 5.4867726935050082915 5.3609801225293152172 2.5792670164745667982 14.429377891021486491 6.3792159257973679587 8.3283336875817361999 -15.565541262921650301;3.0501541606241300286 1.8998244838563682624 -1.8954678394225576987 -2.0269103370739891545 5.287041593855639654 0.5887499549938686183 0.34144775368670671201 0.14628125724454929535 -6.3562816623345330669 1.7339597357395348798 0.138314655786681312 -5.4229678859624419474 4.0198230357914885857 0.62534453030752368274 1.9144104138077464583 2.187416255856163616 -0.47969828842267409064 -7.3070577771797688982 -0.19901299114005699398 2.4733838894720130419 -3.0593662619986048412 -3.2210936823355349468 -3.7255779773244523945 0.91158332040525735795 -4.6975782785812922882;1.8480045060337622065 9.3873148576793443709 -7.7145618612488604526 -0.08714139894387480656 -2.1898229661734087692 1.3568487579986383551 14.781653714322427362 2.7848306374202396185 1.0602005801638003657 -1.1568139809136788188 -1.7062391125466995767 -6.4822157339626507166 0.61185989644168814472 -14.471390768569529328 -8.4284158121387857676 0.55747241399825031571 0.037940562688527613988 5.2126769823915415358 1.1796744413880559943 -4.3923115837027424391 0.24426259550567608336 -13.467054273835259437 5.7021526200665775264 -11.998400718259006226 -1.0906898967178251247;-6.3114571547836559873 -3.8696619133749328157 -5.2570977206104787527 0.2387522517054398441 -3.2991299542138392553 -0.96207909869794050461 -3.2169976648293263644 -21.167204434678367875 0.49904817272570345388 0.94520975238415005926 10.48762223553748818 -4.9074799524753496272 -2.3405859437712295268 5.827589854482604359 1.1657647140782949613 -3.0998405375401061335 5.6372816702661951993 0.19611204089314701604 7.0607220692492091985 -4.4381331279581397453 6.4456771660224507769 -4.3095559290095062366 0.67669317366694525795 1.9755670105162190797 -10.001314745995399491;0.77991364231198978541 1.0707512305217496973 6.7358624564684674851 11.804602674250462613 0.8004755629437061426 -0.093788271160862357601 1.712063013164295544 18.019951218374611557 5.455865064689626287 -7.1454372988524230692 7.8833235863077604222 -8.7981769760699748417 -4.86553764963678681 6.832601340204146112 2.3171788888953281926 -1.7845028999569061412 -12.706479979110385514 0.77257010589783792209 -11.892777145910333303 -0.29615661151772287862 -10.334470128943156197 -6.9575487275365759388 -4.3983186833672105465 16.713813016475935314 2.3253498714940894487;9.6909637979974316835 -2.0573002014734615983 -11.893867192019071055 1.4941719073380321792 0.40798304050730171966 -0.22772028009941375215 11.65588291823295819 -8.9231764252237510959 4.7688681078712757966 -3.2949004226461888756 -10.598083872584256326 9.7897229044566351774 7.237598592250721552 -22.202660562535264432 0.41823337315592512553 9.0693159784136021528 -11.236280718112549692 -6.5489927576734752179 8.9312610742726974422 -9.0056409322352699576 14.386674561964758468 -15.498333887461610558 2.1935075557170371496 -4.8172763438651271528 -11.210427707071355385;-14.004284039980872123 3.0313596224705423587 8.2406283218356790599 -4.1670011705247498313 -4.6806340932778134345 -1.5994065830255490557 -16.375085575258019333 -6.1813287208329121825 2.8168923217950214699 4.2755074089282603111 7.77919413136251503 7.4778088721973379194 -10.103801382800673636 11.080003947149764443 -3.1592092258889699075 9.9777159566421520509 4.9985261567619332368 -2.2772705265760611049 -5.7703498121331504578 5.4301262441956188454 -6.8845069918795518404 1.3663451691088817075 -10.461795931936912041 14.643513390390449302 16.078709071576692935;-6.2381355553118105206 -7.8645608689026262894 3.9168448489816882763 3.1979611073788309916 6.7178193464249416067 2.8928758810085573039 -8.2107939787971577772 0.098878336534324748341 6.9553095533633184644 5.9646001086158291216 -14.860339837647753569 1.8777215468581804902 -4.0254869299727804588 15.043484141683725497 9.2728441151641760598 3.5166878018740432843 -3.8204640881612514924 -5.7403366933076194201 3.419476489863419566 -0.5046728811256887548 -13.788380011122303159 -2.8975080780696171878 5.2240270807125668284 0.85151571760736510797 -0.32114249263885596974;-2.4033274104967792262 11.82228963286546275 -1.5959655717634435579 0.014385168444642227811 1.7565148031143684726 -0.29480204346131494431 0.033946686944628587379 -4.6217438545112869974 0.5254110750533103058 0.45512366972293205603 -0.64597343613997748601 15.693008129995913436 -0.16794549116264639865 4.9017799658971288324 1.2006505325633727299 -0.37261523778897659431 -2.7372449087033365878 0.86982837270803881058 2.46160693912947659 -0.65836899926319880638 -6.0109715017245566315 -0.30520802529169238637 1.0816183080577057574 0.64540998336560684301 -2.3328771821549167953];

% Layer 3
b3 = 0.54213042289939383789;
LW3_2 = [0.023950942872978008696 -0.0052890042988096063262 -0.094408038134950439391 0.082333491154153212732 0.022234628185623096092 -0.091151238033931589189 -0.12422388053821367682 0.2190000400760176158 0.28455151867626510276 -0.52734708729223134238 0.0056687973810258911189 -0.0018789001747846214262 0.0078943242147091041011 0.19397912349424281664 2.1245534682402791304];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.00492610837438424;
y1_step1.xoffset = 2075.6;

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