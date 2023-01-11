function [Y,Xf,Af] = ESPER_talk_2_Atl_4(X,~,~)
%ESPER_TALK_2_ATL_4 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 12-Jan-2021 17:03:44.
% 
% [Y] = ESPER_talk_2_Atl_4(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999999390765;-40;0;11.29;-1.99149636107278;-0.2178;-0.12];
x1_step1.gain = [1.0000000017546;1.00031699124092;0.0153846153846154;0.000310077519379845;0.0705602618858232;0.0604107502911417;0.0479365703301392;0.0145655815308426];
x1_step1.ymin = -1;

% Layer 1
b1 = [0.40122540972694748;-4.4069439963351105;0.43713833748461767;-3.553661035948021;-2.3698911304959198;4.2178886336351065;0.45407792405907277;-1.6032115170393828;1.9806059616951603;-1.4959478905971157;-5.3316088288471732;4.0429322854689218;0.58198371695398921;2.4881471953122838;-4.2305137958909205;-0.073095761295140282;0.59520945001001646;-0.60965639813633832;-3.1822053853057106;2.8256156281034093;3.1781226881399101;0.30129675533761685;-1.18823563758622;2.3574547406292448;-8.6834910253297171;-1.0793768616799291;-2.7074439594731143;0.1885764068393106;-5.2446231361871156;-0.90555248518534137];
IW1_1 = [0.71755997912515368 -0.53542183827595402 -0.945752698736953 1.0886389898640991 3.4071667993989712 0.37331044452213363 -0.43350394250692892 0.97892028440087342;2.26057871557017 3.023817527084415 -1.0639143807293869 -1.2845650216272728 -1.6258711388258507 2.6630447015458749 1.7120616840099756 -0.40986381790205706;0.31956480255831538 -0.5333527633073577 0.69565755771325877 0.42259868481103535 0.63759786874194491 -0.45463377909334013 -0.40821061431433192 0.54418897823467449;-0.092197821617077347 -0.23869491514903105 -0.25736990823054501 -0.14589857439567955 3.3359775619730989 -0.82399090355268712 -0.07386547120700207 -0.21761881199255118;0.35541763175340668 -1.1403108632846641 -1.0785288867264797 0.23942501020969886 3.2008001144011242 0.67230627497397477 -0.19107321017283344 0.06337859398881486;-2.7023838563851266 0.77105664183883194 -0.81046275175312987 1.62537248304418 3.2151866907616173 -0.32474492185476678 -0.56346884859023982 1.1172163519894855;-0.13074107459457693 -0.27318263162934819 0.14980638179348552 -0.024251476838009309 -0.22794705499078183 1.0033468661349034 0.012952544861082115 0.50762990151748077;0.93938782655104303 0.31217673711960936 0.15909011797002243 -2.225155236994679 -5.8623507803635579 -0.6335295701997995 0.40236484204161305 -1.6718872375531137;-1.1644133310133142 -0.46362563792961387 1.1173188684649222 0.34573031857290626 -1.9394133233336046 -0.23237142082134349 0.83921520011182971 0.45879767156100237;0.063352540561934428 -0.44223863689205256 0.34326099460528586 -0.068855731001039061 1.8765127391895899 -0.17093356887368896 -0.045456170729497293 0.67685817164793971;7.1975502191821228 0.026147120311144585 -0.98534171798047154 0.4536940745105123 -1.8681624074325405 2.0938892518340282 -1.9063697776715218 0.89414789016613905;0.23129144542395652 0.054638996613511868 0.77144176551401089 1.0524650841269019 -3.588291998265055 0.71370828137793818 0.65938791954880671 1.9206359400193076;0.80043076166388338 0.89061195748957878 0.23364151546834647 -0.34781735278846804 -2.9778291357985411 -0.57395446645442172 0.38999789028491161 -0.68768897401176454;0.020062269186950855 0.66356332793715866 0.18064625780082763 -0.5036172606062258 -3.782922837348452 0.74940640569462458 1.0252542518438514 -0.088000813255171104;0.26081081871933787 -0.57049941158301543 0.75293939512110164 -0.38390060871891696 4.8023686090150752 1.0990598551112072 2.1707757745496536 -3.9879749538204901;0.25886019708014119 0.63961511642730484 -0.47926428072617833 -0.51132014745243048 -0.96114489917584633 -0.43252263829371423 0.34463651237132598 -0.79104856667037637;-0.21634303722089937 0.65990168110074077 -0.40431342689866134 -0.55547660705074098 -2.073740312857784 0.30363937240147532 0.65664573959230788 -0.6453643727854429;4.9635670735531372 -3.9439819229896571 -2.8790803198649328 1.0353752867395531 2.7462434997485179 0.98429560620150047 -0.57904629185556411 4.3425931720833129;-0.62651996955565958 -0.83465773317890812 0.053547704895017895 -0.12207430415791895 0.098433258027237597 -0.76439097281413582 1.2251511974124396 -2.0463618301685385;-0.20403201494110162 -0.20777125765248919 0.028269441868293426 0.051115322625626763 -4.1981082497976496 -0.028997434842586996 0.01327580899599274 -0.11739393141889037;-0.13583630955271722 -1.8471697656298744 -2.1370790925011502 2.892449999534179 -0.14902628852625419 3.4115377924275001 -1.6087070313913758 -0.18215920099379701;-0.69107203535172246 0.33015197629186216 0.159656320300057 0.81079944492353606 3.1685858895936034 -0.48292207844663082 -0.61587937476014054 0.5562170725545782;-0.77420646991716979 -0.54308014205143385 -1.1625348399903437 -0.15549381104422999 5.160012932238879 1.1861142102145124 -1.7314157859849573 0.95554259790166252;-0.47864142647040403 0.78513758271374989 -0.8768803531190843 -0.81673936842287265 -5.5806081841056594 1.1946891614221895 -0.12963156655424182 -1.5878181247300112;-6.7760597651279415 16.249194696358565 7.5767112428660361 7.8502000365705591 -3.581643882176607 10.416775267129687 8.9847290779764943 6.5716548858458372;0.43987500932263951 0.21266321672090996 0.96066280606220722 -0.96486683963293796 3.6920894273007243 -1.5764744591719886 0.50407975731356869 0.28057761382905044;0.087430066299196724 0.19361549589302754 -0.068408926699979236 0.030062366238198676 3.8569481450712626 0.27984913738550315 -0.057866214841511845 0.13705548983156407;-0.1601294152636181 -0.64011326756422826 0.15151201036455542 -0.36223809646391009 -1.0810706871131306 -0.49519459680917921 -0.38232860603485458 0.098415403147678257;-0.10679836263376574 1.34908732811639 -1.5710010359518238 -4.2454707649394159 -2.4097796233688253 0.8170904963142972 -0.65566926467033271 -1.8906637294858333;-3.0515951252384768 10.582455404081749 -2.5312276124073811 -16.173682881838886 -3.593621043436626 -3.0338964056295841 0.46276019558933296 0.66840635899320922];

% Layer 2
b2 = [-3.254424831509001;0.66830939265162959;-1.3666512993996451;3.094515557124212;-3.4980977852962334;-1.6310049510210267;-6.6658455163881154;1.3994738135906442;2.8858913143594216;1.1179925936682347];
LW2_1 = [-3.760383921395706 -1.6958537878032616 -3.4433812639399992 -3.9634419300454047 -3.1600334933310492 -0.015915954646605779 -7.0672837676692373 -12.113210215926225 1.5463986359125359 7.2112206441861417 11.89747442718166 1.9854995915368021 9.2797299729511309 -2.590511192717186 4.7304175323700246 -1.3540649204447976 4.3506001280190478 6.9393737465946552 2.7428617454125597 9.7286869907892797 0.25449473818927532 -0.09547943336326481 8.2088591574124905 1.0654611965251362 -3.6295003627533684 -2.3942128840710741 -5.1369905228842567 -0.62261898836565222 -6.6362710424530373 -12.01935724513322;5.3519383755575278 -0.91973124116202187 7.0481680783220728 6.1682257774889981 0.80577970723494707 3.0884765837155346 6.2520909390294399 -1.179507430196733 3.0377740790467427 1.1397543974864519 10.35882563607457 1.7767091817830012 1.7604109812283293 4.4433900782451925 -7.4661526419495781 9.0988384282636705 -3.0330193459752812 1.0221990172707525 -1.2392564521619123 -5.7658848461348766 3.9220400358365466 6.5053720940961313 -0.7324393381060833 -2.3743308720837417 -5.7054604577928361 -6.7903378809076091 -1.7068953167516656 2.3398699301836539 -0.47779486625810841 2.284689520942103;8.9401753061741172 -0.0088738342388909128 9.678907390453956 3.8961198736156772 -1.7365125881156271 2.3145948512651997 4.6501864655857537 7.9723596492362381 0.94115904986730925 -2.7575227230068311 0.16601519573533599 0.49878507322125715 -2.9394509161902573 -1.5653405749110538 0.024473378345964705 1.396143531882464 8.1096385068599588 0.1268596811629929 -0.68353657492770492 -8.8053045512322132 -0.30455594212984777 3.0330433641584245 0.59968815514156515 1.3692720616340186 4.9661086359470881 -3.3710048001572614 -6.9891732993248796 -1.6115900442321343 -0.3247828715024309 -0.0059505431531943391;-0.14122480264970302 -0.22126995684034434 -0.66238021067208752 -0.55704136690294959 -0.79882868582941513 -2.9633185414449721 0.87244608194073958 0.43170545066638816 0.50909720627210475 1.8079375564979474 0.17145906968664765 -1.1004973573656009 -1.1600644362513075 -1.5522568733734621 -0.68794865498170177 1.7137215143695428 0.72917482151451651 0.3687051500824004 0.016800440276460354 2.0395650162515033 0.44793515453685345 -0.050892970008636192 -0.80785475745385404 -0.39881652991269895 -4.7411980954405388 -1.3720081695812438 5.1222588063930532 -1.7868670595706879 -1.4418737704290276 0.064608895145298412;0.33376825287272693 1.4925390929134914 0.2573576961638317 -3.0279300664041173 0.63775066834817173 0.90252804813732956 2.4979333150043481 0.2762005793115011 0.14294107112707458 3.5568288284551599 0.06318443061806589 0.38728024150761958 2.5259318869179856 0.90942184678076254 -0.2942231418986771 3.5545721698582389 -3.0329660741847944 0.26145633749198804 2.5240757096380273 0.71983790821179927 -0.0099276720101007243 2.6963983760272345 1.278093668694773 0.86358615933248184 -0.045106886969443184 -1.3998188095647071 -0.58774537800981164 3.4465068472191356 0.85772146540682348 -0.27366176433807671;-3.4782927894948354 -0.12342035192037298 -5.3349904139800062 0.21578989470329774 2.6585253875764367 4.5203884584959848 3.4411769289128031 -1.4628869724698359 -2.9584551495771847 4.5642239303940393 -3.5382276358147293 2.5312752936294363 3.454241360115176 -0.85386096803727407 0.00052222666442094945 -2.5036361772853621 -4.1780386746316314 -1.3695724957499387 3.0589278246265033 -0.18017248436210159 -2.1267247814711903 0.69161475865539679 3.0610037017147564 -0.95348896818241347 -6.7013024269595363 1.5158517621565437 1.6553243710522805 -3.05463141985352 4.0221497739945162 -0.24492652689215788;-1.2801314262886694 -1.1724116652170948 -1.3291755473669333 -3.5673694398567162 0.74256517719588 7.2723399167598135 -0.6600012644425054 8.4086907278545553 1.3043112686433544 5.0858847269900247 2.1115197212694081 -0.4225736006178078 0.84627917670185193 -4.7874701157986523 3.4227943305447632 6.4042218798204678 0.52618778406615796 1.8277416410121894 -3.4279656013989137 2.1951009533264996 1.2367635601650533 -1.24882776881783 1.1054184954282795 -0.59076923317628383 4.954446218261447 -1.2972312681244342 -0.59596183524733193 3.0137315026097524 -0.39968214031834071 3.0619034862282897;1.09523959194168 0.99516642511702824 3.6117132896809414 0.39109277971620304 1.8846294812187254 -0.81051534286815119 4.5272725372092459 1.726315105883627 1.3650867277328551 -6.9992950356069521 0.12840217432442233 0.49242271166399254 1.8091905347928161 -0.35707685873217593 0.59511739478188708 2.6636471331499267 -1.4661366248767018 0.20818888327717092 0.80387244658682544 -2.4230039204045415 -0.35353507674745321 0.19819929125506386 0.54099588409233734 -0.56416231171815445 3.9844599140928278 -1.9483622398265332 0.40329732208859459 1.7931951738584155 -0.67778170844807573 -0.018105878294056337;-3.3182114106893486 6.5334918073900123 0.48713845136823919 -1.324478226187551 2.7709784905572445 -1.7400995235456576 -12.089781374339189 1.1914772058220191 1.5874822124758119 0.86688429710615678 -0.84279500194654045 1.3222103195642527 -3.6788665680022601 -4.1615579165312102 0.75698118349715904 -0.16698551105644338 -3.1486395308066193 1.5550924321055268 -2.9350228379584742 4.4547755267907885 -3.4734659547027733 -2.0009007246655255 -2.6182762359189007 16.858499187821646 -5.5927574309269668 3.2500272247731088 2.2913939293041095 -5.721047660172216 0.2799895971463815 3.7473660917726863;-2.3599822245753042 4.8522590438737039 -2.391933434020256 1.0131862115415757 -7.0086087185283716 0.77195413918469646 2.8573315755979389 -2.6041005090725453 4.5173649250323953 -2.1239116004419896 0.64587608531445762 1.9740285709523158 2.7447496112326935 2.2876940565216723 -1.048849044564923 5.056562742596272 0.98448218364335438 0.18410625272936737 -1.0648418859411775 1.8227496814981214 4.1501541388837033 -3.1382060935994631 -1.4528215954861841 0.99414525638256102 2.7494308678711334 2.9768527097217015 0.84526187857935076 1.642882516485362 2.3645631130529483 10.315849315859579];

% Layer 3
b3 = -0.48350971139807059;
LW3_2 = [-0.00089530559414679946 0.00153242328799493 0.44662923559334572 0.10998037238964174 -0.80653430509388546 0.0077803109769900534 -0.0052012399333597808 -0.42970953679755869 -0.0034283493272640356 -0.0041348333127444042];

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
