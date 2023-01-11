function [Y,Xf,Af] = ESPER_pH_15_Other_2(X,~,~)
%ESPER_PH_15_OTHER_2 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:25.
% 
% [Y] = ESPER_pH_15_Other_2(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999540494077341;-75.751;0;29.715;-133.803853032615];
x1_step1.gain = [1.0000770825255;1.00022983623563;0.0140762793578401;0.000302617642608564;0.251952632905014;0.00432440760873659];
x1_step1.ymin = -1;

% Layer 1
b1 = [2.1996038398816106785;1.9465331802902454683;1.6138120500392660883;2.0469010602667920828;3.9314028194757364432;0.61237380630569315798;-1.3248733951449012114;-3.1462669542871832817;-0.50851805082221535059;0.4302988629572163104;0.93202840063634828649;-1.2993945220909455252;-0.54544811728165698561;-0.48338467777866500619;0.33134432904136168485;-1.9785720609359320754;1.0043350729009148381;-3.2619188466679935878;2.3549190069953196414;1.5961865399564449941];
IW1_1 = [-0.22845307042056239855 -0.25113090573431129382 -0.68938007731478045947 1.6969485524897329753 -1.3844167014824304562 -0.27669105489278705834;-1.5068718819151230015 -0.86572622406420718377 2.0239977357735936536 0.97248240571698651991 2.2006002110084335044 0.98102506468498440473;-2.3052798380867076666 -1.436225680787076131 1.6386748550395071611 -0.12174374618048733332 0.53637758776827382246 -1.8657925881960986114;-0.56648611529004599952 2.2010449225266683904 -1.3580347850920835295 -0.26932283468115003622 -2.1911817136907356662 -0.095841332801224587978;0.42886670237499840974 -0.041334240546651933312 -0.7599535589997137297 1.6173922507591413922 -0.47086967540559532264 2.6076539178089648274;-0.48679100341202008462 -0.25867792895356328531 -0.45714895919216158582 -0.8130596141356902784 1.4318666508328670162 -1.2755151082617894875;1.1871023681729278021 0.27063844556425487031 1.3796814184470871378 0.40075718361204737406 2.017328323089685238 1.2651299371962780871;1.2353191006168169341 -0.041856081036871717327 0.35297038952360210162 -2.5747067721574783405 0.2252416100276845301 -3.12560662027546865;0.49949752291189158937 0.5684584408352849616 -0.79655482211352335309 -0.48933040125706894763 -1.1074470175788149184 0.56966571314965597228;-1.254110232814231507 -0.14128800442030062245 -2.7848225294955866183 -0.90865376920400209304 -0.29101744394215256495 -0.1011527949224242523;1.1688478942080573475 0.26514786182368715872 0.82073477244358006288 -0.21911269849507006291 0.87240201538832551353 1.6231654558713064507;-0.0070241116415155338201 0.0079451962622573129441 1.5656338165051360001 -2.2246374460504996229 1.5411131881805248955 -1.6528747984156115791;-0.22320347351332117669 -0.42886498789340415927 -0.78458222029938429198 -0.17519738498251660741 -0.34815489804035965049 1.3322143126181724693;-0.36882218989931714326 -0.50802472519410823981 -0.17155016494451513775 0.15062614536155766176 0.65695746871484239637 0.46128718909260113534;-1.0478787164223257733 1.1579491173714229202 0.95546547024983374818 -0.85697087590623655107 0.71931183540000975807 0.7402227820996643004;0.10312392016892492352 0.98266230126091036468 -2.1017899265948241094 -0.87429365390625768395 -0.71296835122073121127 -0.29923708037394763926;-0.053618075112706109486 0.48373984427342320291 -2.0079778956455682781 0.44405619013469777023 2.1516570316980900834 0.15378842172881532213;-1.1189701211265583591 0.032551984732682927304 1.0127091256935434327 -1.249803024047940303 -0.37616866508499868793 -0.73523825091847516688;-0.18219382700617997584 -0.19941955211916684299 2.3460128756662177452 0.64931767772460768917 -2.2806489324740382152 1.8430655439588505828;0.35183960307626366282 0.23486324091188368057 -0.50904916127748089938 -0.51991887357744626108 -1.1839716478487356532 1.8809319000367221708];

% Layer 2
b2 = [3.0795808219279305895;-1.1478578342402288026;-2.2558080116177650076;-5.1191572466530717023;0.10043096900892785339;0.79927655276642062265;-0.44762436896777274553;0.5232604912465007585;0.58046000773096040337;-0.14269522229249206879;1.4760847269854144859;1.1188137166123448552;0.21826903382795054265;0.7348354572126947204;-0.41579537139633138132;-0.75890847280529094565;-3.4361472934321506401;5.3298823157469588097;1.3748385866206764838;2.4808544180294695991];
LW2_1 = [-0.85700703806482037006 -0.71349286407238665664 -0.56002595902474516798 0.092684782764036646729 -1.2087526990814596584 -0.30630978397105756983 -0.39919725984433235277 0.021831999498487485167 -1.2744538371413438593 0.095069622048803370529 0.75712274683588209889 -0.4156108263148883819 -0.88567394327767634543 0.22670250099448585601 -0.078047928472665145661 1.8864863351924980872 -2.0726646455950605841 -0.214887097462313853 1.4343815597047613331 0.41959200214140629503;-0.40514089867783759091 -0.10324017894204018386 -0.67271642044715473929 0.84981810366219889818 -2.1760277087135087015 -1.0282474738750366772 0.74712481329456392753 0.91132256630357433291 0.50835440558295252877 1.2540019516851215187 0.09568351010999967543 -0.88150534664719970568 -1.0303648309413544482 1.5329414305463890944 0.34655469821635387362 0.099612917693117128404 -0.77080954496361075989 -2.5615910978135536524 -0.47816283188322761255 -0.81390991810662793871;0.1054504661055610526 2.6155814691305350195 -2.0316767752853017726 -1.346372074622334658 -0.055734015087389982002 4.056277160748609667 0.66502620294291492353 -0.38436964769449766255 4.8965036380813673489 -2.0171649392051578253 0.37283070991538092009 -1.5206893984877720705 -0.0076861029790530263628 -2.2421832106766235171 -1.6443450864020390245 0.27964375834981902713 -4.2966605328250482287 3.5722980704043072464 -0.70804091016344761211 1.7903180631234183462;0.98534762727039859609 0.55118474553006158168 0.9462998125453477849 -0.42042548786977906872 1.6273269347284899933 0.47549276520632210996 -0.12389859110191592251 -0.6136182236102223353 0.96475561241961937942 -1.0280214905231068734 -0.49444170351349886472 0.85585723307438343888 0.90284115935706432676 -0.16502029906298992801 -0.46345982334611240949 -1.3483278536838576578 2.5761992745638879754 -0.3684324339604692633 -1.405341287993434829 0.41631751680029088769;2.6268439752620484917 1.1460474944174012801 0.44441657973697273887 0.35517543292519093701 2.970184775380190878 -0.30566454684996069435 0.54882328498220145363 1.0151341664604014436 0.078609508532180821794 -0.77988981574065086022 0.66437880535596471887 1.5468757189564799059 0.91093479104355412623 0.96216960127603268216 -0.12357013069418296136 0.56897029474646865221 0.5286620866518647599 3.1392484263510831433 -0.26352034098863436062 -2.1557806932968630598;-0.94066245717540009075 0.45107354194147769588 -0.47416600008406506461 -0.49727629694909847302 -1.7000153804541748848 -0.52947182913468271437 0.01629310856543208022 0.1139293106771301356 0.15045622754616908012 0.89713651548334816255 0.12597243626945853201 -0.68278438519664241824 -0.65279059709619546847 -0.53089335504292578527 0.098864613475346696525 -0.30625066110137044717 -0.28993822172386352065 -2.5683888728829589176 0.18037854742052861212 0.36011150974480898856;0.57774671335920579907 0.77602589278024436403 0.016790364984228194067 -0.068148489023167876932 -0.71360151321840492589 0.36135242796088529182 0.047538671276841887314 -0.29078986795344780036 -0.17772533130920722955 -0.9214502601671235027 -0.79086415234904305915 0.11072878167737656985 0.73006933543909913187 0.093286316805737262703 0.076644771650654613948 -0.51291451596549808567 0.58968706454816466334 -0.050124311654804153127 -1.4623004768413554455 1.5620947820133463946;0.27502485196315840099 1.0466307176140998703 0.01633536065389966721 0.38673257418911866079 -0.98322503977906805517 1.3741492406288984274 -0.16432471979633836434 -0.82900147320409578455 -0.52399930885522238988 -2.1389484484026448641 -1.8078276755557398658 -0.31627297181982011454 0.099321736362651832919 0.36247924174128720587 -0.41946665690916207225 0.10491235351753946947 0.090222518996686540405 0.1727364084008334022 -1.7807092345586694737 2.5065909988905055172;-0.56337264494536476445 1.7807595974472469891 -0.67129001548097999574 -0.34594982844353189266 0.13295703303216768743 1.0137505874623786717 -2.603842718490009478 0.0063840492144035219935 0.90580191161041823111 -0.81437601969755168962 -0.23679740816306760198 0.2772350093250109393 -1.2455726127547595006 -2.2457612705040284062 0.68790402735312938809 0.3915344091739200727 -1.6010776342375854231 -0.83455044011136836346 0.81025365582580899382 -0.66875461796521595836;-0.60886841621492582277 0.63246403788793414957 -0.062291161602914063566 -0.45537657268248493692 1.9509555030270724352 -0.34042113741384621983 -0.89525105335504939763 -0.031218562406719094215 0.13482375323700257486 -0.56178617647765050069 0.077778464727981277971 -0.058929903986549410366 1.0781645093236149613 -0.55660323931081956417 0.50568099798819587942 -0.5601397846197728736 0.36498815501307912967 1.8966651319735017278 0.042047894639021157759 -2.1495021671254117557;-0.78597969957327984414 -0.41196886363448176338 0.2140806733264192796 0.15150009541841710781 -2.1025747188408576704 0.14107234232454779477 1.0734079276324726138 0.41299442647767004022 1.7539645021211718312 1.2705747333290438039 -0.22014148889215315319 -0.66577969377741919033 -2.2758333637163090657 1.7336414216842839231 -0.21416784756540147017 0.031486015252021722077 -0.19522149110782560655 -0.77068402299727101212 1.3740016261692808985 0.40029968896101603137;1.8281879120028650121 -0.008970903707493930318 0.28018398045688669429 0.28644217449087527871 -1.6398981241494412497 -0.70595720087596802728 0.52964917226070895584 0.034041420735663720232 0.98238039731978898939 -0.25619398722000574775 0.46533483007803194687 1.4361228249347166575 -0.77108871976666670811 2.4183950580801871766 -0.2402649443923544903 0.34997919580778469317 -0.37684079326311226144 -0.82202699907240173172 0.11005793325554312501 -0.45520934591092676591;-0.54751367801194572582 -0.86348666051828015799 -0.017399461568470692141 -0.026214010390494670716 1.7860973907520496518 -0.52717072581856949665 -0.0087605690771105203946 0.42287279449651166141 0.32317173491504919136 1.2337025150895910564 0.96067938477914249784 -0.039280274463340890401 -0.51022402081870099622 -0.35695941259620961317 0.06582181521045854522 0.30698042974306688491 -0.57147281901396917991 0.50588173221171961469 1.5014898716742781648 -2.3192683615955269616;-0.93169287302262038253 0.034591135813468808224 -0.10478074603811822485 -0.56634994525186010961 0.70140688150017282965 -0.77582484006970986457 -1.2600282591292382151 0.069474668729175376169 1.4833132840406788588 -1.0077677747716102807 -0.32361194360635248124 -0.5730821881559960973 0.034336663653420654296 0.95862120867201294772 0.33702140655491663601 -0.35449591716477185654 0.36840335770539328619 -0.049541464196226840666 -0.15089104568854883048 -0.72405195602371774655;-1.1382682595663826675 0.17086392699369246095 -0.22632998775198065466 -0.53619350732912141044 0.25217042859129584764 -0.60796292992807610123 -1.2701699245407087613 -0.015692687531116741168 1.0923167074857340619 -1.0785483186427147206 -0.39125054039348455825 -0.66317367103788282812 0.16892476636574499382 0.60553832058983969588 0.24551676887180789999 -0.10954236822962379327 0.27909952811074734136 -1.745972462711595119 0.073810957246519787178 -0.91854336274366232828;0.6331111290471712616 -0.1656804703462801398 0.13741544400162089556 0.42363393989548348006 -1.540590985102284316 0.6067130193409076 1.2902743130000855487 0.19414631224391362685 0.36854293831896495215 1.1031964063869454939 -0.50591921086102209504 -0.080750164033943069852 -1.2362906942433908863 1.402224960996548031 -0.71766900320854110085 -0.054403724944844300548 -0.58536677275962201161 -1.4896635254655916647 -0.23960562607096683174 2.8433696402179418072;1.2974461482914412258 1.0496936070580034084 0.39484383980805215497 -0.14893620939643656409 1.8480311001172107233 0.37097865682190039882 0.25657671907223639218 1.2178349582869070211 -0.75387967427505497842 0.42277621251627467425 0.15549461280535181729 0.55730847924604154375 1.9572612190327185466 -1.3815265558048266392 0.40341811635758773047 -0.23534033257677550166 1.407307321194885219 1.9509287342758145911 -0.53858511268503972413 -0.26023458038137942117;1.2305223104493923092 -0.12245024286750717801 0.35626495432653565443 0.31523195581307084101 -0.3169065186911713683 -3.4026060389516512572 0.33216870256849767573 -1.2437162385113273988 -0.86667106227484358438 -0.28648447298664164729 1.8848412656795507925 0.68630059921813479917 0.54027818428903773285 -0.76432416867010510853 -0.77782417328778952204 0.82779279118032333962 0.11026648025159029032 -1.3268808017269178023 0.56561551451513480959 -1.9994926643107915343;-0.34403341006737625518 0.33775653823840617873 -0.41503495682620256124 -0.10261926063594647984 -0.76300402555155333317 0.98847497466225631157 0.7314764850302176713 0.12719148609792860216 -1.4993928548800952161 0.97191642152171753288 0.39916412132179424077 -0.46287968859408551436 -0.11509604255289371255 -1.0575936659861751021 -0.21696869276496519308 0.48084865306813856511 -0.22766557985156365973 -0.45946351348063152464 0.31020804726864908973 -1.5224169448641016178;0.73115573528758204258 0.38652111633066388308 0.48506303578854481628 0.25909101414578400435 0.16440702592350797673 -1.0807980879772653982 -0.81581683811822647279 -1.2988576185718740952 0.87943217379902260955 -1.1953146701716232947 0.45342558059596776987 0.92721777075511691013 -0.15902870680918351765 1.3420720359562332558 -0.41258649390875373753 -0.21555306907910579128 1.3326923891233561914 -0.50777635541092569493 0.42659739134710727981 -0.94303305850005225697];

% Layer 3
b3 = -2.35568642714864529;
LW3_2 = [0.62873000675496648615 -2.4408935931939641506 0.34063412070229642792 0.39702592478531739673 -0.81187455080858872414 -1.2211823251320408357 2.8151692388538531731 0.37849739089312739448 0.11517158632868697288 -1.4179017794819532838 0.46462375548536249958 0.69512950466311285158 2.5025442455219422122 2.8930082507797045821 -2.7983652471662381345 -0.80789956186885858447 0.59626650920961110636 2.4719165111633625287 1.0579296665899362573 -2.27946947646708864];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 1.90239167587468;
y1_step1.xoffset = 7.36300175998845;

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
