function [Y,Xf,Af] = ESPER_phosphate_8_Atl_3(X,~,~)
%ESPER_PHOSPHATE_8_ATL_3 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:29.
% 
% [Y] = ESPER_phosphate_8_Atl_3(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999994516886;-40;0;9.1619;-1.99149636107278];
x1_step1.gain = [1.0000000017546;1.00003626070655;0.0153846153846154;0.000309885342423303;0.0656325902039552;0.0604107502911417];
x1_step1.ymin = -1;

% Layer 1
b1 = [-2.8898181134487788135;2.1010925455429436681;-2.9643496047748616462;-3.1806480807413470835;-4.0992017265543916338;-0.11188202191436813493;-2.2767161994380860079;1.0745603982992677494;1.8582779992242399736;2.7434280034695666117;-3.8001642100078543152;3.3547634075658439023;-10.072445892650183197;0.38943605895630473146;0.52643528808223660764;-0.3792192879947032802;0.80223842452426574745;-3.0405151583047151931;0.87964183286360786429;4.521554982833134062;2.6357382276433862778;-0.021811813182745061479;2.5030764535627669787;0.430875373149194274;-3.774070280524109311];
IW1_1 = [0.90261150941440693174 1.3929985653225998021 0.68704194369744930881 1.373070518255036454 0.22276720815537137677 -1.2611383925795838223;-0.94625189414592547621 1.3772203519641130764 4.7039700866800302492 -0.63797773457496875782 2.4677415274067651119 0.78837793395100397298;2.488083102301074323 -1.6266713728634567548 -2.2940108139465205106 0.4482870994302266543 1.3661872664743306416 -0.61792394709153641852;0.26779994070636947701 0.410293730910901433 0.30497750652359933987 -2.3186064591078778818 0.51885605461887696066 -0.69138335668609118834;-0.245427093493424342 -0.26607298693572045689 2.2741267905528026461 0.126264214250986051 4.3860705709521274898 1.2653748575766228246;-0.02455254405433970194 0.20173768760405605405 -1.0504554276031081717 -2.0919283404247077662 -3.0743619844193963431 1.3405582228983883297;-0.58246688819295322492 0.88161311699778299733 -3.427799466643591586 -4.375030885542739334 -0.84745481348710405278 -0.55276454120045026563;0.022703618068380244238 -0.70619377078105805001 -2.1506297236043936394 0.4423556041361398905 -1.9209366455865513057 1.1179468270014192299;-0.057875291508445954325 -1.4631195225744995092 1.3623719295696712361 -0.7446687292576811279 -2.930118632954558322 -1.2373423188276249096;-0.95271901432759897865 -1.2232002319556050107 1.3485728024696981286 -1.4058382768938839646 -4.8970408278079204933 3.5339606710444404492;1.0461945088742883314 0.41981698123955912827 1.8623040085133208965 -0.92639129117526686219 2.6137437721003387026 0.050127921066925713867;-0.094105135588117311274 -0.70086249880784268207 0.25876935135063550142 -0.28304528594490308002 -4.2131026789938115229 0.91357160431093409336;-0.35037990921684819545 -0.05620914119825445171 1.8070805290297968781 -10.958107244448969553 -1.1845858415933905317 2.3040302610996969968;-0.52560579889758707317 -0.069551708537029896529 -2.4559141820966234881 0.074494920201457939801 1.5486021499333788576 -0.2362776415861486401;-0.30021536223957695499 0.66640484753992301759 1.7277183390778356742 -0.25738442629819330287 -2.092311487699490602 0.04643400836790764763;4.4982790232947191811 0.49744781385888314773 1.1603673786653792011 -0.84816952237964515771 1.9403009242878894725 -0.86804671659513321647;1.1035031018465459862 -0.41629145004473105951 0.90397038412023333898 1.1695028795539854638 -0.85737480069061777765 1.2220843977848532269;-0.67314029562035826437 -1.3864263444964997074 0.96458153775431143551 -0.21468480978096696199 2.2471245867118834205 -0.30474361954141510056;0.94164079508135467567 -2.0244559747449724085 -2.0560757050909383992 -0.15356307768850188666 1.6932743544091151122 -1.610091485156554425;0.25856743675928850479 -0.37459678751880226555 -0.1105074003242291214 -0.045713475499595875029 -7.412562104925059181 0.027196252948885921846;1.1369778408194599617 0.67742165736060233616 -1.7759377114163317213 0.28401176177971515502 -0.35968052468252814657 1.9559478663941332055;0.34805406987617698578 -0.37372484947923589882 0.20986769836994689786 0.80250127540024662753 1.1740333465301744376 1.0614950287020115116;0.67370961039834964534 -0.65342694409619628271 0.031176740109455515493 0.49898247990389732465 -2.0022203620553837844 -0.95496788911285490631;-0.59354666055481586362 0.049358044636380445092 -0.18181170054837192396 -1.1028365165372111534 1.2295374999613171152 1.2943520052815749022;-0.99953951393839779449 -0.42501861191217521485 -1.3611836310495608604 -0.22249110384415712871 3.6363239077315152592 -1.773328278586373985];

% Layer 2
b2 = [0.53098112063813041051;0.032520882680570122047;-1.298694791607694965;0.64745958850656748496;0.26516817420532518002;1.0243752867583233446;-1.7830375697518581202;-3.1280295527060695271;-1.482023786821901723;1.0387614647090119835;-1.4333917048496995772;5.3927997129191176384;3.3333213145988369597;-2.9306596193677787454;0.067890593323398348846];
LW2_1 = [0.043568298470256811017 -1.00327618352245862 0.59624005956667547323 0.35241249295853749324 0.024928867795749279357 -0.45563063123745439764 0.39413744175787357449 -1.9894140110161069313 -0.56430420649324919413 -0.4958352329298553185 0.62199483232307928304 0.10593749283373145309 0.45668117354445464384 -0.14251697125386936227 0.011759065671024668959 -0.13923815847116349786 0.90584945819555118263 0.53355407164277157683 0.11749610513950425883 0.64293239323292761256 0.06355993918408475718 -0.56736562210231500369 -0.38866854639876957656 -0.10067812168626438796 1.0093911261682746439;-1.1289607211417513888 0.8372484145793904764 -0.15151157197623812967 -1.9209975147327182921 -0.46350636036151132702 1.4811020015419542606 0.050420290554267416239 1.8821779410241927799 0.38796847980328613703 -0.45612435534065459697 -0.71711508306038129756 -0.17402534235509953264 0.2004322588468743227 -1.5019744875023930408 0.51998887324931108989 -0.11325267939215767632 -0.8848829081174025557 -0.40987446385294162843 0.032354706656208973015 -1.0957305295265147294 0.027065206050347800221 -0.079231646255440657534 0.044636497293466803615 0.7935687360357756015 -1.2331372654029846281;0.14568670169909195278 -1.3692424328407408751 -0.86102252033073833903 0.17401948161211683308 -0.040123030620434617377 1.0544102657639813359 -0.66085333016933800465 -0.84864201890801560513 -0.74323077837845252791 0.07975195215662861381 0.10404681056631298608 -0.52572531550060297434 -0.28153962811323185456 1.6216717731348924136 1.0943594445396498394 -0.14404004006690293571 -0.25060174263842627251 -0.59006517120270984833 2.1045311989926522855 1.2768087492227728763 0.1627859110842685697 2.9389832085290956698 1.1375338642322825322 0.094405710118458288638 0.21304081494833768118;-1.4455644735355328567 -0.092141979463524814609 0.724900238670998176 -1.0446147808200070273 0.39109937479441692076 1.1195227039636526101 -0.331511174713301493 -0.056741897741136519684 -0.034854768365359292681 0.54092230721015910344 0.94833992359365881786 -0.12781412152417800931 -0.67458727632235016003 1.8142765222056624896 3.3725412269839503132 -0.65189982480412056454 -1.7466515625178216098 0.014967866575968897386 -0.095240309074541201761 -0.040820164958974697322 -0.16644614613849997853 1.1413612180010548336 0.88975044306777140157 -0.66017296099150435396 0.26023066193154764658;0.61502648805574022806 0.62400100591118301541 -0.70617572261909311582 -0.26386543901874659479 0.32350969720310229771 2.4679097831591390033 0.61205541461175161366 0.17923489923453631012 0.75279977180094137523 0.51386483862805254486 -0.75838851001967622079 -1.9052710031030715143 -0.20983712773316753886 -1.7931593319544583132 -1.9728390474987749226 -0.29831526494637522484 0.039815916765104993058 -1.1310876789611423732 1.1903395414239601457 1.5659130393725866348 -0.77320064077369943245 1.0299668534106527495 0.28751488062830710124 1.6421102524851400872 -0.10380070079372819003;-0.86988903979754028661 -0.42971211261908626167 -0.2668602374500086416 0.79294071113465691703 -0.0016276412906544154499 0.030780278136225931035 0.099723671442833144996 -0.046376932789285986447 0.44013961436573534325 0.20637476082167738323 0.18511707325043474937 -0.1483667873092308076 -0.30760320965885640243 1.052237963938991383 1.7680188248328436718 0.33146557028401218981 -0.56868193721198834201 -0.73643391558217785597 0.62017185763601390125 -0.05612976045206638509 -0.62328731493010547382 0.9352184608237986696 -0.53093920391680493598 0.11868124187679500448 -0.25316777663638001528;1.5386831693589664294 0.87217549679752659486 2.6158714340073867533 1.7666676793221478547 1.1199809456276847985 0.017800341887862409368 -0.15010369128020298324 2.6698726201892273302 0.10518249717257291409 -0.69410202985239044526 -3.4563023067328972893 -1.1835584396996676748 0.44909303092431285265 0.32739308469445860261 0.6607906288904742409 -1.3983139706838270122 1.0368090904515701034 0.28145684183877561768 0.10917888667835679017 0.97095813807963959885 -2.5509318238830949177 -0.041152081495988507287 0.7415212409997813392 1.6390295828149532564 0.50117513326521334616;-0.81609674295987144621 -0.10472794268517621186 -0.10333115927936431522 0.43053956966131662343 1.1636846222886050306 -0.77370431214031010203 0.93695541059125808658 -0.029768221151148684805 0.76426159129425919048 -0.0010933747011610366051 -2.3802515692905874545 -2.5565970459450366725 -0.31988222129170768016 -1.4289205040402064295 -0.82486752559775256 0.10202004522128993447 -0.55840806652023333712 -2.1353377351576554588 1.489743439051253393 0.53318295057318432217 0.69569105883433957072 1.0129834253467739913 -0.13299878036376139301 0.58214857861619218049 1.3037824451942972015;0.58348819365865167708 -0.68839940300489699254 0.92277049640804331698 -1.8823228971479157856 -1.0291587996119637705 1.0666926088255850313 -0.71214002267672660285 -0.23547380919116742137 -0.34832463455842549438 0.13272343030968530853 0.25856003400182242524 0.34743968704742300702 -1.100718420988063162 -1.1054590186532426976 0.011156485696626826376 -0.98410779246040980262 0.67727175600072797579 1.3290023507173847062 1.4912141646891632352 0.66657337272231453884 0.84194869391353166499 -0.0019978103081874747093 -0.84303341579003721051 -0.028037832389065178024 0.061462395814032906582;0.75738256395964398404 0.17718237731841268712 0.67166090485450091752 0.57973587965358275298 -0.53165273484478459132 -1.1763036964647792804 -0.3739251824890525433 -0.43528738181751935876 0.54867469364447307711 -0.028213542650671899481 -1.7972427566595841064 -1.6542717820986332367 0.18509338983874148998 -3.1380818929113023508 -1.8672714317875698242 0.56685755475858723429 0.92954295229665595546 0.14810609106946945324 -0.60974027499955762988 -0.97613822044332820838 -0.3458382497213493556 0.26285873783699975625 -2.2000038425336811976 1.3507208014960077058 0.92791582388687576799;0.72921961699251325317 1.6122779095717922893 0.012641402812881787424 0.16266025944483644694 0.07059127103605875031 -1.9464104098493510975 0.12343736415906204906 1.0211241643594701767 -0.57531711642077776325 -0.5483160804700155655 -0.35491958908808673989 0.20056600381530756416 0.75969904466126236819 -1.7367188188895323897 -1.7435295407004913759 0.2211197913534135584 1.1364912156653967923 0.33015658111514922801 0.05924349007218185359 0.19632699373107703744 -0.025019451540663543293 -0.69581479064354456732 -0.39849285664399369811 0.43291381479594137627 0.11516270732846049796;0.044073792859038403436 -0.037966753566031616329 -0.32989780462599771482 -0.32759896840841168375 -1.2919794639828126037 2.5304287149077921981 0.43045874826523988776 -1.1434244055092492864 0.011752279603539106737 0.14851583254152467273 1.5845617537662297636 -1.5734763961990685388 0.30571867367257332671 1.7801452583926853368 1.9646387354445367102 -0.24381367125783839911 0.38105037499463068063 2.0792466121912887012 0.51083837711300883466 3.337479187671243519 0.83603892183522510173 2.1548262408413796365 -1.4791524250773881999 0.29660117677958591775 0.67546315736707529531;0.98422237589755645626 0.80336474993162576563 -1.3779501590531966215 1.3513937563940263686 0.34847569117145121842 -0.48296658708154621698 0.33097321570093996002 0.26545111811884664066 -0.25844922480152293653 0.01274358870448257311 0.60469078978822332626 0.84112972623778858861 0.27815917766984460391 0.8314969338871099902 -1.0031331756321901771 0.88518863692790605846 -0.10535162391490154876 -0.21704930267870103933 -1.168285411045383082 -0.48056394301622679111 -1.2530575640485457356 -0.31970390693802802629 0.92506993589686226098 -1.125156229555828169 0.14445146094139582615;-0.77545633723086793143 0.086183874285053863051 -0.014666157690463251656 -0.035800588979643603593 0.37833809065123225057 -1.0259057841148693235 -0.0064605202156753612289 0.73149178790526681038 -0.46372353927365267401 -0.018714358827906046184 0.57190686400009105306 1.4401647783445195028 -0.69838700839191925063 1.0250704596044835082 1.2011627353594416245 -0.095600750972296213992 -0.51593580490161672358 -0.42662444636796303943 0.98865578775935980538 -1.400323635507689346 -0.057700404168441243868 0.49574523061843245308 0.088576830276573495504 -1.1103856637510747252 -0.021376203321190254819;-0.24141029066733357689 0.70991201168082940232 -0.33733107707513754203 -0.11244453382252160567 0.82235185167895141767 1.7262379555061921899 -0.39081068555823000876 0.046408391121264339074 -0.08716014545707077843 0.14458069688549810006 1.2203977667177901001 0.11156120708428966803 0.36831563965010166717 2.4831273085285969771 1.4393100301476668967 0.01629487708358255682 -1.3890399412439606852 -0.51446186169691132672 -1.1071258036079345555 1.5432900063136911761 -0.014785196235622223571 0.45936710820891241136 2.4298056308083055832 -0.69622058941899622742 -0.27351675037093925447];

% Layer 3
b3 = -0.58710796502061535218;
LW3_2 = [-1.9457325366301150282 -1.811253902124920101 0.59873307542753906318 0.9509202707528839893 -0.86567517169015506528 4.1206496466346544949 4.272603400578077526 0.52289877501200454635 1.0172030027513336847 -1.6855111412431460494 1.4510203365104710205 -0.99746287451466419949 2.0774876501833974274 -1.8630643482579347658 -1.2173955449870585976];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.604229607250755;
y1_step1.xoffset = -0.03;

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
