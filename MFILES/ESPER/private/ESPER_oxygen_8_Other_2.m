function [Y,Xf,Af] = ESPER_oxygen_8_Other_2(X,~,~)
%ESPER_OXYGEN_8_OTHER_2 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:46.
% 
% [Y] = ESPER_oxygen_8_Other_2(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999998793562952;-77.961;0;29.715;-2.18999082656054];
x1_step1.gain = [1.00003298043486;1.00000060321889;0.0138606862425759;0.000266418009857466;0.193068829037552;0.0582073296585987];
x1_step1.ymin = -1;

% Layer 1
b1 = [14.163133803910556097;-2.2064071397582880607;3.3551907391486377463;2.1143535200840797295;2.3029062633486985234;4.0055060431653970454;5.7599802775649511233;0.44607925345843874032;-0.32185234292653780441;1.0309132394269118205;0.67709333416064332312;2.138723000393668805;-5.1711447656950415563;0.46782080235465928197;-2.0941635894270591578;0.71270512154668741456;1.5712836505752953986;-14.834653786725320757;-2.6999872964266491415;0.978320332204120513];
IW1_1 = [-0.6621517601635982242 -4.2149687811517315694 -0.060768065257366485343 10.15413117737285198 -2.2360225118479535311 3.2561180639575710849;-0.49312114292529923842 0.24758505974710701136 0.21462914592307744122 0.37472365215370240232 2.296220108425457429 -2.0261832797981322862;-0.67106686766804923838 -0.86297713654760932478 2.0709717899299051602 -0.87562424750530476469 0.090941108332309575357 2.8775748411784740632;0.10266488988290978224 0.8200832805991169483 -3.8011085846898233775 -0.86806584148908971788 -1.3983651226755093866 1.0438616859240921197;-1.0146441128676015886 -0.52517946541107674019 -1.5484905224882166141 0.1265491940893019307 -5.3849764259359025331 3.1314136326058359394;-3.0782630006133109291 1.3506191329198977691 4.0156718055139855039 -0.30445400647574405451 -1.5636931225078283259 1.9367044651168423552;-2.9281986758696718631 0.95361541367025071114 -1.2674719012993125311 3.7130393039194879989 3.637921972353356459 0.92895342168894945445;-0.22426701897964185228 0.10674058155424710759 -0.13950114862372209323 0.22100200593596550092 -0.87825743073921558324 1.303426939379748406;0.050563577974636508827 0.17152652229204104484 2.7934449351999601774 0.017139514011737820576 -1.6961304269545600665 2.2082097189443481433;-0.89903421962651031141 0.56035293106727934731 0.91606014687778269945 3.8462325252393361197 12.717036421905051569 -6.0970304089537155079;-0.25794588833812465145 0.83212540086506281423 3.170122553167263213 0.10407677792534532346 -1.7472786417295076422 2.484773182829602689;-1.2859425974577156637 0.68557096092193781089 2.131593661683406804 0.6836264873000407194 -2.7795601548535984193 -2.6995951828556461471;-3.2606895651168250616 -0.11608228056105976411 -1.0637454094833069806 -4.7312545864230548531 -3.2540170165250743572 2.606756956382999757;-0.059884345277407512176 0.23647854583027846931 0.23014202055643689926 -0.079109782860457472653 -0.95620658710704975469 1.4892285640532652025;-0.32722876100684150646 -0.073144539922480461547 -0.47613023173345830896 0.21646514357767540626 -4.9344509044588829028 -2.2552791158906670965;-0.37051902072955017475 0.41912446549892579117 -4.3328798077465009442 1.0322451607969607146 -1.0681941791645799267 2.116157877830138645;1.0291084941248911377 -0.49656083492966657245 -3.1195286339930672348 -0.36637134013571859281 1.1532325824552012161 -0.51206298373219694131;-0.039575866960826489027 -0.51380613947475861636 -0.059906279807513375391 -14.949946211015655351 -2.2869463208541782606 2.2656258816601009265;-0.74243199839272711049 0.42705653566554008238 2.1044954266952045074 0.28725754079571019783 -2.8254508600357288728 -1.3315352477630484351;0.02766442233679476323 0.11813945808690198835 0.62998084864708281394 -0.071476685177731383303 -0.23197789501900853426 0.53430372095399780719];

% Layer 2
b2 = [1.2083697476987531605;2.4959585715885137525;1.9681777230634180942;1.3232935936043115355;2.0696586885153802626;-2.2763776400856423621;-1.3785825486662199868;-1.0327106606811693457;-1.0288109732692072296;1.4225939946242729395;2.6085663972772699104;2.0559614662910719041;1.3105252213985809639;1.0692694079610103319;-0.53016249177955676597;5.5467958340808092288;0.27384537939082681524;-0.78879605371261252156;3.8457595560102375565;-0.079714075538023634881];
LW2_1 = [-0.15096419627464796598 1.7865856398418591411 0.054926749593941372185 0.65780465706201440845 -0.26759261271044226804 -0.12803607327700053387 0.41355632748769499685 -0.12308862009487120392 0.82500451018418041738 0.27632439864538688168 -0.30019613039629189366 0.24418274089874181243 -0.14600978548126422596 -2.2141184517773604412 -1.2151076881799909302 -0.43312580317613730285 -1.5942401763378732937 0.23581957870872280458 -1.3309919784277115529 -0.42684230250498333969;0.034726287613286208056 -2.5875238022866366983 0.60577736089852729506 -3.4011611939920460834 -0.45901976028338986602 0.064683359522741096237 -0.77636948833838759221 5.8009226563672857679 -2.7703160633313692607 0.26088448150350224131 0.36151290743587877197 0.18524059284094770983 -0.027914973085835435157 -2.502166948826545223 0.33014649944335638931 0.014430120932465248518 -0.88855761251150178204 0.84385042581209801504 -0.04106763280985325365 1.8473641186889429022;0.13201926794415463906 -0.29474734457455914471 -0.46127203399896754554 0.97001554812561563423 -0.27389898802876216921 0.19267646413398781258 -0.50613445189331596197 0.39214665669345138976 -0.53508929969903840096 0.32279507011487090784 -1.2571103056528953879 -0.060124104987383006349 0.25939950051908888806 1.1043349199881116185 -1.1705471733980052917 0.2341328137321561087 -1.9505762311018444155 -0.16922898295665267376 -0.025941426480337124311 -1.2455424541282005269;-0.10006696574207987749 -0.44046278602740901453 -0.16896087138439661035 0.48043354016928146999 0.070071966433684546027 -0.0016686913035819529033 -0.54597110146982930079 1.0993366216700517501 -0.048217255432683750471 0.1365494411108279238 -0.69228361672825222684 0.17454439503804700107 0.19273576693485092171 -0.27097883078163748083 -0.86577349544291304362 0.072299706276197181487 -1.1260935649131418312 0.22036591211417957958 -0.4903953192359819524 0.57261817315829599551;-0.083988326551150738197 0.0055698271602737242647 -0.15108551238562387131 0.87981587361589819185 0.2217352949710185428 0.82187843355661360789 0.16329181421881838743 -0.55796109726123166439 0.99340544233650074357 -0.69403345771787605667 0.51921099137489501452 0.47376900406719846259 0.0098834850384821731906 -1.0677044969733580082 -0.54145908758428806795 -0.32436613650097989581 0.62116468912412092784 0.58108018910789838163 1.899227625744614123 -3.6060622458396336398;0.097540761950712720485 -2.0799571329109300954 0.1617633417597078016 0.34887530663834603883 0.16295039010628362131 -0.12863121297460158865 -0.66318788422939178862 1.736911163229777344 -0.034852923613906551525 0.51397305918502889099 0.34814205610620119202 -0.47305023620073810564 -0.10477708876110328196 -0.34360570147346813918 0.58949134104472211693 -0.25338844660135734221 1.3660741926190429485 -0.58703752638179929058 1.0779112498122371822 3.5480991919536615242;0.49105733613841068275 1.8508043238513287587 0.70317409906177696932 1.2010612635386173075 -0.36587811516552576352 -0.56137276671953473262 0.033997812579194812421 2.1931155044487131711 0.85569672942448871122 0.14235594635957871112 -0.77754651521546624959 -0.018772337675260666878 0.22186375237883598777 -0.52787600626350095201 -0.90798001037332498608 -0.62701166849622003152 -1.9278204581995708189 -0.40579565382407406515 -1.3012608005260357658 1.9078620090140168131;-0.11461540295632247077 1.116732109562836639 0.63718313561029660264 -0.041935136957116869727 -0.012707967670369938948 -0.45664259690895109278 0.77593211184203580988 0.82947393330293484048 0.7958893425999050697 -0.31529886412109464633 -0.51366234594342008446 0.19359368836837675953 0.21631897256302184962 0.40850944870472816195 0.48161265512132084288 -0.15970162560796008488 0.97081105544612389391 -0.54700908564714711968 0.62876142813041846491 2.0504352039840072486;0.28489340728416473469 -1.4387076078069689355 0.42173778222591684628 0.69730690374575265977 0.28059598992160128539 0.059523784259189603529 0.16779844000853569308 1.445106849945586136 -0.060520278086992611299 -0.17469498662082863505 -1.5648925141611955691 0.58708555492887981586 -0.18247347188436321974 0.3868227843586978798 -0.52856029273190141549 0.028106156586573029926 -3.9745845239729145071 1.0969166103268637169 -2.3912573183043130953 -0.62394132053807915828;0.74288971762539990351 0.33788992619029573072 -1.8426159669697204091 -0.36960878790225121238 0.98679825780247953926 0.59043302402646113958 -0.11665425141214515181 2.2570904604459132159 -0.43873006286697846523 0.40251117593660112881 -0.0030488681941719291205 0.1916680698480270062 -0.78557082705804159861 -0.73429502924791123508 0.18979493095285138171 -1.9468401007423128402 0.6885638362130362955 -0.4235536405034543983 0.54042310089635825499 -3.1965144382520218613;0.32650151818387740388 -0.48429369081468848224 0.0045767513924386202251 -1.2165266205456877202 0.2991920032501595994 0.10762772980535703637 0.034687505854429333263 -1.220161490631334944 -1.1386788488576222722 0.0012202053303441103872 -1.5745360879236214302 0.43046801021798025699 0.30752358256254735958 5.4097690924676555824 0.035873992879305705017 -0.3124286593018315572 -2.3852727303626903854 0.36702855636867942613 -0.30840076890426931788 -2.9589720579403557821;0.0094656859530322552676 1.1551715829396129731 -0.17111097033646655019 0.77206667482021318349 0.15248732178353338607 0.30106129733641184831 0.02014460681535377648 0.39902543684414976122 0.98337266491556474257 0.15173154613228412302 0.68931642141589499762 0.38316166595677442208 0.063531073480506197737 -0.61915641009269706263 -0.64055537333605949701 -0.20800532693846171539 0.73245890884687936406 0.3105321885789055214 1.591133246789527611 -2.8977670212363157987;-0.024242037664706728162 -0.8702999921357337243 -0.62910119686577281595 0.75226590812206750591 -0.0050858651425768237886 -0.070063047628938765476 -1.6173625251568286565 2.7334391135798177253 -0.72417328569234251656 0.15556418076113154303 1.0143947337795899966 -1.3749876258128208573 0.065767148136207168041 -5.6396108812155327783 -2.5095793767016671971 -0.2909307342760429349 -1.2941817296441866514 -0.11984214214717703317 -0.61188543811700824016 2.0262864813241763251;-0.071251084870252820447 2.1230777882047755512 -0.37366107054498809426 1.029906688531070591 0.25918903298164214233 -0.27786619960889508141 -0.50087525038008950951 -3.3885761204832500759 0.20414423379703694472 -0.11415100629057960135 0.30636427901999974122 0.18150018376343035609 0.24504439261956451057 1.2466422416305968923 -0.4084273201903498296 -0.18356034288913808394 -0.60413146259308192487 0.76743014459143443773 -1.1826683496677665719 -0.12478108927397477124;0.086285047964455477421 0.56052725074961573526 -0.1120668152825509134 1.4716177069805342903 -0.6845631966938444668 -0.91025770358361102197 -0.70637516153818480369 1.5323904959752090438 -0.68066530214537424026 0.2825670579245764813 0.26921288695327660889 1.0161184884391698535 0.42410602508126615806 -0.95027003031703438651 -1.7043571090288782699 0.38434905618200199129 3.0765308993761784073 -0.26107921726725885625 2.8516911519135392794 2.0042236089241844965;0.2817946418287160526 3.3308756795498024594 -0.3937290133857611063 0.95721591062735666977 0.083496562682977948366 0.62450794430394052092 0.94279010665599494168 2.7631428118441827557 0.065289297355474978923 0.078866830138336355138 -0.087403441697727138338 0.52359628411656589098 -0.89483190362090614389 -1.9121906098224399351 -0.63277597352325343216 -1.7124276001022820104 -0.35440789975830727654 0.86425806271788152557 -0.026323019866626919117 -4.4992478689623567334;-0.10751320479721336587 0.49520195374769593322 0.25826945813538948471 -0.12498822547268192973 -0.5162953258375506671 0.065998167888006831205 0.0083624291933311039465 0.85910631545888949123 0.32750585025242601178 -0.35581862363850774722 0.45368278220161262215 0.40447210168035685474 -0.12327158397142629065 -1.3326808753601686597 0.14198840866912010306 0.091912179858197579074 0.99198247649562942829 -0.52758194289291771018 1.064656383300852216 -1.6852554112876911852;0.24022984598032784409 -0.4528792189057294082 -0.32513608556003853467 -0.053547931398700981953 0.73247877716134335291 -0.056772412548814404909 0.047141945273779267689 -1.6965868733684459357 -0.7438452856526514001 0.47973893155590169313 -0.60385210823194646057 -0.26950481362743844427 0.10039571281445296147 2.2630060593326000173 -0.074609218770262566989 -0.021247563098275110832 -0.67792597066119397731 0.75818696634830329284 -0.62546804608512052148 2.0782992796135690838;0.83787479587864543973 0.88133913467728786095 -3.0382628431140852321 0.64752910824244136467 2.1258820815480583555 0.7109196989133975686 -0.19863439155307507922 5.1372389924813495199 1.3819240510949695189 0.42795295590588777745 0.73779516989354099188 -0.5787113395708470831 -0.91189671291192642144 -1.9691297689673166182 -0.12549385082303268768 -3.3875852145480789268 -1.3485476205521191861 -1.1452687387937303676 -2.6234325275283407741 -4.9135258828763408445;-0.00073559327908113633643 -0.51776732836086925005 0.32768762810004670172 -1.3276474017315236242 -0.71887487483618661699 1.01696772470840191 -0.012866821548837559241 -1.5133784436360413483 -0.028896072702792534315 -0.49376748935732667212 -0.091363996165100991309 2.0777861358294731353 -0.28103778463188527903 1.4662522340899379891 0.82882637549235149699 0.090729939643694926255 -0.63329390196017121273 0.062268404531366859656 -0.80147947103831118643 -4.1475026999338195921];

% Layer 3
b3 = 1.0524776832012072703;
LW3_2 = [-1.0495704446870046667 0.61945234044757313541 1.3612669205545164797 -4.1017547480996325859 -1.0085745885308383674 -0.97724134751021729528 -0.59061895604455372233 1.1673205300881532853 0.593360583631281302 0.56503180099321714813 -0.79598733179077774125 0.99141389076139896375 0.41294913943553696045 0.99731660240296871045 -0.47125876433061159032 0.40357807217400798461 -7.0716710741577495369 -4.1538150417160029093 -0.28559466983858400324 0.78072409565384548458];

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
