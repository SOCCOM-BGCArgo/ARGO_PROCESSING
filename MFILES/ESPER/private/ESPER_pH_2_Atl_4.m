function [Y,Xf,Af] = ESPER_pH_2_Atl_4(X,~,~)
%ESPER_PH_2_ATL_4 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:21.
% 
% [Y] = ESPER_pH_2_Atl_4(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999025224415;-0.999999994516886;-40;0;11.29;-1.99149636107278;-0.2178;-0.12];
x1_step1.gain = [1.00000048738803;1.0003169936794;0.0153846153846154;0.000310077519379845;0.0705602618858232;0.0604107502911417;0.0479365703301392;0.0145655815308426];
x1_step1.ymin = -1;

% Layer 1
b1 = [-2.5531508942166034615;-2.1368963130571008158;-1.9267392004741772382;-1.4926120604680965798;1.3245173977091690976;1.1926040558655301993;1.1656243414555340632;-1.5244003169274475962;-1.0153861187159536339;0.61940076171844959951;0.45778713992553671241;0.66066765984502595632;0.41563945518098549536;0.13352033022514375737;0.7694955537048147054;-0.22030694121366040106;0.73971316270040210483;-0.023963515390706853492;-0.6230302745216629301;-0.58811382455399241476;0.80505025905430371846;-1.1918340078907903035;1.1221525445769753304;1.3669393318016038741;1.4148999064309650819;1.3383167388714074519;2.2220451735806014959;1.6843687660801938133;2.4964686379436260566;1.734045934845744652];
IW1_1 = [0.18415441452835262792 -0.24237796308988426852 -0.62292876329240221001 0.63692603993396157414 0.523715986386557808 0.57566042440425035931 0.88958942804111651448 0.26921179889692026244;0.37152921198998922758 0.59708746168945958654 -0.8148176210471950176 -0.8876110010011550866 0.052885570003571504005 -0.57190277265824396302 0.41898492928227143439 -1.4088703465288521244;-0.23910685657031252593 -0.43334566493128895415 1.8116770224742702045 -1.0984910292989611413 0.31794968351262775652 1.1938247965236061798 0.85039999884652039963 0.093761894496625627915;0.86480978030098731857 0.50300564905557565432 0.039496411075440510108 -1.4795176472815692037 -1.1474368126190834793 1.1717363710036796665 0.98965672812224136923 -0.53599306731201057374;-0.3083984515716317043 -0.603659077644289499 0.64029855008586333831 -0.45954185358919447113 -1.2859545502714428977 0.19754440568019340785 -1.2784237879068420796 -1.0754377874668170634;-1.0795055702693314714 -0.53067215138328660551 0.26420974374917982619 -0.69029811001431717887 -1.3040000670775617753 -0.33144213959994262009 0.17514293821454571254 -1.0192260216932949035;-1.0835233043878134129 0.41369355942279062432 -0.94095891996360170317 -1.4607781119709553064 0.45900200816183406305 1.0855669361341453438 -0.43824283654561685397 0.24704859332452266951;0.9859574954900357513 0.22064811762678915308 0.16126680778451221721 -0.12550671974170657386 0.2563717137092039744 0.20039576474418271945 -0.96542723934500973559 -0.39804102228033433342;0.94323079665901310431 -0.94707781211305364089 0.35759562981513665081 0.11631101954752005867 -0.044421783358503101613 -1.3657542375650166822 -0.74335449577754264538 -0.79597601680961127535;0.037534622573636967913 1.2725018882284244714 1.3750427959037081482 0.9282530094586706193 0.68591417508176077966 -1.7636969988898432149 -0.43004016153667357081 0.4425141229285450839;-0.66629037403568813769 0.46860884421785409248 0.50538568169289699306 0.70955994850556736164 0.71515299641573593092 1.9445114173238553512 -1.1221020456699091206 -0.26592836574948353068;0.1263956632717334605 -1.2581888534914451583 0.47070308560379425789 -0.37028664087796547832 0.3475918756807016563 1.4721494811311142126 -0.95840060467281307233 -0.50043431995379539234;-0.66057456220077026554 -0.71082048908139716659 -1.1275715671850115029 -0.75667161260199700124 -0.87365335977972868875 -0.71003047962833465423 -1.6933420358886879242 1.1779960165807823458;-0.27593268915526109053 -0.29207827688627330254 -0.46623078694413516709 -1.3044019721978892257 -0.70626260111003402731 1.4478988935836987206 -0.77514608639397641898 1.6238321281779917182;-0.96043558376455062664 1.1829953136156790006 1.5330203465418494879 -0.055425982875796940552 -0.27074876592755686699 0.067939559038497895904 0.1368198104897038192 1.6767868242037973214;0.35044659908406794457 0.12801914019488844221 -0.56649289325770979531 -0.03497104367021129645 0.89165760756614820259 -0.18901005209749763281 1.1680906012650491554 -1.1762642404234160409;0.67600567176545245385 -0.65192876482032491303 -1.1389315088102383022 0.643158584889169882 -0.77910304562296539554 0.86064497059986455607 -0.1220400864530968299 -1.1031056473242901728;-1.590876998426681066 -0.1569998438431305654 1.1499438502842400656 0.18102336718046621589 -0.27573625358851078238 0.24482172126467827145 -1.152759661164177718 0.16710109458280392181;-0.9570279049396664961 -0.83773225578051402351 -0.5803228514717692077 -0.50161947761710046212 0.63905131614422494568 -0.36383743872129131747 -0.37977497602956666878 -0.95803638877680974417;-1.0874605358576434622 0.3497273236664681928 1.8871540309536094515 -0.72195372544627545164 -1.1905724229648035895 -0.60650636985295447445 -0.069223677707440228746 -0.30615302158412727174;0.11545274635986189915 0.83332375325445284808 -0.49593671682816609092 -0.74001934030929872321 -1.0489896472496238466 0.25281267399104345728 -0.95540635905340287604 -0.6641354095402157709;-1.2903440718874943371 0.45178268697006984667 -0.36522915188527105279 -0.34264406248618256567 0.87211347288230933739 0.062969006655331838385 -1.0984953320935590959 0.0031109020154866532173;0.73094558894918704173 -0.62855632928903182322 -1.2784445803265838126 -0.777632988302151551 1.3376207301017646323 0.14144561853225726122 0.16389633025738303829 -0.46863911148054238076;0.59190074611252840864 0.1307986674182973208 -2.3803530587459076351 0.037645588459755525979 1.2622791887097986052 -1.0730268134798748303 0.3581470793351507198 1.4340513966774817156;1.0244409445048454632 1.013139627579728641 0.27985629211578372244 0.54552156881767233099 0.32897070231284464281 0.026682913679945210322 0.59473540209753439356 1.0061011758476343658;1.5201924514470841476 -1.0964803111500254218 -0.35052129946394300353 0.028866958505414687386 -0.9707925844626662526 -0.86205404857648060268 0.47020026955502797206 -0.080667898284076264281;0.57803956086150032334 -0.91178284927189978681 0.36801872946250163698 0.96076129869356619029 0.060638189908223630309 0.16512177548784190417 2.037914873042083741 -0.40860381633716813532;0.44930385508918829762 1.2922295778650056963 0.75929484115907996955 -0.37204522448838700521 0.3676413002464916091 0.74298461115816627487 -1.3951874426936077267 -0.23552944674725911911;0.88449098941306103505 -1.0491814948603792601 0.6858051407986313075 0.075955252410147625075 -0.22017082276351998549 -0.94825390486332516993 -0.57273947079664921933 -0.0076858319315852002135;0.88608386091440172194 0.035409570048008182086 -0.87923719365954267069 1.6262641427972703845 -0.39511669484426076293 -0.22321088340749375978 0.10494501525827146238 -1.2781942599837410324];

% Layer 2
b2 = [1.3392509147951441317;-1.1934588612617027525;0.73315094652265855313;-0.50585268152203755054;0.18519630170059520058;-0.064409820717733948237;0.15505432811945885074;0.54905109500734605721;-0.93941537815495246289;-0.98018301088903581064];
LW2_1 = [-0.37512067792714715031 0.23144765738788738996 0.30677278915829925277 -0.58259836559039868131 0.16264531903397783097 0.09745033885089346426 0.19166036143010942161 0.085672700209568963814 0.075067200689635904531 -0.44043388642553971124 -0.14324412739044223297 0.12126109417047882533 0.063555509368576124096 0.27322171657749344531 -0.34605688027506892324 -0.58836897413740163465 -0.50281907582411700108 0.22947247227216593934 -0.10222127169721029194 0.36783220585969045313 -0.16079538841562374119 -0.9492405946285598084 0.18451408720848741507 0.47710951325643374821 0.40062431166826428797 -0.25800870924281016361 -0.016346580374023360338 -0.34647254541266092609 -0.17290481653617029156 -0.97923483826650747908;0.38905838982096896395 -0.85982711037602732596 0.47203750709447039879 -0.5640646801734052751 0.15308465812014176999 0.0042600315208663431729 -0.58429050004519977879 -0.094097320062394726459 0.22374257361735752636 0.34482926335247299221 -0.44394727714859583134 0.51569654597128999907 1.0597473663946925804 0.77405433923786659545 -0.52967033403383745327 -0.87057732530629128309 0.15451166164930488223 -0.22037687712495562797 -0.18321481963954275307 0.28530396114160172516 0.21612444447137518777 0.50796891802679311478 -0.46311399324749641382 0.43780908610199009701 -0.23037095028581092149 -0.14269672218823223342 0.34859872077290865455 -0.49633533067742741096 -0.50528183538414450648 -0.55049334982286068918;-0.21229497338681740937 0.78190077892083054145 -0.070834965298937666023 -0.71447599926473037701 -0.40069212976808193938 0.34049593293219820866 0.41798317893739111994 0.21297897566072865327 0.10521365572959337498 0.35831344521484848054 0.73402443128118777249 -0.054709247144974826327 -0.47623439755769919612 0.26995160996959777222 0.0034057287872799393846 -0.91677010310699347251 0.81352550378520993579 -0.54835439226006976732 -0.26193514626787639532 0.65412140604584734227 -0.41796994468435094694 0.1004542275225473319 0.63843971800334908817 0.738336176891173257 -0.16062240049398870956 -0.19984756815543600061 0.86891515510857630922 -0.41586103061830614891 -0.2531817137899707304 0.055864144552919707543;0.42659287481632673922 0.1080408696945859004 0.17923137547803288938 -0.50917320628097795243 -0.081057321872513229954 -0.14323851530812803379 -0.38694327406962819582 0.27630983468092218391 -0.1094459955540707935 -0.1641889322965678355 -0.35172375685944712309 -0.43172849653010902937 0.19860045921862276574 1.0169340731312817017 0.69630225415658619781 -0.075294749995895332573 -0.066851147418713258164 -0.24066215239745941679 -0.13156557223293280368 -0.39099768959337849461 0.64224388658036613009 0.41036321331703173865 0.15370759278187484664 -0.52587669969108596923 -0.52447196859129807667 0.066064905781421198738 -0.28099562824267115824 0.45576944433170585302 -0.25614735260610749235 -0.70512889111124876784;-0.067846958812902050107 -1.0498033900325278722 0.72681344811852710563 -0.092648968871754000531 -0.46196359337285558544 0.28059743818063292986 0.075957767373997800431 -0.14223836451436133732 -0.8473643752839378207 -0.14199816846380386437 0.47987983827735064546 0.15222491895210576995 -0.14907572184822223682 -0.7370334802772423588 0.21927846326334932137 -0.23996326718326463356 0.068018262090618814009 -0.63037679616914776837 0.01551014140348545671 -0.64092444778226365099 0.67786413554631319567 0.1601651231617887694 0.11847555201840723194 0.44852308593812173498 0.058125646136403398556 -0.40506965646266773762 0.57214820376813979586 0.89615495885834939305 0.029951073147423439708 0.95745952672018475038;0.071850266671108359295 0.16503576488182522541 -0.11194693976258197898 0.33361725493221405792 0.68036680243094838083 -0.64972861931016467896 -0.092559972798754219792 -0.6148641716745268937 -0.16053184880656137801 0.3537842479698766196 0.75243918152049382719 -0.14430131793872255641 -0.50236630639576929092 0.3669134163711173513 -0.28998373536390054106 0.24377107293670052801 -0.69638342116563023421 -0.45664346279917295002 0.072348430855062692557 0.14465613639141228219 0.28482808081496741037 -0.12169754727225240432 -0.021587550287073643124 -0.05420594219848132167 -0.56017100383562778987 -0.24695591558713908542 0.47009162461025355562 -0.011193842698230946006 -0.34599369487121628408 0.49327763729358864442;-0.11525087705166023322 -0.5609141047465149299 -0.48076645411776869343 0.74056086235756812819 -0.070826459404711480938 0.37667960539352152161 -0.1249452389259479157 -0.39397552181754807998 0.043249509893400137073 0.4056879128616979191 -0.29771013920092515015 0.19315548544317048818 0.33431136722356979973 0.041044647154964836178 -0.13862113043144397295 0.24850702804409574576 0.19677019910108423106 -0.22752935237695967596 0.015661688659423286807 -0.11276959104082591445 0.50112782422567114793 0.3523533712265222162 -0.75893180610860155877 -0.62923627784629709847 0.047125413545273900362 0.11970043481774018945 -0.35661757719439574066 0.0045598903698978784232 -0.3435617873652321852 -0.14308580543968421606;0.37143866692370508531 -0.31066972334304643333 0.74441374204660815117 0.15585024244554818162 -0.30788793126091795749 -0.060909555547161819966 -0.16989065254218449486 -0.15426604641782171168 -0.35436634630624247455 0.27506959072211617512 -0.2947257455499505685 -0.37781796698918157507 -0.025567625597210319466 0.023186743010207983007 -0.34884863690804529535 -0.800269007288635148 -0.17649827146484367568 0.30352783460957405426 0.17828959743913250224 0.00104095447523775992 -0.32191197981738539635 -0.66187501869592735293 -0.22208931700327472214 0.4347114513647223899 -0.55238121138661622211 0.7463242433237755602 -0.23994622312121235974 0.62838293426699187716 0.07247993742072582346 -0.31654569000557453151;-0.2525481541687080167 -0.28575275009510392454 -0.31497288331611905665 0.25347854714974377144 0.066465438020309339229 0.3504294705796647813 -0.23366486050890830639 -0.48348307768257475781 -0.27905236366278662885 0.27099729886696088954 0.21813463511768615977 0.18698164918935544776 0.48915290275509981655 -0.24371437121306846474 0.90107594592445050541 -0.33464732724738294189 0.070739330499529184948 -0.12217695709452651431 -0.8046380736270448697 0.74125635097712061494 -0.059600149053721310222 -0.34401788307587477522 -0.20824621545124732958 0.36640192161763901479 -0.032945992943229564232 -0.2002904444104366688 0.68288327239365220755 -0.099010035775201440944 0.6684565476653477889 0.80542918677824659746;-0.1691630421315672983 0.19249401069322244484 0.30214345225712668785 0.43699586678193313993 -0.26330687959002746501 0.095240340051862590331 0.14174371727336057547 -0.81915534992018335103 -0.40500944332317484475 0.67515968674371784619 0.32587316624877610272 0.57530927616140348224 -0.49048431258037139546 0.94393657491341242949 -0.21219131790462059106 -0.75219820112826307756 0.62467953084693472743 -0.41185453670991634123 0.40749775806622756669 0.60152857506574752833 0.31498832407530574784 -0.66819884932365114505 0.21827740644013360338 0.13231975448460234213 -0.33657551203630681202 -0.17061907038887197396 0.05727993483800089547 0.02575736161656390355 0.29621657973585291401 0.25767625886507389277];

% Layer 3
b3 = -0.20026384718501460958;
LW3_2 = [0.87708528891640369363 -1.5667592441378053181 1.8724165202828761956 1.3426878062177525219 -1.4257782011465263139 0.98868630619353492506 1.5509803835088520341 0.71427990496112681029 -0.81758707468930758022 -1.0985737795924905846];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 2.03830178109291;
y1_step1.xoffset = 7.47181919175367;

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