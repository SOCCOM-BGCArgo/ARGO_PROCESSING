function [Y,Xf,Af] = ESPER_TA_8_Atl_3(X,~,~)
%ESPER_TA_8_ATL_3 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:02.
% 
% [Y] = ESPER_TA_8_Atl_3(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999999390765;-40;0;11.29;-1.99149636107278];
x1_step1.gain = [1.0000000017546;1.00031699124092;0.0153846153846154;0.000310077519379845;0.0705602618858232;0.0604107502911417];
x1_step1.ymin = -1;

% Layer 1
b1 = [1.1039158086542377735;4.0252326467863763781;-0.94866568098873449433;1.728635778285600999;1.5660025817470522025;-1.1264356437944651468;-1.9219819167127256687;1.9856213641655364732;-0.53985705519367155336;2.465321750932462308;-0.34329703504408370929;-1.0165185737303357172;-0.40187288664459996923;1.4370231408368971948;1.4474395352993763009;-0.05480235302483365406;0.17388079106473861013;-0.13790794633971131811;-0.76247839103079706469;0.55456887015433342647;-0.72707201538300614274;2.2830588778162224095;-2.0344020029843710162;2.0907825937554189721;1.3749243058351168223];
IW1_1 = [-1.6253708234307640002 -1.7452372525329391539 -1.0248625748918842593 -0.089340653876464559691 1.6218632219330677557 -0.63622406077297810612;-1.7784178638278120843 -0.79658332597575121081 -3.7977737606323280772 0.067810267694463799426 -0.41900928173995477533 0.61156488278785881008;-0.77274731647687744385 -0.98903482480014726086 1.2371151279782757193 0.2122486461284026904 -0.38006201410515078054 1.8474648170383043411;-1.1668005090730404483 2.0382404143144756681 -0.91902308943670163455 0.24375780715484018923 0.31614000900502930991 0.13121912665198626202;-0.71753070147244846222 0.35599252620798294 1.1439483259676865856 1.5751116313970254801 4.8642326764823335949 0.073249290823445278864;0.62939543111072138171 0.37706934765408911536 0.71813208284075791266 -0.97971669805090577654 1.4384405323123399256 -0.073312670422251105529;-0.010942329223531555016 -1.3906367539108031028 0.42087062492106297462 -0.37791463046356410516 1.135353500033783547 0.69776889811640430139;-0.11899759124490601192 -0.77514611835762758929 1.0731979734472270049 -0.19667900399845369974 -2.3385232755260569881 0.63051740384957444263;-1.1884034039332251886 2.179497369819210828 0.75321762309051942808 -0.094596486746803709611 2.2666136198687931724 0.7749387070698771085;-0.59300074142103686992 -0.63955581322575560943 -0.32709779293368201891 -0.65262171320668593832 -0.98994626370409288807 -2.2179398638449132619;1.0052061921800035105 1.9076853249807108615 1.1127956410742869675 0.24823493768113794045 -0.27064757543680573182 -0.020928662726182628573;-0.11617898324912748209 1.3112451895755645737 0.99114796406567551301 0.26042032672723824938 1.460852509212418715 -0.3633954191299708647;0.23370365009195565031 -0.55942130708931492311 -1.5220368225280249685 0.97425473442347076958 1.8442523368687657737 0.34758920512701424466;-0.42592497656246219861 -0.46456740870456170178 0.77491698026760835738 0.30264453219922243798 -3.6264679045856835771 -0.55575113137343445935;0.49074453694532121917 0.58456013087194613753 -1.1570031284429440621 -1.0359693568821022325 1.9182517217585401959 -2.2627206981301459798;-1.6519968234724899947 -0.13943851777115812141 0.018348925303692475985 -1.3468887540007234982 0.70285910685333441439 2.1765966339816253949;-2.4073677379456044356 1.1075281301014541224 -1.4915160302618375443 0.47040564895070430262 -0.22815464152256137442 0.7336576604037412519;1.0477541979212625645 0.64276576294630416442 -1.0379554874723286861 0.17857951162150978441 1.305481102148338568 -0.31428627993690894105;-2.3440650395729232258 0.6686176538425578908 -1.5120907686651219315 0.55386490542648247892 -1.6397186958510312493 0.75336451841941654006;0.82514184357809716719 1.3972178876465033426 0.49244033488559857403 -0.39037364618872794075 -0.40129343484415802568 -0.30261198971995739715;-0.4521580979688356039 -0.68409069407923672745 2.0178041747979613696 0.37325886880452852079 -1.8452772654352860293 -0.67021161023147046709;0.70653646088525301483 0.74218693371474109277 -0.42585819783445544306 -0.8818602460202586979 1.6138616187220979903 2.6410875964294571183;-1.1055266895250237535 0.93821002504542116984 -1.4168283317990442161 0.37593386559290481452 1.8395904772437992047 -1.7232389190345425156;0.73615737835696271052 -0.12447664856522623611 -0.39328956030687067802 -1.9658121089712430329 0.84388716712800460051 -0.43951859058693787397;0.28403865006120110337 1.307066601965602981 0.74248536720773761211 0.099986785125247681916 1.0099150506178793396 -0.83843477477668093556];

% Layer 2
b2 = [-0.85139180542036219101;-0.54199686062575724854;-1.035356531404224345;-0.94356489040528856993;-0.93944232898906598894;-0.21890285870888182007;0.090901770474175164205;0.11428204742337629007;0.26143377677212581434;-0.39431553573154859427;0.30513236000851590157;1.3091953048120894554;-0.92121687592114276377;1.9449119183898897489;1.5500181372167238791];
LW2_1 = [-0.21714546517764665712 -0.88751317953003383554 -0.12986935019581852258 0.61302029906631738676 0.43316364175154692395 -1.1219929318060761769 0.30150810948199568973 0.30999406180446048564 -0.59934178535288196166 0.81521419072823997354 -0.4384052841296154357 1.1391890860944442476 -0.53393190830032255523 0.07525372227918339374 0.31982262328489741154 0.036972612850735969703 0.40930813118534631334 0.61318537693347718864 -0.30906776539495567357 -0.096352076740473427408 0.60028827730053768708 -0.2274776262496367718 -0.24117022743631555404 0.81972180415190443181 0.42721329511775318899;0.64903493333052297842 0.31883196867634883143 0.15852455242359381193 0.72269489379227824788 -0.10456248200499450307 -0.10419321953159646488 1.3345154795301275019 -1.9527829090491470865 0.074054122792499579853 0.61882145842337155095 0.63026343029980425037 -1.1594960801254956539 0.84056583011315333387 1.3446361283196845982 -0.23924160693787455134 -0.097851319486208940268 -0.40331287990275266342 0.29223623291879513575 1.07989111844785457 1.5134073274280457344 -0.59680617138599734162 -0.33672718120554123722 -0.63288564931930901825 1.0956126261687451073 0.36257909252783793175;0.24647120151361578966 0.36436626565149277113 -0.017740905048936463778 -1.0838346031073040621 -1.9951389328124224942 -0.87676412046571683412 -0.16435177967658951004 0.24957325834114885144 0.57503891204440038543 0.34936741747734867403 0.67929725525779782735 -0.053096275032709708275 -0.08475151560948081253 -0.87261892242077321313 -0.16890638587327622133 -0.055364310048978494783 0.10081667974359193041 0.11709031427303760731 -0.92493612774735112225 0.20911785330582544074 0.01655698434192359389 -0.057443060175055259131 -1.7847025258509086054 0.54930609837207611257 -0.67098896450357259358;0.078791097391434722352 0.49416449830389069797 0.1869281882414658702 1.1215048930306268371 1.0237286428662435434 0.73844514930380267259 0.17972785524930537382 -0.19289619652168657482 -0.2984375135481889818 -0.032177902368035776748 -0.18940816005735186534 0.19606898986408138885 -0.3411466831889036766 0.34328569948628262498 -0.30590577464518858841 -0.13104709391404972818 0.086191062803240264811 -0.11949554801782250046 -1.4990400033426316551 -0.53709937822572273447 0.31211280631315041489 -0.17953743517955822462 1.5804475691737065457 0.35244866011100062231 0.65513353432972665225;0.010460693575569965028 0.35374136142528234128 0.94283556216943209538 0.078127312895411055171 1.0341259558879671854 -0.46465977826395921113 0.67601548974034042505 0.019371066208861146413 -0.17524282248922162797 0.3657822369943909413 0.53202458892794224354 0.95690966491425888396 0.47084136169018414586 0.59899211712854694145 -0.70599515720917715988 0.30670184043783488281 0.32140773319746973602 -0.14047182561753085439 0.65455492388588587538 -0.12973975121912628694 0.79843032169312189161 0.46429160333264379856 0.38094756965834469797 0.13502988660349501959 -0.81557495986077754857;0.010453275828443985007 -0.64716110157400486003 -0.96242665064302435063 1.4949326307075576636 0.60510539054346568211 0.37161558277097117564 1.7231727795317868246 -1.2160364475462890343 -0.11727981445146461181 -1.317943914689296836 0.023427579073358205386 -2.1492749550741998554 0.29824217078873083153 0.24013887434911138086 1.2665979790515555781 -0.2055698178728193326 0.71013180633510974982 -0.1694002801564730154 -0.87563945159905487436 0.62815797501555603599 -0.82677684568394660403 -1.4117743250797605281 -0.13473032918947902559 0.66973703151873231754 -0.31245905564685055866;0.61649479817268959536 -2.6803627713801381205 -0.42742581533927281878 -1.3726863571878342718 -0.76901618005510596898 -0.83366913760071159611 -1.0166116820421196021 1.2336145089155130528 0.11456752563532932032 0.035155224787570500677 -0.10903996311106912975 -1.7951850220423681836 -0.26211581697417357661 0.18021926570404642187 0.30271883601501980987 0.79661802747642795897 -0.38538289107911677389 0.67075251919795020328 0.39312910939663181509 -0.70048523213097080298 0.18267466392917006135 0.46314981042985758908 0.65909380751754442773 0.55193446287430691921 0.016671992781187989358;0.12459624445532695158 -0.15618099914575409026 -1.2179375150029498887 -0.51556559749378028101 0.695067116217702341 1.3231158028878093891 1.5509346791730169723 -0.37430374987532449849 -0.20128958004838010121 0.91262828102474435887 0.7748566770329525788 0.47238685211565700772 -0.24610220926250894746 0.37569935992423098003 -1.7434799709382760025 0.25675957174158975382 -0.029188123206761194817 -0.015297388998823576237 -1.0604669502411769866 0.6107851440675117205 -0.63234994451486370437 1.321622076951777025 -0.70076735358400621401 0.017232915260363947108 -0.19623890255987538889;0.96612039771509694575 0.53810081692211941817 -0.13517175990804303209 0.73659416424902979248 -0.74109604161277831125 -0.92981895193363561791 -0.37252502805256393348 -0.1075395008124659102 -0.26002539982164535326 0.54474101151972031953 0.20554353500821367917 0.3864122378301891314 -1.2093102410637288635 -2.0275084822603246515 -1.1491391288817331962 -0.54574389562602187986 0.28763839785477451905 -0.79853457033311181501 0.0059422516301702614158 -0.43801253881465945783 0.84477904488474697686 -0.0033301186258269960955 -1.1289596549720584129 0.17146447537798673033 0.27232090825566185321;-0.078870544891542410104 0.16405537715659790354 0.36194806515382332845 -0.85395131684421254903 0.50853906862293085656 0.37865409905565228232 -0.3950109494972787183 0.33260955573469513924 -0.47620650479098886043 -1.0814020426700272637 -0.38117920055702053661 0.044591914199999367618 0.05377897852901192266 -0.6220191898237056316 -0.1551523228401251453 -0.35592289221507650288 0.61281294076858427999 0.67285807652704032744 -0.26554537757515073171 -0.1922892152954909506 0.32180324840298624167 -0.048775622173113229729 0.66038395327566268289 -0.44847960566627398382 0.60913325722859512634;0.019214520608377350791 0.1675657048213543554 -0.023074694312843372496 -0.17077972339022495785 -0.83042104394696458769 -0.018016251059977818338 -0.83125287166840466924 1.5174736424371608745 -0.57285540394230705097 -1.0611477185911868393 0.13571395944897926555 1.2053129016254948702 -1.1693514672138460853 -1.7424553286804334018 0.16686228327319854725 0.079061792816346912471 0.13079538157394832165 -0.043009615696447428546 -0.87426666204989011977 -1.6276533647608748545 0.96609916812645457718 -0.23878933349202671743 0.043797256978266953242 -0.26895915696666616856 0.27993944117501035063;-0.25489840710482014607 0.061769567972043658288 -0.31676051150708350868 0.70570043334966603421 -0.96687269614138837692 -0.93588572233067934825 0.58198360432682150201 -0.51689240277127879164 0.37060439367926456544 0.83845825124851347532 0.058620592036188974572 0.19846159044587941001 -0.518685919356487668 -0.27888880919947034442 -0.51240477153887720174 0.2865438306831818438 -1.2964687762132405258 -0.59537005345672244783 -0.87736582409913732583 0.46276032046442516776 0.43415773985282873904 -0.8423099972403332858 -0.4813616938425329117 0.68928868526679087925 -0.71379308654079076657;0.87312101762518756765 -1.0146075492898036785 0.03572298382200092548 0.016913303332598546985 -0.23783293814293249291 0.4026205591671837003 0.43009787643114888445 -1.502998527846267196 0.092850818231968612615 -0.4029262042370695962 -0.51750673575649175362 -0.58737869450420188144 0.61719309219582063619 0.56504798115116727786 0.23722224963642471462 -0.00074242443912391560719 0.36126832191670493755 0.72792703626825194796 -0.038159717288646359534 0.64760212439101361159 0.14497352532376020506 0.70831934500088833406 -0.35036080557714249784 -0.31500517560503865644 0.50818347065211388891;0.47734642514703934379 -0.63811977522882934455 0.4282703601131571669 -0.26434399877679537472 0.72353656442357783263 0.1996486065306034563 -0.46687626850875690865 -0.17344588829537743724 -0.45761686686249974931 0.17211391874288151094 -0.89700873745848108687 1.728470481930547642 -0.42042139661624450264 -1.0876555445911697539 -0.071571253501542947406 -0.26562278471465383678 0.55084329980291379947 -0.6364314401799228138 0.11849271673298425456 -0.91108903366938254198 -1.5941665019177155838 0.39108476849200457526 -0.39233186375102813148 0.34364451949094787508 -0.27939034128776202293;-1.3832992016593315121 0.4320270524630821507 -0.28123116586173224629 0.11118983727863988864 1.0432286343368140091 -0.13841217481987355686 -0.42130639833770311409 1.6603059753526685505 0.13818949323198989609 0.40800754109924036728 -1.1577929525260206578 1.5629651358298308761 -0.72921261470896547863 -0.28553898008150513377 -0.67306114944862660554 -0.13169967244232688719 -0.55677914990030508147 -0.74632798959515456172 -0.15485500539559535205 -0.65002713883167939635 0.030532559328386089936 -0.49392082728815267512 0.8198253368809039765 -0.28219942216172028138 -0.55385151984384650703];

% Layer 3
b3 = -0.45471984132076420648;
LW3_2 = [0.30112436044904461463 0.69108828246500653147 0.76824283569927298565 0.72128264380843010173 0.19693551265372513659 0.18283700247666898053 -0.053250203006932388117 0.0039140806427374998447 -0.2721943637794451698 -0.41288183143900886174 0.70800509857867444108 -0.43848484919155794604 0.5375677341127156339 1.2465217238054433491 0.50748047175545918375];

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
