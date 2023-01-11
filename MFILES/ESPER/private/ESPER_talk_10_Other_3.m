function [Y,Xf,Af] = ESPER_talk_10_Other_3(X,~,~)
%ESPER_TALK_10_OTHER_3 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 12-Jan-2021 17:03:48.
% 
% [Y] = ESPER_talk_10_Other_3(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999540494077341;-75.751;0;29.715;-0.9;-0.1386];
x1_step1.gain = [1.00000001035699;1.00022983623563;0.0140762793578401;0.000302617642608564;0.251952632905014;0.0422119037568594;0.00826143816767911];
x1_step1.ymin = -1;

% Layer 1
b1 = [1.8569360703099742;-2.4359187374253355;0.54018055448852909;-1.8569575048743017;-1.8001784670881285;-0.0010587436840623658;0.49637695681608035;-0.32643663674345891;-0.0014208159208149943;-0.83393030080741171;-0.81308733259857768;-2.5297352935534212;-0.22356407849746038;-1.1982577491750495;0.28860483566499684;-1.3899611585891785;1.7976013939272972;1.2005855630738944;-0.80340643483012786;1.7749676723784271;-0.050391902183460145;-1.1119909644575214;-2.4977629970903048;0.71066333218900835;0.42268357100648318];
IW1_1 = [-0.005121517584827661 -0.44407389931541585 0.86790476646715886 0.048856462537868309 1.4671679800276405 -2.0603216420960551 -0.73188899849647937;0.035711719086032871 0.081325127943580142 0.033579577759477744 0.57762810733422865 -0.42651068220614002 -0.20375468804645577 -2.0222422382981255;0.042608012341231433 0.32279366798598813 0.63341293954830524 0.45641471922775295 2.1972242944047062 -1.3798800610821365 1.3118979131873167;-0.0021564185527458203 0.50931904967521213 -0.96561010062070762 0.35428426756153508 -0.99806638086011745 2.2604482688839549 0.68329032889615937;1.2999452846674413 1.6306320666575218 -1.0426855982384906 -0.062192125013356178 -0.54706208098737297 -0.01375770370613995 0.35811493979573211;0.22172951257946305 0.15982880200276392 0.8554479478900342 -0.059042798426681253 1.54078887327215 -1.0456504757966103 0.85629437981931211;1.2566815397469446 -0.14502097250034032 -0.11995166852734969 0.88632225801440889 1.0284955191582836 0.46946477216592331 0.57671960164066349;-0.26936267728048813 -0.64918194821871766 0.21078846265489587 -0.079639487619268132 -0.71822894792674874 0.31516095674998901 0.28049954809481292;-0.93131715271262161 -0.38174397105529867 -0.18114241696463479 -0.25444214371569257 0.038326132023669919 -0.055346034623740846 -0.31932527545684436;0.34163432288481704 -0.33278220556617072 1.1498202139393461 -0.37575786264075639 0.50350111908525341 -0.7782870766426323 -0.0039999975435692694;-0.46547323614556863 -1.0953905179279306 -0.26776302621540976 -0.27127889100806152 -0.65989804921210538 -0.7666785143447632 -0.67393825387713013;0.45030596891704266 -1.2081471247440796 0.62388112827952868 0.17247407984258981 0.13355055052732176 -0.67485396510118889 -1.0082772028233933;-1.0163797881651866 0.69052438334747401 0.1274961317664576 -0.5823661880261044 0.094332537646629591 0.47892449241153295 0.066184652624566118;-1.1936635976636467 -1.5501612481915394 0.57548077190882685 0.098757561955611806 2.8666471115590135 0.21723360135495523 -0.88537677128395209;-0.34559642661620105 0.15524097551926555 0.059849958868417409 -0.085638799253822534 -0.52256419519775488 0.47283705939405285 0.18487226723179678;-1.6360010871776878 1.1774966884507068 -0.60690154891899761 0.63026812746246907 -1.2446936770953991 0.54871730224249515 -0.1351730336545745;1.1417828925622964 -0.93340235568154095 -0.93071291337895456 -0.41291691405909592 0.25615832047076642 0.1921174791390865 -0.1591130806774555;0.065709693764369428 -0.098639070741091475 -0.38735112712280045 0.07542914397595607 -1.7319480105568421 0.013213157591687454 0.61391867131967526;0.040698168685129314 0.0073110402133582299 0.50694529526150856 0.020658954060995653 1.3065817341991413 -0.13409415609648082 -0.14472711987658821;0.21249047820919048 -1.0974177241336904 1.0015701801280039 0.9318642452956406 -0.50672120778010099 -1.6022774304059173 -0.87341159862135764;0.39618697139358428 -0.28532715040130885 -1.0629506923773941 -0.56376400941242943 2.6107206270333032 0.40360040857131518 0.069617444476066576;0.37951051845987427 -0.21776846004042488 1.1049532261609083 0.16945496441582603 1.1636672252730125 -0.58040977499044188 -0.61899217594535461;-0.37322244741918487 0.64713600237616875 -3.209433830891276 1.1794232799046329 1.4123898207651482 -1.3974200447513261 -0.43180422340767288;0.23932272010372216 -0.33153344214716662 0.68609378308277214 -0.47022079061976002 0.99866400507071607 0.45828169912917871 -0.045711416988478533;-0.083604644447698259 -0.36911323462136009 -1.437395168145958 0.14387171033310026 0.5491066889911782 -0.48782873746789451 0.62717171939968175];

% Layer 2
b2 = [6.4485465389243579;22.558044850869035;-3.0516086399047593;14.092062218552034;-6.492419288754494;-55.773651578745614;-8.3345385657462394;-7.6507894326459347;-0.046998340257400135;14.105669803991141;-2.0178193383030374;-3.0508242994656602;3.2587279297022929;10.067845417908307;9.9667351407410685];
LW2_1 = [2.0788350841473573 5.0136114543230947 -0.62729096991944544 2.93622752336656 -0.23316424617946591 0.98760715569256585 0.12317550830098468 0.80100967021736313 0.44810268106525147 1.2328892357172863 -0.96029574151876207 1.0179078437419626 -0.98226054222925008 -0.11764397701468562 3.4226241618064401 -0.53428190164592726 0.93164391432242133 -3.0314758752431787 1.3281682327072828 0.11471287293372261 -0.70476898086254547 -1.139859962208817 -1.2263608370560328 -1.8986576982753938 1.2673801045697355;-4.0459167470143411 24.254388239088915 -1.6250060130442168 -6.8488789509697447 -10.795746287022471 0.67920697949669784 7.5211950202865836 0.53606155343956863 -6.1180779146036866 1.7215173749082655 1.0365417497398723 1.8748879547390773 4.5080799608450084 -3.8254717811542442 7.8491311878708698 -5.6297061065070952 -9.9170613726801466 3.3669788252852197 3.5396370755234625 -0.96193848047867458 4.2773971643893507 -0.77144034635352066 0.3004399047480249 0.49026185620721113 2.2969246945593813;-0.46519055095734563 1.0462063910573365 4.0431235467879034 7.5672679919190013 -9.3326881708395426 1.0685761934423526 1.8124006075491563 8.1680355314576012 1.0741621459391208 -4.9342977421930305 2.0207615849040232 5.38413627910223 -2.643502551993071 -0.92492729861018197 0.89104303899555382 -8.7228570808452552 -2.3369818704520995 4.098616667243693 -17.015551872757712 0.71851589734595611 -1.8800039158481288 0.10386319667571429 4.1831779710409558 0.3181310739330413 -6.8578580176682777;-77.065320935178249 -53.684994745729554 4.171946592108573 -58.088979825795683 -16.278011511314176 -34.368456013372402 21.92456362243874 -7.4791735836474436 20.686973416847326 197.23400336261966 4.633856919037183 -57.223397505379111 6.4139424796444118 -6.2195972175119154 54.53401425453707 -14.570654199694326 3.5730749981089995 0.33469831907224606 -10.090805852659541 -38.49540393146868 -43.205176125287721 -37.390964983886185 2.5383837605727835 -30.396403785935401 77.028381522787228;-18.082303447882488 -17.798117400392972 11.792460321021334 0.80033702829015296 2.9959560052408327 -9.3178133752533991 -6.1718693760887096 4.8014360360793926 2.2350975314818387 1.6635336614005327 1.3195356218629275 7.0044110752083002 9.733533500156474 1.0622690883223986 -19.104990017336387 0.19720017977921819 1.7767108060375274 -18.817869833629832 -16.559972046642102 1.5248239836096196 -7.1316329180302684 -16.064601394146301 0.13230646146519223 5.8471034163529412 -3.2302572970279684;5.1645614699731901 -31.410651489170711 6.0825193266972644 -12.32917955046614 -0.74632588346540629 4.8306397830866779 0.83825290184866441 2.2844931499247632 -1.9318093846561155 -4.3054147229478419 27.196115727669742 6.3397894430939399 -13.052755716655534 -3.5778655256197309 -13.320494404752129 9.2767114893075746 3.1861280690115406 -14.081822346626755 -30.966287543194724 -2.3391657158863612 3.9993311630220392 -13.048615123673367 -13.840772653258076 4.0222965830537332 -6.9793946517565093;2.8380073249413176 -2.2020187851783346 -3.0846597166197451 3.6017055182046587 -0.4256785707246496 5.1490462649810445 1.0549957398709486 -2.0011412633308434 1.6969284642464977 -0.72331238766993544 -0.35140950192698539 -3.6810151279010328 -0.11121977726476859 0.65621387622432137 -0.060718937918629624 -0.33500506502306016 1.8549826364544091 0.20199772209830844 -7.2822808260623608 -0.50633413107775549 -1.9813616150234508 4.0356035919483366 1.5912530586353542 4.1056871403118009 -0.025036436428298323;7.0869356749632519 -4.7314024556520904 -2.8802175408649768 5.8149880013451876 1.5310310736069317 2.0570379215788521 2.4090245538368595 -0.2053700214096342 7.7503911022413181 6.5075191439009519 -1.0546302558718799 -9.4170605420366531 3.7570016936275374 -0.056327767703596945 -0.26677681792932489 0.49911056006594662 5.930985784811174 6.4013245058357215 2.7540250848990921 -1.0024402045464389 -3.6697099624131995 5.777484745102301 0.32177708711986763 -5.3679195935793995 -1.4184094301422732;-5.6866221377826998 -5.9494635992201932 -2.8621388339555529 -5.7255201954681985 1.051922172753915 5.667044218962106 1.3729332094245628 1.1582537825313619 2.1648122714331448 -3.8987735667016188 0.45986246101006145 1.4106485954652639 1.4923292901692415 0.10721981226971118 -4.8526100810854214 5.0611226525114059 4.8659692217336037 -2.1782914235014537 -0.57676838918317375 0.0093753668796161226 0.76726828004109315 1.0319905319322928 1.4520118555384591 -1.5525721327159403 -2.4363977732331104;0.42474378761117776 12.983090582095842 -2.1690072589800238 -3.3502835565072062 -10.892073123930563 7.4088448932928506 1.3828994038129452 4.4734644542813058 -0.5314017194294729 -3.5168190220352868 -0.053534253997686895 -3.3549896339413898 -0.33574322417571234 1.8984281406249641 0.70127568385157646 1.1272798627724978 4.3316350109249049 -12.799688556074399 -23.168474285523537 -6.6286673469520805 1.8644835046006889 4.5683539420851442 0.71982459538683929 -22.209880779978885 2.7353297786726505;-0.15080348077600697 -0.3129166825987208 0.18748650133925943 -0.15296386304402809 -0.31001747100712407 -0.36951191764184893 -0.53872920522481282 -1.0973523862662402 -0.26410914122531387 0.14125961738095918 -0.3961386374272689 -0.42936838537346483 -0.8627941380170252 0.13662840676402777 2.934484904456566 -0.99989744969218697 -0.36833064006843652 -0.5414208715756258 -1.4616196566403725 0.18737919595705532 -0.11929804364151547 1.0757355605721912 -0.18580443448784545 -0.69193042178377817 0.57143548222944351;0.1323443424261791 -0.37606906292271275 0.13695808454402425 0.11346903325508556 -0.36192613226245773 -0.3080657048005695 -0.59060224436148812 -1.4449516976076879 -0.14936487292091133 0.05994839868330501 -0.46840904369947428 -0.89467523683949779 -0.97628955334350676 0.16728686953107536 3.3487165674006558 -1.1863592580556883 -0.44004308210402582 -0.031797330183896695 -1.1645069444957006 0.18332983327817992 -0.090250339165815346 1.3407022632697587 -0.1523389865909294 -0.45307124163336643 0.68331849686255508;-2.7126386496106316 2.9810289633563776 3.3860560237861885 0.42390625376566687 -0.22081086833555613 -3.8830188905558631 -2.3175822201878002 -3.2936789446686507 0.13031508423340629 -11.708761181157032 0.63055132790305257 3.1313534278446071 0.088585082713572419 0.14545925687526542 -0.95391762302704075 0.79832197963018259 7.4621126245711817 10.519059739276296 15.98125148349078 1.0988074102795804 -7.9290600798824915 6.7077132326569568 -1.5727983999571273 2.7201900059249287 2.7702672050041515;3.1443850603127492 5.5257432948184251 4.7385402754647403 1.8746808066866723 8.9645837480373753 -10.411197098552341 1.9476084773586624 8.8032282062105445 6.8023787979003094 3.2764927738708298 6.7498637044457643 -7.0969923876269796 -0.69798306132348897 -7.2679196026883979 -4.5467529246602023 5.0545169539958126 -1.0884033676003351 -4.2606307364747487 11.979692490592313 -0.17813245681479217 -4.075931410296894 -6.0669773130179321 -4.1979881381626685 2.388084255220166 2.9043637640475279;63.712184301608488 -7.0386686593906953 -4.1829291015873196 -5.8988349431890335 25.820874478073506 1.4796255571991348 -60.841032358958664 -15.149155167741187 -12.085940191213703 21.029838484884227 -10.891248148195881 14.423065591327495 -21.018440229513363 -19.701218624902719 6.8037849916496613 7.1046117246064018 -17.290783428158338 19.399558811312847 -4.3346891441729829 25.377812076754275 -36.229986976169634 11.095499903342896 15.858622255384057 7.3910756703241436 -3.5333204640452816];

% Layer 3
b3 = 0.1870224896643122;
LW3_2 = [0.15873383122705154 -0.015106888310850965 0.0075783651025080508 -0.002434504862757723 -0.0072594542281042233 -0.0037900161963029446 -0.12459403157189611 -0.023618258300121653 0.048948958936404295 0.013094506259731908 -4.5512107522440211 3.84784020456943 0.020635126771801831 -0.015829217472152431 0.0024440761074528789];

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
