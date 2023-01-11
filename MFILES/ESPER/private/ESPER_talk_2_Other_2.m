function [Y,Xf,Af] = ESPER_talk_2_Other_2(X,~,~)
%ESPER_TALK_2_OTHER_2 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 12-Jan-2021 17:03:44.
% 
% [Y] = ESPER_talk_2_Other_2(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999540494077341;-75.751;0;29.715;-2.05132667284751;-0.9;-0.1386];
x1_step1.gain = [1.00000001035699;1.00022983623563;0.0140762793578401;0.000302617642608564;0.251952632905014;0.0584431850247228;0.0422119037568594;0.00826143816767911];
x1_step1.ymin = -1;

% Layer 1
b1 = [1.2854859797072324;-1.1347358878588238;-0.065571192041260884;-0.9277183906944988;0.83207206053008986;0.5899099433466124;-0.2434368258954579;1.0383405419927669;0.73474301957891441;0.26091900320582562;1.4622015892346272;-0.38033597407375075;3.0961320565243309;-0.10243630382272219;-0.14727054136525941;0.069592124700848079;0.79711204333026797;-0.10560819526810108;0.018631118441760276;-0.90168255130573927];
IW1_1 = [-0.61489067211326165 0.83340985849295723 0.35890264985809611 0.34397541433700884 -0.012887371245302539 1.013394457962197 0.79264959141823266 0.46959112560637922;0.42576781850460743 0.21862629849028592 0.85041257515766711 -0.42592691995261067 -0.048929910695594506 -0.88898853208902351 -0.34587112387179131 -0.011180979422195196;-0.085377506558552002 0.02432230179932323 -0.34923205743782265 0.16033285719190338 0.85448643227772303 -1.2511060734722201 -0.79445229864176514 -0.3135724165576913;0.14737979728241629 0.64858731459705787 -0.29577591843474255 -0.31636481045873355 -1.4177880244785561 -0.38920958646222303 0.12541957126653053 -0.46393526458492246;-0.063938312407798917 0.064559233320318404 1.1809669282105459 0.044621781108982642 -0.85724355369015981 0.23542523139407992 -0.30042306862512691 0.1552270150128138;-1.2536990060880673 0.7389663786284254 -1.1355154818441342 -0.48505660199671413 0.37045824767254698 0.47645306467800452 1.1642796564105553 0.37849179201427446;-0.12819740914752584 -0.038267003283622855 -0.60129494376016501 0.05850207297864541 0.28520309472856209 -1.1403306824874078 -0.93053984640193121 -0.51659217154877524;-0.29929157720346811 0.042306847730451451 -0.54747425808510852 0.15097365958660214 -1.6289482137243068 0.2032592771834166 0.065264451307403562 0.75562741643242803;-0.56426187023807128 -0.073202718270652847 3.2710886883863934 -0.65170068790797142 1.9144835819777277 -0.00030687677492434798 1.3933278767075599 -1.0340132344024726;0.040733422122454899 0.13935445647196945 -0.043945937306288041 -0.021923001633454051 0.65759021700853448 0.40057682764978314 -0.086965337058722753 -0.38679406588343535;-0.24346340557709339 -0.28799014817816304 -0.40644763842239845 -0.87334734320318286 -3.389946149435922 0.92321119914294558 -0.33807772007326442 0.22688640519010433;0.97622767208749239 -0.6222371911318777 1.0295775259497897 0.12776108699706509 -0.11144026940931098 -0.021591691094465001 -0.7223724485905324 0.13055640943847105;-0.38011228452812762 0.14324015905339332 -2.1954009220192918 -0.48837428258189636 -7.2688430043967784 1.8905665546130677 -0.7755088031074383 0.21308713785072689;0.56837660710818583 0.30979401902810699 0.040040135148789208 -0.48379307834860802 1.3737494475770506 -0.41178526620189521 -0.18256797868977889 -1.0691432451645713;0.40017150154243103 0.40885665604968396 -0.35742248293733891 -0.68045819810263808 -0.50702942101658677 -0.97875659210619081 0.019135768760153989 -0.63006651678584857;0.093609276838583777 -0.012164646299945875 1.9782371523200437 0.12887372487472343 0.10089893754293894 0.4954248346141697 0.55991593715866261 -0.064215799046163882;0.00053967454704384737 0.12755940524726378 0.16081336064765875 0.030478530818404535 1.080399629377073 0.21787108704480795 -0.14875155139612112 -0.1664294261895054;-0.10806984632527715 -0.67059829295569329 0.91055712961179969 0.61027492776716563 0.52523472622601242 0.77028486106157124 -0.11976833558630667 -0.4792076144735184;0.41007854791366888 0.3354550328067864 -0.21907983970099573 -0.61228265978768126 0.39043949138116169 -0.53877501490322999 0.15373028398307678 -0.2136726911044062;-0.066303638064183448 0.68859155486666834 -0.69179099299913904 -0.47153039420038073 -1.8222446645590302 0.083801544166806533 0.37312003337222971 -0.092714987691098955];

% Layer 2
b2 = [-3.3553958632105636;5.1573784902862014;7.9713039392104754;-0.39978859693645985;-4.4119817274626252;3.6987112470368535;-0.78694016557488033;-1.4555804504225851;0.019218445679018741;-4.3283216326270928;9.5711172458439755;3.0880078844071499;4.6756782905606791;9.1850459213392437;-6.9131381624556836;5.5764296553278072;-0.94802115636404749;4.0226981573478904;-0.96046274242255825;5.3649116699477082];
LW2_1 = [-0.15802497918209524 -1.1890776483345489 -1.2411223157963633 -5.7027734176924856 -1.8821322414509996 -1.9801749702918805 2.1228182562854934 0.87039439145961828 1.3399118158677674 1.5088839938704777 -0.048067743777941797 -0.076486853079112571 0.57314201378192731 2.7041395723150417 -5.4877743693925209 2.1960126902225761 -0.012829856421139719 -2.0980992113278165 6.5385003993878037 3.5750236490883589;2.2293254883737847 3.5016081565774959 -4.3900090058765429 -1.8126407396327915 -2.4363831229395281 -2.5345086154323844 5.8965300869931845 -1.632096451246364 0.50178503248804329 -5.3988390850250481 2.2812614365943422 -2.4916494797687214 -1.6477794412951003 4.6672288150407404 -2.5727980942829696 0.78513556216322089 -3.6672339380427856 -2.6110606302654835 -3.7789077094666101 0.79356807019435993;-0.61028849549414499 1.2910914794122141 -5.5725822952042821 4.670327064427239 -1.8078736281375267 -0.73540866967076834 7.8411908536745116 2.5772370145367947 1.2440612427488691 -1.5389819250073851 0.38412910084933521 -2.2876268737524388 -0.48144190273953635 -0.50743882838953702 6.9342088065985736 2.3515477152721949 0.97276063490637044 4.2793545117988323 -11.307956271508546 -3.6064073776240932;1.8835741407308677 8.3355825276977527 -5.3415817873002016 -4.020056478394995 -5.5394050267339159 -2.4463325109965481 4.8432606248659571 5.3575235432317783 2.255904876322985 6.5178214093654594 5.4093628714200115 1.9228748313332957 -3.5202104190663066 0.81002199163729194 2.7005496480520215 -0.26747323766675307 4.6744602731109453 -5.0341343711165525 -7.924332009941768 1.9731484635145935;-0.97836595454237241 0.49155832167591662 3.963630282631649 2.3545054882520251 -0.013929510657020238 -0.33918234693308952 -3.7451064312001012 4.5004592760636122 0.12051961572488601 6.057766833000886 0.84090815207832237 -1.3294867170396822 0.0134915584042757 -2.6750143309541574 -1.6956105641804209 0.25463256449004235 0.47430919958257012 2.3985394650787542 5.5648042672443054 -2.716132540704276;1.8993742983142241 9.0063425035046087 6.2481633325227932 -0.11357951806491623 5.4645549156472644 -2.0338799247253476 -7.9549335275240143 2.5972847922561013 0.41185456583721886 9.6295119080514233 -2.6131168628718346 -6.0135888306768734 0.97812194872100777 4.2117439386160465 -4.4192578781552854 -6.6693625889830281 -14.270092505486895 2.0704093180580818 5.7418249763860345 -3.7639915915486015;-0.55723936466108226 -0.95437197122125794 0.35587513588221792 1.987061951691335 -0.088034458494449019 0.72958331043768243 0.86236742919297982 0.38386573833521054 -0.19336758440869117 1.4746234698115972 -0.11712029632973713 0.1115356312501929 -1.0271536210429131 -1.8254744011456567 0.90147132113966877 0.65693259403696125 2.3510051817723139 0.19886860059990122 -0.75433174644654832 -1.2283117996699529;-0.4908857288494044 -0.41765043125107276 -0.12669360433891541 0.82454502482767777 -0.65448289625498779 0.42617704686272262 1.3021924986367492 -0.034403832696200516 0.074112871332626035 -0.41294769968980799 -0.29929157528498168 -0.11331990096833885 -0.64539040048969132 -1.0473720068742178 0.65747021951022466 1.0136980880471995 3.7894816465416312 0.073782973977246868 -0.60952429711448586 -0.34567100136672591;-0.37838091994807194 1.5233896188347291 2.3637303004629477 1.8549795372298827 -1.0321270576281274 0.63518598637722656 -2.2418305040901076 3.729475202993755 0.53125676229211538 8.0116807365832212 0.55301517285585133 -0.077116819601216122 -0.050363568759736202 -1.9754311298123202 0.90045149232085542 0.79924152037611451 -5.588074109216163 2.2700278475274098 1.9344705477951558 -2.4405864070340986;-1.6543549158033839 -5.4879939040445098 5.186287233703915 -0.20710236586464803 4.1127455043046863 8.5677501558514191 -5.517313806032452 -0.99671485908575752 2.8259534426882835 3.0933883033698351 -4.1480394542757502 11.384122388153932 2.1047614528679923 0.28104527882843899 5.6614397850187625 -2.0147097509563059 -6.6065963242633634 -0.19903999676216957 -1.0456543996556764 1.8546268795811427;8.3256576356589562 10.668553010692325 -5.491567552401972 -3.2493815756684921 -5.3314578000733253 0.11730133923315879 0.74291175381030394 -1.7562419178244231 -0.3476161667728474 3.8277328260296355 10.437027083711813 6.201114331598224 -6.6863679474385833 -2.4683880178396977 -7.1502855935593752 -7.5165393060377861 -2.8997403401589059 13.338253062253933 1.6900617168473151 1.8186279288599239;-4.1954135944882172 -4.2348425372812954 -8.4866020250371061 -4.0436798339493452 4.0723240400971337 -4.6774489650730064 10.729052339059711 1.317799174567341 0.89714286463485593 -3.4835280105673005 -1.1448679853165522 -5.2481912847311873 0.77200211964897913 1.6465644742185928 -13.56120962707519 2.2187337028265248 -8.6908437180191882 -4.9001835758808392 16.993929700094949 -0.0031181740416631382;-0.28284209575370056 6.5179704941524648 1.2432249704787586 -9.1643479507653911 1.7061651902192969 2.5265862510286108 -5.5900006005516323 0.089802270404777881 -0.25367970145504493 7.5302051845758129 1.1041684849065476 6.2515611282671282 1.820106326759479 0.9709701814603372 -1.0396977748861418 -3.1813574600871757 -3.7121106442265912 -0.82431628478983576 2.3026627976688134 7.4306971002576327;2.60428295695176 -0.1971130428920117 -0.80752978087309635 -9.487569125767898 1.7687268118550858 -4.0729225902143584 2.9700710098326546 -8.4517262114979701 -2.722284026091625 -10.637579392219262 -0.64274182185895146 -6.2094419890621744 2.9740357014628795 0.72627214151531538 5.9925159778024293 3.8641415395122967 -7.9722005473861426 0.58982258154770362 -2.2160791127356707 9.1322298861244295;-0.94080160969702986 -1.6104239682432679 -0.77710514212354798 1.5526949440703066 -0.70993727360462788 -1.9212333382812485 1.9925131712806829 -0.77647476378130498 0.64337758022505898 -6.3559447660618771 0.93755215029417382 -3.5846053350182352 -0.54053558872255147 -1.0237772804755485 0.88597708118481344 1.6902267032557734 11.45449144488321 0.71733688060409995 0.39890256448575856 -0.24328164757775503;0.96692317695391794 13.51814820120798 -1.1144278954678368 -1.975573003164391 -5.2764899099871112 -8.5057200379361362 -1.42231351850004 -2.1974709282430989 -2.3429093547376878 -6.2380353558309203 -1.0927856529346094 -11.428123264696795 10.426624393662218 5.2993215776760136 -4.7881273931283364 -4.434247524636068 -11.385895327675833 6.3463699704419048 -1.4587872987134254 -3.26898355288307;2.3586347775723002 0.43000447859049484 -5.1626562786129488 0.96380814929394476 -0.6462060651775815 2.7803415568027656 7.6413834893450447 -0.33255989329734309 -2.6726738680864539 -4.3341008015430234 2.1330386939456987 7.2098666716089435 -2.6086168466864823 3.3237445583270229 4.5602597749364362 -0.68276540819107867 2.606662703926979 -4.6381853286240418 -7.3172336182068554 -2.3144120832863719;4.2030738232558935 -5.1851491356775705 -4.7304215580067419 4.1606259132433046 -3.8213896129238494 -0.21729879011092787 7.4771439493517304 -5.5661039717428302 1.4125300262769556 -15.292480683685561 1.3768811306153335 13.047647905610534 -0.89429827183110733 0.761274417474613 -1.0219579361741953 0.12451349633392322 8.0314237796362367 -2.5242683223949518 2.3966699649102621 0.16892647047085269;0.70698943927546209 -5.223332454812053 3.0491683902830098 0.11738149756179815 0.56768209298289729 -0.3635200829538075 0.27774098614493559 -0.28209505459755491 -1.0778441848063947 1.2588044802664871 3.5504104274430834 2.4007189557838196 -3.1572733385857732 5.4465156694401751 3.2415194844600106 1.2421828663576735 -5.3862185134821425 -6.3669437896118977 -12.876723257985923 -2.7434924265783831;0.27215157950918484 0.64633113684188048 9.8757455715351128 -1.8700256125126622 3.2488397916825011 2.4174087812657374 -10.679029966339961 0.15464589563461809 0.53864746075566394 11.061915827963547 -0.29251459475124153 1.9666229895968983 -0.058867466574608715 -1.4192082513628035 2.6173987498196367 -3.8650845891098169 -16.327938243965459 0.62373507066208833 -2.3540980581415192 0.701193743085606];

% Layer 3
b3 = 0.1318340868779323;
LW3_2 = [0.056513582652450753 0.079210561362681731 -0.063837782280746047 -0.015609392761487794 0.33171875042607235 -0.029713972044320799 -0.49030845865908479 0.77184802387144547 -0.33068319511625327 -0.064520410912202578 0.0057290846012256515 0.027523940017655887 -0.01777629242961571 0.027671220517575403 -0.20522299479648676 -0.021436929833662934 -0.039934464659410085 0.02059389458154072 -0.027972486955673143 0.10618143101896919];

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
