function [Y,Xf,Af] = ESPER_tco2_14_Other_2(X,~,~)
%ESPER_TCO2_14_OTHER_2 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 12-Jan-2021 17:03:59.
% 
% [Y] = ESPER_tco2_14_Other_2(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-1;-0.999998793562952;-77.961;0;29.715;-0.9];
x1_step1.gain = [1.00000000274156;1.00000063368067;0.0138606862425759;0.0002963841138115;0.248478071810163;0.0419287211740042];
x1_step1.ymin = -1;

% Layer 1
b1 = [-2.8727911530134507;-2.0727086162699635;0.73153501187735237;2.036880546277744;-0.13058385665627814;2.2834948876718819;-0.69355356068222473;0.400804139721543;-0.070816371628108596;-1.3021456331738026;-0.4848696235057901;-0.47575974218123085;0.6448566388607484;-2.1239404198197636;-2.8133629209555933;2.4309443338581933;0.97970272168728212;2.763238691711793;-5.0166493647491857;0.74532856167460848];
IW1_1 = [1.7185462587928217 -0.84478802056366376 -1.5106831043333624 -0.058514881298972594 2.7308657878976557 0.52400407613294708;0.50434163119536035 -0.40343906854522754 3.2105505019982337 -0.70898303801007168 -0.4951231969137555 -0.19973288932904179;-0.42384582093780687 -0.29086976986983787 1.418050523757072 -0.010754357816449064 0.19109018430129232 0.16038305215493789;-0.9357379896290331 -0.083984841405611493 -0.76617976797535914 0.16061886951204285 -2.0847406160758815 -1.192846183160992;0.70561777469219311 -0.29040015627997895 3.7052210174689248 0.12951929164315093 0.16214084024971159 -0.082807707845110226;-0.017651554326861673 -0.18054401445874155 -1.0799493358630032 2.5163232440659629 1.5780240657381408 0.87391207653998648;-0.081513278722725679 0.039182240881760401 -1.7586412008737249 0.36369180708845766 -0.59514483199770307 -0.32907069157363972;0.37666630266835643 -0.12721379396274687 1.7722567476563049 -0.56398527334431769 1.7069545723937776 1.284174814355308;0.82551746788350855 0.1696616771732643 -0.61858297961439601 0.50186173288655322 -4.3852369262972335 -0.64722426165091884;-0.12328019489745966 -0.48690786298873029 -3.8711479498989019 0.55318722148807564 1.4132326465629321 -0.28355369282491816;0.14722868601186745 -0.67629659671499864 0.59917696282761557 0.16551436779304152 -0.95829009879929505 0.55738088847680556;0.93191853887634613 0.25981215954346509 0.43736270624145235 -0.37652450520642089 -0.98792713400610255 -0.23770818020870621;1.3895508248193704 0.72017865991858698 1.7861442489530441 -0.59194453584785511 1.4584617854844053 0.17121189224810349;-0.41313108936273774 -0.99596653733919316 0.6466311954375622 0.31653870514626103 1.8100228524512023 0.7466734622727117;-3.185296471674059 1.4496978273897607 2.6525994892895666 -0.56634684457588402 1.4040903244350991 0.20235680454945035;-1.2427130169857745 -0.1444273825409105 -0.50330858471791218 -0.12108011350271722 -1.4162541631300694 -1.6445249754358737;0.93992159373286555 0.74593297369238265 0.38986570865963521 0.17248824869329352 3.0152397331719549 -0.17087884527074929;0.78627085620764392 -0.54868244039772185 -0.25747528824020455 1.1618436816547164 -0.15316019373751363 -0.12439192507253016;0.22526898520113756 -0.12103397435208944 1.4883504749872807 -3.889370206866249 0.30049953026987225 -0.55735041340910707;-0.050227515291347072 0.51138227764285882 0.23683997583036914 0.37335322982813957 0.50123256989656317 -1.3018724270940492];

% Layer 2
b2 = [3.8271076414844103;-0.69643661037520366;9.6992723040482307;-4.1393907205557348;3.3785018431642593;-2.3078085415662613;1.0324921361910349;-5.5171976278607024;-0.50484269544262883;0.03912804681960972;-5.9525869743925144;1.9304447262810558;-3.0786525164928449;-1.3155447318546949;3.4070757399984881;-0.20694936502479883;-5.7609486051529872;0.75280434330977686;-4.0417725812021299;-3.5440445469365103];
LW2_1 = [-0.073744541561250501 2.3161966510581746 -1.01533757792663 -1.3852472917657057 -0.1283352950404737 2.3831419121375395 8.7715597890958286 4.2467456806317774 3.3858166539935897 0.15901025374841693 4.5130629597292362 -0.66527926312517105 5.1443349253996971 -2.3717301228889696 -0.60827241729643666 2.6278141972945335 -1.9559763417571632 -2.4633188206634227 -1.1317519216345713 1.894148982739571;-0.28138906195927404 1.1578082806617589 -1.989018564383497 -2.1853084716448987 0.25129015485403328 0.46762471356830465 -0.71165952679367361 0.064740345210492128 0.78567372784931577 0.16442158571054782 0.67124339927038423 -2.1522370435536069 -0.079472854819729363 -1.3918063335022755 -0.0556683320930563 0.62644661489250997 0.55978981377454406 1.4557111911743972 0.27474938413715716 -0.50531156929067367;0.33732527631064785 -5.3930087418592052 2.4692673447415796 -0.91983983778787692 0.19245715020266987 4.3127163716441492 -3.2921976377328441 -5.4294217433465279 -0.94015438395746154 0.14319295332664331 0.3413003045461771 2.2187500402596205 -0.53385507950715916 3.0553252242366296 -0.68247633134518748 1.1897895705360868 -0.19508919728570723 -5.1886140114588981 10.26647416297293 -1.2014530295698285;-0.62659480244181021 -0.18816379298625052 -4.1732021612447259 0.89997769702544006 2.1483171427822261 -0.82841276169492073 -4.1370262954126549 -0.3310375261471849 -0.57493471871532931 -0.083539541024731045 -1.1073880765032229 -1.8521695490249148 -1.0704278652908066 -1.896653956831267 -1.2858730767734867 -0.79802843021096626 0.22326624086586808 1.8160208108152425 1.0890448596901487 -1.8203975613015169;-0.43191985244336051 -1.061584982938254 0.59646548852761927 -3.1380855620945689 -0.82154843601506422 -0.50370459499132214 0.99078268597640096 0.26442437220001308 0.040599094135061145 0.21527236778500614 3.660931335342029 -0.13903694599784966 -0.13038369964299262 -2.5338083156968985 0.20028641584023907 1.1021627709610369 0.31775721547069552 -0.39268678582275574 1.8516497850258959 1.5221658871492925;-0.67507300851587504 0.57726667980949242 2.1237243721769055 0.41960820737111981 -1.7994920936812928 0.057575996068007257 1.8692702452810077 1.5100915942365922 -0.85668931884054877 0.24604120401840565 2.6318350168472424 0.39641732265605822 1.3660740566503917 -3.058650470766981 -0.4508865938321997 -2.0356627078886431 -0.58474223548556681 -0.32087728730043147 -0.75662635792714583 0.93203895464551634;-0.28778542436612814 -0.027789284510527264 -3.0185201261545882 -3.255914070607099 0.11612062827181686 0.42193966740739841 -1.8277123397824711 -0.6207085410372587 0.45659878671261445 0.36580977686878735 3.8846345669771933 -1.6172023437874439 -1.2555909262126572 -2.8260029136099916 0.48754745923144266 1.7089063998868046 0.95781333222142451 0.26234843366554539 0.51951145819168354 0.6042072599698165;-1.6441891046789776 3.6487626418614285 -1.3309054782679306 5.0262826637111795 1.029125414025081 1.1696808685282407 -5.4152333163160549 -0.37176867879328196 0.61302287849224646 1.7932898142942544 -1.5191274978612177 2.7464784735471608 -2.6246215205806589 1.0308106041604854 0.098743121755343 -3.1491135626770808 2.7664423793550093 2.7243748890691228 -2.3373408591203528 -0.62783121426959665;-0.93703115878401155 6.4137984918347435 3.0362254202989178 -4.4388596789730892 0.7206939206641334 -0.79164679530665094 -1.5432738126588157 0.062978583534773799 0.13935103318856076 0.92464896427097254 -1.6354710656804656 1.4899390658867051 -0.23452331960035666 0.1411580230174691 -1.2630199829232078 3.1690520338926569 -0.41043703884870025 -1.9190591527111476 -8.5393585335805948 0.73598041864019315;-1.1872021379964668 -3.5855628137978033 -0.39731029810943969 -0.58419172252562457 0.41107785994242951 0.88549019771756621 0.082265387308093588 0.02954179464280491 -0.46090534975172648 0.84540069376147697 6.1315348089561512 0.54023859346928582 -0.20175206086524011 -3.4062679694171702 -0.39810616488625339 -0.84364659864951774 0.65548431439561139 -0.059938854989731755 2.9874151300849761 1.4976012934452505;-0.54191076745382472 1.5589758649186016 2.4481936510622297 -3.030113100912077 -2.4595361608377257 -0.45694957365847066 -0.93614706050311003 0.1030283915820959 0.074739056466500017 1.4491269716202899 1.4945575266318518 -0.08954695837418257 1.884823698656982 -1.3556041393340588 1.2559859684432269 0.92264087347433787 -0.030250156989368497 3.5970097817503417 -3.9241971408771477 0.89320199216411111;-0.68423467908917046 -2.8667016068282209 0.11948171179389576 -1.6564226576336789 -0.53511731418455344 0.27882482576063899 2.9591975444170027 1.4673570607585056 0.1900163876255343 -0.084744154099424401 2.4337137810969591 -0.77650423956107528 -0.39327028623646071 -1.1493948758069916 -0.026091017742155408 0.70967046236651365 0.56552812108219952 0.10914833607840796 3.3625187819687041 0.7223411639349433;-0.84482415536535449 0.39697176423543296 -2.4814121894810621 0.19406879565172241 0.7538003012849962 1.2016631506754398 -2.5585608792619201 -1.6252458633440305 0.50571271251816796 -0.88484548407063623 1.4922376849425456 -1.9583522506727684 -1.0479294598231563 -3.1596753935417596 -0.5564288587861782 -0.96852724488374631 1.5684470322815489 1.8211679637323213 2.1816257802695116 -1.8511975723938923;0.34430814590745951 1.401256619877673 -3.7752561069929218 3.7502157249597179 2.6434858743240381 0.20843347483636288 1.8167525110720475 0.53495790847227442 -1.6319319489413939 -2.4547350317061101 0.61909932926515054 1.6933006313582695 -3.2312581165075729 0.47154613933599676 -1.0496820727066889 -1.0573659330793825 1.3822753847981992 -5.0018896398432098 -7.5929989402442031 -0.8067056451955239;-1.0225869128174396 -2.7252171545801738 5.7692430488145758 -7.7716918714607477 0.53943973441254467 0.9150888071654093 6.3188750213014613 0.6682281751539334 -0.82057764023055069 -0.47916891435819942 0.3911374141121583 -0.19204319333804132 0.91199061565519102 -0.72711895453320396 -1.9558166059753304 5.6402906764785303 0.8393410418285312 0.45478658882482487 9.2563194771439044 -0.02324945073692836;-0.47936352696529194 4.397261461068851 -1.1138551096188842 5.5157506215648535 -0.5974996697487126 1.8715755805083334 -2.2394746639221377 -0.48640703059431584 -1.8058481194453007 0.78868744499482679 2.0349360656839934 1.9131492102784107 0.36376742738398116 -1.6149325052314041 0.98053544707078388 -3.7845883401322209 -1.8898811009580019 -1.400057864358853 -4.0839567989070371 0.61002167831356735;0.20517933088211068 1.9366183473930316 3.8532456095431851 4.4684792058520388 -0.018809835820308714 1.4660724515177905 2.1344359082049547 1.8605648572729208 -1.0893675390297972 -0.25611401942073286 -2.0179376705773615 1.9876323966534921 -0.5120300318431863 1.5195655846808016 -0.25294985192089325 -3.6445413601837764 2.0106557628659765 2.0479154439025522 -0.74938290044574385 -1.1247524607221109;-0.16074824719427885 0.12071628491610693 0.59047488713361451 0.335277388539925 0.21509116273158854 0.49870530221001275 0.44821880369902189 0.081311490752841276 -0.19728555368896547 0.67466734944344842 2.9506228046791674 1.1621650245990036 -0.25751227242346492 -2.3767788182821774 -0.1939281793005109 -0.02028159373450341 0.72601586048545164 -1.9389544733209882 -0.8361455785918025 0.40721327748136366;1.1518140248902675 0.20994483715874951 1.9219713647711145 1.7678758521118727 0.42409323617739564 0.43134814939313315 -3.4308294709201617 -1.2810320783090834 0.35943096663622504 0.24449191022229583 -1.0278426032915722 -0.4070858173268227 -0.61392809338684751 0.79429209960914204 0.068365734355661401 -1.4807714429590277 0.15532733108917332 -0.033930069085690504 -0.49705152060296842 -0.95958910408945419;-0.44926364595524826 -7.6735556061602619 0.66613028155704135 -2.7460222541888366 -0.65438262147350557 1.7831764152816969 2.6145280813665193 0.38476348970738938 0.042954556838951025 -0.87776905714586728 0.52678006141440448 -1.0762262938510905 0.047172946308904701 -0.20602952509644465 0.27811670549758538 0.57940816195214218 0.15328718962383545 0.84232728977938764 6.1178752615716911 -0.20366153488005348];

% Layer 3
b3 = 2.0039218061600463;
LW3_2 = [-0.044858060426235274 0.61359488014531882 -0.078571281481383343 0.14444897477689211 0.89556217310329134 -0.2104595249527878 -0.44520059502624254 -0.082270106652839667 0.13451070426473447 -0.25998944789376272 0.17202856405020395 0.4603249199431792 -0.29956643471305411 0.089106768598350594 0.05676313304420956 0.16028951654182608 0.29739849486940056 0.64595681508400138 2.6569721572923011 0.37142529996709289];

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