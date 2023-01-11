function [Y,Xf,Af] = ESPER_talk_8_Other_3(X,~,~)
%ESPER_TALK_8_OTHER_3 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 12-Jan-2021 17:03:47.
% 
% [Y] = ESPER_talk_8_Other_3(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999540494077341;-75.751;0;29.715;-2.05132667284751];
x1_step1.gain = [1.00000001035699;1.00022983623563;0.0140762793578401;0.000302617642608564;0.251952632905014;0.0584431850247228];
x1_step1.ymin = -1;

% Layer 1
b1 = [2.0342618710967102;-0.86242743814777256;2.0054575538995127;1.6412236585065767;1.6013060342962624;-2.6327066075292413;2.6631802179003574;0.042650711989917808;0.051778115085482336;-1.6438174318975296;-0.51781826672893605;-0.5280102677104398;-1.6144437430498018;-1.7272671135017179;-0.30630101774437962;1.465058742285539;-0.77370549507302178;-0.50635704923514691;1.5206520777544192;-0.52320995667268522;-3.0567869510288292;-2.3896402378498482;-2.1203231597446073;-2.6125512132454403;0.18870440769004071];
IW1_1 = [-1.0790231167799262 -0.29885698109570863 1.6333043751657537 0.90109917256271366 1.7343905921591962 0.4468408616395263;0.42152862699838473 0.24412360222050111 -0.25041927013388227 -1.082901000859418 0.54641832442860039 0.4061495491712741;-1.56682457920746 -1.679404505937846 2.7514087411247852 -1.5517370359524885 -5.186959652158559 -0.74927995529785374;0.038565257336850844 -0.38894123451116791 1.5779233606842149 1.4808185008078558 -1.0246103276392695 -0.25371332138497421;-1.2351021047293593 -0.14660858157559156 -0.46319615578737744 0.85257561406403581 1.1143319079527811 1.5498167109486483;0.51777540846434722 -0.95166219340315938 -1.3864376947931911 -0.95397373410349828 -0.52611513863064074 -1.1331376342647355;-0.19534118657633223 -0.54292700556880247 -2.4995237671645296 -0.54458254797548622 -3.1840751689539375 -0.10621841805728156;0.011527033149172909 0.17544603152489549 0.84273785162725345 1.6722822580785597 -1.9380462393494471 -2.1923919068503426;-0.26303797244635324 0.10136181499295287 2.8641344634997474 0.39596030514599567 -0.47655104927917269 -1.8092746556636812;1.6472376964244844 1.1743749803352712 -0.14702343330052364 -1.334472204246403 1.0370479133884087 0.39343774113120938;-0.14700404245128751 0.24688246129227961 0.5347121562905921 -0.79556557132338379 -1.9996097371291539 0.43899503095372439;-0.15034498042682529 -0.02281721699898287 -0.72193149021151137 -0.31298568291729767 2.5448895630863766 -0.8608763906289888;-0.36982668938991381 0.006323004218821082 0.87706143730289099 0.63607309221189212 2.6964681503365169 -1.0068180329371459;0.14507015689822902 0.11120855301151894 0.98691809521013407 -0.57026584570896544 -0.19823120964438809 -2.5117548893299841;-0.46781775492267375 -0.11742006606120708 -0.66046945623480746 0.668376508553855 -1.8127519674425385 0.05180938721413153;0.037438757398644606 0.26130132473184364 2.4960380747837596 0.28604146959609272 -1.215383947878991 1.5776545130511952;0.43267856583631548 0.98353870805242871 1.1495511269597356 0.90638776379436248 0.74407218848009704 1.8020446217689732;-0.48001448338020752 0.87947238797754657 -1.0125752175274221 -1.5547494114352789 -0.79360784900450332 -0.46085879168495852;-0.035099410896002374 -0.028697628627582086 0.65335201541964705 -0.54483962806592379 -0.08436298630119525 2.2297376159852056;1.0152877965629092 -0.8546861679580936 -0.78633714492015294 0.57347067372106575 -0.50292144415511775 -0.93437684302383728;-0.46674571555648714 0.14395402380945338 1.9959162557391417 -1.1474516922804912 -0.9241840272213836 -1.8351951934763795;-1.954079131850629 -0.22038479050073348 -0.49562550620958001 -1.7881370098228275 2.5380871034549584 0.79967142073810127;-0.0162614121947607 -0.088582044326407317 -1.746879861465934 -0.83927559712914401 -4.4878840004093448 -0.51164697258232461;-0.37404885602410615 0.17067095409188435 1.4259658584188615 -2.2553491366870886 -2.2083464290140005 1.34730555582976;-0.34107317283963223 -0.049368903037907502 -0.21855709577837457 -0.72149189177581963 0.18486686573525341 -0.42437719398500889];

% Layer 2
b2 = [0.089504706471897316;-1.231257171941285;3.6607976309152308;-3.3815440397904277;-0.77958341828989508;2.6289858336811949;1.0369873796211999;-1.5152166666570346;-0.4625733873186581;0.20712504307667087;0.045705902750497587;1.4708603201425434;-1.7467138474769428;-2.9291686268513732;1.4999314397554897];
LW2_1 = [0.49083352585416551 0.87872755160471872 -0.089745035227156947 0.0053250681730917966 -1.3284761918003769 0.35061204037355587 -0.44590242368699118 0.52922193677232621 -1.1316807907741271 -0.24160826774874453 -1.6759927249505453 -2.0088104143499077 0.33171670581813767 -1.6481123788279963 -0.99569796967655122 -0.28404859409226407 0.14189673012094633 -0.17230844536376669 -1.5167855175381346 -0.29495245672111975 2.1612497923606289 0.39577553275574556 0.76592156159780156 -0.67420874614385429 1.8380370939724207;-0.16240837008595629 -0.65666682449380476 -0.24787389411644814 -0.21832163260534712 -0.58424929740039389 0.46003359495885104 2.2262836833421811 0.16353652257708748 -2.0576792979841243 -0.013897504867650755 -1.3135382523053867 -2.4595952194827975 1.0965239527048489 -0.50273836608307376 -1.637297543287485 -1.4026085726892583 0.50469449325998172 -0.86218335144871927 -0.54977723697834024 -0.26496729717548173 2.759844994839217 0.30086651568386774 1.0142343643253378 0.056659566009517347 1.9339319051421047;0.10764272688982618 -2.4128184115738418 0.74550944804542341 1.9655030865640268 -5.9450140908985025 -4.0001954306448333 0.56840128477938734 -4.3780077024589481 -3.5410536760608857 0.82638668427644335 -1.4905566967371762 -0.76433208072449932 1.4970349620829924 -0.31770670812117308 0.093405975410708747 0.56509723458751582 1.8206702511081179 -1.6309567840566703 -4.438699574758461 2.8993349994684166 1.4824892052682663 1.3345450848544207 -2.5589859630469047 0.97878838000861301 3.2897537745219569;-0.0020633645738455525 -0.76256778508097367 0.39194430912842498 0.040951772341628126 0.096474926378673151 -1.6085704944676342 -0.17666720533717359 -1.1523556837997542 -0.22655071637109625 0.13904127867014843 -2.3231300275788835 -1.3115345489230368 1.0037308160531666 2.0361951348891791 -1.7976029846834007 -0.91350940330316166 0.94058633402939806 0.038205877627932969 1.1972284521245402 0.64032403725556142 -0.57497001035272866 -0.2111262071814837 -1.0822978311812457 0.065649956853350386 0.75031044542052538;1.9450760535482281 0.077011041368595767 0.16595604014678966 0.29032219409301735 -1.314353253592216 -0.62092171558438058 0.6876587397936087 -1.6208367409130733 0.16858512532357819 0.69789704994957358 0.16617556515846846 0.642297537143243 1.3951583292306478 -0.86111724937287437 -1.1671906411567781 0.49497090428630236 -0.24254623552513713 0.04200346341661429 -2.6372157726697854 -0.17591253878973134 0.83591367670820904 -0.14188433893002694 1.0203371986824994 0.87932118229066791 -0.03719141684584943;-1.096332605473388 1.4689623522989474 1.4043900034161472 1.7136647048426457 -1.2707612854751624 -0.16786499925212206 -0.61960063355227479 1.0192721363273598 -1.7132193750533242 0.20215156811133636 -0.52499750549604218 -0.57594313695389721 1.7667218669950417 -0.45232444707714881 4.9440748833434478 -1.103113224713645 0.15229676951546345 0.36801580652759347 3.8828958457395681 -1.2010759861228932 0.53311292800622512 1.7372943247662482 -0.36047106787477018 -1.6559888599134247 0.46244283345684445;-2.4493936651370221 -2.6425940157420706 -0.53028254247468809 -0.10192597676907574 2.6407416319058252 0.41390757712201692 -1.7132836922665062 -0.92646931324948156 2.7400414923328422 0.68876639054503652 1.5051801958441853 4.5803338375072968 -4.4817837978858481 1.211867305269589 2.4381113194344404 2.7509775998749926 0.15777013489815733 -0.039405654068227478 2.2657282995459629 0.93602030719700524 -1.9263301436573641 0.31977919970305002 0.052282490382843348 0.9438536653132189 0.069830834237519263;-0.084362125678453442 0.4784235265917165 0.058878480004426637 0.19505489049088381 0.11303831278331972 0.022604053428679476 0.39153551382024238 0.087384700036621216 0.093383659773157768 -0.23328325395881158 -0.70655102929512537 0.037275622763445357 0.81854135774632653 0.62863460811660088 -0.25865696662571613 -0.27428105740596304 0.35721141953126767 -0.26374413388678219 0.1179908799753131 -0.28297062144738178 -0.37423294582781458 -0.075595925321962118 -1.5090198612640553 0.18942050496618612 -0.043779141441950567;1.9845610502291642 2.341571455788666 -0.56552938417648546 -1.0679992080189846 0.76475279006099639 -1.5323228868560226 1.5339802330243624 0.21145827746587093 1.8261935016555493 0.38606706315618189 -0.72613138282466694 -0.39118567356671513 2.5875396475423038 -0.90568657460442337 -2.1516553023316374 -1.1396042680206402 0.17152039458357823 0.20519767676147746 -2.2696301592962107 2.4058767666555858 0.98684176000113566 -0.43843174540195051 1.6229272378022481 0.81682679855501028 1.7142890852764276;-0.77118519795362905 0.68282087684203474 -0.75615183079741422 -0.40216173792338167 1.3116171197076913 0.26889180513589217 1.0341906208569585 1.0176738792333071 2.1894395565249676 -0.34191385901812132 -2.4423347019039303 -1.7144474681826667 -1.7384601284456622 -1.7506421177262106 -0.0072892008243525958 0.45393892825073023 -0.1000594441038842 -0.733310059157262 -0.0190036842319307 0.77144046428900037 1.0682138345179142 0.18973887389525601 0.93607980407751323 -1.9248422303675914 -0.73738283463183973;0.14600720238478715 0.60139016228446585 0.10586173372808687 0.42207746865772372 -0.058110590034879749 -0.37637702897321063 -0.3790577626668396 -0.99608957996425695 -0.13549590053330415 0.17747051150438797 -1.2749674703745624 -0.77884942269050739 -0.21651339735534489 -1.9466235094615689 -0.081919590057577416 -0.46841182965761191 0.26759823309186043 1.2360119854794072 -1.7778872899510525 0.48298216963208784 1.6799328784096657 0.25747815055259504 0.34507350933470216 -0.74275272229011913 -0.044240382119621467;-1.4931578503773597 0.42903299554658886 -0.52812669203820806 -1.1026372678143321 1.8486668376365729 0.61889705461154576 -0.84530634693392626 2.346623320180973 0.40812774869979679 -0.47963172630743406 -1.1726115729757687 -1.4752502949339283 -1.8031283695924603 0.26002674468789649 2.0880251751496166 0.22320157040475314 -0.63042270008401591 0.032146355885907293 2.2212413227705592 0.91432112595511184 -0.84575119087285167 0.13854925178426705 -0.60924608782564693 -0.18735305711791589 0.53829617118797768;2.2345175837029547 -0.68877498028818229 0.068581960646208834 -0.54227238436595016 -1.1940160240524815 -0.85957954365436273 0.10515819135742327 -1.4967612214324615 0.37422457046938579 0.73596388083438846 -0.5269620925023607 0.19503371341643266 0.16333445204459371 -0.87318905958488979 -1.9883548621646545 0.9715701959842229 -0.68123286130936878 0.269599137610049 -2.9019576616541656 0.20754513780375833 0.70233919161843428 -0.16470401512717225 1.1791474516748979 1.2631760981906572 -0.47791155559744347;-1.7917462942082834 1.9086221484009667 0.21734307463996022 1.9054693009835562 1.299725847435244 -0.068456598079385492 0.79721784820398445 -0.99962628860148672 0.26074331971124803 -0.66494178270364879 -1.4743265905392253 -0.63248133772020032 0.66939549159424205 -0.33102063826685679 -0.59356078813223456 -0.30827481439240934 -0.044472362512430745 -0.14397091058147501 -0.86254863220493827 -0.44775890379326344 0.72091379968826841 -0.68258271112350943 -0.017825499073260566 0.21506007525345031 -0.37435698380925292;1.9448410938548442 1.9331396149398277 0.56544937648354421 0.28958664361344005 -1.2279688068078944 -1.1819492450462026 -0.6902961913757002 0.79973301106034311 -1.025538124039125 -0.15567359080818902 1.8078749810575929 -0.60335042081589141 1.7006029407137204 -1.3743884451314148 -0.035064164273765158 -1.8069820322210455 -0.18615502288485758 0.86598148775023898 -1.5201344505443262 0.066700509037230904 0.86347665347531832 0.40933066027992548 -0.88611730433953995 0.88286982706925821 -0.22690061991493343];

% Layer 3
b3 = 2.9267808152199337;
LW3_2 = [0.7264131293102174 -0.60431318445909976 -0.044415900408530551 -0.2531799189498859 -1.7230354429453372 0.040361685961444592 0.28930771151154616 0.82041722060318878 0.10712054906505114 -0.16116447342361215 0.24398272515871894 -0.58400254645016936 1.2974658853674583 2.5755547678066257 -0.15158806063034361];

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