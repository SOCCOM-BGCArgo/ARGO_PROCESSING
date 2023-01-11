function [Y,Xf,Af] = ESPER_phosphate_4_Atl_1(X,~,~)
%ESPER_PHOSPHATE_4_ATL_1 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:27.
% 
% [Y] = ESPER_phosphate_4_Atl_1(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999994516886;-40;0;9.1619;-1.99149636107278;-0.52];
x1_step1.gain = [1.0000000017546;1.00003626070655;0.0153846153846154;0.000309885342423303;0.0656325902039552;0.0604107502911417;0.0145232735458572];
x1_step1.ymin = -1;

% Layer 1
b1 = [-12.692202715234428112;-2.8894761197622838012;-7.015987049682488319;1.8678528660223441271;6.5644591909673710717;-0.23808565291554847909;-0.39949435668365634866;2.4413962930435033627;-5.5637646079115450348;9.2173860694564861973;-1.244470941744289183;-9.8413577522317225998;2.1751323662138060122;-1.5291287967198856634;-3.173221290013847451;5.5541358247702827811;0.59265913447069551445;-1.1426620393932069852;-0.59107321310313298834;2.5554786758772012334;-1.7512559787958510427;-2.1355590177702574728;1.4871534117619822357;-0.56461612723209253506;1.5153787707833530352;1.6021990665384269192;-9.1738732365919872791;-0.79524811942884798821;-1.7047104962383488047;-0.10441903089128770632;7.680675330730426964;-3.8703708902685498572;-3.0907400001966713887;7.1389417590996764318;4.4634365777055702296;-4.1553659071520963053;-2.2353543327125340312;7.758364490061360641;6.0551013656781771388;0.3033909700542641863];
IW1_1 = [13.289051104750759791 -9.4882306215246696723 0.015023337314375611803 -2.0196524823210912913 -7.3647966182565838267 -0.74780708769870007924 -2.6877586807586242301;1.2945439780903991167 -1.4302205545948025733 3.7325489709971853003 0.78514871782281647761 3.1193368577524278606 4.4188176959275260458 1.2789463014989581602;0.38376494254131127137 -0.68735428180008506516 -0.085915129043843047052 0.61156129159470340451 8.3973691004963963991 -1.3181866739316785431 0.73597279993448794233;-0.48682661963752055234 -0.096150086530021133813 -1.1011724089481949029 0.071996462758900314771 -1.1301337709636678497 0.091964054700709507717 0.95351878267884926377;0.44319803076831659761 0.30786747656739277268 -1.2711762439600908792 13.888222974085067563 -0.10509362132551383628 -1.4862927496399256366 -8.6687092401549374898;-0.77260110715428975681 -0.73806535453042820816 -2.6373487218966382706 0.79273192224256805449 -2.1754597156257760027 -0.58235975808116491326 -0.015647941035304806234;0.51145064167192610949 1.4105611868392931019 1.5229348699344704698 -0.7349191270146236965 1.9580740816797415071 -1.1828944050881629035 0.46582578263955060072;-2.1862224542117023596 -3.2761359256473947532 2.4784983056367306808 2.5295822562206926598 -3.617477824018518362 0.087972155794429407472 -3.9281729343563678292;3.2382060696833496749 0.14047739137080900163 0.83476352258423092145 1.5140434571968022226 2.9311179291136131297 1.8871012954895522995 -0.82512831441040146441;-0.14998811604136280673 0.024427997066750539029 -0.10289143593705868396 -6.2751496804152644415 -17.303504534274797066 6.2689889254434376653 2.3015588402850863936;0.14001069189529546932 0.057699731822357561151 0.30140707073813743921 -0.12653200515175910823 -1.5658358864269279476 0.84875698162169010974 -1.9648256961532486464;0.21250902415883946017 -0.61373619846609861472 1.2101225003634343835 1.7885330933682395127 10.24673366354043047 -2.2774390222587239663 -3.4204502423520404264;0.58764866738781462718 -0.95044579842304199868 2.0360617139381318452 0.36338046562028997188 -3.2946328280847541059 -0.006538864487464053632 0.72255384372410280669;0.58958280515587402348 0.082236294600470299931 -0.11352389187556094496 0.54829407515528860717 1.9901086620514112369 -2.8610981671506228174 2.7356919948535556308;0.5038402874693997946 1.876702865761706196 1.0075035433459369472 -1.5759195375707923947 4.5875642322326015687 -1.1746644369395802254 -0.83855705372156963406;-1.0812752634003222063 0.75159215766110831769 0.97370443024691388789 -0.54218806156454135703 -8.2308362443376985595 1.0944505464281562013 -3.0564657145358387069;-0.50996714713000634855 -0.13683830361883789206 -0.57336548387464403298 0.065764096970619717597 -0.2129087841025355532 -0.56444773630169720047 0.76725517366911699657;0.65926757843382510771 -0.88524765562241491956 -0.62601851606973590858 -0.76930971394136982155 1.6594569723768157754 1.1566211174615061985 0.37517488145659083632;0.07122904687278508451 0.54364024462393623072 -0.05434728139755846954 -1.137550842444809307 -1.6814779242746749066 0.4578619780719957455 -2.9605669597094266265;-0.53168280162163683578 -0.44792904550882056691 -2.0009081924811566999 0.42866528598981790132 -6.5043020580839359113 -1.3879357557744871343 -0.12965789739549846149;-0.57231459200427414746 3.0099676251059337417 -0.21304425161324022642 -0.0043867186075967233225 2.5913853969188296311 -1.1703926709459566702 0.024049120939849601142;1.4600346221273603486 -0.5835938049237420211 1.5078685899845456486 -3.7780658279707459712 -8.6224569823707390981 -3.9632723342237112796 -5.409301720069269237;0.20942842392067587132 -0.27626352677128013235 -1.2411273916174465626 -0.62352501781141356219 0.86244726008951788643 1.1628418564652303058 1.7729339816605489411;-0.205946340672992112 0.64349550247691578164 1.3231379315496492044 0.3663842102694933045 0.11523231181941025192 -1.7352011778693474664 -0.0050065793953659931334;-0.38621774887787296837 -0.59267558808900810963 -1.6295710287922919868 0.50425895877180559879 -6.1157898714438525545 -1.4471441845054713315 -1.044169655062345381;-0.028870349014457091785 0.32455422130353722698 -0.37199281318847993161 -0.37947979075468835086 -0.51444560700283947785 1.8579246115107923565 1.3837371467118433355;-1.1404011035291572451 -0.15865811864323373093 -0.47165909714473874814 1.0137614636706202553 9.5151586959569876001 -3.1498097637061639453 0.011278553756790123813;0.085505557517539201862 0.62164382260678618941 -0.066623798016019600277 -1.1243077265309882051 -1.8368871259810157781 0.0653780198660129791 -3.0158558217126509859;-0.087974887727772413393 0.21334311806677769274 1.3105586235273249951 0.48327266151709968867 -0.47730360778087882467 -1.0868239108969464102 -1.6134519878327067133;0.70480972875428604674 1.8199827616736714564 1.4711243149464925395 -0.69813859259048627326 2.0034272732968783082 -1.0334091493149173413 0.63288903566085463748;0.89826684792390110168 0.033248798584098085418 1.3136962984744289251 0.11105993140794219454 -10.009127371500245118 3.3193666959716039422 -0.59983322474167100058;0.050461398667801082851 0.19324501358279855268 -1.2132353795323591594 0.14462198201168577349 1.838763347535762982 0.53742657878274724403 -1.2766011815342079849;0.14459568239365933051 -0.71979499561514503903 -0.20556098886070303666 -0.055199416574282032799 1.5152105718045947302 0.11621922379565834127 -0.068043759680442489346;0.60422382486590020534 0.32254569491191331654 -1.9285935093119950245 -0.40166635636031366996 2.3756682816177923101 -1.4913141950613835185 7.5982279725422525729;-2.6823711356146402451 2.1345623645191613704 1.3991806218079672064 2.3365563801011557388 8.0216755132980175347 -1.0261481732834736302 6.7465972699513159228;0.083484730539902018975 0.32440757751687165911 -0.91030750070350141545 0.15101140347750061088 1.389046189078719129 0.66532364641187557996 -1.8949621178129498755;0.26823856679674890069 -2.699462241064656709 1.7977327787999899389 1.2695158737498437596 -1.0985296919463296561 0.33457436252906064267 -0.043819418107159419096;1.0275174303933913844 0.24248193337014323223 0.29624405830446348675 -1.0391964235670858585 -7.8582113120373335846 2.7522123628896753011 0.015632673206692496159;-1.0185938320217142738 0.88555331256683189611 5.2067012264971728186 -1.8754258267157204454 3.653208000946448486 0.47944743016212248987 4.4104684019152964325;0.93983275491087792108 2.0954980813715637922 0.91264814083638123865 -0.23871526856037253439 2.3335286778300456945 -0.12640421349452643507 0.81233285299675728197];

% Layer 2
b2 = 5.0991466096913322659;
LW2_1 = [-0.019473747996828116164 -0.14071577025283402151 -0.79213681097785448859 2.7593216627348362557 0.089600831240212966411 0.56361654830305418784 1.6562630403358040709 0.070425304614277928206 0.37556068222637178922 -0.10237020958451543973 -2.411561425203879061 -0.19434530789501214265 0.45687733999658186601 -0.45264197636695630012 -0.19402844530339485862 -0.2918852205369342756 -2.5230529424234391378 -1.4977645140529813883 -2.8164618317225600386 -0.93642981110009515522 -0.24147027820449062929 0.046741411885604003396 6.5983302802679046906 -1.0209570415287341216 0.96723376818151540757 1.0329705487922802298 2.6378710007290786521 2.664018361930949208 7.5416338537919607177 -1.2430473085898181562 -0.42269351824104511861 -8.357306155555036753 5.8314473050392985343 -0.41074125395977856279 0.10536284325838950193 8.8598853716865022534 0.33057256595342038308 2.8038502360578503314 -0.33778596442366282115 0.59054705337424495681];

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
    a2 = repmat(b2,1,Q) + LW2_1*a1;
    
    % Output 1
    Y{1,ts} = mapminmax_reverse(a2,y1_step1);
end

% Final Delay States
Xf = cell(1,0);
Af = cell(2,0);

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
