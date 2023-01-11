function [Y,Xf,Af] = ESPER_TA_1_Atl_2(X,~,~)
%ESPER_TA_1_ATL_2 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:29:59.
% 
% [Y] = ESPER_TA_1_Atl_2(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timesteps
%   Each X{1,ts} = 9xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 1xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-0.999999996490807;-0.999999999390765;-40;0;11.29;-1.99149636107278;-0.2178;-143.241523403244;-0.12];
x1_step1.gain = [1.0000000017546;1.00031699124092;0.0153846153846154;0.000310077519379845;0.0705602618858232;0.0604107502911417;0.0479365703301392;0.00458552991639557;0.0145655815308426];
x1_step1.ymin = -1;

% Layer 1
b1 = [2.2402647306240015546;0.54189121585852162433;-0.89751948053560659524;6.5293654836180730072;0.81890413060013644841;-1.3489283441698394217;-0.94942062891410849979;-1.4971153811496151764;-2.7165890716632468127;-0.42412889348037435822;0.014955679396662585959;-2.5607889138631989567;-1.3102007382343916486;2.2934099900231728064;-0.94832648871220492914;-2.1307972952648732701;-1.1824508160511242316;1.4703760271016230732;3.3550202667430726322;-1.0165518117515395247];
IW1_1 = [-1.6145449097129271898 1.1783878156647558111 0.72549426339189593271 -0.35708760557621721521 1.2125700153315996932 0.62437228699882207916 -1.4400897987786496213 2.1276139264389648531 -0.11730995315519131272;-1.2501742783032856465 1.3418138915914121778 -0.31259879211672675359 0.29309819573511708413 3.0513113932854620813 0.013861966913817830305 0.44561740978419384085 0.096043169382489149144 0.79496658945735387825;-0.093779176322430449209 0.28296264454630171148 -0.3236774279876880045 -0.22366957488248992436 1.407964972355735922 -0.24018722892826274706 0.084142372010384131054 0.37142596189632631765 -0.21175538757214484864;-0.79577951367989563369 1.2725672665056517729 1.3621026376045655137 -1.746385635434558381 -6.9899154651791128856 1.2398697300063425963 0.03689594729045451299 -0.72574420652164284284 1.2634132635907542319;-0.60299115843587913766 -0.2890896091889145203 -0.15943047367322904506 -0.028613413972728110601 1.3675398784318828493 1.8639403379834165175 -0.04471237121901910222 -0.28424820324891564294 -0.013059658314414034769;-0.82660566307789984375 -0.73247216992805097568 -0.19217788121668849532 0.71716598047523527626 3.2236696305584482403 -0.89203208587624471804 1.4717520810610675852 0.78264054693231599469 -0.27008509129746716138;1.2467672209476239331 0.21251992241578179788 0.70918115500317679967 -0.34062797258706795089 0.59385909548150650483 -0.89065622159874779928 0.084776594191854501292 -0.51378743156609085752 -0.67305283704231588793;1.0155699300682246022 -0.94431187156036811015 0.34316076545465601377 0.22341462895013350365 0.35972345950492451427 0.20377560170619024582 -0.4561281246859123506 -0.49512128783766440332 -0.95644134360574606379;-0.41806900893661674834 -0.82130140352815217408 0.55721227503725512253 -1.1046800456926497702 3.5401383177746672537 -0.52325350333444287099 0.43405713186741057275 -1.4771449348663447942 0.12517705591813185761;-0.71271173839263346483 0.51092659740152868153 -0.5549555713342594343 0.50829149311890264862 1.4726041632509347945 -0.4224442566826775991 0.52587705911818638516 0.78255524844561397124 -0.92291219608053831625;-0.078349521927048512704 -0.084341314086029850472 0.30681868109217069751 1.3321002372041628181 0.33108506385495034552 0.91734074844964608175 0.78240163158817332967 1.1398785952596734994 -0.32629912506972641451;-2.0594251659493854056 0.5490604470015555405 1.6166548119283767893 1.2779558077733068622 5.0683742324059402051 0.56061441622833041443 0.73665938824061094614 1.3476649248396588288 1.5393954926940656147;-0.37916739787235481662 0.7764498887876098987 -0.071455144847028667643 -0.56832043314431812675 1.2578063539401933557 0.34331583561106621127 0.27533207039068297028 2.133008628812425922 0.13513137080911666321;0.78866321446415210694 0.28144878332935102661 0.68568694522957318593 -0.30723712704471001222 -3.2843721214376180306 0.29680065957418555467 0.83840460117514969962 -0.98214674220145981032 -1.0365125325860955474;-0.031690039945032853763 0.43220173227087943513 -1.2291254166294551808 -0.16826563363291915931 0.56031770962914295708 0.95484839497207063541 -0.95985402908065553529 -0.85565875966413185161 -0.22872721910069881979;0.08346340821699184398 -0.54994581624563210198 -0.050845621332308771345 -0.34006479557460106422 1.6458102836646659739 -0.85350529893551874139 -0.14556602201138210684 -0.96984739078788206879 -0.058817455922632419507;0.75415879401586927155 0.066743576226690803921 -0.83741680685000929873 -0.1926118119216211344 -3.1201332920365141277 -0.69823533640189194838 -1.0327628378334139381 -0.84504128116952570426 -1.524902275546992314;0.28879921588860490989 0.10382227689343181964 -0.098114608291171245824 0.91373586979676657549 0.83675838136153013469 0.31557415368201502393 -1.3657228316561897952 1.8340133717710123129 1.1616350010218685629;2.1557719174045786303 -2.3906074145680475951 1.9044832443520531839 -0.7714937960407618478 3.5333179469691904373 2.8515449027478401156 1.836338940947497278 -3.2121290534440847253 -1.1081970050861813526;0.19948635627774499723 -0.51430335438284346949 -0.19608907317439802642 -1.0328899970991030077 0.0051643366385317563891 -0.98403890595468690439 0.3834631139825066537 -1.2658624635783819468 -0.37358164006999117213];

% Layer 2
b2 = [-7.3330183204444132627;-0.34120208806275625912;-61.115006821343854426;-22.923170054502691073;6.7890536733859958574;-6.6680257985581015134;271.77418544775702003;386.91581052296663756;3.496654946183661572;124.43569596562889501;-108.35275911409827643;9.2995830031797588333;-27.939341242968342982;-12.079730127584500821;1.7834361871126038768;12.188096394903084274;47.852071519855982729;-4.2865611530463558765;11.780067817478835579;-1.5663581431113560161];
LW2_1 = [2.5476971686452749566 3.8135723297968917578 19.157108098088556858 -6.132245337288333431 -2.3606502004659604133 -9.5401245106523653305 -34.518641205921468895 8.4981258707978533096 -16.280324966509734708 -2.4292055072859501941 3.8130398008038568491 1.2439209475729351428 -2.3201425042659722386 16.419329377785899737 -19.839571225266627863 -9.0777214009027584751 -14.963187152231503418 -20.175216037212901909 9.3917613429754904786 -5.0408032265791646864;-0.016191783748569086449 0.13573398539953002162 0.27674742737188640085 -0.072266188817651955123 0.057123220577816991028 0.01408498236846239017 -0.047074658212490350095 0.030236584258797348512 0.029703327938445951012 -0.20238556758174086614 -0.020317258062507604494 -0.013750944385531649328 -0.012406866259071383171 0.096788345015952514938 -0.038897693548897938065 0.16739476850589060652 -0.0016173078050245484442 0.021287620705077267558 0.20906026926300028124 -0.21826627167415388331;16.042114110079701561 0.29146950975825181018 -21.218454258563227199 -16.967444363179705391 0.55278679772217687027 -11.129180866544725603 12.01519719658974239 -29.77891108179201396 12.127057819446454445 12.360503914418110583 0.6836088272029130497 -21.599515925049992404 1.5445817430025394756 -4.5889656700445042858 23.101890247505334486 4.9440367426002103102 -25.038204748261051691 -22.739351370721191614 37.526889196209218369 -16.877675301395136387;-3.5875642227686079799 -2.5108787689136002896 -4.4811933733991535789 3.3993930154178451097 -0.57712732498748597898 2.6662919534553206624 1.6230394943897252791 4.410257313082011521 0.022073767468434698363 2.6831626012484801613 5.0735820005717560122 4.3511622567414764973 -0.55248379075988707232 -3.1993726016249648403 2.4684588074472491037 -12.277178996230414754 3.89057927786322999 2.3845876354249710971 28.352010479285340239 10.685333327953287608;8.1080020162611265988 14.692068712093124816 21.427379981053562119 2.3232039745205841363 -10.622610365073237304 5.5781259744232398745 13.408292284699321328 -10.419513640464616699 14.317391835488962215 0.12290643635075472473 7.108947712101843841 5.4938479488914495619 -9.6569223107286337182 -4.7182799000087580765 12.864828477137232099 -5.2975472311818823101 16.228176813643905518 9.0500595947726036883 -16.770145899028914727 0.68734249191959295722;0.31120673553851663229 0.26495220743555430998 4.2616362693629268321 0.69383328619630746381 3.3522059918940341205 -1.6391586830045774548 0.77238228194856128894 -1.2154569457562627299 -1.7693422765858144441 -0.76528067140515254607 3.2673133713366815378 -0.49701473707130305879 -1.6769814005807537338 -3.7398304894181149116 -2.7535299501351913065 3.3442340255321534315 -1.827919416466034841 0.30811085537369481413 2.6895382517915344955 4.2292931517340761616;276.82500812726976847 -232.84196753600670604 313.10176562858856641 24.305115714085143708 120.64025448127944173 0.50279908897715830385 -203.02139100118691317 -94.126237917051710724 -139.85517410294235674 2.217379630741672436 105.75282555974081333 -139.59224240338213008 -72.120555170497212316 -144.98182679234506054 -8.7343985419707355078 113.10930513950822274 -22.329311521924658734 -268.66713137854554816 -145.48402270515597934 -85.130128614037261059;0.062323206956516864763 7.4524474035335961375 9.6060936599265822622 0.56618441875297320642 2.173511361204815806 -3.0134263185830669762 4.9906967821473759983 0.25467848683008914579 0.27561511550345751953 0.46484490859180394029 2.676648768267300138 3.2942067728439203833 5.2746329972782701745 -0.22612414679075512614 -2.422558603381149922 18.795339098048703619 1.4467171826344060737 1.9294595852823495097 -388.34880904614720976 -0.98360476241719529966;1.3761276635032044791 -7.2259004136170625898 -0.50503518209640962944 0.35630978841258642564 -3.1239506203922795535 -9.6665466521428200508 -1.1154591151772978996 4.1148599865664987618 -0.44124027508439589429 17.599222113209357587 -7.5304312741452550384 -5.804000754154385433 11.612146685426562698 -11.881441789865929692 -4.1403193023762998592 5.5368650296793688881 -1.1980941750702540638 -0.52675565077657393598 0.65028035973165110306 -6.0783061900234782726;-54.063031674547573857 -30.793597859907499981 85.016800187828877711 5.0900189370158068769 -43.958909139504811492 -11.859898736219847137 -39.263363637578187593 25.882318695766947769 -34.079921025236529886 -14.88280608296318519 19.966041441313329585 10.285996855216245294 18.84617790577365426 28.842029667255999925 22.360953877833534165 116.67354733128459543 -114.68323123550725029 -33.740244780884289355 -57.102745592364911431 -25.567800766703360438;4.3557255602287634488 -1.983239336095855565 2.277152065146597959 6.0951546289102926934 0.9653409093427827159 -5.6366044976708629122 -2.7955250274538996358 1.4733272295794317408 -9.3705741692739099591 3.2671531681764953703 4.826177005116287333 -0.53194809585872149249 5.5292530763349887479 -2.692200262717124648 -0.28603076980820929531 11.109213927395680699 -19.578446722907219879 -0.96619709782783336216 100.7558990652286468 6.8914922076969435949;-17.196752558263650457 5.8315601582429694716 -11.811145930776214996 2.9366644288693373177 -17.424753989270175936 -4.8442129214410316607 -11.042303039125254571 16.583256105732509411 18.543655961191827686 0.3372325855695600394 -4.5327532807237229662 -0.91421590661794782662 3.7435955334879591483 -10.341535844970602298 -5.5634595242707947094 -22.905320019016894406 1.7527527355316636637 10.569820011871399856 6.0638102877179296257 2.8263774013565017285;2.2467546592339218137 5.4840021105767435827 16.468278743065241088 28.503064265614625583 3.6654865990469587267 11.879276328554906783 33.050030356383118146 -3.8080389189692542828 4.3982895738568563715 -0.22938135457051758714 16.211671445753491838 -14.226369964612000985 -9.9743291779632556882 -32.803219668043823276 12.397873846587666335 -4.0072176827276750899 -16.257467849016549621 -5.7859890861149816388 3.562247018774206353 -0.08380848499783544181;-4.8375081629669702821 -11.068525324892874551 -23.971912637740764751 -4.6038792119586799956 -7.4464225446502334549 -6.0481191859465157634 6.0470986579537528982 14.341150052659259018 14.416195999108822434 26.849469049559619549 -15.274658992303074001 3.2359562675139650345 -10.741374891961589455 -21.823667392535302412 -0.055197603705229830728 -45.401740841339588428 13.06087291757894242 -9.3505949817768243548 24.117363168279062791 0.36317316370976060469;1.1652587692864837532 -2.1042221355822681161 4.2338660543679678128 -1.4006082078205248642 1.9419323961757679609 0.035998498420454848323 -2.5101291160096486088 -2.1521061145694035943 -2.1839158641822162821 -2.4627929932141592317 -0.6193922539492999757 -0.66061724924979792206 -2.1991941048902026168 1.0128060695086382736 -2.0883706585147834289 2.3065032950564994252 -0.51558235468077084906 -1.1783762601336331688 -0.55044612894640532197 -5.0533880187391995165;0.78997157544395868456 -4.3940011237851877013 31.246659859136581616 -10.156346765088706263 -4.7461953551628024428 -57.13659041127807825 9.9586328775628523147 11.570926723775249556 13.555698981405388537 26.894676400250947523 2.1790960505618430787 30.872174886015915973 -28.035617091617723418 -15.645744118403111855 3.1850307811520104373 -19.288653221708734975 -40.513310953933242331 -11.110721641447844021 -5.0276681715323086053 16.093414328770506927;-13.570480321832798509 -1.5670151697948571101 20.054872268107484246 1.7236408741694080327 9.0963994696257071126 -2.5984866104834405398 -3.1024085926127979995 -1.7935189623297307371 -0.51559883814531803914 3.972905170018703469 9.2690890605622620058 7.0830739704867848872 -2.8055756621915506344 13.948354942539252121 -2.9341755169829539796 30.708700251368433953 2.889121355128905666 8.3134161712301697378 -44.787045320756270428 -6.4612059423424259563;0.8046128414911738469 -5.1492683261968084096 -7.4902808103192652922 0.45164147990602387939 -1.4374687255614002623 -2.9328296190231544216 2.657043162528756941 -0.3088489630496997207 2.2075088651539944706 11.014333510171564257 -3.5353076006166426204 -1.5985326410059370605 -0.86232380846309408806 -2.0359953285426435343 0.46234708588417733344 -5.083198931868148307 -5.7269591063276301668 -0.73960800792714753715 -4.9258042578757521568 -0.61902797722667091129;19.123476704464803078 -34.853639671500985742 -69.129572399354159984 71.168892502051633642 7.6793500560623160212 102.50537153078298047 51.783109562928501646 108.87281461035904329 16.196728764711281201 2.0085116918343475234 -102.39198508616325967 -16.880345147125936478 153.06314510160123632 -173.35565875549914949 44.807622751674479389 -89.597732588051030689 66.822400779440329188 21.896356601645816653 153.46802771087703832 -66.6176006821890212;-2.7911371337772408197 2.3409253736477566576 -0.89040978906565648732 0.1843088834696507694 -0.32578499994453707922 -0.4548300268189252038 -0.35360412458174012018 1.777825378791699773 1.3475102747686296123 0.64640610903389505992 1.1888599881705868011 -0.57801797374065133628 -3.5863990972139574431 -1.4673841909268616757 -1.0705017624772024121 -3.8141678290474421154 2.1643886980821571697 0.68327471965349562133 0.4882399651815980679 2.353200174631555619];

% Layer 3
b3 = 0.53904706436588967211;
LW3_2 = [-0.0027364762497659711991 1.0996639445410336666 0.00083862397695242600452 0.0079396377789116689278 0.0025190846812410214867 0.011321792688845847344 -0.002274100591642778859 -0.0056114724734415127178 0.0051785255855735946168 0.0017997856510160331091 -0.0013430315666628164523 0.0030818283520481376016 0.0017457767989745005908 -0.0036629007212721120985 -0.025714809155657557954 0.0034613477095804638661 -0.0032587782676394179415 0.0019371822539216365048 0.0011944639441946905238 -0.02257947343245904695];

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
