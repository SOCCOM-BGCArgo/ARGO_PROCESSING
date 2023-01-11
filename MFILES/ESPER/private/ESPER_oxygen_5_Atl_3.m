function [Y,Xf,Af] = ESPER_oxygen_5_Atl_3(X,~,~)
%ESPER_OXYGEN_5_ATL_3 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:45.
% 
% [Y] = ESPER_oxygen_5_Atl_3(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999994516886;-40;0;9.1619;-1.99149636107278;-0.03;-0.28];
x1_step1.gain = [1.0000000017546;1.00003626070655;0.0153846153846154;0.000309885342423303;0.0656325902039552;0.0604107502911417;0.604229607250755;0.0470499670650231];
x1_step1.ymin = -1;

% Layer 1
b1 = [-1.1897188932013209595;1.4497921175933266014;5.4693064101152835121;4.7502693425901698632;3.4819594498632424617;-2.4482532802637924085;-0.84740739279324650557;-6.602724349545810334;0.71607955158498026726;-5.0281008994063087769;-4.7940777553335607308;-0.030958908877368636337;-2.8338821512096679101;-0.20945374515644524238;0.74878288133012460559;-1.9339807225459624274;-1.475380411336961517;-3.2872949199414054355;-0.17121807987123113803;2.6344657157243442924;4.0238656943805803934;-1.8443070681862345772;3.0027022151033677666;-1.0298863597201457498;2.1441732517966207716];
IW1_1 = [0.64476743391386270865 -0.42477295634902073607 0.29581791087649428329 -0.21302267320223661207 -0.43745777634640359288 1.0026098055116179353 -1.1114502953323270873 1.2009695963704232202;-1.0754171670662786564 -0.20744310149478342731 0.20125481686192361885 0.46914052452530757842 -0.69309305888824868092 -1.1016464969480872504 -0.54558130267027482585 0.23502894383729922012;-0.010427956309867522644 -0.12521717737855575514 -1.3103558894532785217 0.5062791353222225732 -3.2059555868093205788 -0.17857002400547852816 0.26828736711031009321 1.8521477750993851075;-0.4511244586523365796 0.66524502186088474431 1.0978742969501893434 1.5079500370381615504 -2.0217378638795815426 2.3142186551947649953 -1.3606748417687810626 1.6555915279092436876;-0.15131270348077574539 1.0197093004869508537 -1.5352872345782053731 0.49356729159277384911 0.4763794353816792615 1.5924992256942078761 -0.2847869431660834949 -1.4979090479430521299;1.3656793892213658914 1.069489968753033482 3.8517951991449375093 0.25414389631848288698 -0.34931919850202430844 1.337193693330652744 0.11464418059245973236 0.40307879125286089073;0.83432335332165674657 -0.84899197910995505278 0.081654819591770580467 0.96772855809619640066 -1.5479869138222341896 -1.1601101148085892767 0.0045623346791140370615 1.2086770395192756311;-0.34127584709091368653 -0.16539268092704889623 1.3082689967845682055 0.24937796765067912763 5.7620840848110939092 -2.9165270021077804863 0.014874954754593047029 -1.2569945842354197652;-0.56893329668994874204 0.73552781197141625658 -1.6833036012239499524 -0.18168164924825547235 -0.73920193280423307503 1.8588238246107180451 -0.29882312952020562147 -0.022288002279356025626;1.6989735686818392413 -0.30722790692361162668 2.5176031873584672205 -2.9850500435941218491 0.2526264327393877851 -0.38721848349057563787 3.0671861092284364858 -2.4859538266302387832;-0.31850025900468459383 -1.1010025442416400132 0.94151184347092076177 -1.058810297686520574 2.7692558179553228292 -1.4206179092762645233 -0.49328482980658389101 -0.0081837454229081955243;0.51375598188581406145 -0.21952983941589379602 0.20956945787828365479 0.22193418847012788553 -2.1935915375356485058 1.341610461278446742 -0.2199588847296692351 1.6236333896587222814;-0.40982276923872940078 0.37235726789267442616 -1.0053010770037964505 -0.48352741622740480487 3.3024768095875787211 0.15331339402963789498 -0.79264483003999552935 0.4382147983118352963;-0.53970951450077941924 -2.1110717195460124707 1.5266666760555236415 -0.55763113719395340784 -3.3410269610309835109 -0.18131829943333740984 0.06427374300205310198 -0.71009417468665070317;-0.19759708056997915548 0.12697590001476230515 -0.12673778497289400757 0.30616590165581814009 0.22432009886794213038 0.030216299860885836859 0.53850727591669445005 0.0052306017783554986308;-0.79707346394285938374 0.88983279715273877475 0.66671543631281604547 0.21946495751431857602 0.98689920867601432786 -0.8554660315022745154 2.0453999942029930104 -1.4403855356718062719;0.45678270257570491308 -0.80747779100623717774 1.1217983301353868786 -0.029560520408500986361 3.5420238740811558209 0.75631213049299528262 -0.99292993110612959295 0.20424613175965214418;-0.37917955951993098651 0.12196175342124279672 1.0836987769983119634 2.4572700372137346037 4.7355630314347907373 0.4226266401093017322 -1.6266503525664128826 -0.91576714838434358779;0.24805255078468968488 0.027216730302822234044 0.225604371988691349 -0.95578661510379325783 -2.4469149885216170404 -0.44126248569299503233 -0.70963875755023908543 -1.3259847845141305989;0.69438850534483798072 -0.63920023476991039146 -2.1977490108499164734 0.28451693743589645536 1.3869087968786157905 1.673042863331284158 -0.8483934603659659679 1.1644247407890568891;1.3918916280486539616 -0.81159328185498469388 0.88422093740743445967 1.7476798875123971122 -1.5863768859210463091 0.49354218141027222932 0.31836596486844304099 0.61563795280870370341;0.0018555686850003569321 1.5888371764614648285 0.26516555954735115952 0.29161552957741743164 -1.4675952286502795818 0.92585968772348070654 0.22238538185519771129 -0.42253260754630173057;-1.3220200582608263584 0.37023023844034569363 -0.76607493553672001418 -1.7463614654783339297 -2.0473209251427868516 -0.0051763771361390377782 1.7252670105343437257 0.32623742030759261823;-0.24192114006411460925 -0.71129276231374738426 -1.0300765513605349977 0.045804964836359117475 -0.2258176427900163985 0.2776816979626333115 1.7969212225861501686 -0.80156677720478430871;-0.026858637776191629509 -0.73613528240377734679 0.90839498576188304479 0.21451080288529478235 -0.7324320649829855201 0.60688291600571231399 -0.99238102313000764454 -0.2618573094683475655];

% Layer 2
b2 = [1.6333774292815501639;0.98096795997144459189;1.5397920880008908462;1.1920288999214609937;-0.77051135928882119597;1.9788613090201436773;-0.25890700575862374277;0.68287170801715890711;0.57787553646114109274;0.093842138316858197955;0.24669677570329137661;0.68218470601129510555;0.51167554197536235794;1.1701092325425055396;-2.1893107580439732374];
LW2_1 = [2.4544128237179130281 -0.63233638395559566003 1.9167617010206301487 3.573742763238360709 3.3761688978296895414 -0.67828541250902774085 0.12607438957714772743 1.7564374038202004158 0.70940783902512705161 -0.010976025382538435321 0.6718636971206178865 -1.5304586603165413194 1.6038661303215500098 -0.0028389465270426059162 0.90592578599435757969 0.67071913816308303957 1.521350912331212335 1.6411835491296382639 0.44580241219677246844 3.3166858590360392967 -0.61703307575141397301 -0.54376953676023964679 -0.32968803538044638701 0.62748035996212625509 -0.40142892145102454515;-0.03061065975643869827 -0.82950732382229286355 -1.374236069394274562 0.27115394926539498899 -0.079390407823277886479 -0.022053188917154610549 -0.39150082197224284064 -0.39731682527864103927 -0.21953519917352515289 -0.66827453817375515044 -0.13059694592267515056 1.0458642280133296776 0.32211290716764201525 -0.030772138461777273633 1.8478336946606126467 -0.54365851977173129672 -0.13867231783475045259 -0.50770044900298194612 0.15229786289520277909 -1.9731244187450178718 0.39759757550407598403 0.38930387623702128463 -0.35918328823237405389 0.65963145074743756435 0.53438355737864495776;0.017265047158310068559 0.54646531104171724991 -0.617262889694719874 0.092752854814829560626 0.27773075301926369862 -0.4919880342538804352 -0.19066958594403021765 0.75182915588984722266 0.86373025456747087869 -0.032300392253146820165 0.22210653492321594249 1.4039181716210302753 -1.3725610784959583199 -0.29354793075385865464 -1.1103266744130051435 -1.6422553463322615563 -0.0031857119370667982133 0.1059296117315375596 -0.1256220530831398674 -0.73879002702869012165 -1.4498213458529027875 -1.098569706999758866 0.28711874586062358095 -0.082055535091632950118 0.15216209731440577535;-2.4110926123472373739 0.78491927918793324981 -1.8126351442567258943 -0.80584115787065302428 -1.4457262680130167354 1.4292788176623807406 2.3228239669419976643 0.66717335855846859882 0.99668789668792345804 -0.70054665734246024389 -0.12323869420138583441 0.43704676926165519912 0.98482564425932994645 -0.13953126709137031458 0.42390282450019184379 0.73434845162282680331 -2.0296344302916700109 -1.2654406172213228388 0.97903921164814622369 2.7686464516460937979 1.4277862458679515534 -0.73280339362021351679 0.97844795995280675616 -0.66533631790375880044 1.6522692977849637863;0.66241392054495751296 -0.53455995800939748452 -0.97478804085455639239 -0.33712716892309746441 -0.94416045842615614614 -0.45238618225621235469 0.35719740738135269398 2.6697838006503959996 -0.16629645362041620138 0.0744811797351420668 3.52572413249865102 -0.38070827875415269581 0.50561267605394133451 -2.9898826268892428892 1.3825995105257626872 0.24618650681888426535 2.2346561290981941106 -0.67133738697514488969 -3.5531178163421004434 2.3559878499979247657 -1.142337153306618891 0.4103994994339415503 -0.022439150728762763398 0.9507730069761919367 -1.035923976155973536;-0.44776864380045577896 0.8131286490438537351 -0.038713112728597828816 -0.95158826216295533396 0.82161985910559454105 1.6560462549975867308 1.001254307988404868 -0.036699701523773692324 0.91195836120280937553 -0.89557338880735004594 0.8175058918029149968 -0.18715304348545280799 -0.87983458386529200812 -0.0306803471369720554 -1.6029714084921466988 0.79367271803079297854 -1.9422421592258605205 2.2239320135820404012 0.29362947451350929207 1.5049409329854062101 1.3028310146374124656 -1.8080563002130873862 1.2814873313687125123 -0.6483341148791015307 0.75485781434227139375;0.84886552571982853976 0.25301181455219473593 -2.1095455025646026037 -0.75791866786242201304 -0.88319730460203194067 0.23521772492890455353 0.27911565375112962206 -0.45288772941566241004 0.44906663900867932471 -0.77946025494779669351 0.08356388375044772332 1.100548142071168245 -0.021957211481063339253 -0.038095083445993840887 5.6042322471220753499 -0.86972266515061913417 1.1409244216613525502 -0.71251927577947815617 0.64947243484302219496 -0.024764832271438282085 0.25291726743635806196 -0.34502528532699300312 -0.63679275032572824067 -0.94645978576603939647 0.58308231436793944713;1.7312637641344543749 1.1095415041749281393 -1.0020766031441108446 -0.62058910496340180529 2.8839993191452752974 -0.39844717896934006562 1.3762474854310173811 -0.57037513714479470828 2.2691395341670506092 -0.12093044799279842327 0.47611206312929132967 -0.7425239101601317504 -2.3052349027352869726 -0.42518279196359576311 -1.3141196232884693185 -0.90900459678393408236 0.070394053400497896833 0.14991488147947684695 -1.6642664206513473335 -0.40001075803812163878 -1.1489886755525602346 -0.26560066774630369801 0.42513795593020370367 -0.60724578973332610232 -0.2646663795081877768;3.3588755036000899068 3.1135578737839173513 0.24308885832264948146 0.87991094750378917588 3.0293049561550247084 -0.17931041141743866985 -0.79316404198408430037 -0.7584003982835624269 -1.2761069770934703183 -0.66023602651624391235 -1.3382614934024292808 -1.2055078519472988141 2.0063455047529448194 0.97695124867082328013 -1.0288728176422994309 0.68752556112674745314 -0.83851603357110349179 0.62050587962773473549 -0.40127844440394483394 -3.2307897628833126547 1.1688446239371854674 -0.852046086028542482 1.7214595869447006926 -0.52844660793298170454 0.7168996884543232806;0.37048909819913178332 -0.86731485510167893871 2.9671791932044202511 2.6188873362756912577 -0.77539045868517508797 3.311547278040125164 0.66987802050322187153 -0.046035018388204992923 -3.4568978546957964681 -1.0232873313351911193 0.46551920075255365239 -0.075052862899086969506 1.8202744207793180031 3.0384552348039934877 -0.33721266684130091207 0.097896183803534847501 -1.2769432180425386303 -0.19462733708781329578 1.8478524224117558461 0.39002973861175277426 -0.91999567410192450989 -0.37621035130389129941 0.44612805557694734659 -0.7856008905462209535 0.25274526631379268249;3.0207336067041206107 -0.96018638953425872185 0.38469150153187114949 -1.1004172063543327109 0.76757546313731650489 0.97603849972445244632 -0.16502112879450342442 1.1101553814440936474 -0.59377895416591697231 -0.042568031566339162297 0.35556852945774136687 -2.0866756246577748968 0.010028522690849750454 1.374181797761077295 -2.7856826818460342921 -0.60549514172397223 -1.4502047955055858974 0.29485861819871633793 -2.3126299276450317244 0.31530171824708841388 -1.3132900167937826552 -0.40092461115691718776 0.61600389411588540867 -0.16468167951441237129 2.7390916520567509806;0.46499729747455254003 1.0801177959651111493 1.3240228726017471139 -0.33146045588907763202 0.57740465004682417494 -0.059230432079569121651 0.2596486391122665438 0.56709384928645678592 0.11650835786642228797 0.67285319289062250903 0.2640832766773129614 0.26146033615917052551 -0.49286528736634588332 -0.17466393518369635607 -0.74835296043127896848 0.092801315870249839746 -0.15798053856124522065 0.3008381238174424599 -0.17334398362842232122 2.051849875510057597 -0.17042564050092448569 -0.11399745599638702598 0.25712978587737345126 -0.19428356459554974034 -1.0594761979650542116;1.7169063741881660512 0.70669514891500584408 -0.32712394708709785496 0.20949403565986463027 -0.02893358446529643857 -0.6115802285877388611 -0.62901368906149557336 -0.039544319944827821534 0.83443528312682635484 0.1012518812039536753 -0.49843248843804649484 -0.20594921218216269421 -1.3374864800971280054 -0.17614654306130295147 0.74464980155900195591 -0.71041826875122793261 0.83649715084565268164 0.27891307513708357702 -0.018584324950703501322 -0.3264861446830997771 -0.57954462119580074653 0.87306504320813782538 0.33836048253813544306 -0.1804296424667246046 -0.9802275004345015752;1.5782358745334834627 1.0492000736759135027 -3.3388320900345740405 -0.42428944110768190834 -1.7479214787290062727 -0.65524197316295307303 0.90156500748813761525 0.29947015905001261871 0.1811195553283100057 -0.87933506538587791201 -1.289854468434861845 0.41190021900450340953 -0.02350319214927168443 0.0022803942316297013543 4.5115408233575058361 -1.6570297555238207021 -0.57368038430282075613 -0.90643951811878531455 0.62251524925786749787 -1.3659314421930743855 1.1818606489211815891 0.120236185518278485 0.31660006356884884404 -0.34621168099564003651 1.1518535393721842652;0.21779197415164677776 -0.67991112073063542098 -0.10351134329894588759 0.22389874733591377809 -0.53659023327152377547 0.3452688444699136272 -0.578462147270774385 -0.21663661585343044202 -0.4107930409260811544 -0.20494655465659583626 -0.8145258370830932293 0.80467524777600196728 0.50940221758963977194 0.06371735922326417223 1.1815777738573640665 0.45715579292426689051 -0.34053525766318948342 -0.39835620225456103149 -0.070642893126398029202 -0.043749221053504365464 1.9128193078812056882 0.30876420544253457834 -0.39079062917711704372 0.52934953235964943286 0.65776612541910250975];

% Layer 3
b3 = 0.10863686684975956687;
LW3_2 = [0.07381632727744723288 -1.4096420301449001933 0.84267794003168128292 -0.099472381182444710723 -0.13506297248662660349 1.5627920427521375668 -0.51247195885818275229 0.37150556410293433807 -1.0998054256944327367 0.032582302225616491254 -0.13337549241322996241 -0.95766067831255796694 -0.70488281548187425862 0.39835401839797418555 2.1311705974666739216];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.00368391969055075;
y1_step1.xoffset = 13.5;

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
