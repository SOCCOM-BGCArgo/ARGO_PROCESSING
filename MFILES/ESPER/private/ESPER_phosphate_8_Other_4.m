function [Y,Xf,Af] = ESPER_phosphate_8_Other_4(X,~,~)
%ESPER_PHOSPHATE_8_OTHER_4 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:29.
% 
% [Y] = ESPER_phosphate_8_Other_4(X,~,~) takes these arguments:
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
b1 = [-6.3736981889533312895;5.1146976768432095284;-2.551530923964120845;1.7326897582369304995;2.3117601130453993008;-2.9380221420725125547;0.24218600430504833021;0.14353624564359068749;3.5378298281291344907;-1.8191544533768699843;-0.63158814204389102098;-2.7684966518734106522;0.155900571961438833;-1.6147263766958435127;-0.25246352864389093318;0.90429393314347872312;0.40871662061233560159;-0.83524980562198714029;-0.86541842441123562057;0.29392050048522244676;-2.4361477433765066358;1.3363795896972241461;1.6093404101184425325;1.0143427565878990748;0.64334431847240225455;-1.1052029316483069632;-0.88675461462807614321;-2.6671762512858920857;-2.1914263134942468803;4.3079057014135937465];
IW1_1 = [-0.35760908821125247936 -0.5816851364558608406 -0.33177577431706595235 -6.3819013377684878918 -0.28010054840439968515 1.9856093084538490867;-1.0348741964358267875 0.5007782062349928065 -3.8929187999626462258 -0.25744305528725042187 0.82391261795690928782 5.2504095951876674775;0.045297754725794073249 0.026210000761739078723 -0.94793002630295086774 -0.22229809865836619442 4.7990156111163271291 -3.7615787978007251269;0.34440998704561825994 -0.15589494890163968499 1.2843226443945729986 0.3575979673900994138 -2.640874192647887142 2.6814439038463495635;-0.034618856982625352547 -0.072611357816115565389 0.6692040449784012468 -0.31593198420104556812 2.4246686785163298339 2.6822584795104011945;-0.033864769368313923925 1.6770817409594414027 0.85662840623340752444 -1.9414301892908643055 -0.18729777537012912147 1.9801405688279325457;-0.07776110540984607078 -0.26007423818361091739 2.0313210528602034799 0.050573077957523221437 -0.91853697103487352837 -0.53750815058752077213;-0.17741936489126830523 -0.12803672619359707241 1.6881146314758832538 0.18424426389453726904 1.42082218795874482 -1.9103962604334130404;-4.1817707048637808853 -3.2700044523194136126 5.1377444033280612601 0.97395054985503703993 1.8762047863707385087 0.39684679449588722289;0.077359909477488136842 -0.12541158826503270851 -0.94803529156265042843 0.96874275334683246363 3.5034372079100957187 -1.6555421683632889351;0.51202677244861893868 0.10961394709923631685 -1.6536178761832946993 0.99578906323766769315 1.5201417796209120148 0.27615354725912333311;0.68691753104893071846 0.82654164186511092094 0.36413522832492578596 -4.11312880987396845 -0.5852880974530486613 2.6239435242227373557;0.90025341057452223925 0.23979919622915352173 0.90030950369985507376 -0.095918818233641942039 -3.0782092971451051611 2.4348404666538483454;-0.76565014585581769868 -0.54507035875548670578 1.4309547003038567414 0.64798386295616372532 -8.4368932415698427718 1.4914324030213810257;2.1478346964752565995 -0.86617422951115685681 -0.53780399419626956359 0.51397691930599831434 5.2346935704479742313 -0.16288502293623893791;-0.49582866701378625063 -0.083776295231884595482 -2.0817212157729487743 0.44092513911703951068 -3.0248140505737888262 0.22473927048870073753;-0.93171527555601174608 0.28559887064901601228 -0.3314365306085562124 1.5536687739660854923 2.1505045437344834092 0.59617586515543030679;0.61402375426565625194 0.4754214382426156793 -1.2829383224578481126 -2.7178020314837585403 -1.2623116079114082844 2.2473700626392920654;-0.61199759744860737776 -0.13349211097823407424 1.1340435320385764406 1.670865665021848967 -2.8837488550768530793 -1.5239979603166804445;-0.60824636645220653275 0.23733974681292260844 2.1052850286703819194 1.4194732177783389027 -0.81685050644773959938 1.7004195062940046235;-0.30415872350282907588 -0.51111917965542508213 1.9069307476614381081 -0.71591676088567768321 -2.8364447703723651983 -1.6011089693530307088;0.0028430251226957653325 0.14127915713683616983 0.53309188721352085416 -1.3146283373879086298 -1.2725271043832104834 3.5289017888762330699;1.5405508617963865614 0.9052739709305225313 0.093399774699014767587 0.34272071077904298075 2.4842563574593681253 -1.9823678571729823616;1.7955538956597880684 0.71951848190456668686 1.950623525218840193 0.028862993404787807422 1.5050768549005313357 0.17961490880269892489;-0.14348534896750947398 -0.36112720562302247362 -1.885565469088977153 0.28507290374635085861 -4.3344072441455399058 0.69905754995598157198;0.25312719084666668312 -0.16266684196767666748 2.0275223584878334648 -0.96107477766109006279 -1.5428236898541416444 -0.59475397564193666078;0.0035981794536900446949 -0.084533130210935794602 1.98066370726921126 -1.2628167862105212915 -0.081556040515856298434 -0.75800730988587949888;0.083766952476088826773 0.16733172177499669453 -0.33065642689542790178 -0.058585058119323497883 -5.2690434759028237366 -2.3624836328061022961;-0.50672011045208376068 -0.52449746591097279325 -1.6153441936093007048 0.34229364161707365177 -0.16168354207163299074 -0.25412489887349209994;-0.44140566815916049137 2.5506004074735568743 -1.5960202559160623359 2.0375320117270230291 3.1555606883341198454 -1.3818196381818286245];

% Layer 2
b2 = [-0.92368593041301427782;-0.82804053404024668961;0.33472321517099157262;-0.29418417864014850815;1.7921238575926237147;-0.82034293964863214121;-0.9425880732168421261;-0.066306206906672465751;2.2808533629655229902;4.5844137771294146688];
LW2_1 = [0.31539951930731086049 0.23672754174038299069 -1.1884734659455273853 0.1471760426353965201 0.030598881414337163187 -0.46059427567050537089 -1.7508837195559496802 0.030536188683587706105 -0.13345479552667663392 0.46065823614244677531 -0.047391657026975037015 -0.80348718593813184263 0.47300897021779864104 -0.0071024787559501858722 -0.00031891914159173603057 1.2696476508140279726 -0.35722222514666041038 0.23483512097562037058 0.49617746386181110774 -0.32203973040110106529 -0.19428557233746854283 -0.94493779020849100192 0.31016166670028449959 -0.11593155298267804332 -1.2584339124022192458 0.3883793400208573332 0.75946437644833064695 -0.54294625103770188801 -0.31139486741249966517 -0.019006057552275999051;-0.62651407277905446502 -0.46067661290698941334 0.9958462902554149343 0.45099709344603994277 -0.55571647135981294685 -0.54187951596664685194 -2.4658981255470986405 2.4165579284838574203 -0.048365518749290022094 -0.69077754656139223055 1.2967726827386989452 0.19448700309575195355 -0.56332089333000534559 -1.1823557690454988833 -0.15278844999562871121 2.1507362053530609813 -0.79912175288061360057 0.095254364352681433759 1.5819343990600651662 0.18563413458673583278 -0.88384115039170974004 0.6227864291392706253 -0.40554592917394349794 -0.26495755995329955956 -2.047409266517536075 2.7062918523955432626 -1.3830992320445634203 -1.1133244723650885799 0.34327026938192373162 -0.85223421229837759139;-0.16294654218520759903 0.34356490921908888758 -1.0202384846416074904 0.095555038130014424702 0.85861475859105429453 -0.18758777142848873098 -0.39891574287230519014 0.58458719631439137743 -0.10815260529030426562 -0.50187297231174199386 1.5527104965521558722 -0.54834296776274205065 0.36955552658033624658 0.093688423039573490159 -0.23934293590439401078 0.23227290149730356328 -0.12642139688970976952 0.60074113340884693901 0.20815795501482350915 -0.30743584458068917709 0.091957659788917747723 -0.99673510883562144169 0.2403612413999875741 -0.21858173406723419019 -0.33305419945998115949 -0.14833555869544476358 0.5779378798168034459 -0.1554180548744187873 0.27487187784326222451 -0.12578221514961490657;-0.6130807964595286208 -0.50593806418141573644 -0.64355065480826589841 -0.7166517031049917108 -1.6441022045639848148 -0.32230243015497583858 -0.76192052143731603397 0.57355650412056402576 0.56181055312023064907 1.3409156886775803752 -1.1388640372407539125 1.3025885125950764643 0.4519294633001973982 0.06202593557423644044 -0.12372373486759220829 0.35082426859350740456 -0.67819778272812936137 -0.43211453440478880816 0.017023530320957179496 1.0917617533910761729 0.038296487758630225195 1.0983511495628348609 -0.063717439768687753054 -0.26427407933931612583 -0.54399347681916554453 -0.75770268223134351615 1.5047296898870932136 0.014864494057676855276 0.85172930694109461047 3.4404828284348475467e-05;0.18370737201780362735 -1.002250643521505058 -2.6306022870865013452 -1.9060064809531773733 -5.3111417222417829365 -0.92152989926124107001 3.3403056575084146829 -0.7554987552865397582 0.030713660659074909498 1.3581253947112643665 -0.79832356857752861945 -0.77594489665947041779 1.7054474540829509532 0.14429016221780430196 -0.12907335506389919222 -0.3096001344531623678 0.84409296266742328552 0.56150705599882844155 -1.6545841781972987494 0.45395024096166619598 2.0810364864792125772 0.64078762036160785698 1.1930401179884653207 0.24119893063701675073 0.30406146656926413829 -1.2384555113353790112 0.61494359283713828646 -2.9497653012151965335 -0.87525083794375158064 0.062633229787716099279;1.2885772402365964773 -0.23863530687697975807 1.58411763567014785 0.68415172932949064322 0.44369230030065948123 0.017998962839971478084 0.21213841037760072794 1.7436427025750864672 -0.39714548805669225295 -0.80500259664578988072 -0.079443826363333816531 -1.3160640624522426112 -0.19228873009163394148 0.24001470418354173275 0.4650931129495335381 0.24378522935949120676 0.95093849703319333067 0.17135744577913292197 -2.4443099438622684083 -0.18873460067564662435 0.57321396280553138514 -0.61131193628514945981 1.4749400748372083481 0.1824516453447930231 0.9278662025220106857 0.18224922354635086919 -1.0065997192586961795 0.9620867162410503548 -2.1117701524218350784 0.18231017115968034537;0.5506249795556407145 -0.13252521871411915289 -1.0212425825095623466 -0.12026896725176189351 -0.25778146682775576748 -0.55780670005685106805 -1.832655192698721347 0.37879403995136873284 -0.1335745909645557894 0.46253957716218224716 -0.89681161960456823312 -0.48880072616048098055 0.57985031369719020056 0.052101807839665739697 0.18334823713801706724 1.3749633305775597059 -0.16202198746733148527 0.037781042885133848952 0.19203390293892994012 -0.21434140588591157894 -0.0042304271774273741538 -0.3449706121242496204 0.50927075735828675196 0.032418724474558985316 -1.1535996562189927417 0.56632223466344833174 0.16860449770591551255 -0.098595364841473270845 0.062204572177699340996 0.43560774027298521105;0.60370035385392084581 -1.6855619980354956677 -1.4219729000190450829 -3.4997173772217005094 2.5962586716042679313 -0.086846596614739507425 1.526035415692863717 3.4554795087695895361 -0.24761930979536447839 1.431831683753373774 -0.84566064275991525978 1.5915292523462480911 2.7334501149168404233 -3.1801791880101664134 0.076632021964356325516 -0.5420353506588164727 1.301286840099452613 1.159860222399522911 -2.4241394059896896884 -1.1907048524764667352 3.1535303121101305734 -2.1423922219876043727 0.88322652678592183939 0.25486970281542248307 -0.65954118842646880516 -1.8532556278463512278 -2.221416978236441242 1.8161005771097455952 1.7531385131964034585 0.20555501845719445497;-0.4120608621417318207 -0.52338147642387278236 0.16930907458822597178 -1.4819140861434347833 -1.2817025580797527518 -0.9739394506587776279 -0.92114347230671922429 0.4355446948739392643 0.25012975399069248272 1.6800220800603016258 -0.29133083214374133618 1.5552849493959521876 0.58080065026524274607 0.10621724484413168355 -0.35647569131640821638 -0.58056728136396651418 -0.5009371792408788826 -0.20238114956368413733 0.21710557499102817181 0.88495060089952937155 -1.4814901854365194467 1.343257668487833234 0.047976819330845278666 0.22501514955898696546 0.45392519736379544071 4.1280668591926463762 -1.8426507369517168833 0.16185374161908122614 -1.3329359580976447575 1.5101881386186926104;-0.45388397704499217067 -0.6416323465945894533 0.78589767627249162985 -2.3735071328244274902 0.82430530559433501381 -1.4572218938993595216 3.4997717963147159104 -3.8583117091416352373 0.43224417991021135688 1.9542599995317417871 -0.68441510214983702909 -0.11013379444775867055 0.99111099119093837295 0.16028180489908702699 -0.27508779652475157107 0.28551132349481594952 1.1951923572269718132 0.95933247691055090201 -0.78133335949726367975 0.61035407222286552908 1.8310128583946558845 0.60312524126553157533 -0.60219541645504892635 1.3543543060494998009 -0.85705219509945929079 0.44711240350838327151 -3.3678802455357041978 2.2052335064010462062 -1.038455734510197459 1.5188131346383784148];

% Layer 3
b3 = -0.16136261179697902746;
LW3_2 = [-2.9516724842743897739 0.9773159877134564022 1.7511688385742756324 -1.2545827563352487388 0.33570376826331232989 -0.72473997817594715265 2.3214228731957136098 -0.42810161413213781989 -0.7901008732451691996 -0.63211057810370940846];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.568038194888224;
y1_step1.xoffset = -0.04979;

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
