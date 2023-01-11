function [Y,Xf,Af] = ESPER_phosphate_2_Other_4(X,~,~)
%ESPER_PHOSPHATE_2_OTHER_4 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:27.
% 
% [Y] = ESPER_phosphate_2_Other_4(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999998793562952;-77.961;0;29.715;-2.18999082656054;-0.9;-2.02];
x1_step1.gain = [1.00003298043486;1.00000060321889;0.0138606862425759;0.000266418009857466;0.193068829037552;0.0582073296585987;0.0415843642790311;0.00819772922900357];
x1_step1.ymin = -1;

% Layer 1
b1 = [0.93501244223617852036;7.8469162949353696845;-0.75780944057109955736;1.0989162967124566084;-0.96746330912379463474;-0.71103035559927763032;-0.74435499558628326344;1.609882322792993703;1.406225297678767916;0.87289371315905406323;-0.21997123953913788608;-1.0276160386666839397;0.52711914971084661286;0.29850008462216087723;7.5035488322851016108;-0.62868310227114609923;-0.027973224957888748882;0.46751894663363341431;-0.72345134142508227271;1.1338554179980264713;33.482286432546516153;-3.0423850939296341878;-1.6499550134045071204;-2.1252913922348284359;-0.53704463208884978265;9.6544666711824689287;-4.4957394414977756014;-2.243867584845604668;-3.352047398618550833;-1.1562413030872653597];
IW1_1 = [-0.6172143383269395045 -0.23283148886562263224 0.061472585680133071118 0.14168451030898981946 0.14706422435577617969 0.24820604076374755409 -0.42943385478275669653 -0.67165193051035954053;-0.68948397141615935801 -0.90639040975489026941 -0.76852711402621010084 5.6674370555551671913 -1.2935891876665834843 -1.2460312900741672415 2.3820478602387122891 -0.65090754716980836481;0.13183173426813207807 0.56095834722463755462 -0.27902166187491500571 0.075315232375746998805 -0.1686384625168030893 -0.43689982685080003266 -0.38495004766096352355 -0.70137983665633141239;-0.12616749451776890223 -0.23716933464812031684 -0.13205898837023363757 -0.00051400499716360644083 0.12455500555009557662 0.071707934453724794799 -0.46913017531411860261 0.39499855104083142487;0.31698838420107322555 0.29146624403784776591 0.40213436691148346513 -0.11772793541483195234 2.2816879621312122772 -1.1307308768526456166 0.37636045884901753178 -0.17570886589443510117;0.20004380838581192759 -0.1834872778408528915 -0.54043847096930330842 -0.1295030754480406554 0.057775467114106374145 1.2138022499593144055 1.2527890697389310048 -1.0017387033444993261;0.2605241192878905987 -0.0099278563513271762125 0.0061936338353372330165 0.12993419038774320362 2.5923862925436353954 -1.7666565293981772999 0.17982987230982383497 0.29792884322109713935;-0.47281370165867026012 0.38875026760990455399 3.1677445098850776084 -0.12693247126594228713 -1.1470583856088112373 0.34809078282698679985 -0.38483837337111886923 0.51525716120535269837;-2.1796996411913114322 2.4557111161912823505 -2.96099887301880349 0.49567330395046949043 1.9711060739062660208 2.2189228183538789452 1.403184747554216294 0.59998793941448802514;-0.3352572087687075153 -0.3616410474933556829 -0.5324815307999120062 0.17412530873618145555 -0.60481212776860948654 0.79955000827550803955 -0.24795316826321509018 -0.54807905466228157287;-0.65735060811965906336 -0.84954411386048411803 0.37481993875118557735 0.069602592613308048031 -0.15928140453505210128 -0.14787556735184859136 -0.1999453898122448281 1.3881332664099110819;-0.0020721151864692881686 -0.18113552008521993653 1.506916006549357423 -0.33083824841296260999 0.29491720050821329524 -0.58199481614566483945 0.57382623718299585747 -0.68651371883974765531;0.049774757067319964909 -0.058331287678609046454 -1.4280074516504179805 0.33445515593911756547 1.2169455864688507596 -0.29936833772137733112 0.057267163711166552131 -0.36568912211358295394;-0.13579336702927047309 0.30691098741325317789 0.85578716102348184602 -0.018059743694237340855 -0.93225930674410684151 0.22444550820434996941 0.29552270076482795025 0.53765554929898862202;-0.49452885148203851173 -2.0895809490934809105 -1.2287051425920245684 1.7088158073951817961 0.61566891926527012036 0.28849299608548156781 0.23661844599570314163 2.706869339229606819;0.34041280884134994533 0.13834968587052098865 -1.4217856139232933899 0.3373612246343258847 1.3325684618701763018 1.4188569694949630939 0.35088293409348492169 -0.68876837589258688155;0.12534337244997018068 0.12723548559071820585 -1.1895542554335958307 0.28956286566121386405 0.8560716081811016176 -0.77044626609394195604 0.092716952883825101694 -0.50907306385270856275;-0.21188596669328302324 -0.30321800643867879987 0.51973017564222934173 -0.074128872361880426833 -0.3510108794730712467 0.59007280164093456332 0.28243920581953124715 -0.012142358471426450625;-0.89174500345110829791 0.96032908722004772795 0.3952988926100812761 -0.61780729494593245033 -2.4084702328212572731 0.62937918124144287546 -0.83136007698468583627 -0.15808336548935100763;1.3868702141096631575 0.038417057034017064454 0.7396944168159181876 -0.65761048897564611249 4.1530527686148310806 2.0141788087234249183 0.86643993737119007204 -0.9190782620799989866;-0.78043096332408434801 -0.72538542450154130936 6.7146040023138482411 29.220948357032415998 2.2055829294207627278 4.2325239091865736185 3.3544013976979165648 -4.3475744646469296484;-0.2542394083525033599 -0.15346710065174043458 -1.2516540476923720959 0.35589374960159875227 -0.0095377467071379919833 -0.11794422663749690172 -0.24322919498619754997 -1.5535333692479440515;-0.013667904834628397462 -0.1979940332559522731 -0.48378656291427768465 -0.37508888129785322718 0.2957619734666069311 0.1915482317072868812 -1.7377130312977751014 1.4310252359118205145;0.1431271089293042198 -0.27020727383456599746 -0.74744822448740988463 -0.075484361774147787361 1.2610946760719914383 -0.76157290417302270225 1.149815362622718995 -1.3388820007832431536;0.21309468381662999281 0.32319676643629763335 -0.20745485307306768541 0.01724506729926914006 -0.61392294919491108285 -0.18044200199186591349 -0.35466454970589977869 -0.90794953676400635967;1.4310975858397876337 -2.1225465157068881261 -13.382366862739209168 -1.893745596287127908 -9.0614735372873855823 2.0176196528067347202 0.44753251646012642828 6.2667291979863808393;-0.38084339718409210285 -0.47052711251560436212 -0.53065268167621604256 -0.27519060256260463593 2.5966740054200427856 0.1465846437430942395 0.71811614640666121545 -4.2876643553166964296;-4.0678280798801944584 -5.8597256200525134062 3.4156121065721616326 -0.42727152361261078184 -1.0362797856532919294 2.9103601964174408145 1.1295013717656920971 0.33826653721242055939;0.12067856965942382053 -0.98455343303913733788 0.9595145300294484203 -2.1421966764143283157 0.79276159916897337432 -0.38177281026785914175 -0.027474531371282546327 -0.21176452086645772277;-0.26246427829492691153 -0.18689353292190258315 1.6505017334239155957 -0.23582421855997925952 -2.5768624992657356287 -0.64065927707221237419 0.17364780025593795632 0.47977223711385141858];

% Layer 2
b2 = [-2.252488076511170334;1.0651207931614934399;-0.0076387298811781724484;2.5912074607702195728;4.5876218147948852888;3.3692231575277888922;1.985411502470512124;4.1605077908278005339;-12.736938914054151795;-16.098744253099763313];
LW2_1 = [3.5091808350062287936 1.563330296367564376 4.5592362322758930304 -7.6498523378586051891 -9.5325414572644611155 -1.4520577825484279888 0.66741980440463100077 1.9774914676685475534 -0.79705121469668371859 -6.522298218421704874 0.14323178503891068059 0.36461782457355940146 -1.5636511514490929464 -0.95486038528414352022 2.5854579107604012123 1.2395841851097095798 2.3721066349715189681 -2.2112856938548874908 0.98593859845822029442 1.049549804836119149 -0.49840342383144159433 -9.1882396471693539297 3.644847847067496982 1.3937832241201431049 -4.3823740001480206985 0.22850507193619470914 1.5407553889302625372 0.58982418644711565037 0.5846442003899032791 -1.8176415389173272708;0.59436688668847936778 -0.22315722422765976729 0.94649598830657166726 -1.5377670423419014245 -0.24133955313395233366 1.2159656841876311795 0.59989882435117958259 0.30659737478592985482 0.015762346644570702364 0.042075574398113921304 0.019257717353882791711 -0.44139700603186576844 0.62156962932254700149 -1.3257055618273823061 -0.21236039110492174742 -0.80018706756144875047 -1.316289355895879476 -0.42099725138193849761 0.25530826699095687804 0.11864620069050381412 -0.20943799043611321942 0.53601366501499603245 -0.44122357090948782776 -0.92435065790085169102 -0.67738180413133286528 -0.17404396142476907805 -0.18873449828573118925 0.12519078276448800113 -0.31800749654103444009 0.5957267107014687646;-1.4316667498098951672 0.31067511336891723284 4.8379655755749713109 -0.63307522966711515799 -2.7462653601649886781 1.3555646310646189967 0.68521141393986895185 4.3353543318283849217 0.065812974895356993188 -1.5717103043935141127 -0.25151798368471645073 1.8202633572898705427 1.7587418299644062003 1.0352447058840013039 -1.0683792190825733837 -2.210214702078781901 -0.53876698217430396554 3.6110350226078011815 0.088998824891537126058 -0.62924776633100776291 1.3071291493993137944 2.5935692319160730257 5.0264883682315844027 -2.2806010968133048067 -2.98265440633935075 0.30154459479829820445 -0.87944636068232795445 -0.040526690625064816509 2.6578108078225377575 -0.95918051816123217623;-5.2261803533923618659 -0.2463986109403678848 10.17449730236093508 -1.5533543004298402401 -0.9272696531292211608 4.4094513632655996105 7.4158809425486804656 1.3437323319330736204 0.19029838848437966026 5.8976052360327342328 -1.8375649612592852566 -2.1150393003139322801 -9.1543249001434325862 -9.6890015568775833543 4.6404985630329296242 -3.1558413172989889794 6.4017226726281108995 10.515565190811823726 1.4227151892911795716 -1.3947699440892400435 -0.11004545447675541381 -1.1038430351634391791 3.2005001565803548047 -4.0243121953482754449 -7.4143037373675895196 -0.22066887894565229145 1.9144413414382375205 -1.2854435510589437008 -1.3855035761009086848 3.18638484888429252;-1.2370436078497371035 -0.30167385361942933297 -0.84175363508680578484 -3.3941238388904451106 1.8932187855483626482 -1.9922575057568951973 -1.6545107289166647302 -0.44608123420342410403 -0.12285154548561789956 2.1281711502395412872 0.4892719493895382743 -2.1942233280446039601 -1.8112037568675338761 1.8175409172753709264 -4.6147807539449727798 0.13808351331502508041 -0.12319390541950721307 2.6282353872719279586 0.23006113615861162325 0.028547902034521437792 -0.019555473711573533036 -2.1983124079313833832 -0.41068716603712701607 1.4180623935388392276 1.505806402764327423 -0.10344309864330261672 0.20373906056816787347 0.11249248554440903491 -0.71707702719142740744 -1.1581754860372706073;0.30940512238364048514 -0.15850554721552567727 -1.9463643025959940225 -5.0093367558001977713 -1.5874140200047213156 0.5744947981635386558 0.28481448184391289224 0.55311645278266996506 -0.079154345330128883518 -0.075406527336195364186 0.24250853995897461712 0.010220682000204320966 -1.0724542594605648294 -1.2031956171159803315 -0.48876282326966524749 -0.015204076158276163405 -0.42789853802067012722 -2.5103934062338075783 0.47210383915623621576 0.10599779184278429467 -0.19816075616489772782 1.1369316952829424672 -0.45058389302065010584 -0.61117719870917353031 1.08608293905372566 -0.14878720540469422828 0.0021894926128897193585 0.013696473356287161907 -0.15546908887143540134 -0.79588012016935372195;-0.084523709217295947549 -0.012827734567123045065 -0.62258715211406767676 -0.95926470486042714203 -0.059727757326946240291 -0.16866510384908178777 -0.13290684139481304116 -0.059917936618461536957 -0.011606952586268938152 0.011157544287241770442 0.11780371971122248487 0.090665685844184246145 0.16202267404706813259 0.29821988248080721906 0.061584473631692393525 0.099605718296598907968 -0.188693477749138111 -0.52210549752120183342 0.05973149566566023777 -0.0030809549295838156366 0.044172918683743150892 0.23696151822778477314 -0.19238788257903513546 0.17092127584074354174 0.31335398884576515233 0.024013401463786283874 0.034808326739773687297 -0.0068907482410847627555 0.096685415984723493965 -0.18210904803716482236;-0.34358245531780762239 -0.32771154824756987267 -1.0047110206064140403 -3.378758299426017242 1.6163087281630172942 -1.0120441945160036923 -1.0514019355956452095 -0.28773024979360622977 -0.048410538882150980222 1.3530836069238214581 0.3685574348260940436 -1.6119436406781240123 -0.51316478136784815511 0.69187697081688726009 -3.335781998399663717 -0.1845943942647587932 -0.85345907051101899121 1.0736644557768282304 0.31826730493085736429 0.065411132229172608277 -0.059198526430689769884 -0.84480175994031581066 -0.63118610054694190836 0.53778429256223614363 1.3452359943670895159 -0.11726813848436930254 0.10796833878046796507 0.13468054666578910261 -0.53330120859929996957 -0.34046199002318472404;-2.0393930349327740181 -0.10275119491559933649 2.9375973400659960433 2.8124781348698366656 1.7121952658335588549 -0.074981876837465155639 -2.6349358446590929184 0.23004754446669481438 0.023641712049951885133 -0.62787681527667860149 0.19669511996616606542 0.16852520689308989055 0.63356723754514654434 -0.27248613615446554137 -0.8251102053370690248 -0.42500708141051574129 1.9502370745820769571 6.9635672035027793569 -0.26141803232435195659 -0.058851632115763378994 5.8614466463557715059 -3.8925671012931859494 1.252998827372855617 -0.66955205140335405378 -1.4032418852705774714 0.038710659176380661506 -0.2211652128883226931 0.053657201806759480256 0.69125686533946184564 0.26123083593159041049;-1.973947764600539001 -0.052540405556991957858 5.3970116467855380193 -0.68867945569956867669 -0.61517280392031892244 0.88262432996678874542 -2.2391703463886427627 0.13570015866081711065 0.055513698308284939453 -3.2591304480656546616 1.7132734989542550696 -0.72990202337180909442 -0.060678301904606685158 -1.8760480568196171713 -0.7655254300938013623 -1.0246172656035765769 4.8596112017718215981 11.842409725852339974 -0.5182567038479108934 0.07745096569975186418 -0.73459890842581199166 -15.513643716302636832 2.7574928633190194383 -2.3908429004115099303 -0.82352868557917890957 -0.23259398461381605716 0.26596032694740995295 0.026781834253160025083 0.018024955849045620759 0.4166759301903604773];

% Layer 3
b3 = 1.5757189203095456342;
LW3_2 = [-0.064364078708738600709 1.3154588465490464344 0.31228994314430263968 -0.46794637687592621766 1.129834634796575088 -0.70495929256107225935 9.6463089409937747831 -1.7341006309086326453 -2.680924621872174729 11.785210796854237003];

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
