function [Y,Xf,Af] = ESPER_phosphate_16_Atl_3(X,~,~)
%ESPER_PHOSPHATE_16_ATL_3 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:31.
% 
% [Y] = ESPER_phosphate_16_Atl_3(X,~,~) takes these arguments:
% 
%   X = 1xTS cell, 1 inputs over TS timesteps
%   Each X{1,ts} = 5xQ matrix, input #1 at timestep ts.
% 
% and returns:
%   Y = 1xTS cell of 1 outputs over TS timesteps.
%   Each Y{1,ts} = 1xQ matrix, output #1 at timestep ts.
% 
% where Q is number of samples (or series) and TS is the number of timesteps.

%#ok<*RPMT0>

% ===== NEURAL NETWORK CONSTANTS =====

% Input 1
x1_step1.xoffset = [-0.999999996490807;-0.999999994516886;-40;0;9.1619];
x1_step1.gain = [1.0000000017546;1.00003626070655;0.0153846153846154;0.000309885342423303;0.0656325902039552];
x1_step1.ymin = -1;

% Layer 1
b1 = [1.2909081080119384488;-3.3727644935199325715;-1.196474917875081756;-1.9021173406571063591;7.5845026075754242711;9.8060717195015065784;-4.0211928489549917387;-3.473836893254724778;0.64701786254820403688;1.8988883045873135025;2.3255598759075577497;-5.4853957878719112884;-4.6618363121201147337;0.89754439445701172939;-2.2699286490834675867;-0.22619341102942358557;0.27865178630994835718;2.2481812905535805136;7.4910955385579658383;-2.9776780874303288016;11.704745917504499886;5.2651204548516101767;2.0561364425169328385;3.227901349660391439;-8.2060518695116773102];
IW1_1 = [-2.2426178620812335929 3.7604012757185021343 -0.84685981412561039949 -1.1956866377007353552 3.7543046615839150704;-0.59935793800549808097 -0.72153760390265686553 -9.3456115306714373503 1.6790825339218755641 -1.5025226287731172636;0.94226501947196006359 -2.5115188699232500191 0.99101256337717247291 1.9893013469827187656 -1.4188502686087394178;0.52908986957714032595 0.71496585231770426461 -1.0732392988845780302 -0.47601042062585735559 2.0752285463419490874;-1.382248814613261878 -1.3800411694288892139 -3.2652791184508394018 -0.5333420061814809765 -8.9452543787607847747;-0.0066829971936347518785 -3.5546127636473827671 -2.8631353638700223208 2.9518494412913907077 -5.5626841233722332802;-0.2055554587008586731 -1.1445647049846743304 -0.28931978766451105223 -4.3459452823210131811 -1.209059146549303021;1.4363783294843464056 0.53728076155325921182 4.5200806110821698525 -0.48147263935912065769 1.5287443674957956752;0.2285242044806620465 0.958945165396426491 0.82063239160833478802 5.5688378605689949907 6.6328633864473536264;-1.2793890251350459319 1.9690026263032835985 -6.4782733005461770048 0.12885201678069152509 0.5561867178009920476;-0.049689967642340814036 1.0929013956759932036 3.3369959679694427557 -2.2090312296469303277 -3.495824809295575708;0.0080879196810948340152 -2.1804445308686966065 -3.5611886303989930269 2.0154467541283076493 6.7228458650671072405;-0.43904045551652981505 -0.23199748919848087181 0.61168209882130519794 0.29885300327695779599 6.8666743849024820534;-1.1045852447798620855 -1.4628223403940407987 -2.483432955598592784 0.0486749285160754247 0.74624507181927679511;0.10063769555283738621 -1.4897224478788531066 0.5753239313141632616 1.3479431524387308627 0.86717329242635987629;3.0005468075223249436 1.4590655457051513455 -2.1685003288050150871 0.51109245559711380391 3.5707570253810598793;0.030580314373936302863 -0.13122904314997588249 -0.004728579015078900917 0.52381092493353620387 0.13603209691872658604;-0.75056526119305755351 -0.64552338690830013412 2.8342233955101692189 0.85205953787134458555 -4.4568514459655457216;2.2405452138332653256 1.4290111937882186144 4.8151866951075508183 9.8873714886795784906 2.2667167597371737919;0.54796011783707698761 1.0192475180730460327 -0.90285813267608783228 0.7550349642874745193 7.5824097333247877373;0.37484168453660188236 0.44668844438642046901 0.95756218632233669474 6.795921887420501406 -9.2290434066956894554;0.12318241760212567115 0.50345933935156439443 0.088688494673390197431 6.9689273549135952734 2.8690484898327874852;7.6637627776036021388 1.4671364217079494452 2.5673124374639439971 1.6802413489904528543 4.8981951085313166772;-0.098336345979993994293 -0.25612239874366082848 0.55529403542462540955 3.0396761339208380193 0.5709602063973570063;-0.14121835532469087515 0.003903712686454781193 -1.9242892622567719929 -8.4873558578422585441 -0.65440964719169136821];

% Layer 2
b2 = [-5.1785978395160947585;-2.991429479099493971;-1.5592227460411756113;4.1803968711949481829;-4.0534710717432815841;1.5350727273182429578;0.53966100137896377831;-3.0497582430356939831;1.7291970818961124312;-0.41966368717699020774;1.3738462513317701852;5.8722797905395296425;-3.8408854564641274898;2.3900655891813666898;1.9524339972659436881];
LW2_1 = [0.08445644296104061477 -0.21729946775519318192 -0.75178307784428355376 0.89184655693217029793 1.2025532220493380109 0.40587288536174792597 0.77830039756544344698 2.8127939638338435024 1.2585928035059994468 0.37074996953331951222 -5.3368450131882658738 -3.775177022944377736 -0.92848660270217209245 -0.37558663334973868553 0.81009890006512730309 0.49767568259220923421 -7.277132067375256419 -0.36477843349110294913 1.0297536086055860682 -0.89506285543323049758 -1.6796697217680705005 -0.29517753195073709316 -0.52615521228492034833 4.1148429274457454241 1.9544119487066773377;-0.19617514871135996457 0.072074105893845413773 -0.77412972096894050456 -1.0644224798350632177 1.6052238857451588405 -1.4285060935778717273 -0.060398997022448989991 -2.3973753806674880096 0.73260653520789120741 0.57649894009777324921 -0.37146434472600753374 -0.12951143337647291753 -3.5073783750633382716 -2.7034222941770651438 1.4645404894520965833 -0.60610089115744925081 1.0861170848254353416 0.14013793762797627229 0.54504933718034909162 2.7969823094796351093 0.24242594393693259192 0.86505434129251945752 -1.8191936778661983443 0.30329625927543218955 1.9488787763705319733;0.98260912079164897115 1.0124463932021789958 0.12795859143381810297 -1.7447293925540441339 -0.15871590770575982887 -0.13234326404545720957 1.7936630931090660201 0.58815659173719769814 0.93565432730867159616 1.0833169791008649963 -2.3450079086529256323 -1.9003422924004280681 1.2232931792735186605 -0.60080035754083394028 1.8933917718586701362 0.15177634065902911864 -2.2797837302854762065 -2.0185060151516487359 0.71797371004188415 -1.1318387935113707332 -0.31663285674568375505 -0.71872426880179351993 0.1306268447922950271 0.3616572297238846434 -3.812883630662226242;1.1831866403080775196 0.19064674952261745355 0.52703571987568154089 -0.40510340587275617263 -0.49339866609270749409 -0.30893731063378532253 2.7560403130528365878 -0.78002026594058659903 1.6959347113553269182 0.54995964879270919923 1.6119321344659769402 0.38354948866845661737 -1.2935870186317750541 -0.059742034463625318108 0.98778333352068559137 0.45805911278617383875 1.4292812443742395345 1.0158180265444241019 -0.38836088793483730663 -1.3426643800925255334 -0.10319581089524215711 2.4833111517597239448 -0.93985206420970746066 -2.648634103741833723 -0.92911300819917241878;-0.68889218693194009813 -1.5579599995609227925 -2.7774239426084852411 -0.25777398962825143203 1.1213963717860930824 0.75153865274196129942 -1.0137302696169780969 -0.84571857136917305553 1.1237657832664098567 0.62585553729371556919 -0.60684325722688459148 -0.73030809479901159786 -3.7477302905372682496 -2.1655918546081163178 4.529675415247511161 -0.9269700812158896408 -0.37036734450373198779 -0.5713957012229833321 -0.41693266271333900841 0.71685223927030294888 -1.8646469031113945825 -0.68936834422490944618 0.10219871164231195126 2.8532744938345868668 -1.8637696999873847314;-1.4770831605562786404 0.49327246547093184237 -6.7034024532525435802 3.0577516168560601173 0.27016639291468147466 0.14280422422519151171 -1.2521183452673725345 1.9838179819636303503 -0.54561470770619480142 0.24607871685976784315 1.3381058551018780722 -0.19161452930436398412 -1.5887969844447247247 -0.19808168150394966922 8.9972678390722116148 -0.74780704328722780172 0.066473838959467582144 0.65697443931640142623 -1.1240835061739349321 1.2887252716686985465 -1.1585109395100645635 -3.03116741124014899 -0.72161602380845502047 1.9221855180927949203 -0.1915662734711959525;0.66544418513840219997 1.2306914957648136255 -0.22074832873464794725 -1.9646164192671655702 -1.0794873441519849955 -0.64000680685991229169 2.7180135055187735205 -0.60750087365562432939 1.2296986321368146289 0.76803513185947502873 -2.4090640512014633856 -2.020971798965356836 0.80647694476271758468 -0.94524670549899392213 1.5833781706767926512 0.31059005603348044389 0.11293182889540419322 -2.1761697742936001454 0.82398583016553139657 -2.8320104104409589141 -0.47414401179462323244 -0.44205558284462109775 0.24548748895653144109 -2.0862587640996408389 -5.1810631928152535153;-2.8074735457480812784 -0.18627179650108557918 -2.6392012897873726729 3.0999280659003276561 0.56858829802596977387 -1.065142935426617754 -3.1289466159998884365 -0.12314419426918865597 -4.4032367814797712668 -0.83611835149794255795 -2.068582888918578444 0.38046825192821592809 4.0564307581938168212 -0.9926580006634979636 3.2204905369413907579 -1.3459456727345702198 -3.1480546880929094478 -5.9694936394405875291 1.3960369179417233543 -1.9197404426830113611 1.3633911756084897693 -4.9126704609391511624 1.5240878247581934879 7.2225828014006729916 -0.40658535997888906044;-0.68683071711407284177 -0.57212039477761456041 0.15282371935928332207 0.88677504824443720288 -0.013891321273785695381 0.037267876157979169804 -1.8343960205688571996 -0.47654850676788956809 -0.87574529259036282269 -0.8361836539832920856 -3.2081322950274713257 -1.5715046155179777099 0.024918986618021467944 -1.264695248505085523 -1.8243476927355710693 0.071102815207445629886 -0.24173131480523121772 -0.72835516034984992828 0.25094119170389700413 -1.8804262577080341767 -0.34099726500820365382 -1.9183318427952085194 1.3561812386688198995 -2.0291589555863565053 0.35434873595616356301;-1.5043008534665847087 0.0093938388957020262787 1.1733648027851528717 -0.60921768329692294053 -4.7144148019064919808 -0.71220778042691590404 -2.8014789394958383717 -0.70844041114322608976 -5.3402463878461503555 -1.2512255282091462316 1.2703279251776091741 -0.91138345926751596071 -1.8227153364749990594 4.9944064827304730514 0.88517488274741518151 -3.6389694498957738666 -0.56683929838370106236 -3.283383821782945855 1.5676767168626613014 1.7645961218477308474 0.40710092565093020101 3.9068809220514553182 -2.3132571412141431644 -0.38707949922356876815 3.9330858913123178944;1.548173074078877276 -0.0078191022026017509977 1.0449110440685527301 -1.8338551872231034778 -2.6765389851614602712 -0.11921381563788387148 -0.46941313925132643892 0.079566687272643876749 -1.1907226775382320216 0.70062862896642530419 -1.4560092257390535586 -1.4973388262392568482 -2.1084807093744495354 3.2940739415122517464 -0.52477691032060991283 -0.26551313352004474178 3.0419521712871251395 -0.90077276457626753725 0.37103742746844459166 1.4517057361267271531 0.74420887329342322314 0.065698030369069571055 -1.9800215104452785209 -5.881894284660657668 -1.3242549927182063474;0.54373895998985188882 0.061162580040656173985 1.0649974045532144551 -0.32359837166935656594 0.65849893364067246804 -2.5732200193291108192 0.5040879948553642631 -0.59802689420761956374 -0.53567664310588969379 -0.18298408642475014196 0.297074044495876044 0.2964067612509574956 0.6725799156668574863 -1.8277659351597101978 -0.96162839125374111049 0.29364919210805362493 -2.703336855199205413 0.93777417467877888235 -0.45883339744504159086 -0.92003989072065772703 -0.41874740170786822757 -2.2716225144355695065 1.3286702361219662194 1.1772113547543168188 -1.6458105276370560421;-1.966487656893472602 -0.24882978124259996777 0.79341627401230951921 -2.8918443580661494785 -0.021073006788363518593 -1.5077891163142995712 3.2977446161752776987 0.18876576140868472908 2.3921250216055809368 -1.821066158791333045 -0.59787474691290876105 0.85060940048196875818 -2.3633145682793004205 0.99758128691447733871 -1.6797384075465673359 -3.8378284826665414187 6.0385142095790893535 -1.2438774853872509585 0.1447573509027757499 -0.53858704458726380793 -1.8626836792076142402 2.6284307572779446183 -1.3274670806663195766 1.0911336475808210977 1.8534176333709533324;2.1675257882284268973 0.31989376072468428402 2.1268861237475036319 -0.46394417297176931969 2.38017372420131057 -0.92742978187999369766 -0.036619816620251675887 0.92317534079102847588 2.0987599695577414671 -1.7949938673113345899 4.878901784363133487 3.796018414363904192 -0.047424622960999800181 0.26620265896111000581 1.1061517181279529876 0.73267713270128687064 -1.6542901530894726925 -1.8396233634192027395 0.1221784100405844814 -1.1523336024710484349 -1.304530945514104312 -1.5129403899905859987 -0.96586056869242664646 -0.52051805803523454319 3.8623231176734651449;2.0789734844149245774 0.31346760798136474868 2.0587976304956496776 -0.55836818868837090335 2.3759046636103144934 -0.93810392117728813588 0.19646299256399130106 0.89945864248465690061 2.2444291574281569623 -1.7604934344497626952 4.5871356270066163674 3.6365727942278196494 -0.07680956415549695937 0.16202322277690309194 -0.1484308464182182985 0.68767673221123970517 -0.12242379519522894604 -1.7586539603995023917 0.12566115005355735224 -1.1574859540204833763 -1.3443424788292988037 -1.4101298241731419214 -0.94999168179343074225 -1.545124279748930185 3.4814651082298713192];

% Layer 3
b3 = 3.0517245939391339782;
LW3_2 = [3.7231133530842162038 1.6425571224378201496 -0.75483207309213384839 -0.51429039466794679303 0.22265152574400573493 -0.82194278455456926125 0.65312764539727685698 -0.17432330499847925598 -0.29835084863706901048 -0.23562116780645483982 0.36328253315091446618 1.4814801602142118053 0.3858363143308918275 6.2575292528987676022 -6.5593033732195111796];

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
