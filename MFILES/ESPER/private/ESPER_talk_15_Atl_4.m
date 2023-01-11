function [Y,Xf,Af] = ESPER_talk_15_Atl_4(X,~,~)
%ESPER_TALK_15_ATL_4 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 12-Jan-2021 17:03:51.
% 
% [Y] = ESPER_talk_15_Atl_4(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999999390765;-40;0;11.29;-143.241523403244];
x1_step1.gain = [1.0000000017546;1.00031699124092;0.0153846153846154;0.000310077519379845;0.0705602618858232;0.00458552991639557];
x1_step1.ymin = -1;

% Layer 1
b1 = [-8.1436973109601176;15.666592212891626;-17.408917289891043;-5.2923705743055276;12.894298128159996;-2.0365014290074699;-14.986924569107421;7.8682675595882596;6.1117300659056513;11.073715361426538;-2.1004714570637999;-43.031904227769616;-12.92762027762117;15.194511102242959;-12.752216125997181;11.507508912170056;7.4407437195320547;-4.7156861368964584;-65.028186912966888;0.3423280626318001;-1.1678956755684708;18.198605642109431;0.63390747548193249;-4.0581602572801527;5.903062562457051;-3.576262014024052;-4.2604620613054456;0.59808423950993539;3.3677165820330743;7.8675092989667732];
IW1_1 = [9.893044306607603 -6.4812877753242741 0.68106869758683186 -0.33428769893214999 0.18209835713195155 -0.62492662991996628;-40.385452175425641 -0.077735505124667995 -17.393073578159207 -23.895992583187507 52.758360885961011 -37.200789127845823;12.505131517012984 14.197347680750855 -0.3128846491294166 -1.7874764580671803 25.583199550978996 2.231296090938339;3.3281835153348425 0.57179629745525506 1.6094002017269957 2.7825921368476516 8.5996873096327615 -2.7794167731956438;5.9282575651905418 7.2062308656148675 -2.6357553153052358 3.6887711645862651 -10.407558490908958 -9.1555712764447907;11.546327319010281 -10.82124470880278 -2.5022404067352997 -2.7136055354708839 -17.73614673556774 0.76060293316722449;2.1882513047687686 -3.7401436894059317 6.566972040223602 -2.9524805452293621 18.591398243617157 -2.1938699528895613;-9.1487387561677842 -2.6231092704095764 5.8721064391907714 1.4658233570412729 -28.798688469154321 -3.7373933846169693;2.2416844254645323 -2.266746927665408 1.1887971652328899 1.6002702502651183 -8.6281463802874931 2.2105711040315859;7.3237250858370455 -0.22304997938459262 3.9739172534734073 9.2520697555516733 -25.63444717900158 -19.960564812191787;-3.189737973192913 -1.3314370948640704 -3.1152998353196151 -4.9032227666695736 -5.2055392045255262 -2.4198528394605541;9.7593852482301156 -21.118676331751804 -1.7421597383682885 7.2075225336081781 29.004963347579356 -5.3536780500589911;3.8064723205465714 0.69017001824465951 -1.11447093670652 -9.8997914671954526 21.740706584541503 -3.944758328902501;-0.29090781175248742 -5.792100092768635 5.04112793121433 1.0906713573129372 -26.978744266107892 5.6715978711883794;6.3577256531203927 -7.1642161627906837 5.9095784251340993 -1.7339873845620357 20.675194067259579 8.5918241994706133;11.936719016414182 -19.863063571598097 -5.0357613832135835 2.4101117509635799 -57.895197865666539 1.0316951460462775;1.9583479225015099 11.84141776944338 -7.8543672646621303 -2.7754784403751533 2.4745579210036008 7.0842169456288762;-7.4070249019777217 4.9171740460764264 1.102133354842493 -0.42220059265187571 17.284712857522045 -1.3343535594775737;-6.7140231001467221 -13.758267317362542 -1.369550247585966 5.360881120315705 90.133427881878774 18.974632194178657;2.2196563894953329 3.0551445120378493 8.3292859811355058 4.3720952508257085 5.8362912700289789 -0.77819310152665999;0.86666083056580945 -3.8240628846372506 -4.6367017630253393 0.81998432121115594 2.9154637367236482 -8.449283604910411;10.03094673768822 1.551641713737554 1.2828700698889346 2.7611883259405694 -24.482607760884882 -5.0250267981361354;0.018180175115445962 0.096668566551288318 0.11330292754648071 0.58006340275414248 -2.8555378485452341 -0.27485592031354344;0.26443265517669123 0.19472495252448205 -0.11811787875211109 0.9141287762217124 5.6673690409558892 0.63775124423203444;13.88045939900285 -3.1528209621416767 -6.8381643561262928 -1.8287851268494608 -24.422461226983142 -4.6236896721192586;1.156349404891106 -4.5283007837170688 0.14969162115017262 -0.17307826299319543 -3.4618923072028061 -0.77613020246035525;-0.62545149783873788 -2.3793555312621426 1.2036616282260726 -0.13968314236164209 5.7600246380591669 0.22196377601331332;0.74825233916783473 -0.29823290937432512 -1.6019665933257641 -1.5101081229880153 -3.2309632356276072 -2.8568884836348407;-4.2184752155686196 7.1370280381547859 -3.6787096920700648 -2.5338037101621311 -4.7107375824013591 -1.6949213672157299;4.6488937249623996 4.1436296421627157 -0.28926562973258979 -1.5872894034592069 -7.815144897836003 -7.1021228587753908];

% Layer 2
b2 = [-2.6325466675950335;170.16633711416549;24.755849602088006;32.189575732244151;8.3902969056350312;82.884032478788058;-5.0414046293406987;-7.6089188907821894;0.88753042309992736;27.543110080038499];
LW2_1 = [-1.4128612171901824 -6.089720599545335 -29.595945081583377 -8.8856698018194074 -21.304713205612579 5.4328328169652469 24.884885656206361 36.538498406310865 -14.955024585550861 25.251315059469789 3.2384390495955486 45.823640753053091 22.497033217049616 -18.414373625084149 -31.536745577167977 -8.7448001152904951 -6.9216563668443483 -28.398090166280252 -21.566177325145187 -20.957830123163586 -15.213953687514891 -21.575426199388836 -3.2707036602522339 -13.589882892918414 36.553350885278682 -21.896506886338205 -9.1874784613607954 11.963431496566786 -4.3289859192246176 -40.568470083988068;-72.579231368261247 -153.04598285191724 -0.82202798322622828 6.8133852093458911 -10.295248088264932 -45.714624573520616 13.299656372881552 163.70084568564968 15.241924685393206 -6.4090847259059389 8.2069279386079135 33.142710890636955 85.974846965837756 3.3392261359470741 -120.83824785580423 -33.717574560691745 50.197994264464405 -136.28484323249634 55.864583328059879 -8.9716689555303972 -16.168584364201259 75.685188510613713 63.13235770108534 40.509765604358897 27.882861468204165 21.02089915749885 -16.091783419802791 17.632368995395957 -8.2888166069938674 2.6196603831327496;-8.0108611706693029 -7.8976915603522624 -1.6940243187787405 0.20073236651363491 -1.9583866493652791 5.3950530794110101 2.6108120365080687 2.5089676150942912 -5.9948681458472342 1.7264879942581934 2.3395177090753205 -7.1529801361836167 -27.933112699518258 -7.8871901705741303 11.295075677644668 1.4816283887908233 2.0025524226228879 -10.121706425505435 2.4675474642489701 0.59538708101226101 -1.7629193128455911 -0.7862938215682459 5.4314940619689462 -4.2777942153264705 -5.6213264419243769 21.541783055106603 2.2834762256505394 -8.7635331187315089 -8.4177226227676609 1.099597222836292;44.516097531257103 29.640458820618743 15.948831177529572 -38.273364230114993 68.39368362708494 -9.046359034676275 -24.115466590034547 158.08684848763937 43.3278468914624 31.598668518881226 -5.1671487317926186 98.069203009652156 -75.485345269464375 40.440062178184768 83.492253279429789 41.201576967957941 22.359715800677254 48.069807431747329 -48.485070168519812 -97.946516548678375 -16.935792673416394 -83.511843768095574 19.857376004398763 -89.624926163059428 46.015001163823101 3.2177572916930242 -38.986378066089053 -23.578980585353232 -144.24947243510422 46.811805645982197;-7.4942804377635959 21.564108487809882 53.04386758430303 -2.4228615619206901 2.4670545586328534 0.80560466715669021 -4.8876317321903597 11.713162998071944 -13.908074435876472 0.25096296936680806 -14.905947105556967 -7.91792503844054 53.983229808406513 -1.8834240399662652 -4.184537889970243 -9.7754146029108373 -3.695037806474919 -9.5394243061842392 -1.6374418641694435 2.8211920157536299 -1.5339311445205919 4.9111764738923753 21.328210840658546 -7.4933737801116198 -4.7467495282292038 7.857633080374157 14.092739492001119 -11.45977968167796 26.994325858646029 1.2904251658593895;-0.2013066421636564 18.736810182492725 1.4203052160920133 -13.700495784321058 1.2414322400400282 13.030689862285934 3.1072351299166532 14.613842653427774 19.5624633149658 -3.3815941247954742 -8.7297859337373254 -11.854522699303113 -56.443446779055158 5.5645277965411868 9.9759320630219861 6.7932036764210588 6.8808026489058722 4.9949358219697171 -0.23614467213537244 -1.2623853550338786 -6.3973112140639419 -2.1491846807771435 23.160336181323409 29.902070526843783 0.39321410937480239 45.959178658596763 -6.3558116846280823 29.204215829929332 -4.2236619751435676 -15.867365250252632;54.655695844651866 -0.038001150914470962 15.483244745039364 -7.2862575006642212 6.8974129334188659 -17.896830705897923 -15.512770408953356 39.293958794417719 23.434497020436638 2.9633700039731701 6.2501989282458323 3.8920382382149938 6.5107774971679877 1.2853703484377761 6.1712480471119653 26.150246927381332 15.345295579986981 8.1932428156615646 -18.209067311182181 20.701263971404408 -4.0240496197822431 7.4355471242107658 7.7981698900026482 0.84658610070591434 -12.304506721523556 37.41527352684394 11.359927392311473 5.4878965548306704 26.073972701567005 1.3378330370749727;1.1315430328756566 8.3512054582184714 0.70636244934369108 -0.76026308668330322 0.19372550771303287 -1.2081354456062217 -0.11580148788769269 -0.17001564810220698 -1.5095609405010493 -0.11954714270381328 -1.8205897001446625 -0.073997301985783137 -1.1656682880621916 2.3467206732816037 -2.098242066110779 0.32911215654202197 -1.0445441394005914 -0.50540629885122856 0.1290302964189306 0.52653033722784348 0.57739542500326613 0.31172624408974814 -1.9665016189361855 -0.73946902695706829 -0.64073219410542726 1.5505486277373577 -1.5438443328521112 0.14097305776936325 1.6352112436791839 -0.58772956679787192;-0.0084449530743151684 -0.008309841784786675 -0.00072256208588961156 0.0020340855509996618 -0.0020774670159112063 0.0085055295762195941 -0.0015681844467559946 0.0037385380572857916 0.019517892566027501 0.0026362939149139259 0.0028607998169994024 0.0018321003971878344 0.0035459904772589985 -3.6111158858410347e-05 0.0019575221635117169 -0.0036427486998833894 -3.4253716915669075e-06 -0.0077172978174576204 5.3079182946491748e-05 -0.00579998196354001 -0.0076172359259401206 -0.0016386242530073588 0.26611759737145796 -0.080794314959528074 0.0051153329940365807 0.046038400512742467 -0.011224130865713116 -0.010926502649623991 0.0026286021634265833 0.0024062814755915459;-20.27742748091141 -9.259736392926273 6.9784064819066476 -32.499184188138408 -9.8448953605742915 26.538732323662018 -39.890794901984364 -87.085389149178766 70.881389105906351 41.415369409300084 -26.084904974342855 -51.157790853470537 -171.88832611956957 16.575685935285509 63.226030171868089 -24.651129952953045 8.2395299029400224 53.603695050052153 -39.388318725113187 71.899497995377729 -8.0098457258614619 -9.4668095617047978 39.653331780509617 -53.137723944843799 -30.598383574243371 6.241786303747487 -8.2786346833990851 54.993379811804331 49.241645077932915 37.352231544816121];

% Layer 3
b3 = 1.7554494356790535;
LW3_2 = [-0.0022258722178146477 0.0023484799374994085 -0.00061202305669598347 0.0010654828388500666 -0.00024759863852080581 -0.000868230028343503 -0.0070192735930410689 -0.016245384070620268 -2.3709468765115331 -0.0009561503212233009];

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
