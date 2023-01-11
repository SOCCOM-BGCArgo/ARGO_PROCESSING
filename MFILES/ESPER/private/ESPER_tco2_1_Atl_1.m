function [Y,Xf,Af] = ESPER_tco2_1_Atl_1(X,~,~)
%ESPER_TCO2_1_ATL_1 neural network simulation function.
%
% Generated by Neural Network Toolbox function genFunction, 12-Jan-2021 17:03:52.
% 
% [Y] = ESPER_tco2_1_Atl_1(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999999390765;-40;0;9.1619;-1.99149636107278;-0.2178;-177.232103709389;-0.12];
x1_step1.gain = [1.0000000017546;1.00004398834998;0.0153846153846154;0.000309885342423303;0.0656325902039552;0.061098595529566;0.0461352500991908;0.00425400504498441;0.0145655815308426];
x1_step1.ymin = -1;

% Layer 1
b1 = [-2.5188884519329817;-0.28623718047062385;22.623248056212599;-8.4876071341927197;16.864719672406711;-4.0260942758081111;1.4661490481271979;-8.0589614160978016;2.8497875388775613;-0.58483274021999065;0.98406173904680783;4.8899527402668408;-6.8852774594789183;-3.0678679765384214;5.43716035501261;12.882457426666816;18.534345391967928;6.7514450855218113;1.2821690005402446;5.1204111763542839;5.0174604026780791;15.284862523431773;-10.493331220597115;1.4650977060060388;-1.3162846085617279;0.94346078897605223;2.588722954864398;9.9145672301511514;15.748439809532531;-16.489336104457585;-2.084009102916712;0.2810905802497613;-1.9395338152499164;8.4418189074341932;-1.0399155598332777;-2.5134079473853412;-13.457413033489354;4.890387361488238;-12.505725543165649;0.21748064971266381];
IW1_1 = [16.697966979591428 6.9920392404409482 -5.5560663702992388 3.8013711211050878 0.71404845074095535 -8.700280117547198 6.1625216394640816 2.4571531145657932 2.4855420220016891;-0.31496619935457354 0.65066899735017747 0.18402341462690452 0.034226388951443337 -0.37198170331664121 -0.48719104953194087 0.28114900438033441 -0.054411086239899639 -0.61127971891166877;-1.7938089659420724 5.3655109204192595 -1.9504535999893551 2.2652253637624655 -18.611427759985158 13.409680083782741 0.60163416499346822 0.24360505183130435 0.45940600660169467;0.72179941655792257 0.22040994292182858 9.6654661491262601 -1.3910991783421855 0.50503179067152737 -5.966089643828437 -2.779026078953962 -1.4937735034139852 11.546846237578531;-4.4869174978052575 -1.3483738331200696 1.2294518444920917 5.1354011553142875 -3.2995179090569264 10.789477150749883 11.005477779149425 -15.966526699503071 -0.29862535652389499;-0.31497376117422721 -0.30158238781193047 -0.33639666777349603 -0.10224833102441143 4.8546518760679183 -1.6142956144452092 -0.91564112327392755 -0.14913655198148534 -0.55824657451159032;0.9101231330405356 -3.2687308788198335 -0.83551367271127319 -0.36155029548398843 -6.441918006146194 -0.3666842115970772 -0.55586381151079456 0.54299031078624915 -0.62833835684931816;0.14795478539534238 0.18396119160078248 0.057705242579555803 -5.7834140997093941 -1.345859334592894 -0.33359711628359423 -1.186790508798522 -0.83440638899283903 0.30102932729072945;0.91831588944019049 1.0613657591683028 0.85909051521437374 -0.10678476292345465 -5.6096880465749681 -0.073380060313913864 0.73138589346911009 -2.2357246732796812 -0.19166003247561569;0.46039286986332556 -1.1124375638882777 -0.22534725291540925 -0.050087202386779668 0.50742921707570765 0.019843684223142196 0.28694282609453275 -0.2630119969992839 -0.021233246573455059;0.047499661923298152 1.4340576133799419 0.20535096834887859 0.037600032142448642 -0.50790527382697181 0.93283650760908421 1.2693650940410879 -0.52843170819583785 -0.17998430668189633;2.5463811682633883 6.2792373417062937 -0.76486758093654861 1.1057795579992156 14.953083115120382 5.3888882315200783 -2.1524756287106102 0.89962259513432563 7.5780308858142531;8.1148508813829761 -5.5329584124872726 10.97278684126616 -0.41622400078654287 5.5399831047754571 -3.0524829654326355 -1.7498639185682496 10.608193564507712 -0.43217213147538786;0.31436474276215465 1.2099717633208202 0.70439812149250525 0.31839395049902119 4.9763306652122399 0.40806724893119101 0.35514286354656366 -0.40506887585848345 0.21250959461397501;-0.086461569451484657 2.099995956199221 0.26188293698524001 1.7724232166949463 -0.86747475080447689 3.4396895451820648 1.9025003743305615 -1.9442327198994844 2.4670045403702727;-7.1935183160426686 11.445227881547865 3.1416855632783975 0.31567760318541244 -18.208766879011794 -0.64432448557357147 -10.46860960740279 16.30556365323249 -7.4892874690744984;-1.1766571301731785 -0.80237735928540399 1.3104590955046167 2.7175906693331084 -34.175755371899378 4.5878432035142129 -2.4996330623585505 3.0837737003225718 -8.1859338810670934;-0.17339571650152408 1.1952154314103718 0.067373645114368663 -0.57981483055835603 -5.7256258378741229 0.39976727351610658 -0.55075961196055945 -0.74529951832036323 2.0236210160386192;2.329032591221885 0.43853403107390343 3.0119593738791939 0.78434240125526333 2.8843518793462537 0.61911715600131223 1.2959589743028861 -3.5806209493897203 2.7738668393141492;-0.10202682898544538 2.4349657494186823 0.83147086434743611 0.60523283679304607 -6.1412444660508454 1.7203340219178873 -1.3659413029608802 1.255728467165474 -0.29523618748970548;4.37671943822507 -8.796201316623117 -2.2799002050770718 -2.6053788514552507 -10.448233622937368 -5.3996646173072724 -0.81125925249069775 -4.2026325516113179 6.3470163569962752;3.4185455077215972 -2.395273092864179 -1.5098907114510134 2.4896527104466846 -27.663135257368882 -2.0113410097790125 -6.4565475372400165 13.000490025468194 -3.8407998132191117;0.12933457661138947 0.87775314752520994 -0.42762571127626353 0.059002948271962966 11.261491664732549 -3.1867817760534467 1.6538282838202436 -1.4432005241899299 -1.9584278228949346;0.42882649047961668 -1.2208389876057222 -0.10999082599076611 0.2049400480274946 -1.3746505130413373 0.86645178926989996 0.52467851357033635 -0.59605400032869571 -0.11311985624368175;-1.0988773654968464 -0.52495901291975566 -0.94205017564578497 -5.4696831641502666 -4.8923478751190563 2.8001465715668314 -14.341021455537154 5.5724016599886719 -1.8885186299386849;-0.59875048743083736 -1.3306201411725604 1.3747518346322845 0.11482775692004442 0.017951631857738272 1.5732923438791784 3.9875346718711464 -4.7943073807799923 -1.2326512215430929;0.84909888458823124 0.94320126014884731 0.75933731775375046 -0.079398932989314788 -4.6270394521948486 -0.20625913533435791 0.60910127569358608 -2.0149832641245666 0.25036780485514759;10.626368624568011 -15.0786064280407 14.901404339193958 0.092016999260892171 -1.8691381096815751 12.119328168430542 -2.6432251023970781 5.2677942653865735 1.877595750556591;-9.3459916487813022 -6.9185311489316597 5.0292549488179121 4.493637410327084 -10.011730784337656 4.6231519552151426 -9.289305287947613 -1.1767364556687683 -10.774324787575811;3.3181581390665951 1.777763691702277 -0.55459539839651395 1.142111432186389 26.200460588671294 3.5239443825285464 3.4081483319762511 -2.2382022935220189 -0.73623614008010674;0.28764462205277769 0.31339372030733703 0.90627764725821935 -0.46774620824204705 1.6008791295130329 -0.42919429229002726 0.21399094987626144 -0.67905526425825835 -0.4766026058654611;1.9068368526512538 -0.26856008279550603 0.87217751950664579 0.7206909373767233 0.077844657262430728 2.5001271944780497 6.2537220510111533 -8.5719917062121969 -1.1798487814115939;-0.16459889682068879 0.61000330161992544 0.034173389835125328 -0.13604248783752573 1.4093938083053394 -1.1247408955545946 0.41514782335472911 -0.29272902311628435 -0.57168910914041449;2.6956774650414772 -1.612594946977898 2.6327080752584386 0.48495975899681387 -6.9996385232104723 3.9356991528762459 2.1877969569845779 0.52412869265680462 2.8352645676130708;-3.0024760862859168 4.2463544817770407 0.22021552498455282 -2.2851548195795446 5.8501304655859068 7.0942514097348486 2.954347883199715 0.99965409734265132 0.13652264473436052;0.3372067078346111 -1.4412809537554889 -0.4978266743243584 -0.33624920913790562 2.5080366078786218 -0.90180543003204638 0.53124672808127382 -0.37086084523636781 -0.013290442038677071;4.0728954887919047 -3.1517190117267271 -0.2642264109239203 -2.1250900050818067 3.8580485272765519 -1.3255695003818784 -1.0359106724611524 0.57514471319595128 0.37489444849349218;-0.17025808920372471 -0.2123924683703747 -0.7296023875259795 0.46943313405349985 -4.6110221942322589 0.85389506331343989 -0.23114600569097304 0.39098076239320756 0.98855358904791268;-0.11046036631048779 4.6970724509278918 2.6511813549233372 -0.55470117884586467 11.735262007833983 1.4555498898241745 -1.0336402753594636 -4.4425006752423055 0.50069177577826862;-0.47062284232834417 -2.2342003578609235 -0.67118062010821911 -0.64912705114376279 0.40672743868311312 0.34080528152421485 -0.81531190294481182 -0.42077025617597713 1.5856083463969566];

% Layer 2
b2 = 2.6856850225500759;
LW2_1 = [-0.0011550050624775448 -1.2435333673821141 -0.0069519110848407801 0.024685353378500006 0.004782973012033784 -0.1571127192615942 0.071025038463248633 -1.2180652449754235 0.40826074359272879 -1.897726372574013 -0.17366444580615104 -0.0057806622283049433 -0.0073477814947670445 0.22080360108586217 0.017088009383973261 -0.003585836636601059 0.011424497501522994 -0.14577316595886697 -0.016620114388146721 0.12357096178031068 -0.014126346454068553 0.0062942296792333416 -0.079224079229530966 0.68974649767492158 0.0025783144060611304 -0.018618568361989881 -0.52449022264994383 0.005095569068589165 -0.00092700743474724578 -0.0075230575840620744 0.32864023322196945 0.013759558385493262 0.9651661113868687 0.026773278007138839 0.0090682637495756169 0.82884672975870866 3.0081266432422464 0.32507718907707728 0.019061340080553695 -0.065512272064995811];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.00115315130997576;
y1_step1.xoffset = 678.286531932822;

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
