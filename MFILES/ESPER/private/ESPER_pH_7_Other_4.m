function [Y,Xf,Af] = ESPER_pH_7_Other_4(X,~,~)
%ESPER_PH_7_OTHER_4 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:22.
% 
% [Y] = ESPER_pH_7_Other_4(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999540494077341;-75.751;0;29.715;-2.05132667284751;-133.803853032615];
x1_step1.gain = [1.0000770825255;1.00022983623563;0.0140762793578401;0.000302617642608564;0.251952632905014;0.0588888685539285;0.00432440760873659];
x1_step1.ymin = -1;

% Layer 1
b1 = [3.0139985989923792076;-3.3598730951672584055;-5.4649959443234292777;1.1622356990026592172;0.13739045910334302114;-1.0186716167548128809;2.5833463773573135391;2.375064719159389437;2.0757845622499413452;-1.0723286181068654166;0.29764291629221456548;-0.93593879844238148102;-1.6005601509694478146;0.31257525899014143578;1.0999410459256158124;0.64938844163718789648;1.5876707894361388185;-1.1403047859093922156;-1.5345325856136333709;0.8455653796892685925;-0.050337226098691553733;5.6903976516464780389;0.0030128124529946488652;-0.53623946024700430701;-1.8348001931042015045;-0.43586585575084241251;-3.9082022495648591764;-2.0687613076261217593;2.5140622086029869031;0.96444253175405658673];
IW1_1 = [-0.78154375742757065115 -0.36349606663895689396 0.19384880068258597552 0.19630183761991523705 -0.092050645191806002021 2.8675246354543308946 -0.53981039204196101799;0.51810946332480756205 1.7533278603844706378 1.2434842251611426533 1.20458671022346131 -4.5248871772277157177 -2.0348239195050452643 -3.4383192552471615322;1.8528325606316273788 1.8069959155326347222 1.3331835711894446383 0.80039009413789896641 1.7728559933957055428 -1.2057731782550489275 1.4757664313628491293;-2.0309349846771764803 -0.078965602229843959159 -0.64930297506740841662 -0.41342245740119448749 -0.3145102303052047521 0.58517936237846623904 -0.39547915601725014412;0.18795907191196684072 -0.45853162504218542006 0.63509839983268345787 0.23834486377279207048 0.89808705336763161586 0.73157348551128287273 0.33311019993775542458;-0.094254847664343985314 -0.50062246246823405293 -0.010436560371829432803 0.30607310374271134989 2.817489275450241859 -0.85302937807647294388 0.097500921965014669679;-0.51930463339632415565 0.909720667682732409 -0.54275604048284131764 1.11265760147368864 -1.3232767466857568994 1.7006177063922611747 -0.67973038526839169027;-0.29245412522645530817 0.55246247469959530729 -1.0377584340174701438 0.49260533150925078782 -0.27834241345937188461 -0.21627964115584652482 1.2351474619368019159;-0.47223938437188023309 0.66243567534361802451 -0.27933643372808492078 1.2860921332362229563 -0.12329936302818692018 1.418967310748667332 -0.69389754011417359436;0.60896059288560355771 0.7403084263508156182 2.110427308208739916 -0.35511854003564480298 -1.5009223236257893674 1.4359356084365473283 -0.5975376564461966078;-0.08413580599612523192 0.16775472349200551414 0.31756634198150873738 -0.63023487014569712361 0.76300068024965206348 -0.23920935470316007088 -1.6384614565814676901;0.80682435526021611949 -0.054040218348251137237 2.8241321387337845472 0.49772767192463784891 -1.0561278665153150946 2.1918086432186116852 0.84695951539711555878;0.15149020048316594322 0.86663855717161242254 -0.85281161064087596113 0.17858399224123452931 0.096469810685992793409 -0.38687921852455797733 -0.074456355768502346937;-0.10107299696512965503 -0.8193682227983887234 0.91943346594076236222 -0.082139443275190585081 1.3768684875038055271 0.093341636973482508544 0.031895857884581366171;-0.091888480376569173225 0.31679289647729502732 1.1448173563187202806 1.0721440377928097565 0.89511327386666827977 -0.44853155497033747601 1.3766885696151311969;-0.50787145954572321394 -0.29126002905404940746 1.8585730765702896949 -0.071125463118989915023 -0.50682830366884501494 -0.2206149930399378678 0.4856439711908286383;1.597088550670066498 0.57728052352855918627 -0.65614421650678711462 -0.67382890660194227905 1.3303358578254265154 2.2418499137742764482 0.82289336482651909144;-0.76827634006161649172 1.2742488927456909753 -0.016394712265050754435 -0.038866884919012807675 -1.2804789238959228825 -0.38268910736759553792 -0.29082869907605302107;-1.3121653435121956033 -0.54399857582426658364 -0.51347648562516146065 -0.10406791947071070137 -0.39798410526141836385 -0.63050815628773781985 0.2616359668544772088;0.50637112678444573532 0.45227609028011567238 -0.93280436302701719153 -0.8019585171564631576 -0.74558855287138348089 -0.6514362146622962868 -0.83974296861682096171;-0.88220264847127383145 1.4012145996459395292 0.77277983396583249665 -1.1212783184818537752 0.22060274981013805462 0.20441362056540904324 -1.7858429200379517532;0.7928995558656840581 0.93380544160484380001 1.7370543957656545686 2.8426073296357494691 -0.84250068597624361288 -1.9176610210387439182 2.1407053638321547062;0.11938588048976178646 -0.4583132948923465122 -0.3708949553467055793 -0.50762459117745228898 -1.0344686831386367132 0.62650379548356627257 -1.5052317363789637916;-0.55259739749488157035 0.048773306334223928693 2.8911853050805880194 -0.19189777815224134949 -1.2208508189065088168 -0.46987396905338041364 -0.6388437224241079182;-0.44474550642835852798 0.18330234807328049662 -1.8722343414296354069 -0.18757793944304340661 -2.0598982550057884922 -1.2404732693761246765 -1.1039317526173240136;-0.96742631776453436032 0.71591631812692035908 -1.7915209676045764642 0.094468254213125907315 -0.97218632063583254599 -0.58694851367520195495 -0.51946251564048162574;-0.98758072823640830062 -1.8643180018813485521 -2.5303880305813755136 0.046301880427176797528 -1.0234705280289757567 -0.47656996847212357293 2.5995071713212816178;-0.066519214904966539326 -0.18151786270360265041 -1.2149126853150438166 -0.95374955315368947772 -1.9912747556148997496 -0.43525065649318506189 0.51984992830920706552;1.6830771826315971129 -0.6889424109962883902 1.3607866682078493259 -0.9518885985581988507 -0.077857359024664188341 -0.10359657503562703929 0.63345049226173777424;0.41704758377562839433 -0.66045325556344525175 -0.4560890966827588211 0.038460633879961958481 0.84133256970426451016 0.40881361758348977453 0.38090243596147627114];

% Layer 2
b2 = [-1.1412191962545854551;-1.5101024244753638026;2.3552448676651249748;2.2005002616561450957;0.19911435656915330017;1.3647313459239818201;0.4925560476704137769;3.306445319893480228;1.0587538772442937063;-3.0218735995705943154];
LW2_1 = [1.4450311022582558174 0.29156378756059087021 -0.77692960590394299913 -0.37891662696528222165 1.9417850463588097742 -0.52527442678499169926 -0.5068774503646586993 3.5194867202707560416 -0.038606354199491738077 0.51146607394698584415 0.31612958407200375532 -0.46215567409663288467 2.7624413708317301364 -0.59056075032049504259 0.57434470018052552476 0.68437703394478277552 0.08653225013034475388 1.0174191587281948834 0.65188746750526038731 -0.15907018198627467798 1.0564890979518173353 -0.11667044224720991841 1.0941197725177689293 0.0048855140610588113081 -0.94305864421728269331 -0.37351255835646990366 0.54735858292925432256 1.5008850400198994191 -1.0890418541851500933 -0.69793530670253411952;-0.43460659298680154317 0.55955727366788732269 -0.40721173555302098146 0.53005550449227956111 -3.6338637905729136435 -0.53585136547980938282 -2.3942516448780226668 0.1604225188979202521 4.2604736490200068033 -0.2533834995713805105 -1.2377281344986457157 -0.11698063154101763506 0.49653364042690134594 3.6017278321857157941 0.82538532239788520162 -1.4096431110259208808 -0.15985280337079868507 2.0949624738538563662 -0.78224885007552025229 -0.60773607514945870989 0.36652261457723767091 -0.15766663132108529877 1.9604043944680977773 0.85251658655548334664 1.0561690844670497125 -2.960089185423689262 0.44022116131512356585 2.3532860782366258334 -1.2607275586040231907 3.1815623940158315186;-1.730152785993513076 3.686554156379074243 -1.8518625930228869159 -0.10687771439214183156 1.2984261216036179221 1.2536219026580028846 -0.26123256542946088343 -0.075627161056597669297 -0.16375147920777377575 2.2144573788906112277 1.8702731346371868071 -1.8084515900462529991 1.1507336168038508895 -1.6606886111193077671 -1.7116803466937995371 -0.62667979432765930081 1.7459897183309234237 2.9714626649422046789 1.2551235385246455323 -0.17976322707082709185 0.18298890319448396502 0.66956373358619936109 -3.018979332538625826 2.1640464407914299372 0.42622484996214188291 2.6228542992763053121 -0.34527515424509541697 -2.687551474179645794 3.1312724917773127054 3.1276941675975580281;-0.43598316809894493984 1.5666947516246116923 -0.83519834805709392 0.52967201247199580116 -2.5408356267331186906 -0.16506164751581095529 0.4917673945407745606 -1.3549459945104973535 0.25558117380400385832 0.5522317013333439828 1.336992767792849035 -0.4523837929957447157 0.25856974942426841357 -0.42183243641078138619 -0.014955798965042902543 1.8284780607978070588 0.47804102543362431144 -0.37179313720382212161 0.78146526968146423631 -1.4288920530596995739 -0.65967708071216157872 -0.66464050844704869458 -0.27983551044184956647 -1.6996781687405946304 -1.2877490237979909438 -0.65054522728405383969 0.229836433769331272 -1.0697794004514167288 -0.15448034948969716074 1.7824939136360493652;0.11928018256107776607 -0.30307857993876069491 -0.063489799875168173582 -0.28226411891053626579 -0.53498189110089477971 -0.19050332478309500717 -0.054155222990203207012 -0.083854116391800959107 -0.14888940702961514506 0.35011574767154762799 0.25541504780682994236 -0.1592026934719451281 0.088038478877896342301 0.20197705207128413551 0.10572423085322853187 0.21434705135883336258 -0.018613777365201027159 0.64834798782056601851 0.021496119621679566497 -0.58618707930256708227 0.031906496104518039525 0.02337871132508657912 0.035971427893733490255 0.10389510651451765166 0.10204133976754120594 0.011888700274169257051 0.11705218300468227255 -0.14132042007041054066 -0.0095795537304812598822 1.1758359887344986827;-0.077153046150726531383 0.52733782518311544951 1.0902383438496023871 0.13420381766922762079 2.3886086280320677488 -0.23612888235899989753 1.1738266413039741831 -0.8270094826185236947 -0.65207386420188573339 0.6670232831143410035 0.45109326876093636116 -1.3486963314284496906 -0.29451886309289071475 -0.046067147191149238061 -2.4293917533671831599 0.30377713625383634399 0.015331483559627236829 0.34096193729541673934 0.33365260548689507081 -0.14725602948392102487 0.28842076882197731358 1.4598388456447768391 -0.29395080717991151076 -0.51192819890656726667 -0.26881041821286211224 0.60810974693772956723 -0.24982832207333957597 -1.9038239534702030564 -0.49452526330026363999 1.8397546011210759165;-0.36651342937894576535 1.1917155000052033209 0.14784529290388512268 0.37003920678840707126 1.5943094067792427104 0.051207818954223799368 -0.11211764391099218741 1.6700696270848149894 0.16576111336259188866 0.39102069801967481721 0.11594600145238254107 -0.48017012134734982087 -1.0043791832865032188 -0.95796949511812501399 0.033189459796685784332 -0.49289331697372579244 -0.093530586588433148321 -1.1946063721564923465 -0.041031429349667206119 0.38292008493286472381 -0.16087772863633273412 -0.38704417718263345316 0.16352993064622597941 -0.6126560124184816658 -0.056897376005775622543 0.45138773118086800151 -0.35400451432618368353 0.21318916771301138491 0.53983593517149131014 -2.6482384037357702411;0.83170192389080765949 -1.5742212534048709571 -3.0071059564265762987 -0.97839501017761898449 0.32903056908821093351 -0.52342029256584976515 -2.7961002385121700442 -0.63693257494484611936 1.1263694543759723565 -0.26300854654356004003 -0.86088071544083777376 -1.0236672496223337969 -1.3407187495787067366 -0.90534691320267179915 0.79460436756076302522 -2.6076478367866688934 0.60104644491891834335 0.58119147157771888779 2.7783523276433990112 -0.30330059691795280274 0.41017110245426496862 0.27352962761280324999 -0.65257052326874276371 2.9071811663873958764 0.23273196670158383115 -1.2972888047239354847 -1.3954113828847571899 0.50418087147422452343 -1.9969315520069828018 0.25910932539160491661;-0.091798711756007372875 -0.85211606720727472997 -0.51979631768462197794 0.20502458708669191179 -0.57864110950133285094 1.022331888064767158 -0.21102831731737817234 -3.4327022768060402313 0.56661206389452911658 -0.70871495818517404874 -0.48075581728460353714 0.63883688208119449037 0.73315636243811743711 -0.14118734509107325192 -0.32161580729796207523 0.048717338050538752658 0.09121788680614338618 -0.4309928710729932777 0.035685869773270303962 1.2372823050864547589 -0.01471058478264577335 0.23666449589444060742 -0.58739119219636204416 0.4882054081032494719 -0.3276565848000699277 -0.47796197031457804405 0.020658136176099402975 0.53928638110854587584 -0.3535871622041783513 -0.023075496419455710645;0.10193470863421634498 0.50388271119404459686 0.84668532865105239349 0.5504625310815064454 0.35134561451539536581 0.30899470651725918779 0.64511671494972466778 2.0115119593778447715 -0.4308988218605536713 -0.27872287320331101679 -0.031998768219816529912 0.18531231662532213211 -0.32961982463245298458 -0.027138653679574350347 0.077499396398276562681 -0.58386493692598007144 0.2041838474450898 -0.52074851186371506895 0.030186139459900791782 0.67902219426459498042 0.12543855320418334509 0.0034076538700984139625 0.31660541037793377139 0.37728693640692972933 -0.037682386130147196845 -0.1375354918601963683 -0.12374332658514940186 -0.48953661341777765514 0.17620716210484124731 -0.83678061438969364438];

% Layer 3
b3 = 2.7935344138501951861;
LW3_2 = [-0.61819927607782199086 1.3997758631785222061 -0.055792220404219071417 -1.9374235291433639095 3.0286228460136856988 1.1051986044790869101 1.1834668793182832047 -0.85469617444165491005 0.98606132848340033448 1.7231495546530621876];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 1.90239167587468;
y1_step1.xoffset = 7.36300175998845;

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
