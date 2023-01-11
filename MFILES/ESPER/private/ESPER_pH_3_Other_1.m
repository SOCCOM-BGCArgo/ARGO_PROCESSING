function [Y,Xf,Af] = ESPER_pH_3_Other_1(X,~,~)
%ESPER_PH_3_OTHER_1 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:21.
% 
% [Y] = ESPER_pH_3_Other_1(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999540494077341;-75.751;0;29.715;-2.05132667284751;-133.803853032615;-0.1386];
x1_step1.gain = [1.0000770825255;1.00022983623563;0.0140762793578401;0.000302617642608564;0.251952632905014;0.0588888685539285;0.00432440760873659;0.00826143816767911];
x1_step1.ymin = -1;

% Layer 1
b1 = [-2.6774866138822974726;-3.3562819442062497011;-1.7706854852829869706;-1.9775168385958883377;3.2760000659082484376;1.1727869825611476617;2.7558596554206262752;-0.10053367877238020278;0.34364268619399707116;-0.5462442923411664486;-1.8724413770823673531;3.7620731543272567521;2.0130657463030656729;3.7892189208346240825;-0.2965115702958379118;-1.2946822939034658972;-0.16313760885780906329;0.69577159089820395099;0.3787179916670326274;0.80115937002780401865;1.1967874895138530533;-0.52543812038094228622;-0.15999431099783686938;-1.6760894473938074967;0.23394452096202217328;-0.26097600887979954809;-0.80470097955475894391;-0.20026917514522088459;0.58684276318868977551;-0.33696206825637248228;1.6366864594972887481;-1.7482304061689746266;-6.758746472539154837;-2.6719252997096654134;1.1345439843999760576;-0.38480471731526494361;2.8802485890763600551;1.2957699714562036508;7.8338562301956917722;-2.9705001005963009852];
IW1_1 = [0.9264747001873373522 0.37694117365716839352 0.40962574855585720535 0.86791078197480231449 1.0972127521596872768 0.87770459696621794254 0.52797832921259202354 -2.1975258496257006335;1.5522311361104828009 -3.3425088277162782191 -0.29484211470787402565 -0.034653584605142694064 2.5653221505214487408 0.091444249112496295284 0.24342589182524185243 0.18294154010006052902;0.13103370139021860985 -0.006980478546093996442 -4.0004018542031447581 -0.53836011520467952352 1.7630118862393369117 1.9408992674788330124 -2.5288297615984531852 2.0587083080051225004;1.7380341126722846923 -0.16989397977029696563 -1.9888784443962947979 0.10206527926407750473 0.26990916633095368704 1.6214385459288604974 -1.6483766634290943909 1.5464431710341772774;-4.0403792012973154257 0.37420286608646613491 0.73361742733212975676 0.1203582024340625195 0.22400930370278165915 -0.89557798240636887765 0.30137864710395867451 -0.75297877689947390056;-0.70175897072645743169 -0.066071003200582595349 0.51284307229421521868 0.21038442181366842632 -0.016359599309234071352 0.43790269701088435639 0.36631523871300458728 -0.080513540737507174883;-2.2963749918724856158 0.093700688801220718971 2.1720701479698139025 0.038509573830625191959 -0.089494496718768135057 -1.465003414495409606 1.4672304431186353746 -1.1259906337275316091;0.51447568407745347674 -0.53506872707570041214 0.35935958227781228391 0.10138078293027541488 -0.2539485422105769441 -0.66639062208776989049 0.45411976169006545989 0.25955197464685114683;-2.3983170671423157394 -1.7449752367981434809 -0.54922686446177193531 0.27450707749489061671 1.4888408427784942667 0.61954427890461560846 0.60453523917419038636 -1.034310222635573373;-0.0603517570566081199 -0.025755151660034357464 -1.4864258183671503755 0.0043722661797649846685 0.1186610372948390113 0.45197151019052511645 -2.7692324132848513685 1.7486197403391618987;2.1422158750738833355 1.3683594186057259634 1.4854260478238541232 0.61100561540259812343 1.9983462571946513098 -0.66882636715183485521 -1.3095991553552523978 -1.9076654219845443627;-1.241647357026699483 0.079667469248567726603 -3.7182430863224285034 0.1121508018274099544 -1.9808005542283746436 0.13588981855920326258 0.98231734016162230994 2.7472002463904470204;-0.14350473861342707971 -0.28798790456924261427 -1.4777540611964714579 0.48329320863899355443 0.20335777451624509338 3.0372572953228749171 -0.79069582405927973934 -0.71773473016247968381;-0.1735092807019853467 0.2422693339796251899 -6.2825688190650605947 -0.74514045560855513006 1.5512527926669710254 -4.0202199271761918098 1.3282967566520145208 0.69176879817070646705;0.29393621394706132 0.29109396743163507004 0.27813887893316791011 -0.12068595123220168808 -1.0626945874022899474 -0.27338048033848744156 0.32498398123270177962 0.25387412924204982856;-1.1019050833876242734 -0.71102860291269909254 -0.54106880042918592544 0.39724824961448551575 6.1827843245776294978 -1.6098343720861847395 0.39389596093570433677 1.3491364773561633772;-0.28439788310648284053 0.17755858491759282658 0.41967429174411369219 0.25127042858217857324 0.2111989399337960982 -0.38601237959282302947 0.50045261483020331283 -1.4624743512268882917;0.29634454980087326925 -0.36378986788692385668 -1.0790517273387012498 -0.636761305231387853 0.023158782992546537299 0.33569003738460817532 -0.53733371149630981822 -0.64931404562364669619;0.24665992379037840476 -0.18824719067711567777 -0.57231793257528329466 -0.24087132510611103386 -0.12528013187417724739 0.79706448738760948824 -0.3217726707449118928 1.5318463042348873238;1.4789208707744194715 0.69471849390637863131 1.9735434593778355961 0.4963164668674226454 -3.4051970707478700362 0.66272771428232546409 0.005867569486582364241 -0.21031059075432045891;0.70425589389964982079 0.65614595502301131802 4.0087841656803977486 1.3931995625217759027 -1.1512368194972792956 -4.5120916835921676125 -0.023420811125781287876 -3.8902569764000052288;-0.18603557169713594299 1.140584823727449848 1.2058113758762769585 -0.3770527597641856743 1.0333079509321845268 -1.7197842812153489866 0.25923904771685957193 -0.043042159734588372744;-0.45760691122133562159 0.46315443264789096789 0.68225408778373686225 1.2408923587155471502 -1.1823031285661627088 -0.84289838366794511693 0.48750346999974036866 -0.0089551451373128723027;-0.47559707400125184051 -0.6125239658136311105 -0.42605554454993826985 0.3291682427243928788 1.2086757212283327778 -0.13883830619907860937 1.3371908929654170439 0.37943798648557441888;-0.30561998668997314432 -0.51095724542559195491 0.041063831053600785559 0.2039446716397775583 -0.33723090044933112797 0.38036677346269887057 0.092855102914865109143 0.25277060960402802925;0.17007088855702337526 0.60978772517463097902 -0.3279234310891616655 1.0947463761178455321 -0.91775629588902674083 -0.030541383315376442048 2.0636964307975862454 -1.2249248510261667722;-0.4360554311396678262 -0.12683732889546381339 -1.346014855999018156 0.28054066640847591385 0.5024002306595370948 -0.6706724945315727382 0.62869329311584731812 -0.47460540015937879854;-0.34909272518373590799 -0.053521759943581805052 0.30182378152778577274 -0.21747071863028100136 0.31666659999997476271 -2.5593962743196998488 -0.97704599599311681146 0.48582415971570863666;0.30406182268726028584 -0.37006428908959171364 -1.1275679342369084868 -0.65163366744146045129 -0.48036435766865731578 0.094376516161353979362 -0.44995137986403610508 -0.53955441734268916321;-0.45156376765764766334 0.49627521163606480581 -0.29515572074335139074 0.53855636069006063149 -0.22828425705333821316 0.10561611162451220447 1.7608959592826243234 -1.5319032945710053806;0.48322106446325302498 -0.93669300506242136883 0.060568830645976970417 0.19650240903832882911 -0.081517186101824234967 0.56850460782973455398 0.089324775384039878379 0.32986202818687621052;-0.16059348504559564641 0.066456052666176768096 0.092737570391421705884 -1.2553859843453800682 1.3628154859283152511 0.021080600311733761754 -1.3836497438686747863 -0.48161293625049139688;-0.49167443547754580147 -1.6051524081466883587 -2.1767236422784188576 1.4109459562443495972 -0.847915618314185493 2.5445669547819429468 -0.25947596250679372254 -5.1261240455144045924;-1.1017179735397266693 -0.81326804100181349089 -1.1527326879283665839 0.049570342910926233626 -0.31916033995117309985 0.18253708297639809066 0.51515884455139382414 -0.69736221250994190335;0.086457990054946073921 1.1092163951563684954 -0.74995629450293965768 -0.38497944656830668242 0.94072635334619980618 1.1675481388617305889 0.82870704676678319522 -0.39395073700910959724;-0.47700190735265435293 0.47416022847653932493 0.69588841319700345256 1.2038265152987761564 -0.88604253568066704272 -1.0050028840125166951 0.55684746622010561357 -0.019991335969755285457;0.23998206871586144628 -1.3083665745459902574 -2.6982151120816895151 0.17519036375717153553 -0.77476138965221308386 -0.43289022017647943441 -0.11472046899735341419 0.011109927770870128491;0.24705428680681448816 1.0207172327981570703 2.0113046681804607196 0.56390858076801353427 1.1710437076654951483 -1.8678727659406879624 0.76002795151763902215 -0.61162095523265347374;2.6425595000082595654 -3.5301280471361500268 0.78898218050599211448 0.47189346152133448964 0.44404739426922396062 2.2391838129372803046 0.037746415833013025054 1.4741132723495091739;-0.10443040286337183786 0.46481727000514211356 -0.32412237452458553255 -0.66963136796918099058 -0.4937901492981374707 -1.2429656795254422352 -0.52498289302518019728 -0.57857304845516210889];

% Layer 2
b2 = -1.6240015457790466513;
LW2_1 = [-0.73371928125386209896 -0.11551561660080057747 0.1603285209260706945 0.79844036652477523575 -0.17105651641453745171 1.1172858120724762365 0.86332993305803007456 -0.93169817799222431898 0.080147366644392023516 -0.38399780299802560579 -0.07445844180726671202 -0.14860489571520202334 -0.14623496143942821268 -0.084028983359118167296 -1.4345518172674276425 0.073218824300087470558 -3.5934460651727593294 -2.3695303423414082467 -2.9014840903539877992 -0.091688673243689469694 0.088439126261123074402 -0.40210293546835645939 2.3334303982245687692 0.6068143378590428183 -2.7231092105414385074 -0.21380957552820473166 -0.49484978378681743205 0.20084949877883892522 2.3920800160847970339 0.48033355530413684997 1.8222933921045736394 -0.52312620658901853865 0.50755138395493637482 -1.0842760946317560045 -0.24022138583079050145 -2.45992817605441072 -0.21064173908323935436 -0.41618895168588104783 0.28665983812486389715 1.7028744934497637153];

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
