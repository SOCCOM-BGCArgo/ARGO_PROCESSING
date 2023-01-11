function [Y,Xf,Af] = ESPER_pH_9_Atl_4(X,~,~)
%ESPER_PH_9_ATL_4 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:23.
% 
% [Y] = ESPER_pH_9_Atl_4(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999025224415;-0.999999994516886;-40;0;11.29;-0.2178;-101.945317545701;-0.12];
x1_step1.gain = [1.00000048738803;1.0003169936794;0.0153846153846154;0.000310077519379845;0.0705602618858232;0.0479365703301392;0.0054210015363163;0.0145655815308426];
x1_step1.ymin = -1;

% Layer 1
b1 = [-2.4730740493936886892;-2.3602307437834002535;1.698568374537291481;-1.7151513859124316586;-2.0027782717998663387;1.2020716850579702406;-0.99278809287418201546;0.85241763529701630731;0.12071901204937182661;-0.97005371569024667622;0.53036847939623155135;0.57815873584676724661;-1.0879511281945344336;0.41501827185582956803;-0.091206391268339695455;0.64680590416977290502;0.28793455905055437283;-1.4705210978056755255;-0.24925400589703761156;-0.66307608840392340266;0.066838459086193582892;-1.2768469721867881095;-0.21425856755573252643;1.1828449228579054342;1.5079511981358189043;1.5878146702266699997;-2.0000331076652053675;-2.1172972769262146109;-1.5757832515257346184;-2.4124581044064390412];
IW1_1 = [2.1016745301305124194 -0.26463217902619218247 0.082393617898634322061 -0.80237975147745022042 -0.84675660459846824324 -0.41133232435361599677 -0.71406695251178298189 0.82274982990910960901;1.4151418513380524811 -1.3704613112977301892 0.70358742919835226548 1.0083620246212199323 2.4054368537040855536 -0.34624075636967732716 1.2939641348298409618 -0.17146197867063825426;-0.5099608561316150146 1.1573936667549153956 -0.28108902555706882787 -0.40905913797575421675 -1.439251453785777457 -0.69813927118095642221 -1.0689842121054951019 0.26623076395649919457;0.24450267735884032017 -0.085820773513397397014 0.073634925308333737393 0.38016965154114173986 1.3862803898139988146 1.7444566601763675795 -1.4874130985012916639 -1.3658954604725239701;1.3280858617805271216 -1.3900548378168586794 3.2213520519173064471 0.20623076323710270996 -1.031869352055002409 -0.76331113348848167455 1.0140736024250307512 2.4860621696143927473;-0.15330320897452054396 0.98428076630359395693 -1.0639782261234689287 -1.0424169859614231815 0.41404512082571559572 0.07480760911997751117 -0.57626334483748564708 0.05118334165997629609;1.4219178194403216953 -0.41934315196421312333 -0.5551713778540026123 0.57194846106644470218 0.056263023195483603811 0.86930005053276815374 -0.97921243811152591796 0.34730682663199285987;9.8811010402463476943e-05 -1.2015092211353779472 -2.8356002137576843403 -0.32587454413199118397 -0.89055264560286584441 0.97106615636790782897 0.78352301584081973651 0.47785782993096792559;0.52692030369970876968 -1.2909792899231358199 0.54227259691151019627 0.71257110364950904113 -1.1678890022609096366 0.86325966671909881711 -1.4839225911667974778 0.63096478137054834434;0.68810732460969559199 0.75678465338976630417 -1.2884270016866430364 0.16457635205594686734 0.79995277787502039324 -0.80912729750206457968 1.5387731220092910611 1.1566113903082608783;0.06327545676292059551 1.4055438303037581793 -1.8859783934470957245 -0.43499459242720894814 0.11540613272402905354 -0.85800875890745198493 -1.3655186185484540218 1.3408186332074647495;-0.70195102739351611376 -1.4220817444951403541 -1.8868029069346257121 -0.27339665816146641042 -0.3006030762962639824 0.023099156554960986049 1.104578609941835099 -0.34722776170776692428;0.17176975364532803825 0.78936029120431994155 -0.44324218516520785816 -0.43631166716304420339 -1.3723512949043044085 -1.0340507747316030063 -0.47371638210459687812 -0.12914893212073821904;0.93679041309204957866 -1.0760093531250451981 -0.27677381496304692687 -1.0794593520091195771 1.6191437621927315949 0.14573461712237550936 0.50775439122792986879 0.31777768814007856246;-1.1130338082139556999 0.71863653641485902845 0.19090789231714286323 -0.39550561874107070803 0.18672168858456997254 0.61717126619407680632 0.5886150003848059642 -0.58104182928334047187;0.66990931490888627753 0.087885473397656913019 -0.32231462035740654493 1.0548235519747843902 -0.62050661375158966404 0.5932928458470627664 0.18405764625844689419 -1.0414321281536800257;-0.68534402121120752138 0.23749973535500887301 -0.059506011201665778843 0.45705818201704689896 -0.90307379875628901456 -0.6700510642648110915 1.418344226327097557 0.75845663431661713982;-0.21820468294276826327 0.44045536233465693021 2.2963021116902431196 0.31626785097993681539 1.4524540558092360332 0.16803041730235268592 -0.1233160402623864893 1.4234513841409217605;-0.73181300847708818935 -1.6346516615193942279 1.0840458395720100615 -0.1732129613560324799 -0.16542062139127011733 1.0480486909384754135 -1.0011768824664024624 -0.67289338335205883812;-0.44460990149437151775 1.0975139461054654877 0.10286089995816084797 -0.81459107269723673994 -0.72452151227190619842 -0.5613732950790183418 -0.021924735751057286559 -0.40011191741657115006;0.86116006168211944871 -1.4271409126999805572 -0.35088015440819481094 0.39004085468825755445 0.55244938675513188819 -0.66904073486615089905 -0.72773275125050995715 1.5002255043922105848;-1.0438573819512229068 0.95481610528306681829 0.92827889488679427021 1.2775081076403034253 0.27346425909754323857 0.10010012399833503882 1.0885493407339994487 0.22514367769765228355;0.78983947161574008522 1.3382246066535956608 -1.8122231674672497626 0.38652489550443958066 -0.58851533101426600769 -0.90496655781769252958 -0.25970036997017986735 0.37131534613727140703;0.53820469600441345204 0.62494362154941240917 -0.41904310703053704312 -0.74534424331061799318 0.9106303990859203612 -0.71838607952410205026 -1.3843405219232167713 -0.2355010603819599202;0.49622546541768669526 -0.88568302961123934747 1.4833623064003955694 1.1091953613531682077 -0.34020695808661288728 0.03648977451701507918 -0.40542496249571219336 -0.64660076865836069171;1.1029531178367737443 -1.1684649412250813771 0.78007321281820585224 -0.64765899167412166459 -0.90206217039162206461 0.3684199307339537377 -0.13265306555009312439 1.0279688301635447889;0.61245823640711682945 -0.74447062505846361979 -0.4343048180606932096 -0.97845160473730663497 -1.5099622551788924874 -0.82969360977977002491 -0.079658921725040521267 0.97254360008424067807;-2.3804390099440420059 0.38363623942865909822 1.3402393290908105605 0.1842143478419720759 -0.44305059911816324636 0.50640559376306160111 0.18497617355111758974 0.77336148905339818782;0.09614253399071037709 -1.5806671440839812348 -0.31757373457840165631 0.75919066582119398845 0.98248470203401305767 0.61772661229644709735 -0.70359034107375562073 0.75087915360248669305;-0.9858563863952990447 -0.68116783847109219074 -0.70161798039247491587 0.43019257019856299573 0.25470429523376214576 0.67895330158058808667 -0.39579196985847220347 -0.057362330112227563617];

% Layer 2
b2 = [-1.5633844773212732804;0.9203006271309479569;-1.144701515049530105;-0.52876691904400818789;0.38963691985854970179;-0.087241457797386637307;-0.56072313498109338514;-0.39619151764496229884;1.1837680374225252322;1.6684885831743196949];
LW2_1 = [0.075599085057009685062 0.4380586902432477836 0.50632385995110551491 0.25996838063651106498 0.42971920796392282238 -0.22927544540652100569 0.058481512459443843466 0.021142240077264663062 -0.12795457725613457178 0.24977461280702825674 -0.18690980157453432331 -0.53762591968274586751 -0.44045926355137837716 -0.31311128012358507844 0.079593763053117716488 -0.17612817223315824222 0.43097739111352995067 -0.18683353886172682223 0.16425417533698144501 0.050967511727153441692 0.7436870257947405527 -0.29578579857926301511 -0.36755896168080426367 -0.042285147112179691242 0.47489723908461239921 -0.067131593628559804721 0.20992648210587072644 -0.39835055551016090858 0.33466496412997609466 -0.18182773547889743981;0.63019204204625256516 -0.22849383605898518357 -0.070061492317733287116 0.49152357554047704236 -1.4029486588045907425 -0.46907956653576338146 -0.42109487277979346098 -0.54385847778034346334 -0.23409812380967906509 0.45799515605201734925 0.36785273438846255445 0.036674763042280879299 0.43490328441649450664 -0.21865009017652553558 -0.046605806846795362619 -0.53898604644443703915 -0.74526785115265770099 -0.48593822486769094882 -0.20880511032822540729 -0.36498021047190587485 -0.033775510648877564435 0.12456594379745028478 -0.13278141223698985973 -0.54579312315639927355 0.26978269391965109536 -0.1136011569292354334 0.088496047134643859855 0.37400497859384396193 -0.16273540082808510276 0.76985941152443115865;0.52618656870880475829 -0.54928963411958253715 0.31900365827277649133 0.12594316999423910142 0.17273030333703934258 0.0033495256710558644203 0.30932570531657677115 -0.11109327232400628693 -0.16365916135702962553 -0.39712368525145941289 -0.57372427073807596987 -0.25317867880282429738 -0.45681955635306270569 0.15235039530758098802 -0.4280652260834607592 0.25400626669775705491 0.22889799697518797839 -1.1583514736184206129 -0.84711014799107597906 0.098949180190483321251 -0.010782130479011678426 -0.02011905542851382192 -0.40305064620259839492 -0.18465544301807243488 0.028434725169755658419 0.49408415273438277548 -0.37122129742856119217 0.4116225751274307787 -0.27372637058968912749 0.53329153585337396581;-0.16633211830430119638 0.32991365181255932759 1.0000652698969234589 0.013166379463003904937 0.25686132115042081825 -0.057300949556131589024 0.87418762558575324739 0.026052603574794070546 0.22552265275997573912 0.39565620478953700268 0.15788697678519425516 -0.5670062344239187313 0.77018850745161759885 0.84659265384839876845 1.2172844216313254417 -0.11694247177458685916 1.0266289875610112858 0.41303186372933453629 0.74674501821099636967 -0.074365933812977444428 0.15525867609605059627 -0.33752577574396624982 -0.82553781474548371122 0.84931720995759296677 0.17609731772831885288 -0.90901997804495471289 0.50779319203975559205 -0.83369468833649096329 0.0088731950264089344549 0.078740949625034944992;0.90153361592175085271 0.078033729029917828224 -0.13787604037316836258 0.39923838045558834864 -0.020745674194109234301 -0.093972976986967993174 -0.077664972164923934406 -0.25049902837696452762 -0.047866325723187552743 -0.29350162449391170583 0.13008884471210596656 0.41365612108101390554 -0.69470457832664744835 -0.46929820989444742674 -0.036407358356068973693 -0.7924938868384716617 0.29214971931705852048 -0.0025635630673860629151 -0.15999547673831218675 0.11563952256397000418 -0.35763676879477812598 -0.3346631672427937132 0.091844775703351116958 0.030390208194691271937 -0.16637982015261812418 0.38775158284361244121 0.12673989238225211174 0.54329358504007330133 -0.40220027736465807644 -0.51636694485009315692;0.068412582799267290801 -0.10395028249816054555 0.26554814256684283746 -0.13385951423468045451 -0.59571411579030020267 -0.14838798293194641054 0.2698060262438128909 -0.22514055848898217871 -0.1108454817026010425 -0.03249414444597246987 -0.67885528708276365872 -0.10283752825098131656 0.078866073006013087676 -0.0016728207867456506078 0.492980053034716037 0.72103408164882809572 -0.64831356625743430389 0.24312731737757117823 -0.26755810611578406677 -0.06987913429753687955 -0.043344186790648123186 -0.39006095265967727359 0.20397493707954270081 -0.29041033801749205345 0.24344995981265571405 -0.53426664363666542812 0.22569912414946052825 0.67358849776531581899 -0.3420340729327681073 0.042384262346315072778;-0.26742032353528377708 -0.27541169538515791748 0.45459213500869255498 -0.50801238525219805098 0.20803296991578235087 0.32167296233694808727 0.33702105694968154115 0.043032803977419438468 0.4546355004686749024 -0.6858570647840611656 0.57098990970682395041 0.25518370761478831321 0.21959198845052468885 0.010299904354930364964 -0.0076793886425852659131 0.28548918611526019617 0.073636642734239685093 -0.53156187227236340043 -0.38573913577013752008 0.16012777664628990992 0.12521438051551384008 0.41158455274119676304 -0.31184981716219567494 -0.056640103464305313619 -0.57489275648343773906 0.67488995220425196386 -0.1364863756435283626 0.14354736924722993896 0.11549710815718422507 -0.35045969591442710511;-0.15486495641898362363 -0.76258085288116317368 0.27658606866029866733 -1.0593375703671552568 1.1082196337086069882 -0.45530022340848175233 0.58050125364914428872 0.44682042998205395712 1.023311172739819197 -0.31562390072030227506 0.65996008861640220733 0.43468869999134579496 0.34341155405718620797 0.15472854179535089347 0.082097458015797847297 0.025841369965202566628 -0.85077264277579034424 0.26491121127493760801 -0.13426887464876835487 -0.18618806322822270549 -0.59967414037146193895 -0.54502638081446674878 -0.47380242140330847134 0.27321890066701348676 0.29643611718177265146 0.080092330351551832557 -0.5717089491846016891 -0.24898699885755357086 -0.037439128396246752761 0.013174173044601158664;0.61511007844611009787 -0.24402012378502727086 -0.22529115019387613295 -0.72357428963180969728 0.66529026897332665946 -0.1536422820339818085 -0.073293898242750565863 0.48562125871042788328 0.69631739822602556345 -0.15002729709849607631 0.33186126950588473727 0.26988174278951504448 -0.040848220732804883804 -0.12145446822541727805 0.12800111829881938297 0.1579250940400774017 -0.2458623406153526425 0.39928286558875974999 -0.13692976365238260228 0.4871838623826493242 0.22275615360505107709 -0.6173323255723499603 -0.38084003887265066579 -0.017346957165044491617 0.0013253146541875268118 -0.40108483682641749013 0.49263380084025321093 -0.44148177439754787699 0.059126352105963460304 -0.31709200641993501302;0.173053472193777369 0.61129556923253536294 0.097561854241303619384 -1.0897386206684096521 1.2298036715911535577 -0.0053061255662337891009 0.15118314057109824566 0.46298879777066992425 0.86525704983157492034 -0.76886093246576625582 0.41452048393725376974 0.32811749327917921759 0.58276799641238186123 0.21729397317543613921 0.14969630729030311556 0.061968764385711312681 -0.58074016583097753941 0.5364817310609981682 0.061616317062535529148 0.042147625580121299282 -0.56163782256082750433 0.45546804486962283454 -0.054625629739223274806 0.38640940246341226594 -0.64663752279028074987 0.054939134198950578813 -0.65258212093990108915 0.23829449595429208486 -0.20757491028769348151 -0.072395197666398361513];

% Layer 3
b3 = -0.59829278988268641548;
LW3_2 = [-1.2951673024513254973 0.99328667859470032298 1.2149679349155562136 1.3304964462927593161 -1.5423239084004221056 -2.0696762874837340007 -1.1554726717544321346 0.74684324860297091053 -1.7270596903633190777 1.2400520482966934654];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 2.03830178109291;
y1_step1.xoffset = 7.47181919175367;

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