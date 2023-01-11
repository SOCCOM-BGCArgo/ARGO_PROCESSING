function [Y,Xf,Af] = ESPER_DIC_11_Atl_2(X,~,~)
%ESPER_DIC_11_ATL_2 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:14.
% 
% [Y] = ESPER_DIC_11_Atl_2(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999999390765;-40;0;9.1619;-177.232103709389;-0.12];
x1_step1.gain = [1.0000000017546;1.00004398834998;0.0153846153846154;0.000309885342423303;0.0656325902039552;0.00425400504498441;0.0145655815308426];
x1_step1.ymin = -1;

% Layer 1
b1 = [0.2850839600002229024;-5.4764793003597285193;-3.2159774836543459919;3.3138288576541303776;-1.6102538916133728453;1.0779375806436146679;1.5983216455767088604;-2.2426839763207442324;0.082209318423575331036;-1.5184792067867982635;0.65950544206201910935;-0.5892460404973396404;-1.6572481203439586306;2.6423662419308406868;-2.4558318215929171302;-0.97585279507572442181;0.95734211438662641935;-1.4800235425985563964;-4.4443044978273107404;-3.0973952331012943873];
IW1_1 = [0.1655516001124180403 -1.5096469065921087616 1.3892150670714367156 0.19970166169329739425 -0.53787948147913577923 -0.10265632787641999768 1.9549003093216077431;0.47210425942820943845 -3.2762402506252992218 0.48654261607134829948 0.84960877441861859527 3.088877110729981279 -0.49933846907150547523 1.0589036133183666877;2.1982306784659124688 -1.5171240791596276587 4.1428771493040272844 -3.1981486875133962045 1.6493754888365479427 0.080265365599378038497 3.2543023388832539311;1.5496297972761563511 -2.7300545968985598932 2.0236144567642297964 5.7977831901281540183 0.84720337939868861632 2.219483250210179559 -3.4622858354536640846;-0.48125609542118602535 -0.24550397926986761155 0.82871722761343902253 0.69395804023609641842 2.4773794516816898259 -0.30905179516525033145 -0.3206951200884293196;-0.51632575448623152514 0.2907028284930951445 -0.58330904881620593283 -1.168078738717688525 -2.5932694984987332276 0.24327210349359565389 0.3516648563013711204;-0.21187149534586649513 0.3798681544047490255 -0.14344833723470168385 -0.79475103500053390082 -0.81097134529944592263 0.27420119381860907914 0.92274166711386040873;0.030960337984934957478 -0.027506051487058230665 0.15023884668554604582 -0.11508822088296555286 2.388505306411226492 0.39894665716057536864 0.42051247917644496432;-1.0947955051932878323 1.3964762493102997887 -0.39834593900014586554 0.9762765146743558331 0.011409066668914891185 -1.8345339703606240445 -0.29313521039162615534;-1.0274511216191539731 1.5227149873855152507 1.5195601599906036761 -0.10805463418825073107 0.74439021740577904218 2.1319202790966880201 -1.0526347498169332795;-0.30075456188700122606 0.34216018289222221993 0.9314324565614066076 -0.32683800778419169752 0.23123691122921344321 -1.1343978269131158587 -0.85178675635195777449;0.32019705285109989479 0.72994420444436536322 -0.92404253573475636507 -0.67113020386687327612 -1.3636681593240695953 -1.0945846389767694973 -0.9922787571519491312;-0.27713268202120822092 -0.2632876143927572743 -2.2888399300099977474 2.6253461359589573831 3.2695837621320822386 1.4868114216827215124 -1.1429928186943536517;-0.48682034258822903805 -0.33483664306238747344 -0.21692263974522679737 -1.5619730863618088712 -5.0526359005248293599 -0.42890694428175207831 0.30022052029971035836;0.31090155703587391001 -0.001936742515241406809 0.39350890991445375278 -0.049426938487130860078 3.2119316373847359714 -0.15792908497736232731 -0.39295997161969753453;-0.5522879993976334001 -0.50205718201133120449 -0.053383164685383961701 -0.19632161306828202485 1.7653202294257015481 -0.026879026149830259618 0.4691501205761883786;-0.72042232201545075565 1.2807411859522581654 0.48362730124629310646 -0.020840065319075035122 -1.545044217042811141 0.88951931974985654072 -0.65155354024280320857;-0.14447095286197930708 -0.058963475015533936552 -0.026993427698422032407 -0.4161360856005353881 1.9172916925939040844 -0.21971018398939212712 -0.07335216301198564659;0.46454746834266663624 -0.31591537139240344345 0.33093627440884904045 -1.1093230840564589634 -0.92881118155552233606 -0.18321200472044946173 -1.6283695720380446659;-0.12882911677108108139 -0.083050529722858590365 0.60906567841090164528 -0.50975864600815579131 4.6614155278931042048 0.036555354564501756076 0.29101251624701673926];

% Layer 2
b2 = [75.892405234379324952;45.698532332422622915;-45.441034334629890168;-38.047185668063988828;36.729913807561054284;-14.871696621725067899;12.957843871203559516;-40.868717659820759991;-18.405935213533346939;-8.0719890132859095644;17.238097802070896591;-7.744098383093918514;-8.3105114251472418374;-23.622318635312993251;-18.962077660260927559;16.243348544051386284;-20.203067252988784475;20.138099221164893038;3.6940838106982987199;18.799551603554757406];
LW2_1 = [-30.36249518468235209 -37.844741737012704164 -9.3746503462666019857 28.472560011380934952 13.049195096220211454 -22.589942930487911354 8.346290978581640374 85.17873388635058518 27.60933795207965602 -25.737140627302402152 74.780848526359534389 39.171205579476136904 6.1261775638515088005 -24.079177742576931109 -40.55273107793865961 31.859173796725539063 0.47919571931297805234 -35.311871421427852624 121.59047076691211942 2.3901570343206195446;-3.5867287417172231123 7.5684409733935069653 6.2457362483098748385 -26.491666673412506583 -22.825302441033567646 -45.978047162530600644 -21.557889775338470173 -9.2437837578905046598 -15.005004838937223255 16.948751972235598373 -26.034742680857643649 21.698519115817244796 3.0144902445003158853 18.336358328321413325 -8.9574739082785956867 -16.189134123190221715 -4.103565268763833096 -0.61971242416088234695 5.3566129595381397621 18.25490169198870305;19.252767617839769088 0.21230145058703908667 -11.665116622773322064 31.863871347097681763 3.8334986977135616648 -11.695359664578061754 -1.4359982707973351079 21.243134895322292266 -10.043537915072464628 8.2540264096903044333 -3.0615830505249581606 -1.8839666763578652997 3.3126147489886301756 15.874995364454433044 9.9017184048348934056 -22.79103108291157298 -5.3323738223664314262 16.973657631780223909 -22.91248551232051156 23.743094415374947914;-40.032845381174965382 -5.6691768832542654977 -7.9345389357965094845 44.779017449917574822 9.0747900600932958071 -0.26462685462273838333 -21.066867575525442646 15.753835073508588138 -31.029323062378828979 -3.3781701254310823046 6.4049498683327765747 13.753397670727403934 -13.011996701257025677 13.069564050377108444 2.6653776017044603641 7.9622837922723714854 8.8974768952216169993 12.251840946910093066 16.396452256693514471 -23.318445758640589105;-7.2075499379839769176 8.9057993765250564877 0.3513041591825650789 10.882423640502672413 -2.5584270652488161346 4.0937218801150407899 0.062277466597059814857 -11.503216274943545372 -2.1992370743470552874 4.0393600775492135568 -10.24429211030911091 4.0527118820764806983 3.3381410199808998485 -1.0681697038450050918 -2.1200153844842506423 2.3332445115269191405 0.11662706165356367949 -5.105768299350112116 38.746471305424520892 10.694821166156083692;-0.2137011806763182975 -1.1014748561140201844 -1.8837186078534746692 8.6264315283213210961 -2.0879394894810827665 -0.06473335095637003378 5.5978436887166074243 -4.6145955933053404152 2.9674132985329881862 -6.4626991119692522858 1.9749278383664836678 -1.7174260438157711306 -1.8347708354071157189 -0.8748393561576705002 -2.145093924098548932 2.771772112602811422 6.447067651514338138 -1.9782836342673613128 1.4923620199184999269 6.4396013081262699629;54.278875834811451284 44.68631409139198496 -16.278300786235938347 31.990946335478927409 -33.432310045145221977 55.191137360526433042 -40.087153561859416584 -17.46557457817012704 -42.093503922507778725 25.596335058337881918 34.195155659481550003 -49.545901368126784803 -13.425963351768077203 -16.420421256278416422 -17.20362360562725712 -49.92199428885152912 5.236720980736396136 39.906995756730751168 3.8556308189315058321 -24.959613694944344786;-13.439386982699200956 2.6413030301707034475 -17.616874994756567929 -24.707657467483969072 6.2118886154820645018 -13.654016017098477676 27.144675571414850879 26.045237646035225509 7.5267808759088108417 21.753569115587083616 14.757964675088642181 -14.488096568346715998 -18.288888504494316578 -18.067499546307914216 14.718333932334482839 -12.57721609093352555 -40.163704203701293238 35.931639841629881005 -18.970643401039701814 -23.238121641272936557;10.419620652902018776 11.937384976860506924 -23.043491483602313252 25.089871649176657797 -8.3278896109012592319 -11.423917708572924212 -2.9465069289010283882 25.950217442980061833 -22.895164297486772398 13.066673765040658139 12.208009602613245193 1.16983936335404759 -14.843363377289859883 -2.5887119523829120737 -4.750368508295154335 -27.464248124078114444 1.8033756607816129325 -26.069808355799612087 -6.0816022470963426727 38.057432061058143802;-31.360800786771736881 20.891943859911329184 1.6246341495502221175 64.995175337420505457 4.434104679141381844 -27.896572410592995084 9.9668716221917943443 43.638669185181896637 -17.238124673836754397 -7.6830458230180784795 12.8433292844974698 16.931517137893255409 0.89352572182081346774 8.1656272290067644093 -5.6793092455158564746 19.382631539370450469 4.3753946978971542947 -35.93555573697878458 58.101270117879437294 15.048481390203056307;12.72880167590486522 3.9610632078696124125 -1.3763644139582138504 -2.2306535122746780075 4.5237662189660516887 4.245449748954698066 6.9793283894167332804 16.096526512090431282 -5.6487440713021097594 -5.4216245212478515469 16.330719215610596251 1.8481060030063660715 2.4527649412707193122 -0.50182566210054180278 3.7283272934210973659 -3.7512456263797151479 4.8535950676171282225 -22.133699294738779173 7.1220912783450849304 23.539799822040102129;22.386593148624267968 22.244226901569696508 -29.088808903432127551 10.102815882102275324 -43.432466401842198422 -6.5871815760503347192 -40.480625652185956653 -28.106962821619912773 13.103234275209615589 8.0748762181015774075 -12.380938837962656152 -22.245935953944030672 -43.312023142339569404 23.664476852863980838 -1.8878442091959157789 -38.673031399047943069 -32.152282254901443537 -5.0698835418312029333 6.035769420669810259 -2.0129077032864670649;-4.3626254701170266159 1.419331997979452975 -5.5180246211176813276 -1.7984638879578347392 3.70585788596945509 -7.1984495688491403342 6.6479674182571786289 5.7698851961336021787 -9.2027871042984301653 0.80285262765584131106 0.039845913162985080458 5.3811413639008094378 -3.2622880773656999231 2.5820105219385980355 1.8605870861049378373 6.7060394128012754678 -2.8905176083555188526 -12.279018404958936372 -6.9857323596907114904 5.1912382306038695745;9.558528916678532994 -3.7585210624004847801 -3.173748994888921704 -13.625585869272004658 6.3257716979597073959 10.379055274484686677 1.63987528807358518 -6.8158829886842955403 -1.4272428288291976184 4.3304406423608732268 -5.658581332424733823 -2.6105416364047790978 -8.4998199541333026019 -11.43591216260856136 1.4241247943529848907 3.2180410842389726156 -1.85463434338249189 3.4919730898329257585 -44.348769231800751811 -9.217185223121260762;6.423485963469103055 -11.258924441716237297 -3.2584207287850115442 -20.871400284387561896 -4.7085308043345834861 2.0811007386644972961 -15.26890930266705837 -26.257679137715424389 -1.2783824188363119845 5.5659424529828473993 -23.429584606133722957 6.0915702896987848902 -4.0127699926693072285 -6.3381757991479812375 2.1248112921695643251 6.1569117069835721523 -1.8306397623874202107 -2.8486345357500724162 -38.069846377992789144 25.2684021467066664;14.70843949769237291 -17.197012626750726838 0.12803664112749085779 -3.3739408047553185455 -6.3835219780448673177 -21.726681382729797321 -5.4125922881854480195 -14.232432785150708199 1.2758562876167593991 0.7450322709581144176 -21.839541467628734495 12.688961650115711421 4.9429510397879568018 5.2001722847752649415 -2.1238392938419572431 4.511837035040382915 5.8825568301199506749 8.6787248953144686681 3.1879748161760406155 -2.5131235379986587297;23.594850152777457453 3.5683111130206661699 -11.636898664357758548 0.63613404056973354983 8.7816990241098071834 -5.9934114378442115978 -21.685928449412706698 38.870868940317947704 -32.655707299395452026 0.014041906193202986428 21.54767419978930576 5.3106682758509951015 -10.246014043829537954 -13.703583156533373 -4.7980776362107686595 -41.284951416485789366 -0.43938364096431037353 -11.331276904539627992 1.6517889001437142671 -14.083003896644212816;22.269444415103649249 23.967294752885809572 -11.772574805908991635 -11.839513466974796785 -0.5553566804414533653 -26.192931241624776817 18.161492792276149544 3.0322311020831635986 -6.339908516218314638 8.4434539711447857968 -13.356109348088393673 19.950227668245407386 -5.2081028195124083879 -3.104278074280415467 10.16524134332444973 -4.7712709552004799463 -15.992016006826196062 -11.762790501933801934 -12.600140139911435 12.041962777120774675;-0.24915589666444881711 -0.30357718184458426514 -0.1163319559312987056 -0.38268879709233633024 0.55857266783735093352 1.1340554507206621615 1.8076787768151294511 -3.3118295031236422155 -0.26648037504122301655 0.52301074447562667658 0.24348773138924401094 0.3146913763500179817 -0.27680861223282193784 -0.79514531525085074826 3.1776696636071899427 2.1784155461299801182 -0.99541389212428366129 -2.991526587394148784 8.275305337612547163 -1.0166301995025455174;-0.30581618716908914868 -13.174509904049109821 24.988621524762876902 -2.9209137083379004274 -10.465127831552351267 -4.5025664501874453904 6.2678904266298571812 2.1845040585389319254 16.576596593046403427 33.58494609603967973 9.3811960198202299921 16.186293568304932933 19.878298711491655837 -8.5398386925098836997 -3.4876705653947683849 -26.71819829938902302 -3.454912602377306019 39.842295622298586011 -17.906599876811117866 -6.7594759010279368283];

% Layer 3
b3 = -0.15162956321702733553;
LW3_2 = [-0.00181886699267958191 0.0010433812911831885813 -0.0046320344641960845 0.0020912612111151367318 -0.0082535821421151258498 -0.014530152893278722864 0.003161463422704781541 0.0057513920721562689434 0.0027174587542374366821 -0.0019078404448801984921 0.0048732166298890346873 -0.0074616140830718741261 0.018750354231359914192 -0.0065724618701617698244 0.0055753168387068464604 -0.011552359747334454629 -0.0019329610721393775553 0.0036408563912309829802 -0.97766295673835534874 -0.00079033900342678092311];

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
