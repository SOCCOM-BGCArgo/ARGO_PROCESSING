function [Y,Xf,Af] = ESPER_silicate_12_Other_3(X,~,~)
%ESPER_SILICATE_12_OTHER_3 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:42.
% 
% [Y] = ESPER_silicate_12_Other_3(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999984769129;-0.999998793562952;-77.961;0;29.715;-0.9];
x1_step1.gain = [1.00003298043486;1.00000060321889;0.0138606862425759;0.000266418009857466;0.193068829037552;0.0415843642790311];
x1_step1.ymin = -1;

% Layer 1
b1 = [-18.550190761220523683;0.10164031683280350116;20.825100808731271229;-4.2186385312122851587;1.5855164825583802912;2.9346020398248353622;-1.1276357396754619344;0.012883707772801821986;-0.43204127924738811029;0.57356271514922652077;-1.8613305866194607141;0.59570651545451303033;-0.75368191738600498475;-0.25483795131701669767;-0.35345772318504720122;-0.43573494193295930765;0.56580326812667902647;-0.68342276681058955479;-3.2859691942444677615;0.17103045115927933972;1.3571016279874015265;1.5228664914568559308;-2.4191963796428792577;-0.5645019435617236736;-1.197319356269606283];
IW1_1 = [3.5462245798827205512 -8.3730603099434404157 -7.8314049065389257365 -3.1193080950897980053 -0.85601582192787706305 -0.79990141120671243158;-0.00058396292211011628127 -0.15625242968460406612 -0.0906199270656714434 0.10619300531446874969 -2.8412083143498647608 0.071042971143866900929;13.829593032113772111 -24.498699309934902857 -13.278807960958781464 7.4258792996828262289 1.0092633668242572131 9.7011995460586479822;0.57033196290334320366 -0.019280442692728036458 -4.5250080953769975523 -1.0602783298599101069 1.6245294240060543789 -0.072687249790445093711;0.14182011453041970528 -0.099037001159926515981 0.13335676407518742836 0.50390884821096315616 -1.0467397288026660807 0.18654494918061123587;-2.0159885322123645146 -0.79703065954259311177 -2.0112737018600084227 0.950896361519732225 3.3529649070601741379 0.64954963262470255536;-0.17963719544204012624 -0.44539354768884936719 -0.88932479736192571185 0.11036443529488558501 -2.4037133173595544022 -0.43647338414166397813;-0.1893704591188382913 -0.41404009334397701148 -1.092395964564963684 -0.61504150209112740111 -1.7336249081356427393 -0.6761371769754601635;-0.11791034508083036991 0.087246302708337081233 -0.72990718788858210964 0.8661347343035132873 -1.1112004307085272536 -0.48140674535365268749;0.86505969135542004622 -0.31434091974595035346 0.17597977227563108005 1.886937020799372311 -3.9910952135547903552 -1.3489077199345456126;1.1007791446697903215 0.3559456641571869806 1.3675953173331716339 -1.4161381869904585695 1.6864714134317371208 0.35019146290252117781;-0.6924653305600056985 0.35699882308059094349 0.06521899342863878013 0.71444445770122677608 -3.7148497108065008376 -0.35746413385941572249;-0.32824012949858916155 0.27799455372992848412 -0.41031999594866058745 -0.24223929050776749072 -0.40091086941733133919 0.1744792052845340613;-0.14149004074476062565 -0.73105817881898860211 3.0043530277596763867 -1.5055269589675845943 -2.0006404764518577188 -1.087959978546452211;0.27450785495297558292 -1.0988828811029349719 -0.11539823674825985811 -0.11568092644359273535 -3.9632067812156117625 -0.82558178126347026016;0.19320157698109866051 0.73937857883630697664 1.0282864435341707576 -0.35190149796521164038 3.5585876656494463965 0.64578101585915004268;1.1420742892780726407 -0.16490328360428810339 1.232228556600069691 0.6955638591562341988 -0.96812597060107374869 0.1503676570967071302;-1.3791397091529460361 0.21658897682694552844 -1.3793758081892879552 -0.81502423062206941218 1.6535214803909361603 0.0065219822606739155657;-3.5126462890553411533 0.78233581181133327043 2.0613353828578202531 1.3738744823430157904 0.68816334596377304322 1.7805293160276898234;-0.6152013522342769658 1.5322010164884896088 -2.1015789885288582717 1.845320822195490873 -2.4560485946597379048 -0.7250249148394449783;-0.042302458467078481175 0.79717431672510186136 1.0344171754689959197 1.6257091869591866473 1.6386466730512145862 -0.99047411495619042121;1.0641276290744809963 0.72551412231506628991 0.45627661029572708085 0.32426465621281108653 1.0255563773073965717 -0.34286062445759329931;-1.5886201187321793427 -1.2386143496781876205 0.80840829006908798693 -1.1175517278521520659 2.4292295364394509427 1.4418283491492189796;-0.12576753169616847194 0.20258077835867621674 -0.86424356278375769058 1.0899135296196154066 0.36061603834758110265 -0.46144543165517765138;0.078070148609168957266 -0.048868188262251338605 -0.15778962130592716351 -0.89641532502797183213 -4.259564403142563016 0.21120942858022068167];

% Layer 2
b2 = [7.7015762027599645734;-4.2895676622866432126;9.2651474141398999507;13.677089131147846857;-0.26283438603039782633;4.3329368615314995949;-4.7569938203472137417;6.5391583494614735628;21.067306026535838726;28.201760470252899893;15.020306515773741651;-3.499470432968798228;-15.205649630840653685;-2.14627659416557659;22.749388526734328764];
LW2_1 = [1.9128470912225421507 6.6908272927089118554 -8.1070220268620918347 -0.29271168457080126712 -10.369929178389945434 10.672279008824972379 -5.6090800912499121367 10.713497923192310424 59.363069257606682072 -0.65456590569492034692 1.5074688525506172443 -23.083266095326418821 5.9101799396343146498 -10.35762068765049726 11.195344576714020945 2.4335611599743129929 -10.762585695926523144 -6.4752530993768813161 11.270558889322433771 -14.688025450755318957 5.5230383436665526631 -4.6767665909263982371 -6.054847217016147809 -4.9683017409660807928 5.2441638360072699854;0.31559642784694291739 -3.3161094211609261073 3.8096760012626380565 -13.690252108893210092 -14.140255585610367106 16.169691515506936952 -1.6882824842446624292 9.7187200375671469743 -8.322439141893823944 12.848254590455383095 15.262850369271605544 -9.6525484157044143529 3.2054460889764664344 6.0256660533165593563 2.0319024393805316109 11.003751097816616422 6.1705693582353697835 10.903346972924456537 -4.4664981232585549975 -5.1100081630300966395 -1.3690625889104186541 15.4556948733353412 5.4771411284579682999 5.4824621080842357301 -8.3083831548061031214;0.3061661837478091841 0.78076538800694361697 0.0038549180426293242624 -0.089634591618919859379 -7.7695299711149790056 1.0668120165319407988 -6.086823915269108376 0.89821530966991713907 20.777784872503847424 0.26170922236813876083 0.0035920793856717169769 -2.587090779726052947 2.2597938474896759331 0.75786644011459569636 1.7222644531446571481 2.8721203784657238778 -5.0069770762387584284 -4.3199845097710634079 -0.10455182406161905861 -0.42971713309339898412 1.27572696049840828 -2.1746488793154781227 -0.060601349220979448817 -13.76186514097008029 1.9388003585912831905;-2.936205465170147022 -4.5214803458756156473 -0.83387819099073956508 -5.1348369661486321647 27.607470254792126241 7.9023501315688315927 -7.1950509617422717668 -0.6999601950718919996 -57.405325547037392653 15.3576605749939894 5.518539805756955019 -0.26709830696889647061 33.736847847643375076 6.944672656113812792 -5.004330141546761368 -7.6373778086342385407 -5.5506243506507244589 -9.986788640290335195 26.967855957539491385 5.099399089033922472 -27.096433372171844667 -12.51281099739965974 -0.92457625300099866017 24.314180497081967047 -10.762913625985200738;-1.0308320600362064212 -13.984717876892256427 -9.5990221603859957611 -9.3195662080797010418 -10.758261615224951058 -2.6137604609718327886 -5.4072266826625936531 -14.35989847704708211 -7.6519065661006235857 1.8095481700504942868 1.1207782129405685634 -7.8193270992256636021 -30.447213794614906845 -0.14422646733916069661 -5.6977094466848061671 -18.614566149098386916 -14.9815455410166809 2.3960447563971820628 0.58177994118433629112 6.2483335559067763043 -3.6772431081217362703 2.9706166658678885284 -4.5955176434690363507 1.2400089675861289962 16.557679625136827894;-1.2492553960429466642 -14.68670906039438151 -18.767791719749759238 10.777794582408084523 2.2694899877541772071 -1.1016694133515256837 -10.909610273102133604 -17.099105443678634941 29.913407861368703777 1.6722759395451536335 -12.919522503786181389 -11.996228497604652929 -7.3751585140710496091 12.146471053748010505 -10.511235833889433167 -0.62104354765227731683 1.3414906158222359434 -6.2851700611874985469 -27.160956843983775855 -10.414529149715736622 6.7662399720375958623 -10.317612553773287232 5.3744850321662882564 7.5368825737288549504 15.735475144330587582;-3.3062494364606487274 1.1219364365366553127 15.210834677835467943 -1.8810354730633860321 -1.1968858823049757234 -3.6367753248867904503 -0.025261909930722570988 8.0917719233166085502 -43.679218765954516357 0.30905442160921070549 -2.2302374332779537625 11.524932656543851905 10.74582906210115496 -1.9794820464439606411 -5.9702506285207688208 -1.4992469839907005547 -2.5645888766628668343 -4.567646670567814482 3.914502966566299591 1.0750164548999872061 -4.6279047281568725936 0.91906104379729802734 4.7427202430986508119 12.605206736749668295 4.8702656032174731493;-8.8070698986913455997 3.1686342497928245976 -5.6929592221231555271 12.521626651631230231 -3.1738415349338162663 0.64676008430646891867 1.8328977688354866071 8.7732485937180371138 36.893914932945392593 -8.1892414051803203989 1.1134206081946869116 -2.9422300282343241129 -19.614276004586916002 3.6888127511880197673 -1.3104608591313080534 4.4800628088687268757 2.9625286390655949376 1.7504671862958907091 -7.8507766907735794248 3.9712625308951170489 7.0851923289627150027 -2.028605841621784478 8.4715648521088873224 -8.0239641622160977619 5.6561831917246276902;9.0187062630038337119 -11.519000374203534065 -13.742700634538390858 -19.108673037399519501 -15.446738216561470836 1.5342871539114111279 -13.587923877883364199 -20.287948225340059594 28.15255427841429281 26.63227181037272473 4.4736198889635474529 9.2834819333855840284 24.986940037545242177 1.1934666699700020587 -6.7368592277789183242 -6.0651863616998840101 4.7477291738929272569 16.007325963987781847 -11.769458828703681874 -4.2370396582801372887 -0.50666748351873058809 -1.6785716135050310349 8.6386777583820748561 5.6618465870359049319 6.5056448019948316386;5.0060000110512294924 -16.069225794569003796 -23.81598128255220459 6.5734383011408592523 26.477871453988324646 3.8057486499854658035 3.5777543642942015722 -5.2174715395344097857 11.51284346739197062 3.1832641930386751383 22.854400752407176611 -13.05744548072410538 -74.468829904183664326 38.348596724449357964 12.064716518594545747 -10.301125008909805203 -5.9692214986285305756 -13.951238085280790457 0.86060377642544816368 77.861444294287281309 1.1702617948433882145 8.6048098255250948796 22.145571306999421779 8.7261711383399447328 -2.8816885637546167942;3.0918589485794774063 20.243406242914183224 -0.036765747402455195569 2.729827116796407882 -14.029478024075418219 1.9167360846606955604 0.62796563929379245472 4.4531293351255669677 -10.018017093073034118 0.69728987437844069941 3.615729079215080688 -1.4040726223169406328 8.8968347350628800285 0.60777535299556650905 0.57468184882369666155 3.5346894862062061549 6.2866006306543589233 3.6664845937956362931 0.13351508848509480165 0.7565142712785052348 1.3579351016117258144 0.57202896800957969692 1.1956782186486276665 7.5287275327193174945 -8.5918089319693642381;6.2540501554187075683 -5.4261719032469963864 -0.80112330942856402327 -5.61699855229717393 -9.170337894388939759 2.3157263176035050023 -0.84427920683152934966 28.391648209200951669 18.601041965840881431 -5.5634926732835889851 5.1618464254419205162 6.7789705758473042252 -25.147112047251695799 3.6124118064697596786 0.67363801289936031669 13.801906362081524549 11.880293616599438167 2.3594278709800491356 5.0965894714263981768 7.4458533775911348584 -6.2110653066118164389 -1.6645681213295975276 -5.1776723271532123505 -13.230663114725240703 2.0787426749849582386;8.2761669024125730232 -4.9610437491227568074 9.5368767032873975609 -9.6104074498497666923 32.83651467391872103 -8.787574676668402418 -5.7416549330809161233 -19.954859984418931163 2.3570929607762813873 -12.910157246488401128 -5.0347989944916244909 10.185945026947903713 24.120680177287461277 -18.785953944382832503 14.817332389482769983 8.7575861100518839208 26.327749597147164451 14.793980472704260976 -12.651607075423793702 1.2892266655639319506 28.192387509311291893 10.12102335950228138 3.3523027326579537544 -18.469909971383724923 16.041885297984308067;-0.95306626697177676544 -0.12194109088122061235 18.783207109394957968 2.1564989319617260399 -22.27646339327917957 5.6767957676084659369 3.3841493904804882398 10.886489992637827129 -13.48003601095873627 5.9031235298327411343 0.26121034074511073619 2.886611752378115181 8.5680494224271637194 -5.3623675218314978963 -3.9736613065365289366 0.58390605893513447544 1.9915619097933747739 -5.9286238777099029562 5.3316224270418342002 -0.72521621327929319811 -0.87906851117668205475 -0.19994502855872892111 4.8664205197970984429 -12.528769482978335148 -14.157458119832426746;-23.425702994474441709 22.405048198293339823 -13.296962030743008754 -5.7200603523888675639 -16.124871154598029221 -16.299228189619714868 10.955384284759567493 -11.730688557064368638 11.823973082729940742 21.962474584937197619 -4.079664941872219508 -3.8041830365261706248 40.410159372887555662 -16.627051411635431322 -27.346976020792247652 -0.54819949151153446465 -28.84245575119999927 -25.335777814571091682 -21.546264209668635203 -11.981550883109679972 -12.562504894502071195 -23.143301260297178601 12.88498953911703282 -15.950304882719878563 6.7079432142027979324];

% Layer 3
b3 = -0.29828649991374978079;
LW3_2 = [-0.059425277452048834437 -0.01214073609270005516 -0.46575960239377972583 0.008862742599988378267 -0.035221408341092746985 0.017824742925247504671 0.0078313399723019799098 -0.023429806294904789027 -0.010609948729690841726 0.0012700231380979155969 0.21806853773130321916 -0.028236625257220787644 0.0087657471322404017372 -0.016012977094358869096 -0.0089879316429739095767];

% Output 1
y1_step1.ymin = -1;
y1_step1.gain = 0.00819772922900357;
y1_step1.xoffset = -2.02;

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
