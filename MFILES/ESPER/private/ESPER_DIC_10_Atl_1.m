function [Y,Xf,Af] = ESPER_DIC_10_Atl_1(X,~,~)
%ESPER_DIC_10_ATL_1 neural network simulation function.
%
% Auto-generated by MATLAB, 15-Feb-2021 18:30:13.
% 
% [Y] = ESPER_DIC_10_Atl_1(X,~,~) takes these arguments:
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
x1_step1.xoffset = [-0.999999996490807;-0.999999999390765;-40;0;9.1619;-0.2178;-0.12];
x1_step1.gain = [1.0000000017546;1.00004398834998;0.0153846153846154;0.000309885342423303;0.0656325902039552;0.0461352500991908;0.0145655815308426];
x1_step1.ymin = -1;

% Layer 1
b1 = [-1.3888396789565302658;2.0046621843007943475;-7.2684289248249971394;-4.6690070517061936783;-1.1572525294404067964;0.044757011587681153064;-2.8714209178221539176;1.3042888439185340577;-1.3293763241885190318;7.0099843127344731997;0.15021262140067129232;1.6404118420931232958;-1.4367889593503142542;-12.63634917474935726;-0.29412790703844299278;-3.9625275605300784854;-0.17806575629322882182;1.4672420094866900353;3.7841523699852217533;-3.276188328377795056;-1.1957299216046957113;2.4343981744897740249;8.1173963629641789908;2.9949913024664729377;-4.1220075809304734804;-1.5589845565015392204;-0.48861446100883570987;-7.0160182132449104841;1.0893614342850841403;6.4469399060085512687;-2.223845417557830384;-1.2793040038273593151;4.0664230917794261799;7.4066619122527832531;-0.87491057223059720549;-3.6444747512071931084;-2.9517644544009171703;-2.9836407563382927322;5.172086801448172011;3.5198543954285561597];
IW1_1 = [0.77993022343891149628 0.60529687235955986768 -2.2766756615477881098 -1.1377060419726652274 -1.9388576705781022103 -1.3221316845926043282 -0.11514892114496480535;-0.19280497130953161422 0.20234789371222436372 1.7171654376860467917 0.9742376175518075998 1.3799859805213299246 -0.04148903964290670543 1.3969272171329678578;-1.2814566298502922947 -4.7221012700383218075 3.6498157107764988005 -1.1206445429908711287 0.057747058915415656888 1.9915523645779982242 -3.0921584207765664409;2.283251291704569752 -1.9985555334588698262 0.0053623611745890953015 0.19280590161675667527 0.49637248245251047729 0.54647135314276085172 -1.0562572835402761484;0.89863024248958867091 0.82159756922567761528 -1.7723354071316439473 -0.73563285067956996066 -1.5736532962385687551 -1.0706306235346136901 -0.48105237786604149175;0.071966705714179360354 -0.98182482152069161518 -2.0201738966806672693 -1.2775549769197593886 -0.83588754655570574048 -1.5802695825529113893 0.8158044488595275201;-7.4975340776196031101 4.7292624673516545641 -5.208905435622342317 1.2410835057337710197 10.811086492598798614 -5.3250117577025948279 3.3883211765576293217;-0.33892663268389050524 -0.15310851251811888329 1.21876979070000524 0.57567754352711453247 -0.16064275031595742771 0.57270078696426218734 0.2889095932298068603;0.43273801058576477629 0.24284409066710227476 -1.425613048372760705 -0.68729303018381482371 0.28158065644720881027 -0.72549876096053200758 -0.40741389709575254185;1.202498563412273791 0.16891146229208423013 -2.1382235920986989441 -0.82192834393143043847 -11.604537643474122532 -1.244013992939874802 -0.67991308905891534309;-0.073219309475433758116 0.15648143585055110383 3.0596404505257308593 -0.42187099787387610794 0.42413586566181710191 -0.43041414009502465543 -0.11656407031488365267;3.4778323782140372344 0.65094637145702216596 8.4982592850196461143 1.5946869517884549872 -4.4259001536936342092 9.7928494368890746102 6.5936516586741840129;1.2971674164623854253 -0.052497173277100815658 0.78371576716728552636 1.5577775362804799908 2.497854850417801309 0.51035037201889699165 1.1804198564589289955;-0.38330546733347486965 -0.28156369333010450307 -0.12143059214366690168 -7.8562285209217330717 -0.47514671109392986326 -4.6428627284748893445 0.36790277472073545928;-3.3714274356407174515 -3.8059197550128143206 -1.1539535607844075837 1.7721108194174171047 2.7520024695646250201 1.0329959729384321498 5.2842921648607328322;1.8554742622809312635 0.60365773459566196557 0.78698964398283899335 -0.83001274745043429615 -1.3146046120258958467 0.8271174773288563209 -2.3872143653230732951;-2.9518323410622149261 3.0252755793890870883 -8.9611755467228793037 7.0707103823932815345 23.252931135918899486 -10.524003341326940486 11.242403304917733209;0.090027115538417967766 0.11227617154886533613 -0.59997261942954371872 -0.32816080743143916099 -2.5378108321022843086 -0.29010372133000988937 -0.48574341034242346504;-7.5879409427471831151 9.3881431134053165977 -21.01520924804461643 -8.3357617480068988414 -7.7861234799058687273 -28.368368592846085363 2.5676738425828005141;0.44201652581091266381 -0.17460371787639025754 1.0058628595670608963 1.6713585285384642276 8.1222931670844484842 0.36214961353863228677 1.8786089969900350294;0.6003677507905775812 -0.41350879003167523296 -2.2873610632733845094 0.41178322529436117705 0.45295722910750479961 0.59889586816740114639 0.19967389540079338861;0.75624136676567177417 -0.0034986810845369096037 -0.4179063932908262391 -0.58599657226289536549 -5.0683674381779377782 -0.37795692251576712417 -0.66408865350892820612;0.35024498956698035812 -1.6459795272621808859 6.3003296120285714466 8.5322301279760708326 -10.530678281671043806 8.2293502689477318057 -10.891133623906071648;-0.32464297793318269703 0.35342903828260852084 1.6957926598025583331 0.89771393344909489187 -0.45931092471630990959 0.11798117467449317775 0.65471345353467458317;-3.956195814294976465 1.8484021797060044001 1.177818872203000522 0.36904893689561180103 -0.20890887311306011176 1.7589280184226177006 -2.6866541419187899287;-12.859111722848624026 6.6236179089294866529 -12.159465645595235017 1.581048253673382975 7.8106166997178521072 -7.34418814134173914 8.1842171403908512417;3.1917344347086240575 -3.1255958292248511299 -0.16901815266807443394 -1.4469065293322298515 -9.3728078027986185816 -0.098469269812257217978 -2.736478625792679864;-0.584567253837150691 3.2980660916437036789 1.2348722241416036915 0.18634813156123480882 1.8825958463255134667 1.8293611121401209818 -1.7240188168819030601;-1.61793134618679324 1.8289426049657493412 2.5265150080756071382 0.9794954807118474438 6.9151126997391987317 -1.4492362776638869626 0.0092827586489612201603;-1.5451393097959138512 0.52409050249104338626 -5.9218906670929518654 0.89915495319831273413 -8.4469517433070304691 -0.35484813391527186832 -0.27194884502193772402;0.35622871659229637942 -0.06012354556745810874 -1.5309660548000674218 -0.78070704849912098133 0.29221964830666785629 -0.56600133561230447921 -0.25797983601802093601;-0.48961645114185364802 -0.33688063406343565775 -1.1882166352197505166 -1.1655112362183979169 2.3448188956781370607 -1.2603120643148286728 -0.44816363030324468664;-1.0620717472176008922 -0.43186243497136500569 -0.85849902694717472595 0.85398991617104025309 0.3371066936407987602 -0.9633630746117907151 2.3574225413618412617;1.2903595837197168184 4.7458899828215646011 -3.6914709370056351112 1.1291135283333684214 -0.21740524969988392345 -2.0138251388696524202 3.0915206093606442117;-0.15369540564220854351 0.6917880812974748217 4.6011273500149902915 1.7309597983839335456 2.503021149752283403 -1.0072918250637654669 1.2604044423200420066;-1.6356028382380147779 3.7720739462087600202 -2.3948806555611414915 -0.38680722654641569447 6.3774487208771484248 -0.10493492864244625962 -1.4618547665330752228;0.67755810414554584753 -0.46448857134650267664 0.099819368674388370133 -0.73111827498598191966 2.421140663864072895 0.26140283976911937724 -0.95501062999603325565;-2.8548838937683824746 -6.5890507780547018513 0.28755227121890725162 -0.13812546427649952108 -0.32244358857232990889 -0.010107347284078317326 2.3896955465385172346;2.3744150080390244995 -4.4902129819472111905 6.4132710597361599625 2.0285489338780160473 -1.3474461873730323447 4.1362519899921705857 -5.0861014927889485193;1.5180829196379026502 12.219950087576821218 -10.859283854360469945 -0.99934513479412045989 0.44264106285581583666 -0.94039494707853277244 0.061305004056781478572];

% Layer 2
b2 = 0.84239988120369069247;
LW2_1 = [0.57506323369274092716 0.78887808214251498029 2.0937782107889391092 0.25250453795140409552 -0.74302584030734242226 -0.12086409784279929236 -0.010241224435156379841 -5.8981073306159554903 -2.5182358555058228511 -0.11258246381110112444 -0.21719792571365539868 0.0046019710488485890756 -0.072153912449253856942 -0.20565978358895461997 -0.002839967766484370864 0.59996376438462983849 -0.0014288793851956979077 -1.0395159634998694109 -0.0017977432945569073782 -0.15907755468895340423 -0.2204352476223263857 0.51005446282512800771 0.0067738070173349890918 -1.916207728564182311 -0.089040007770363915895 0.0095230972578096775849 -0.046064309926230888326 0.29725695924958778216 0.051853127701314040421 -0.03331101794533213728 -3.8321314746317831634 0.11616787229968272843 0.76750489344134464886 2.0713861651346912041 -0.025805683454172931302 0.052498748223368439658 -0.27200675321490008773 0.015028008707853564643 0.038727244798986294738 -0.0076546627676070460797];

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
